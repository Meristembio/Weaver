import datetime
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from django.core.files.base import ContentFile

from inventory.custom.assembly_classification import assembly_metadata_from_classification
from inventory.custom.assembly_classification import classify_assembly_record
from inventory.models import Plasmid
from inventory.models import PlasmidType
from inventory.models import RestrictionBuffer
from inventory.models import RestrictionEnzyme
from inventory.models import RestrictionEnzymeBuffer
from inventory.models import Resistance


GENBANK_EXTENSIONS = {".gb", ".gbk", ".genbank"}
SITE_LABELS = {"BSAI", "BSMBI", "SAPI"}
DEFAULT_INTENDED_USE = "Imported from GenBank"
DEFAULT_DESCRIPTION_PREFIX = "Imported from GenBank file"
TYPE_DEFINITIONS = (
    (0, "Insert"),
    (1, "Receiver"),
)
DEFAULT_RESTRICTION_ENZYMES = (
    {
        "name": "BsaI",
        "hf_version": False,
    },
    {
        "name": "BsaI",
        "hf_version": True,
        "link_datasheet": "https://www.neb.com/en/products/r3733-bsai-hf-v2",
        "description": "Default YTK enzyme seeded from the NEB BsaI-HF v2 datasheet.",
        "buffer_activities": {
            "NEB 1.1": 100,
            "NEB 2.1": 100,
            "NEB 3.1": 100,
            "NEB CutSmart": 100,
        },
    },
    {
        "name": "BsmBI",
        "hf_version": False,
        "link_datasheet": "https://www.neb.com/en/products/r0739-bsmbi-v2",
        "description": "Default YTK enzyme seeded from the NEB BsmBI v2 datasheet. NEBuffer r1.1 activity is reported by NEB as <10%.",
        "buffer_activities": {
            "NEB 1.1": 10,
            "NEB 2.1": 50,
            "NEB 3.1": 100,
            "NEB CutSmart": 25,
        },
    },
)
RESISTANCE_DEFINITIONS = {
    "AMP": {
        "name": "Ampicillin",
        "keywords": ("ampr", "ampicillin", "bla"),
    },
    "CLM": {
        "name": "Chloramphenicol",
        "keywords": ("camr", "chloramphenicol", "cat"),
    },
    "HYG": {
        "name": "Hygromycin",
        "keywords": ("hygromycin", "hygr", "hph"),
    },
    "KAN": {
        "name": "Kanamycin",
        "keywords": ("kanr", "kanamycin", "neor", "nptii"),
    },
    "SPE": {
        "name": "Spectinomycin",
        "keywords": ("specr", "spectinomycin", "aadr"),
    },
    "ZEO": {
        "name": "Zeocin",
        "keywords": ("zeocin", "zeocinr", "ble"),
    },
}
SUMMARY_EXCLUDED_LABELS = {
    "COLE1",
    "BSAI",
    "BSMBI",
    "SAPI",
}
SUMMARY_EXCLUDED_SUBSTRINGS = (
    "PROMOTER",
    "TERMINATOR",
    "RESISTANCE",
)


class GenBankImportError(ValueError):
    pass


def is_genbank_filename(filename):
    return Path(str(filename or "")).suffix.lower() in GENBANK_EXTENSIONS


def genbank_paths_from_dir(genbank_dir):
    root = Path(genbank_dir)
    if not root.exists():
        raise GenBankImportError(f"Directory not found: {genbank_dir}")
    if not root.is_dir():
        raise GenBankImportError(f"Not a directory: {genbank_dir}")

    paths = sorted(
        path for path in root.iterdir()
        if path.is_file() and path.suffix.lower() in GENBANK_EXTENSIONS
    )
    if not paths:
        raise GenBankImportError(f"No GenBank files found in: {genbank_dir}")

    return paths


def genbank_sources_from_dir(genbank_dir):
    return [
        {
            "filename": path.name,
            "stem": path.stem,
            "content": path.read_bytes(),
        }
        for path in genbank_paths_from_dir(genbank_dir)
    ]


def genbank_sources_from_uploaded_files(uploaded_files):
    uploads = sorted(
        list(uploaded_files or []),
        key=lambda uploaded_file: Path(str(getattr(uploaded_file, "name", "")).strip()).name.casefold(),
    )
    if not uploads:
        raise GenBankImportError("No GenBank files were provided.")

    sources = []
    for uploaded_file in uploads:
        filename = Path(str(getattr(uploaded_file, "name", "")).strip()).name
        if not filename:
            raise GenBankImportError("One uploaded file has no filename.")
        sources.append({
            "filename": filename,
            "stem": Path(filename).stem,
            "content": uploaded_file.read(),
        })

    return sources


def ensure_minimum_catalog():
    for type_id, type_name in TYPE_DEFINITIONS:
        PlasmidType.objects.get_or_create(id=type_id, defaults={"name": type_name})

    for definition in DEFAULT_RESTRICTION_ENZYMES:
        restriction_enzyme, _created = RestrictionEnzyme.objects.get_or_create(
            name=definition["name"],
            hf_version=definition["hf_version"],
        )
        update_fields = []
        for field_name in ("link_datasheet", "description"):
            value = definition.get(field_name)
            if value and getattr(restriction_enzyme, field_name) != value:
                setattr(restriction_enzyme, field_name, value)
                update_fields.append(field_name)
        if update_fields:
            restriction_enzyme.save(update_fields=update_fields)

        for buffer_name, activity_percent in definition.get("buffer_activities", {}).items():
            restriction_buffer, _created = RestrictionBuffer.objects.get_or_create(name=buffer_name)
            enzyme_buffer, created = RestrictionEnzymeBuffer.objects.get_or_create(
                restriction_enzyme=restriction_enzyme,
                buffer=restriction_buffer,
                defaults={"activity_percent": activity_percent},
            )
            if not created and enzyme_buffer.activity_percent != activity_percent:
                enzyme_buffer.activity_percent = activity_percent
                enzyme_buffer.save(update_fields=["activity_percent"])

    for code, definition in RESISTANCE_DEFINITIONS.items():
        Resistance.objects.get_or_create(
            three_letter_code=code,
            defaults={"name": definition["name"]},
        )


def record_label(feature):
    return (feature.qualifiers.get("label") or [""])[0].strip()


def record_labels(record):
    return [label for label in (record_label(feature) for feature in record.features) if label]


def normalized_site_label(label):
    text = str(label or "").strip().upper().replace("(1)", "")
    return text if text in SITE_LABELS else ""


def infer_resistance_codes(record):
    labels = [label.lower() for label in record_labels(record)]
    matches = set()

    for code, definition in RESISTANCE_DEFINITIONS.items():
        if any(keyword in label for keyword in definition["keywords"] for label in labels):
            matches.add(code)

    return sorted(matches)


def infer_ytk_type_and_level(record):
    sites = {
        normalized_site_label(label)
        for label in record_labels(record)
        if normalized_site_label(label)
    }
    if not sites:
        return None, None
    if sites == {"BSAI"}:
        return "Insert", 0
    if sites == {"BSMBI"}:
        return "Receiver", 0
    return None, None


def infer_assembly_type_and_level(record, project=None):
    classification = classify_assembly_record(
        record,
        standard_id=getattr(project, "assembly_standard", None) or "ytk",
        allow_legacy_fallback=True,
    )
    if not classification:
        return None, None, None
    return classification.model_type_name, classification.model_level, classification


def parse_created_on(record):
    raw_date = record.annotations.get("date")
    if not raw_date:
        return datetime.date.today()

    for fmt in ("%d-%b-%Y", "%d-%B-%Y"):
        try:
            return datetime.datetime.strptime(raw_date, fmt).date()
        except ValueError:
            continue

    return datetime.date.today()


def trimmed_name(name, max_length=50):
    name = str(name or "").strip()
    if not name:
        raise GenBankImportError("Encountered a GenBank record without a usable name.")
    return name[:max_length]


def feature_summary(record, limit=4):
    summary = []
    for label in record_labels(record):
        upper_label = label.upper()
        if upper_label in SUMMARY_EXCLUDED_LABELS:
            continue
        if normalized_site_label(label):
            continue
        if any(keyword in upper_label for keyword in SUMMARY_EXCLUDED_SUBSTRINGS):
            continue
        if any(keyword in upper_label for keyword in ("CAMR", "KANR", "AMPR", "HYG", "SPEC", "ZEO")):
            continue
        if label not in summary:
            summary.append(label)
        if len(summary) >= limit:
            break
    return summary


def description_from_record(record, filename):
    raw_description = str(record.description or "").strip()
    if raw_description and raw_description not in {".", "<unknown description>"}:
        return raw_description[:1000]

    summary = feature_summary(record)
    if summary:
        return f"{DEFAULT_DESCRIPTION_PREFIX} {filename}. Features: {', '.join(summary)}."[:1000]

    return f"{DEFAULT_DESCRIPTION_PREFIX} {filename}."[:1000]


def record_name_from_source(record, source_name, source_stem, name_source="filename"):
    if name_source == "record":
        return trimmed_name(record.name or record.id or source_stem)
    return trimmed_name(source_stem)


def plasmid_type_by_name(type_name):
    if not type_name:
        return None
    return PlasmidType.objects.filter(name__iexact=type_name).first()


def resistance_objects_by_code(codes):
    resistances = []
    for code in codes:
        resistance = Resistance.objects.filter(three_letter_code=code).first()
        if resistance:
            resistances.append(resistance)
    return resistances


def import_plasmids_from_genbank_sources(
        sources,
        project,
        dry_run=False,
        update_existing=False,
        public_visibility=False,
        reference_sequence=True,
        infer_ytk_metadata=True,
        name_source="filename"):
    ensure_minimum_catalog()
    sources = list(sources or [])
    if not sources:
        raise GenBankImportError("No GenBank files were provided.")

    result = {
        "created": 0,
        "updated": 0,
        "skipped": 0,
        "errors": 0,
        "messages": [],
    }

    for source in sources:
        filename = source["filename"]
        if not is_genbank_filename(filename):
            result["errors"] += 1
            result["skipped"] += 1
            result["messages"].append({
                "level": "danger",
                "text": f"Skipping {filename}: unsupported file extension.",
            })
            continue

        try:
            content = source["content"]
            text = content.decode("utf-8-sig")
            record = SeqIO.read(StringIO(text), "genbank")
        except UnicodeDecodeError:
            result["errors"] += 1
            result["skipped"] += 1
            result["messages"].append({
                "level": "danger",
                "text": f"Skipping {filename}: could not read this file as UTF-8 text.",
            })
            continue
        except Exception as error:
            result["errors"] += 1
            result["skipped"] += 1
            result["messages"].append({
                "level": "danger",
                "text": f"Skipping {filename}: could not parse GenBank ({error}).",
            })
            continue

        try:
            name = record_name_from_source(record, filename, source["stem"], name_source=name_source)
        except GenBankImportError as error:
            result["errors"] += 1
            result["skipped"] += 1
            result["messages"].append({
                "level": "danger",
                "text": f"Skipping {filename}: {error}",
            })
            continue

        type_name, level = (None, None)
        classification = None
        if infer_ytk_metadata:
            type_name, level, classification = infer_assembly_type_and_level(record, project=project)
            if type_name is None and level is None:
                type_name, level = infer_ytk_type_and_level(record)

        resistance_codes = infer_resistance_codes(record)
        description = description_from_record(record, filename)
        created_on = parse_created_on(record)
        existing = None
        if project.pk:
            existing = Plasmid.objects.filter(name=name, project=project).first()
        action = "create"

        if existing:
            if not update_existing:
                result["skipped"] += 1
                result["messages"].append({
                    "level": "secondary",
                    "text": f"Skipping existing plasmid: {name}.",
                })
                continue
            action = "update"

        if dry_run:
            result["created" if action == "create" else "updated"] += 1
            text = f"Would {action} plasmid: {name}."
            if classification and classification.part_name:
                text += f" Detected {classification.part_name} ({classification.assembly_level}, {classification.confidence_band})."
            result["messages"].append({
                "level": "info",
                "text": text,
            })
            continue

        plasmid = existing or Plasmid(
            name=name,
            project=project,
            intended_use=DEFAULT_INTENDED_USE,
        )
        plasmid.intended_use = DEFAULT_INTENDED_USE
        plasmid.description = description
        plasmid.type = plasmid_type_by_name(type_name)
        plasmid.level = level
        plasmid.created_on = created_on
        if classification:
            plasmid.assembly_metadata = assembly_metadata_from_classification(
                classification,
                confirmed_type_name=type_name,
                confirmed_level=level,
                confirmed_standard_id=getattr(project, "assembly_standard", None) or classification.standard_id,
            )
        elif not existing:
            plasmid.assembly_metadata = {}
        plasmid.reference_sequence = reference_sequence
        plasmid.public_visibility = public_visibility
        plasmid.computed_size = len(record.seq)
        plasmid.insert_computed_size = None
        plasmid.save()
        plasmid.sequence.save(filename, ContentFile(content), save=False)
        plasmid.save()
        plasmid.selectable_markers.set(resistance_objects_by_code(resistance_codes))

        result["created" if action == "create" else "updated"] += 1
        text = f"{'Created' if action == 'create' else 'Updated'} plasmid: {name}."
        if classification and classification.part_name:
            text += f" Detected {classification.part_name} ({classification.assembly_level}, {classification.confidence_band})."
        result["messages"].append({
            "level": "success",
            "text": text,
        })

    return result


def import_plasmids_from_genbank_dir(
        genbank_dir,
        project,
        dry_run=False,
        update_existing=False,
        public_visibility=False,
        reference_sequence=True,
        infer_ytk_metadata=True,
        name_source="filename"):
    return import_plasmids_from_genbank_sources(
        genbank_sources_from_dir(genbank_dir),
        project,
        dry_run=dry_run,
        update_existing=update_existing,
        public_visibility=public_visibility,
        reference_sequence=reference_sequence,
        infer_ytk_metadata=infer_ytk_metadata,
        name_source=name_source,
    )


def import_plasmids_from_uploaded_genbanks(
        uploaded_files,
        project,
        dry_run=False,
        update_existing=False,
        public_visibility=False,
        reference_sequence=True,
        infer_ytk_metadata=True,
        name_source="filename"):
    return import_plasmids_from_genbank_sources(
        genbank_sources_from_uploaded_files(uploaded_files),
        project,
        dry_run=dry_run,
        update_existing=update_existing,
        public_visibility=public_visibility,
        reference_sequence=reference_sequence,
        infer_ytk_metadata=infer_ytk_metadata,
        name_source=name_source,
    )
