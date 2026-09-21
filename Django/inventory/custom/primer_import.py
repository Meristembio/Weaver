import re

from Bio import SeqIO

from inventory.models import Primer


VALID_DNA = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)
FWD_RE = re.compile(r"(^|[-_\s.])(f|fw|fwd|forward)([-_\s.]|$)", re.IGNORECASE)
REV_RE = re.compile(r"(^|[-_\s.])(r|rev|reverse)([-_\s.]|$)", re.IGNORECASE)
METADATA_DIRECTION_COLUMNS = ("direction", "fwd_or_rev", "orientation")
METADATA_OVERHANG_COLUMNS = ("overhang", "sequence_5", "seq_5", "five_prime_overhang")
METADATA_INTENDED_USE_COLUMNS = ("intended_use", "use", "description", "notes")
FASTA_METADATA_RE = re.compile(r"(?P<key>[A-Za-z_][A-Za-z0-9_-]*)=(?P<value>\"[^\"]*\"|'[^']*'|\S+)")


class PrimerImportError(ValueError):
    pass


def normalize_sequence(sequence):
    return re.sub(r"\s+", "", str(sequence))


def infer_direction(name):
    has_fwd = FWD_RE.search(name)
    has_rev = REV_RE.search(name)

    if has_fwd and not has_rev:
        return "f"
    if has_rev and not has_fwd:
        return "r"
    return ""


def normalize_direction(value):
    value = str(value or "").strip().lower()
    if value in ("f", "fw", "fwd", "forward"):
        return "f"
    if value in ("r", "rev", "reverse"):
        return "r"
    return ""


def first_header_value(row, columns):
    for column in columns:
        value = row.get(column)
        if value not in (None, ""):
            return str(value).strip()
    return ""


def parse_fasta_header_fields(description):
    fields = {}
    for match in FASTA_METADATA_RE.finditer(description or ""):
        key = match.group("key").strip().lower().replace("-", "_")
        value = match.group("value").strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
            value = value[1:-1]
        fields[key] = value

    return {
        "direction": normalize_direction(first_header_value(fields, METADATA_DIRECTION_COLUMNS)),
        "sequence_5": normalize_sequence(first_header_value(fields, METADATA_OVERHANG_COLUMNS)),
        "intended_use": first_header_value(fields, METADATA_INTENDED_USE_COLUMNS),
    }


def strip_declared_overhang(sequence, sequence_5):
    if sequence_5 and sequence.startswith(sequence_5):
        return sequence[len(sequence_5):]
    return sequence


def primer_entries_from_fasta(
        handle,
        name_source="id",
        require_direction=False,
        default_direction=""):
    records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        raise PrimerImportError("No FASTA records found.")

    entries = []
    for record in records:
        name = record.id if name_source == "id" else record.description
        name = name.strip()
        header_fields = parse_fasta_header_fields(record.description)
        sequence = normalize_sequence(record.seq)
        direction = header_fields.get("direction") or infer_direction(name) or default_direction
        sequence_5 = header_fields.get("sequence_5", "")
        sequence_3 = strip_declared_overhang(sequence, sequence_5)
        intended_use = header_fields.get("intended_use") or "Imported from FASTA"
        errors = []

        if not name:
            errors.append("empty name")
        if not sequence_3:
            errors.append("empty sequence")
        if sequence and not VALID_DNA.match(sequence):
            errors.append("sequence contains non-DNA/IUPAC characters")
        if require_direction and not direction:
            errors.append("could not infer F/R direction")
        if sequence_5 and not VALID_DNA.match(sequence_5):
            errors.append("overhang contains non-DNA/IUPAC characters")

        entries.append({
            "name": name,
            "sequence": sequence_3,
            "direction": direction,
            "sequence_5": sequence_5,
            "intended_use": intended_use,
            "errors": errors,
        })
    return entries


def import_primers_from_fasta(
        handle,
        dry_run=False,
        update_existing=False,
        name_source="id",
        require_direction=False,
        default_direction=""):
    entries = primer_entries_from_fasta(
        handle,
        name_source=name_source,
        require_direction=require_direction,
        default_direction=default_direction,
    )
    result = {
        "created": 0,
        "updated": 0,
        "skipped": 0,
        "errors": 0,
        "messages": [],
    }

    for entry in entries:
        if entry["errors"]:
            result["skipped"] += 1
            result["errors"] += 1
            result["messages"].append({
                "level": "danger",
                "text": f"Skipping {entry['name'] or 'unnamed record'}: {', '.join(entry['errors'])}.",
            })
            continue

        primer = Primer.objects.filter(name=entry["name"]).first()

        if primer:
            if not update_existing:
                result["skipped"] += 1
                result["messages"].append({
                    "level": "secondary",
                    "text": f"Skipping existing primer: {entry['name']}.",
                })
                continue

            if not dry_run:
                primer.sequence_3 = entry["sequence"]
                primer.sequence_5 = entry["sequence_5"]
                primer.fwd_or_rev = entry["direction"]
                primer.intended_use = entry["intended_use"]
                primer.save()
            result["updated"] += 1
            result["messages"].append({
                "level": "info",
                "text": f"{'Would update' if dry_run else 'Updated'} primer: {entry['name']}.",
            })
            continue

        if not dry_run:
            Primer.objects.create(
                name=entry["name"],
                sequence_3=entry["sequence"],
                sequence_5=entry["sequence_5"],
                fwd_or_rev=entry["direction"],
                intended_use=entry["intended_use"],
            )
        result["created"] += 1
        result["messages"].append({
            "level": "success",
            "text": f"{'Would create' if dry_run else 'Created'} primer: {entry['name']} ({entry['direction'] or 'unknown'}).",
        })

    return result
