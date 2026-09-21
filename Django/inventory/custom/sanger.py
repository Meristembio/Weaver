import csv
import datetime
import hashlib
import io
import json
import os
import re
from dataclasses import dataclass
from dataclasses import field
from statistics import mean
from statistics import median

from Bio import Align
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


VALID_SEQ_EXTENSIONS = (".ab1", ".phd.1", ".seq")
VALID_IUPAC_DNA = set("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class SangerProcessingParameters:
    max_files: int = 96
    max_file_size: int = 8 * 1024 * 1024
    max_total_size: int = 96 * 1024 * 1024
    quality_threshold: int = 20
    quality_window: int = 12
    minimum_trimmed_length: int = 40
    max_ambiguous_fraction: float = 0.15
    min_good_quality_region: int = 30
    min_identity_for_pass: float = 98.0
    min_identity_for_review: float = 95.0
    min_read_coverage_for_pass: float = 80.0
    min_combined_coverage_for_pass: float = 0.0
    high_quality_variant_phred: int = 25
    max_high_quality_variant_density_for_pass: float = 0.005
    max_high_quality_variant_density_for_review: float = 0.03
    max_high_quality_variants_for_pass: int = 2
    orientation_score_delta: float = 4.0
    match_score: float = 2.0
    mismatch_score: float = -3.0
    open_gap_score: float = -5.0
    extend_gap_score: float = -1.0
    window_size: int = 15
    low_median_phred: int = 15
    good_median_phred: int = 20
    low_q20_fraction: float = 0.60
    critical_fraction: float = 0.30
    recovery_length: int = 20
    merge_gap: int = 10
    minimum_region_length: int = 5
    region_padding: int = 3
    terminal_backtrack_window: int = 100
    terminal_min_low_region_length: int = 50
    terminal_min_start_fraction: float = 0.60
    signal_window: int = 75
    relative_signal_cutoff: float = 0.35
    secondary_peak_ratio_cutoff: float = 0.50

    def as_dict(self):
        return self.__dict__.copy()


@dataclass
class UploadedSangerFile:
    original_name: str
    data: bytes
    size: int
    sha256: str
    format: str
    group_name: str
    metadata: dict = field(default_factory=dict)
    errors: list = field(default_factory=list)


def as_uploaded_sanger_file_list(value):
    """Normalize one uploaded-file object or an iterable of uploaded files."""
    if value is None:
        return []
    if isinstance(value, (UploadedSangerFile, str, bytes)):
        return [value]
    if isinstance(value, (list, tuple)):
        return list(value)
    try:
        return list(value)
    except TypeError:
        return [value]


@dataclass
class ParsedSource:
    format: str
    sequence: str = ""
    qualities: list = field(default_factory=list)
    peak_positions: list = field(default_factory=list)
    chromatogram: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)
    clear_range: tuple = None
    errors: list = field(default_factory=list)


def detect_format(filename):
    lowered = filename.lower()
    if lowered.endswith(".phd.1"):
        return "phd1"
    if lowered.endswith(".ab1"):
        return "ab1"
    if lowered.endswith(".seq"):
        return "seq"
    return ""


def normalized_group_name(filename):
    basename = os.path.basename(filename).strip()
    lowered = basename.lower()
    if lowered.endswith(".phd.1"):
        return basename[:-6]
    if lowered.endswith(".ab1") or lowered.endswith(".seq"):
        return basename[:-4]
    return os.path.splitext(basename)[0]


def preferred_sanger_run(runs):
    """Return the run shown by default for a plasmid.

    The newest manually accepted run is preferred. When there is no manually
    accepted run, the newest run that passed automatic classification is used,
    followed by the newest run needing review. If neither exists, fall back to
    the newest uploaded run. A run's identity and read records are not changed
    by this display choice.
    """
    def sequencing_timestamp(run):
        persisted_datetime = getattr(run, "sequencing_datetime", None)
        if persisted_datetime:
            return persisted_datetime.timestamp()
        return run.created_at.timestamp()

    ordered_runs = sorted(
        list(runs or []),
        key=lambda run: (sequencing_timestamp(run), str(run.id)),
        reverse=True,
    )
    if not ordered_runs:
        return None
    accepted_runs = [run for run in ordered_runs if run.manual_decision == "VERIFIED"]
    if accepted_runs:
        return accepted_runs[0]
    automatically_approved_runs = [run for run in ordered_runs if run.automated_state == "PASS"]
    if automatically_approved_runs:
        return automatically_approved_runs[0]
    review_runs = [run for run in ordered_runs if run.automated_state == "REVIEW"]
    return review_runs[0] if review_runs else ordered_runs[0]


def parse_sanger_batch_mapping(data):
    """Parse CSV rows mapping an AB1 filename to a plasmid identifier."""
    try:
        text = data.decode("utf-8-sig") if isinstance(data, bytes) else str(data or "")
    except UnicodeDecodeError:
        return [], ["The mapping CSV must be UTF-8 encoded."]

    reader = csv.DictReader(io.StringIO(text))
    if not reader.fieldnames:
        return [], ["The mapping CSV is empty or has no header row."]

    normalized_headers = {
        (header or "").strip().lower(): header
        for header in reader.fieldnames
        if header is not None
    }
    filename_header = next(
        (normalized_headers[name] for name in ("ab1_file", "filename", "file_name", "file") if name in normalized_headers),
        None,
    )
    plasmid_header = next(
        (normalized_headers[name] for name in ("plasmid_id", "plasmid", "id") if name in normalized_headers),
        None,
    )
    primer_header = next(
        (normalized_headers[name] for name in ("primer_id", "primer", "primer_uuid") if name in normalized_headers),
        None,
    )
    errors = []
    if not filename_header or not plasmid_header:
        return [], ["The CSV must contain the columns ab1_file and plasmid_id."]

    rows = []
    for line_number, row in enumerate(reader, start=2):
        filename = os.path.basename((row.get(filename_header) or "").strip())
        plasmid_id = (row.get(plasmid_header) or "").strip()
        primer_id = (row.get(primer_header) or "").strip() if primer_header else ""
        if not filename or not plasmid_id:
            errors.append("Row {} must contain ab1_file and plasmid_id.".format(line_number))
            continue
        rows.append({
            "filename": filename,
            "plasmid_id": plasmid_id,
            "primer_id": primer_id,
            "line_number": line_number,
        })
    if not rows and not errors:
        errors.append("The mapping CSV has no data rows.")
    return rows, errors


def sanitize_filename(filename):
    basename = os.path.basename(filename)
    return re.sub(r"[^A-Za-z0-9._-]+", "_", basename)[:255] or "sanger-file"


def uploaded_files_from_request(files, parameters=None):
    parameters = parameters or SangerProcessingParameters()
    uploaded = []
    total_size = 0
    seen_hash_names = set()
    for file_obj in as_uploaded_sanger_file_list(files):
        data = file_obj.read()
        size = len(data)
        total_size += size
        fmt = detect_format(file_obj.name)
        item = UploadedSangerFile(
            original_name=sanitize_filename(file_obj.name),
            data=data,
            size=size,
            sha256=hashlib.sha256(data).hexdigest(),
            format=fmt,
            group_name=normalized_group_name(file_obj.name),
        )
        if not fmt:
            item.errors.append("Unsupported file format")
        if size == 0:
            item.errors.append("Empty file")
        if size > parameters.max_file_size:
            item.errors.append("File exceeds maximum size")
        duplicate_key = (item.sha256, item.original_name)
        if duplicate_key in seen_hash_names:
            item.errors.append("Duplicate file")
        seen_hash_names.add(duplicate_key)
        uploaded.append(item)

    if len(uploaded) > parameters.max_files:
        for item in uploaded:
            item.errors.append("Batch exceeds maximum number of files")
    if total_size > parameters.max_total_size:
        for item in uploaded:
            item.errors.append("Batch exceeds maximum total size")
    return uploaded


def parse_seq(data):
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        return ParsedSource(format="seq", errors=["SEQ file is not valid UTF-8 text"])
    lines = [line.strip() for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n")]
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    sequence = re.sub(r"\s+", "", "".join(lines)).upper()
    bad = sorted(set(sequence) - VALID_IUPAC_DNA)
    if bad:
        return ParsedSource(format="seq", sequence=sequence, errors=["SEQ contains invalid DNA/IUPAC characters: {}".format(", ".join(bad))])
    if not sequence:
        return ParsedSource(format="seq", errors=["SEQ file has no sequence"])
    return ParsedSource(format="seq", sequence=sequence)


def parse_phd1(data):
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = data.decode("latin-1", errors="replace")
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    metadata = {}
    sequence = []
    qualities = []
    peaks = []
    in_comment = False
    in_dna = False
    read_name = ""
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("BEGIN_SEQUENCE"):
            read_name = stripped.split(maxsplit=1)[1] if len(stripped.split(maxsplit=1)) > 1 else ""
        elif stripped == "BEGIN_COMMENT":
            in_comment = True
        elif stripped == "END_COMMENT":
            in_comment = False
        elif stripped == "BEGIN_DNA":
            in_dna = True
        elif stripped == "END_DNA":
            in_dna = False
        elif in_comment and ":" in stripped:
            key, value = stripped.split(":", 1)
            metadata[key.strip()] = value.strip()
        elif in_dna and stripped:
            parts = stripped.split()
            if len(parts) >= 3:
                base, quality, peak = parts[:3]
                sequence.append(base.upper())
                try:
                    qualities.append(int(float(quality)))
                except ValueError:
                    qualities.append(None)
                try:
                    peaks.append(int(float(peak)))
                except ValueError:
                    peaks.append(None)
    joined = "".join(sequence)
    errors = []
    bad = sorted(set(joined) - VALID_IUPAC_DNA)
    if bad:
        errors.append("PHD.1 contains invalid DNA/IUPAC characters: {}".format(", ".join(bad)))
    if not joined:
        errors.append("PHD.1 file has no BEGIN_DNA bases")
    metadata["read_name"] = read_name
    clear_range = clear_range_from_metadata(metadata, len(joined))
    return ParsedSource(format="phd1", sequence=joined, qualities=qualities, peak_positions=peaks, metadata=metadata, clear_range=clear_range, errors=errors)


def _abif_value(raw, key):
    value = raw.get(key)
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return value


def ab1_run_datetime(metadata):
    """Return the instrument run start as a naive local datetime when available."""
    metadata = metadata or {}
    date_value = str(metadata.get("run_start_date") or "").strip()
    time_value = str(metadata.get("run_start_time") or "").strip()
    if not date_value:
        return None
    for value, formats in (
        ("{} {}".format(date_value, time_value), ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M")),
        (date_value, ("%Y-%m-%d", "%y%m%d", "%Y%m%d")),
    ):
        for date_format in formats:
            try:
                return datetime.datetime.strptime(value, date_format)
            except ValueError:
                continue
    return None


def ab1_run_datetime_label(metadata, filename=None):
    """Format the filename run date, falling back to AB1 metadata."""
    run_datetime = sanger_file_run_datetime(metadata, filename)
    return run_datetime.strftime("%d %b %Y %H:%M") if run_datetime else ""


def filename_run_datetime(filename):
    """Extract a sequencer run datetime from a filename."""
    match = re.search(
        r"(?<!\d)(20\d{2})-(\d{2})-(\d{2})-(\d{2})-(\d{2})-(\d{2})(?!\d)",
        os.path.basename(filename or ""),
    )
    if not match:
        return None
    try:
        return datetime.datetime.strptime("-".join(match.groups()), "%Y-%m-%d-%H-%M-%S")
    except ValueError:
        return None


def sanger_file_run_datetime(metadata=None, filename=None):
    """Return the filename date, falling back to the AB1 metadata date."""
    return filename_run_datetime(filename) or ab1_run_datetime(metadata)


def file_registration_number(filename):
    """Return the leading sequencer/file number used to order confirmation pairs."""
    match = re.match(r"\s*(\d+)", os.path.basename(filename or ""))
    return int(match.group(1)) if match else None


def latest_confirmation_pair(read_candidates):
    """Select the latest complete forward/reverse pair for the Services table.

    ``read_candidates`` contains one item per persisted read with ``read``,
    ``orientation``, ``run_date``, ``run_datetime`` and ``registration_number``
    keys. Filename run dates take precedence over AB1 metadata dates. When
    several pairs share a date, the pair with the closest leading file numbers and
    highest registration numbers is preferred. If only one orientation exists,
    only its most recent read is returned. Unknown-orientation reads are not
    included when any reliable orientation is available; they remain
    accessible from the saved-run detail view for troubleshooting.
    """
    selected_candidates = latest_sanger_complete_date_group(read_candidates)
    forward = [item for item in selected_candidates if item.get("orientation") == "forward"]
    reverse = [item for item in selected_candidates if item.get("orientation") == "reverse"]
    if not forward and not reverse:
        return [item["read"] for item in selected_candidates]

    def datetime_key(item):
        run_datetime = item.get("run_datetime")
        if run_datetime is not None:
            return run_datetime
        run_date = item.get("run_date")
        if run_date is not None:
            return datetime.datetime.combine(run_date, datetime.time.min)
        return datetime.datetime.min

    def number_key(item):
        number = item.get("registration_number")
        return number if number is not None else -1

    if not forward or not reverse:
        selected_orientation = forward or reverse
        selected = max(selected_orientation, key=lambda item: (datetime_key(item), number_key(item)))
        return [selected["read"]]

    numbered_pairs = []
    for forward_item in forward:
        for reverse_item in reverse:
            forward_number = forward_item.get("registration_number")
            reverse_number = reverse_item.get("registration_number")
            if forward_number is not None and reverse_number is not None:
                numbered_pairs.append((
                    abs(forward_number - reverse_number),
                    max(forward_number, reverse_number),
                    forward_item,
                    reverse_item,
                ))
    if numbered_pairs:
        _, _, selected_forward, selected_reverse = min(
            numbered_pairs,
            key=lambda pair: (pair[0], -pair[1]),
        )
    else:
        selected_forward = max(forward, key=number_key)
        selected_reverse = max(reverse, key=number_key)
    return [selected_forward["read"], selected_reverse["read"]]


def latest_sanger_complete_date_group(read_candidates):
    """Return the newest date group containing both Forward and Reverse reads.

    Files from one upload can represent more than one sequencing event. A
    complete Forward/Reverse group defines the active event; multiple reads in
    that group remain available for region-based nonredundancy selection.
    """
    candidates = list(read_candidates or [])
    dated_groups = {}
    for item in candidates:
        if item.get("run_date") is not None:
            dated_groups.setdefault(item["run_date"], []).append(item)
    eligible_groups = [
        (run_date, group)
        for run_date, group in dated_groups.items()
        if any(item.get("orientation") == "forward" for item in group)
        and any(item.get("orientation") == "reverse" for item in group)
    ]
    if not eligible_groups:
        return candidates
    return max(eligible_groups, key=lambda entry: entry[0])[1]


def select_nonredundant_read_candidates(read_candidates, overlap_threshold=0.8):
    """Keep distinct regions while collapsing redundant reads per orientation.

    A run may legitimately contain more than one Forward or Reverse read when
    the reads cover different plasmid regions. Reads with at least
    ``overlap_threshold`` of the smaller covered region in common are treated
    as repeated coverage of the same region, and the newest candidate is kept.
    Candidates without coverage coordinates are retained because their regions
    cannot be compared safely.
    """
    candidates = list(read_candidates or [])
    selected_indexes = set(
        index
        for index, candidate in enumerate(candidates)
        if candidate.get("orientation") not in ("forward", "reverse")
    )

    def overlap(left, right):
        left_positions = set(left.get("covered_positions") or [])
        right_positions = set(right.get("covered_positions") or [])
        if not left_positions or not right_positions:
            return 0.0
        return len(left_positions & right_positions) / min(len(left_positions), len(right_positions))

    def recency_key(index):
        candidate = candidates[index]
        run_datetime = candidate.get("run_datetime")
        registration_number = candidate.get("registration_number")
        return (
            run_datetime is not None,
            run_datetime or datetime.datetime.min,
            registration_number is not None,
            registration_number if registration_number is not None else -1,
            -index,
        )

    for orientation in ("forward", "reverse"):
        orientation_indexes = [
            index
            for index, candidate in enumerate(candidates)
            if candidate.get("orientation") == orientation
        ]
        remaining = set(orientation_indexes)
        while remaining:
            component = {min(remaining)}
            expanded = True
            while expanded:
                expanded = False
                for index in list(remaining - component):
                    if any(overlap(candidates[index], candidates[other]) >= overlap_threshold for other in component):
                        component.add(index)
                        expanded = True
            selected_indexes.add(max(component, key=recency_key))
            remaining -= component

    return [candidate for index, candidate in enumerate(candidates) if index in selected_indexes]


def select_sanger_review_candidates(read_candidates):
    """Return the reads that should be shown in a saved-run review.

    Failed or unmapped reads are useful troubleshooting evidence only when a
    run has no usable aligned reads. If usable evidence exists, review the
    nonredundant aligned reads and keep failed records out of the active
    sequencing view.
    """
    candidates = latest_sanger_complete_date_group(read_candidates)
    usable_candidates = [
        candidate for candidate in candidates
        if candidate.get("read", {}).get("is_usable") and candidate.get("read", {}).get("alignment")
    ]
    return select_nonredundant_read_candidates(usable_candidates or candidates)


def clear_range_from_metadata(metadata, sequence_length):
    for key in ("TRIM", "TRIM_RANGE", "CLEAR_RANGE", "CLEAR RANGE", "PHRED_CLEAR_RANGE"):
        value = metadata.get(key)
        if not value:
            continue
        numbers = [int(item) for item in re.findall(r"\d+", str(value))]
        if len(numbers) >= 2:
            start, end = numbers[0], numbers[1]
            if start >= 1 and end >= start:
                start -= 1
            return normalize_clear_range(start, end, sequence_length)
    return None


def clear_range_from_abif(raw, sequence_length):
    candidates = []
    for left_key, right_key in (("CLIP1", "CLIP2"), ("QV20L", "QV20R"), ("SMLt1", "SMRt1")):
        left = raw.get(left_key)
        right = raw.get(right_key)
        if left is None or right is None:
            continue
        if isinstance(left, (list, tuple)):
            left = left[0] if left else None
        if isinstance(right, (list, tuple)):
            right = right[0] if right else None
        try:
            candidates.append((int(left), int(right), left_key, right_key))
        except (TypeError, ValueError):
            continue
    for start, end, left_key, right_key in candidates:
        normalized = normalize_clear_range(start, end, sequence_length)
        if normalized and normalized[1] > normalized[0]:
            return normalized, "{}:{}".format(left_key, right_key)
    return None, ""


def normalize_clear_range(start, end, sequence_length):
    if sequence_length <= 0:
        return None
    start = max(0, min(sequence_length, int(start)))
    end = max(0, min(sequence_length, int(end)))
    if end < start:
        start, end = end, start
    return (start, end)


def parse_ab1(data):
    try:
        record = SeqIO.read(io.BytesIO(data), "abi")
    except Exception as exc:
        return ParsedSource(format="ab1", errors=["AB1 parse error: {}".format(exc)])
    raw = record.annotations.get("abif_raw", {})
    chromatogram = {}
    dye_order = str(_abif_value(raw, "FWO_1") or "GATC").upper()
    data_keys = ("DATA9", "DATA10", "DATA11", "DATA12")
    trace_key_by_base = {"A": "aTrace", "C": "cTrace", "G": "gTrace", "T": "tTrace"}
    for base, raw_key in zip(dye_order, data_keys):
        out_key = trace_key_by_base.get(base)
        if not out_key:
            continue
        value = raw.get(raw_key, [])
        chromatogram[out_key] = list(value) if value is not None else []
    for out_key in trace_key_by_base.values():
        chromatogram.setdefault(out_key, [])
    chromatogram["basePos"] = list(raw.get("PLOC2", []))
    pbas = raw.get("PBAS2", b"")
    chromatogram["baseCalls"] = list(pbas.decode(errors="replace")) if isinstance(pbas, bytes) else list(str(pbas))
    chromatogram["qualNums"] = list(record.letter_annotations.get("phred_quality", []))
    clear_range, clear_range_source = clear_range_from_abif(raw, len(record.seq))
    metadata = {
        "record_id": record.id,
        "dye_order": dye_order,
        "clear_range_source": clear_range_source,
        "run_start_date": str(_abif_value(raw, "RUND1") or ""),
        "run_end_date": str(_abif_value(raw, "RUND2") or ""),
        "run_start_time": str(_abif_value(raw, "RUNT1") or ""),
        "run_end_time": str(_abif_value(raw, "RUNT2") or ""),
        "instrument_model": _abif_value(raw, "MODL1"),
        "machine": _abif_value(raw, "MCHN1"),
    }
    return ParsedSource(
        format="ab1",
        sequence=str(record.seq).upper(),
        qualities=list(record.letter_annotations.get("phred_quality", [])),
        peak_positions=list(raw.get("PLOC2", [])),
        chromatogram=chromatogram,
        metadata=metadata,
        clear_range=clear_range,
    )


def parse_file(uploaded):
    if uploaded.errors:
        return ParsedSource(format=uploaded.format, errors=list(uploaded.errors))
    if uploaded.format == "ab1":
        return parse_ab1(uploaded.data)
    if uploaded.format == "phd1":
        return parse_phd1(uploaded.data)
    if uploaded.format == "seq":
        return parse_seq(uploaded.data)
    return ParsedSource(format=uploaded.format, errors=["Unsupported file format"])


def group_uploaded_files(uploaded_files):
    groups = {}
    for uploaded in as_uploaded_sanger_file_list(uploaded_files):
        groups.setdefault(uploaded.group_name, []).append(uploaded)
    return groups


def choose_source(parsed_sources):
    warnings = []
    by_format = {source.format: source for source in parsed_sources if source.sequence and not source.errors}
    selected = by_format.get("ab1") or by_format.get("phd1") or by_format.get("seq")
    if not selected:
        return None, warnings
    comparable = [source for source in parsed_sources if source.sequence and not source.errors]
    for source in comparable:
        if source is selected:
            continue
        if source.sequence != selected.sequence:
            identity = sequence_identity(selected.sequence, source.sequence)
            warnings.append(
                "Sequence from {} differs from selected {} source ({} bp vs {} bp, {:.1f}% identity)".format(
                    source.format.upper(), selected.format.upper(), len(source.sequence), len(selected.sequence), identity
                )
            )
    return selected, warnings


def sequence_identity(left, right):
    if not left and not right:
        return 100.0
    length = max(len(left), len(right))
    matches = sum(1 for a, b in zip(left, right) if a == b)
    return matches / length * 100.0 if length else 0.0


def clean_quality_values(qualities, length):
    values = []
    for i in range(length):
        quality = qualities[i] if i < len(qualities) else None
        values.append(quality if quality is not None else 0)
    return values


def max_critical_run(sequence, qualities):
    longest = current = 0
    for base, quality in zip(sequence, qualities):
        if quality < 10 or base not in "ACGT":
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def window_quality_stats(sequence, qualities, start, end):
    window_sequence = sequence[start:end]
    window_qualities = qualities[start:end]
    if not window_qualities:
        return {
            "median": 0,
            "low_q20_fraction": 1.0,
            "critical_fraction": 1.0,
            "low_q15_fraction": 1.0,
            "critical_run": 0,
        }
    return {
        "median": median(window_qualities),
        "low_q20_fraction": sum(1 for q in window_qualities if q < 20) / len(window_qualities),
        "critical_fraction": sum(1 for base, q in zip(window_sequence, window_qualities) if q < 10 or base not in "ACGT") / len(window_qualities),
        "low_q15_fraction": sum(1 for q in window_qualities if q < 15) / len(window_qualities),
        "critical_run": max_critical_run(window_sequence, window_qualities),
    }


def region_metrics(sequence, qualities, start, end, signal_metrics=None, reasons=None, level="LOW_CONFIDENCE"):
    region_qualities = qualities[start:end]
    region_signal = signal_metrics or {}
    relative_signal = region_signal.get("relative_signal", [None] * len(sequence))
    secondary_ratios = region_signal.get("secondary_peak_ratio", [None] * len(sequence))
    rel_values = [value for value in relative_signal[start:end] if value is not None]
    ratio_values = [value for value in secondary_ratios[start:end] if value is not None]
    return {
        "level": level,
        "start": start,
        "end": end - 1,
        "start_display": start + 1,
        "end_display": end,
        "length": max(0, end - start),
        "mean_phred": round(mean(region_qualities), 2) if region_qualities else None,
        "median_phred": round(median(region_qualities), 2) if region_qualities else None,
        "q20_low_percent": round(sum(1 for q in region_qualities if q < 20) / len(region_qualities) * 100, 2) if region_qualities else None,
        "mean_relative_signal": round(mean(rel_values), 3) if rel_values else None,
        "median_relative_signal": round(median(rel_values), 3) if rel_values else None,
        "mean_secondary_peak_ratio": round(mean(ratio_values), 3) if ratio_values else None,
        "reasons": reasons or [],
    }


def merge_regions(regions, gap):
    if not regions:
        return []
    ordered = sorted(regions)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start - merged[-1][1] <= gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(item) for item in merged]


def complement_blocks(length, regions):
    blocks = []
    cursor = 0
    for start, end in sorted(regions):
        if cursor < start:
            blocks.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < length:
        blocks.append((cursor, length))
    return blocks


def intersect_blocks(left_blocks, right_block):
    if not right_block:
        return left_blocks
    right_start, right_end = right_block
    blocks = []
    for start, end in left_blocks:
        block_start = max(start, right_start)
        block_end = min(end, right_end)
        if block_end > block_start:
            blocks.append((block_start, block_end))
    return blocks


def contiguous_true_regions(flags, minimum_length=1):
    regions = []
    start = None
    for index, value in enumerate(flags):
        if value and start is None:
            start = index
        elif not value and start is not None:
            if index - start >= minimum_length:
                regions.append((start, index))
            start = None
    if start is not None and len(flags) - start >= minimum_length:
        regions.append((start, len(flags)))
    return regions


def chromatogram_signal_metrics(sequence, chromatogram, params):
    base_positions = chromatogram.get("basePos") or []
    traces = {
        "A": chromatogram.get("aTrace") or [],
        "C": chromatogram.get("cTrace") or [],
        "G": chromatogram.get("gTrace") or [],
        "T": chromatogram.get("tTrace") or [],
    }
    if not sequence or len(base_positions) < len(sequence) or not all(traces[base] for base in "ACGT"):
        return {"available": False, "relative_signal": [], "secondary_peak_ratio": [], "warning": ""}
    main_heights = []
    secondary_ratios = []
    for index, base in enumerate(sequence):
        position = int(base_positions[index])
        heights = {}
        for channel_base, trace in traces.items():
            if 0 <= position < len(trace):
                heights[channel_base] = float(trace[position])
            else:
                heights[channel_base] = 0.0
        main = heights.get(base, max(heights.values()) if heights else 0.0)
        second = sorted(heights.values(), reverse=True)[1] if len(heights) > 1 else 0.0
        main_heights.append(main)
        secondary_ratios.append(second / main if main > 0 else 1.0)
    half = max(1, params.signal_window // 2)
    relative = []
    for index, height in enumerate(main_heights):
        start = max(0, index - half)
        end = min(len(main_heights), index + half + 1)
        baseline_values = [value for value in main_heights[start:end] if value > 0]
        baseline = median(baseline_values) if baseline_values else 0
        relative.append(height / baseline if baseline > 0 else None)
    low_signal = [
        value is not None and value < params.relative_signal_cutoff
        for value in relative
    ]
    secondary_signal = [
        value is not None and value >= params.secondary_peak_ratio_cutoff
        for value in secondary_ratios
    ]
    return {
        "available": True,
        "relative_signal": relative,
        "secondary_peak_ratio": secondary_ratios,
        "low_signal_regions": contiguous_true_regions(low_signal, 10),
        "secondary_peak_regions": contiguous_true_regions(secondary_signal, 10),
    }


def detect_confidence_regions(sequence, qualities, chromatogram=None, clear_range=None, params=None):
    params = params or SangerProcessingParameters()
    length = len(sequence)
    if not length:
        return {
            "quality_available": bool(qualities),
            "low_confidence_regions": [],
            "intermediate_confidence_regions": [],
            "accepted_blocks": [],
            "alignment_blocks": [],
            "weaver_clear_range": None,
            "file_clear_range": clear_range,
            "signal": {"available": False},
        }
    if not qualities:
        ambiguous = sum(1 for base in sequence if base not in "ACGT")
        return {
            "quality_available": False,
            "low_confidence_regions": [],
            "intermediate_confidence_regions": [],
            "accepted_blocks": [(0, length)],
            "alignment_blocks": [(0, length)],
            "weaver_clear_range": (0, length),
            "file_clear_range": clear_range,
            "raw_length": length,
            "trimmed_length": length,
            "ambiguous_bases": ambiguous,
            "ambiguous_fraction": ambiguous / length if length else 1.0,
            "reason": "quality unavailable",
            "signal": {"available": False},
        }
    clean_qualities = clean_quality_values(qualities, length)
    window = max(1, params.window_size)
    low_windows = []
    good_windows = []
    for index in range(length):
        start = max(0, index - window // 2)
        end = min(length, start + window)
        start = max(0, end - window)
        stats = window_quality_stats(sequence, clean_qualities, start, end)
        low_windows.append(
            stats["median"] < params.low_median_phred
            or stats["low_q20_fraction"] >= params.low_q20_fraction
            or stats["critical_fraction"] >= params.critical_fraction
        )
        good_windows.append(
            stats["median"] >= params.good_median_phred
            and stats["low_q15_fraction"] < 0.20
            and stats["critical_run"] <= 2
        )
    rough_regions = []
    in_low = False
    region_start = None
    recovery = 0
    for index, low in enumerate(low_windows):
        if not in_low and low:
            in_low = True
            region_start = index
            recovery = 0
        elif in_low:
            if good_windows[index]:
                recovery += 1
                if recovery >= params.recovery_length:
                    rough_regions.append((region_start, index - params.recovery_length + 1))
                    in_low = False
                    region_start = None
                    recovery = 0
            else:
                recovery = 0
    if in_low and region_start is not None:
        rough_regions.append((region_start, length))

    signal = chromatogram_signal_metrics(sequence, chromatogram or {}, params)
    signal_regions = []
    signal_warnings = []
    if signal.get("available"):
        for start, end in signal.get("low_signal_regions", []):
            region_qualities = clean_qualities[start:end]
            if region_qualities and median(region_qualities) >= 20:
                signal_warnings.append("Low normalized signal disagrees with consistently high Phred scores at bases {}..{}".format(start + 1, end))
            else:
                signal_regions.append((start, end))
        for start, end in signal.get("secondary_peak_regions", []):
            region_qualities = clean_qualities[start:end]
            if region_qualities and median(region_qualities) >= 20:
                signal_warnings.append("High secondary peak ratio disagrees with consistently high Phred scores at bases {}..{}".format(start + 1, end))
            else:
                signal_regions.append((start, end))

    merged = merge_regions(rough_regions + signal_regions, params.merge_gap)
    consolidated = []
    for start, end in merged:
        contains_critical = any(base not in "ACGT" or q < 10 for base, q in zip(sequence[start:end], clean_qualities[start:end]))
        internal = start > 0 and end < length
        if internal and end - start < params.minimum_region_length and not contains_critical:
            continue
        padded_start = max(0, start - params.region_padding)
        padded_end = min(length, end + params.region_padding)
        consolidated.append((padded_start, padded_end))
    consolidated = merge_regions(consolidated, params.merge_gap)
    terminal_backtracked_ranges = set()
    if consolidated:
        last_start, last_end = consolidated[-1]
        terminal_length = last_end - last_start
        terminal_region = (
            last_end >= length
            and terminal_length >= params.terminal_min_low_region_length
            and last_start >= int(length * params.terminal_min_start_fraction)
        )
        if terminal_region:
            backtrack_start = max(0, last_start - params.terminal_backtrack_window)
            warning_positions = [
                index
                for index in range(backtrack_start, last_start)
                if clean_qualities[index] < params.good_median_phred or sequence[index] not in "ACGT"
            ]
            if warning_positions:
                backtracked_region = (warning_positions[0], last_end)
                consolidated[-1] = backtracked_region
                terminal_backtracked_ranges.add(backtracked_region)
                consolidated = merge_regions(consolidated, params.merge_gap)
    accepted_blocks = complement_blocks(length, consolidated)
    file_clear_range = normalize_clear_range(*clear_range, sequence_length=length) if clear_range else None
    alignment_blocks = intersect_blocks(accepted_blocks, file_clear_range)
    weaver_clear_range = (accepted_blocks[0][0], accepted_blocks[-1][1]) if accepted_blocks else None
    intermediate_flags = [
        index not in {pos for start, end in consolidated for pos in range(start, end)}
        and (15 <= clean_qualities[index] < 20 or sequence[index] not in "ACGT")
        for index in range(length)
    ]
    intermediate_regions = [
        region_metrics(sequence, clean_qualities, start, end, signal, ["Phred 15-19 or ambiguous base"], "INTERMEDIATE_CONFIDENCE")
        for start, end in contiguous_true_regions(intermediate_flags, 1)
    ]
    low_regions = []
    for start, end in consolidated:
        reasons = ["Window/hysteresis Phred criteria"]
        if (start, end) in terminal_backtracked_ranges:
            reasons.append("Terminal quality collapse backtracked to first nearby low-confidence warning")
        if any(region_start < end and region_end > start for region_start, region_end in signal_regions):
            reasons.append("Chromatogram signal morphology")
        low_regions.append(region_metrics(sequence, clean_qualities, start, end, signal, reasons, "LOW_CONFIDENCE"))
    metrics = {
        "quality_available": True,
        "raw_length": length,
        "trimmed_length": sum(end - start for start, end in alignment_blocks),
        "mean_quality": round(mean(clean_qualities), 2),
        "median_quality": round(median(clean_qualities), 2),
        "q20_percent": round(sum(1 for q in clean_qualities if q >= 20) / len(clean_qualities) * 100, 2),
        "q30_percent": round(sum(1 for q in clean_qualities if q >= 30) / len(clean_qualities) * 100, 2),
        "ambiguous_bases": sum(1 for base in sequence if base not in "ACGT"),
        "ambiguous_fraction": sum(1 for base in sequence if base not in "ACGT") / length if length else 1.0,
        "trim_start": alignment_blocks[0][0] if alignment_blocks else 0,
        "trim_end": alignment_blocks[-1][1] if alignment_blocks else 0,
        "trimmed_left": alignment_blocks[0][0] if alignment_blocks else length,
        "trimmed_right": length - alignment_blocks[-1][1] if alignment_blocks else length,
        "reason": "confidence-region masking",
        "low_confidence_regions": low_regions,
        "intermediate_confidence_regions": intermediate_regions,
        "accepted_blocks": accepted_blocks,
        "alignment_blocks": alignment_blocks,
        "weaver_clear_range": weaver_clear_range,
        "file_clear_range": file_clear_range,
        "signal": {
            "available": signal.get("available", False),
            "warnings": signal_warnings,
        },
    }
    return metrics


def trim_by_quality(sequence, qualities, params, chromatogram=None, clear_range=None):
    metrics = detect_confidence_regions(sequence, qualities, chromatogram, clear_range, params)
    blocks = metrics.get("alignment_blocks") or []
    if not blocks:
        return "", 0, 0, metrics
    start = blocks[0][0]
    end = blocks[-1][1]
    return sequence[start:end], start, end, metrics


def display_trim_range(sequence_length, quality_metrics):
    """Keep the complete raw read in the visual alignment.

    ``alignment_blocks`` remain the trusted evidence used by metrics and
    automated classification. The browser view must also expose low-quality
    bases from the chromatogram, where confidence bands identify them.
    """
    if sequence_length <= 0:
        return 0, 0
    return 0, sequence_length


def local_quality_mask(qualities, params):
    if not qualities:
        return []
    clean_qualities = [q if q is not None else 0 for q in qualities]
    window = max(1, params.quality_window)
    half_window = window // 2
    mask = []
    for index, quality in enumerate(clean_qualities):
        start = max(0, index - half_window)
        end = min(len(clean_qualities), start + window)
        start = max(0, end - window)
        local_mean = mean(clean_qualities[start:end]) if end > start else 0
        mask.append(quality < params.quality_threshold or local_mean < params.quality_threshold)
    return mask


def build_aligner(params):
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = params.match_score
    aligner.mismatch_score = params.mismatch_score
    aligner.open_gap_score = params.open_gap_score
    aligner.extend_gap_score = params.extend_gap_score
    return aligner


def combine_segment_alignments(reference, segment_alignments, orientation):
    ref_len = len(reference)
    covered = set()
    variants = []
    projection = ["-"] * ref_len
    projection_base_indices = [None] * ref_len
    totals = {
        "matches": 0,
        "mismatches": 0,
        "insertions": 0,
        "deletions": 0,
        "ambiguous_bases": 0,
        "aligned_length": 0,
        "aligned_read_bases": 0,
    }
    starts = []
    ends = []
    score = 0
    for alignment in segment_alignments:
        score += alignment.get("score", 0)
        starts.append(alignment.get("start", 0))
        ends.append(alignment.get("end", 0))
        covered.update(alignment.get("covered_positions", []))
        variants.extend(alignment.get("variants", []))
        totals["matches"] += alignment.get("matches", 0)
        totals["mismatches"] += alignment.get("mismatches", 0)
        totals["insertions"] += alignment.get("insertions", 0)
        totals["deletions"] += alignment.get("deletions", 0)
        totals["ambiguous_bases"] += alignment.get("ambiguous_bases", 0)
        totals["aligned_length"] += alignment.get("aligned_length", 0)
        totals["aligned_read_bases"] += len([base for base in alignment.get("oriented_sequence", "") if base != "-"])
        segment_projection = alignment.get("reference_projection", "")
        segment_indices = alignment.get("reference_projection_base_indices", [])
        for coord, base in enumerate(segment_projection):
            if base and base != "-":
                projection[coord] = base
                if coord < len(segment_indices):
                    projection_base_indices[coord] = segment_indices[coord]
    denominator = totals["matches"] + totals["mismatches"] + totals["deletions"] + totals["ambiguous_bases"]
    identity = totals["matches"] / denominator * 100 if denominator else 0
    start = starts[0] if starts else 0
    end = ends[-1] if ends else start
    return {
        "orientation": orientation,
        "best_orientation": orientation,
        "score": round(float(score), 2),
        "other_orientation_score": None,
        "start": start,
        "end": end,
        "query_start": min((item.get("query_start", 0) for item in segment_alignments), default=0),
        "query_end": max((item.get("query_end", 0) for item in segment_alignments), default=0),
        "start_display": start + 1,
        "end_display": end + 1,
        "crosses_origin": bool(starts and ends and (max(ends) < min(starts))),
        "aligned_length": totals["aligned_length"],
        "identity": round(identity, 2),
        "read_coverage": 100.0,
        "plasmid_coverage": round(len(covered) / ref_len * 100, 2) if ref_len else 0,
        "matches": totals["matches"],
        "mismatches": totals["mismatches"],
        "insertions": totals["insertions"],
        "deletions": totals["deletions"],
        "ambiguous_bases": totals["ambiguous_bases"],
        "covered_positions": sorted(covered),
        "ref_alignment": "",
        "read_alignment": "",
        "reference_projection": "".join(projection),
        "reference_projection_base_indices": projection_base_indices,
        "oriented_sequence": "",
        "variants": variants,
        "segments": [
            {
                "start": item.get("start", 0),
                "end": item.get("end", 0),
                "start_display": item.get("start_display", 0),
                "end_display": item.get("end_display", 0),
                "crosses_origin": item.get("crosses_origin", False),
                "query_start": item.get("query_start", 0),
                "query_end": item.get("query_end", 0),
            }
            for item in segment_alignments
        ],
    }


def align_read(reference, read_sequence, qualities, trim_start, params, trim_end=None, trusted_blocks=None, forced_orientation=None):
    if not reference or not read_sequence:
        return None
    trim_end = trim_end if trim_end is not None else trim_start + len(read_sequence)
    if trusted_blocks:
        trusted_blocks = [(max(trim_start, start), min(trim_end, end)) for start, end in trusted_blocks if min(trim_end, end) > max(trim_start, start)]
        if not trusted_blocks:
            return None
        orientation_probe = align_read(reference, read_sequence, qualities, trim_start, params, trim_end, forced_orientation=forced_orientation)
        if not orientation_probe:
            return None
        orientation = orientation_probe.get("best_orientation") or "forward"
        segment_alignments = []
        for block_start, block_end in trusted_blocks:
            segment = read_sequence[block_start - trim_start:block_end - trim_start]
            if len(segment) < params.minimum_trimmed_length:
                continue
            segment_alignment = align_read(
                reference,
                segment,
                qualities,
                block_start,
                params,
                block_end,
                forced_orientation=orientation,
            )
            if segment_alignment:
                segment_alignments.append(segment_alignment)
        if not segment_alignments:
            return None
        return combine_segment_alignments(reference, segment_alignments, orientation)
    aligner = build_aligner(params)
    doubled = reference.upper() + reference.upper()
    candidates = []
    orientation_sequences = (
        ("forward", read_sequence.upper()),
        ("reverse", str(Seq(read_sequence).reverse_complement()).upper()),
    )
    if forced_orientation in ("forward", "reverse"):
        orientation_sequences = [item for item in orientation_sequences if item[0] == forced_orientation]
    for orientation, sequence in orientation_sequences:
        try:
            alignment = next(aligner.align(doubled, sequence))
        except StopIteration:
            continue
        candidates.append((orientation, sequence, alignment))
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[2].score, reverse=True)
    orientation, oriented_sequence, alignment = candidates[0]
    ambiguous_orientation = False
    if len(candidates) > 1 and abs(candidates[0][2].score - candidates[1][2].score) <= params.orientation_score_delta:
        ambiguous_orientation = True

    ref_aligned = str(alignment[0])
    query_aligned = str(alignment[1])
    ref_start = int(alignment.aligned[0][0][0]) if len(alignment.aligned[0]) else 0
    ref_end = int(alignment.aligned[0][-1][1]) if len(alignment.aligned[0]) else ref_start
    query_start = int(alignment.aligned[1][0][0]) if len(alignment.aligned[1]) else 0
    query_end = int(alignment.aligned[1][-1][1]) if len(alignment.aligned[1]) else query_start
    ref_len = len(reference)
    start = ref_start % ref_len
    end = (ref_end - 1) % ref_len if ref_end > ref_start else start
    crosses_origin = ref_end > ref_len and start > end

    matches = mismatches = insertions = deletions = ambiguous = aligned_read_bases = 0
    variants = []
    covered_positions = set()
    projection = ["-"] * ref_len
    projection_base_indices = [None] * ref_len
    low_quality_local = local_quality_mask(qualities, params)
    ref_cursor = ref_start
    query_cursor = query_start
    for ref_base, query_base in zip(ref_aligned, query_aligned):
        coord = ref_cursor % ref_len if ref_base != "-" else (ref_cursor - 1) % ref_len
        original_base_index = None
        if query_base != "-":
            if orientation == "reverse":
                original_base_index = trim_end - 1 - query_cursor
            else:
                original_base_index = trim_start + query_cursor
        quality = None
        if original_base_index is not None and 0 <= original_base_index < len(qualities):
            quality = qualities[original_base_index]
        local_low_quality = (
            bool(low_quality_local)
            and (
                original_base_index is None
                or original_base_index < 0
                or original_base_index >= len(low_quality_local)
                or low_quality_local[original_base_index]
            )
        )
        trusted_for_alignment = not local_low_quality
        low_quality = quality is None or quality < params.high_quality_variant_phred or local_low_quality
        if ref_base != "-":
            if trusted_for_alignment:
                covered_positions.add(coord)
            if query_base != "-":
                projection[coord] = query_base
                projection_base_indices[coord] = original_base_index
            ref_cursor += 1
        if query_base != "-":
            query_cursor += 1
            if trusted_for_alignment:
                aligned_read_bases += 1
        if ref_base == "-":
            if trusted_for_alignment:
                insertions += 1
            variants.append({
                "coordinate": coord,
                "type": "insertion",
                "expected": "",
                "observed": query_base,
                "quality": quality,
                "low_quality": low_quality,
                "base_index": original_base_index,
            })
        elif query_base == "-":
            if trusted_for_alignment:
                deletions += 1
            deletion_base_index = original_base_index
            if deletion_base_index is None and query_cursor > query_start:
                deletion_base_index = trim_end - query_cursor if orientation == "reverse" else trim_start + query_cursor - 1
            variants.append({
                "coordinate": coord,
                "type": "deletion",
                "expected": ref_base,
                "observed": "",
                "quality": quality,
                "low_quality": low_quality,
                "base_index": deletion_base_index,
            })
        elif query_base not in "ACGT" or ref_base not in "ACGT":
            if trusted_for_alignment:
                ambiguous += 1
            variants.append({"coordinate": coord, "type": "ambiguous", "expected": ref_base, "observed": query_base, "quality": quality, "low_quality": low_quality})
        elif ref_base == query_base:
            if trusted_for_alignment:
                matches += 1
        else:
            if trusted_for_alignment:
                mismatches += 1
            variants.append({"coordinate": coord, "type": "mismatch", "expected": ref_base, "observed": query_base, "quality": quality, "low_quality": low_quality})

    aligned_columns = matches + mismatches + insertions + deletions + ambiguous
    identity = matches / (matches + mismatches + deletions + ambiguous) * 100 if (matches + mismatches + deletions + ambiguous) else 0
    read_coverage = aligned_read_bases / len(read_sequence) * 100 if read_sequence else 0
    plasmid_coverage = len(covered_positions) / ref_len * 100 if ref_len else 0
    return {
        "orientation": "ambiguous" if ambiguous_orientation else orientation,
        "best_orientation": orientation,
        "score": round(float(alignment.score), 2),
        "other_orientation_score": round(float(candidates[1][2].score), 2) if len(candidates) > 1 else None,
        "start": start,
        "end": end,
        "query_start": query_start,
        "query_end": query_end,
        "start_display": start + 1,
        "end_display": end + 1,
        "crosses_origin": crosses_origin,
        "aligned_length": aligned_columns,
        "identity": round(identity, 2),
        "read_coverage": round(read_coverage, 2),
        "plasmid_coverage": round(plasmid_coverage, 2),
        "matches": matches,
        "mismatches": mismatches,
        "insertions": insertions,
        "deletions": deletions,
        "ambiguous_bases": ambiguous,
        "covered_positions": sorted(covered_positions),
        "ref_alignment": ref_aligned,
        "read_alignment": query_aligned,
        "reference_projection": "".join(projection),
        "reference_projection_base_indices": projection_base_indices,
        "oriented_sequence": oriented_sequence,
        "variants": variants,
    }


def include_unaligned_display_flanks(alignment, qualities, trim_start, trim_end, reference_length):
    """Expose locally unaligned read flanks as low-confidence insertions.

    These columns are display-only evidence. Trusted coverage, variants, and
    automated classification continue to use the separate trusted alignment.
    """
    if not alignment or not reference_length:
        return alignment
    oriented_sequence = alignment.get("oriented_sequence", "")
    if not oriented_sequence:
        return alignment
    query_start = max(0, min(len(oriented_sequence), int(alignment.get("query_start", 0))))
    query_end = max(query_start, min(len(oriented_sequence), int(alignment.get("query_end", query_start))))
    if query_start == 0 and query_end == len(oriented_sequence):
        return alignment

    orientation = alignment.get("best_orientation") or alignment.get("orientation") or "forward"

    def original_base_index(oriented_index):
        if orientation == "reverse":
            return trim_end - 1 - oriented_index
        return trim_start + oriented_index

    def build_flank_variants(start, end, coordinate):
        variants = []
        for oriented_index in range(start, end):
            base_index = original_base_index(oriented_index)
            quality = qualities[base_index] if 0 <= base_index < len(qualities) else None
            variants.append({
                "coordinate": coordinate,
                "type": "insertion",
                "expected": "",
                "observed": oriented_sequence[oriented_index],
                "quality": quality,
                "low_quality": True,
                "display_only": True,
                "base_index": base_index,
            })
        return variants

    display_flank_variants = []
    if query_start:
        display_flank_variants.extend(
            build_flank_variants(
                0, query_start, (int(alignment.get("start", 0)) - 1) % reference_length
            )
        )
    if query_end < len(oriented_sequence):
        display_flank_variants.extend(
            build_flank_variants(
                query_end, len(oriented_sequence), int(alignment.get("end", 0)) % reference_length
            )
        )
    if not display_flank_variants:
        return alignment

    display_alignment = alignment.copy()
    display_alignment["variants"] = list(alignment.get("variants", [])) + display_flank_variants
    display_alignment["display_unaligned_flank_count"] = len(display_flank_variants)
    return display_alignment


def read_is_usable(trimmed_sequence, quality_metrics, group_errors, params):
    if group_errors:
        return False, group_errors[0]
    usable_length = quality_metrics.get("trimmed_length", len(trimmed_sequence))
    if usable_length < params.minimum_trimmed_length:
        return False, "The usable sequence is too short for reliable alignment."
    if quality_metrics.get("ambiguous_fraction", 1.0) > params.max_ambiguous_fraction:
        return False, "too many ambiguous bases"
    return True, ""


def process_sanger_files(file_objs, reference_sequence, parameters=None):
    params = parameters or SangerProcessingParameters()
    uploaded = uploaded_files_from_request(file_objs, params)
    groups = group_uploaded_files(uploaded)
    reads = []
    for group_name, files in sorted(groups.items()):
        errors = []
        warnings = []
        formats = {}
        for uploaded in files:
            if uploaded.format in formats:
                errors.append("Multiple {} files in group {}".format(uploaded.format.upper(), group_name))
            formats[uploaded.format] = uploaded
            errors.extend(uploaded.errors)
        parsed = []
        for uploaded in files:
            source = parse_file(uploaded)
            uploaded.metadata = source.metadata or {}
            parsed.append(source)
            errors.extend(source.errors)
        selected, source_warnings = choose_source(parsed)
        warnings.extend(source_warnings)
        raw_sequence = selected.sequence if selected else ""
        qualities = selected.qualities if selected else []
        if selected and not qualities:
            warnings.append("Quality scores unavailable; automatic classification confidence is reduced")
        trimmed, trim_start, trim_end, quality_metrics = trim_by_quality(raw_sequence, qualities, params, selected.chromatogram if selected else {}, selected.clear_range if selected else None) if selected else ("", 0, 0, {})
        if quality_metrics.get("signal", {}).get("warnings"):
            warnings.extend(quality_metrics["signal"]["warnings"])
        if quality_metrics.get("file_clear_range"):
            file_start, file_end = quality_metrics["file_clear_range"]
            weaver_start, weaver_end = quality_metrics.get("weaver_clear_range") or (0, 0)
            warnings.append("AB1/PHD clear range {}..{} intersected with Weaver accepted range {}..{} for alignment".format(
                file_start + 1,
                file_end,
                weaver_start + 1 if weaver_end else 0,
                weaver_end,
            ))
        usable, unusable_reason = read_is_usable(trimmed, quality_metrics, errors, params)
        alignment = align_read(
            reference_sequence,
            raw_sequence,
            qualities,
            0,
            params,
            len(raw_sequence),
            trusted_blocks=quality_metrics.get("alignment_blocks"),
        ) if usable else None
        if usable and not alignment:
            usable = False
            unusable_reason = "read did not align"
        if alignment and alignment["orientation"] == "ambiguous":
            warnings.append("Forward and reverse-complement orientation scores are too similar")
        display_alignment = None
        if alignment:
            display_start, display_end = display_trim_range(len(raw_sequence), quality_metrics)
            display_sequence = raw_sequence[display_start:display_end]
            if display_sequence:
                display_alignment = align_read(
                    reference_sequence,
                    display_sequence,
                    qualities,
                    display_start,
                    params,
                    display_end,
                    forced_orientation=alignment.get("best_orientation") or alignment.get("orientation"),
                )
                display_alignment = include_unaligned_display_flanks(
                    display_alignment,
                    qualities,
                    display_start,
                    display_end,
                    len(reference_sequence),
                )
        reads.append({
            "name": group_name,
            "files": files,
            "formats": sorted(set(uploaded.format for uploaded in files if uploaded.format)),
            "parsed_sources": parsed,
            "selected_source": selected.format if selected else "",
            "raw_sequence": raw_sequence,
            "trimmed_sequence": trimmed,
            "trim_start": trim_start,
            "trim_end": trim_end,
            "quality_metrics": quality_metrics,
            "alignment": alignment,
            "display_alignment": display_alignment,
            "warnings": warnings,
            "errors": errors,
            "is_usable": usable,
            "unusable_reason": unusable_reason,
            "chromatogram": selected.chromatogram if selected and selected.chromatogram else {},
        })
    combined = combined_metrics(reference_sequence, reads, params)
    return {
        "parameters": params.as_dict(),
        "uploaded_files": uploaded,
        "reads": reads,
        "combined": combined,
        "classification": classify_run(combined, reads, params),
    }


def coverage_gaps(reference_length, covered):
    gaps = []
    start = None
    for i in range(reference_length):
        if i not in covered and start is None:
            start = i
        elif i in covered and start is not None:
            gaps.append({"start": start, "end": i - 1, "start_display": start + 1, "end_display": i})
            start = None
    if start is not None:
        gaps.append({"start": start, "end": reference_length - 1, "start_display": start + 1, "end_display": reference_length})
    return gaps


def alignment_covered_positions(alignment):
    """Return covered reference coordinates, including for legacy alignments."""
    alignment = alignment or {}
    covered_positions = alignment.get("covered_positions")
    if covered_positions is not None:
        return list(covered_positions)
    projection_base_indices = alignment.get("reference_projection_base_indices") or []
    return [
        coordinate
        for coordinate, base_index in enumerate(projection_base_indices)
        if base_index is not None
    ]


def combined_metrics(reference_sequence, reads, params):
    ref_len = len(reference_sequence)
    depth = [0] * ref_len
    variants = []
    useful = 0
    discarded = 0
    identities = []
    aligned_bases = 0
    conflicts = 0
    for read in reads:
        alignment = read.get("alignment")
        if not read.get("is_usable") or not alignment:
            discarded += 1
            continue
        useful += 1
        identities.append(alignment.get("identity", 0))
        aligned_bases += alignment.get("aligned_length", 0)
        for coord in alignment_covered_positions(alignment):
            if 0 <= coord < ref_len:
                depth[coord] += 1
        for variant in alignment.get("variants", []):
            variant = variant.copy()
            variant["read"] = read["name"]
            variants.append(variant)
    by_coord = {}
    for variant in variants:
        if variant["type"] in ("mismatch", "deletion", "insertion") and not variant.get("low_quality"):
            by_coord.setdefault((variant["coordinate"], variant["type"]), set()).add(variant.get("observed"))
    conflicts = sum(1 for values in by_coord.values() if len(values) > 1)
    covered = {i for i, value in enumerate(depth) if value > 0}
    high_quality_variants = [variant for variant in variants if variant["type"] in ("mismatch", "deletion", "insertion") and not variant.get("low_quality")]
    high_quality_variant_density = len(high_quality_variants) / aligned_bases if aligned_bases else 0
    return {
        "reference_length": ref_len,
        "read_count": len(reads),
        "useful_reads": useful,
        "discarded_reads": discarded,
        "combined_coverage": round(len(covered) / ref_len * 100, 2) if ref_len else 0,
        "depth": depth,
        "mean_identity": round(mean(identities), 2) if identities else 0,
        "min_depth": min(depth) if depth else 0,
        "max_depth": max(depth) if depth else 0,
        "uncovered_regions": coverage_gaps(ref_len, covered),
        "variant_count": len(variants),
        "high_quality_variant_count": len(high_quality_variants),
        "aligned_bases": aligned_bases,
        "high_quality_variant_density": round(high_quality_variant_density * 100, 3),
        "conflict_count": conflicts,
        "variants": variants,
    }


def classify_run(combined, reads, params):
    reasons = []
    if combined["useful_reads"] == 0:
        return {"state": "NO_DATA", "label": "Sin datos utilizables", "reasons": ["No usable Sanger reads were detected"]}
    hq_variant_density = combined.get("high_quality_variant_density", 0) / 100
    hq_variant_count = combined["high_quality_variant_count"]
    if hq_variant_count:
        reasons.append("{} high-confidence mismatch/indel event(s) in the high-quality aligned region ({:.3f}% density)".format(
            hq_variant_count,
            combined.get("high_quality_variant_density", 0),
        ))
    if combined["mean_identity"] < params.min_identity_for_pass:
        reasons.append("Mean identity across the useful Sanger region is {:.2f}%".format(combined["mean_identity"]))
    if combined["conflict_count"]:
        reasons.append("{} conflicting variant position(s) across reads".format(combined["conflict_count"]))
    for read in reads:
        if read.get("warnings"):
            reasons.append("{}: {}".format(read["name"], read["warnings"][0]))
        if read.get("errors"):
            reasons.append("{}: {}".format(read["name"], read["errors"][0]))

    if combined["mean_identity"] < params.min_identity_for_review:
        return {"state": "FAIL", "label": "No verifica", "reasons": reasons}
    if hq_variant_density > params.max_high_quality_variant_density_for_review:
        return {"state": "FAIL", "label": "No verifica", "reasons": reasons}
    if hq_variant_count > params.max_high_quality_variants_for_pass:
        return {"state": "REVIEW", "label": "Requiere revisión", "reasons": reasons}
    if hq_variant_density > params.max_high_quality_variant_density_for_pass:
        return {"state": "REVIEW", "label": "Requiere revisión", "reasons": reasons}
    if reasons:
        return {"state": "REVIEW", "label": "Requiere revisión", "reasons": reasons}
    return {"state": "PASS", "label": "Verifica", "reasons": ["The high-quality Sanger-aligned region is consistent with the expected plasmid sequence"]}


def alignment_tracks_for_ove(plasmid_name, reference_sequence, reads):
    names = [plasmid_name]
    raw_sequences = [reference_sequence]
    aligned_sequences = [reference_sequence]
    chromatos = [{}]
    for read in reads:
        alignment = read.get("alignment")
        if not read.get("is_usable") or not alignment:
            continue
        names.append(read["name"])
        raw_sequences.append(read["trimmed_sequence"])
        aligned_sequences.append(alignment.get("reference_projection", "-" * len(reference_sequence)))
        chromatos.append(read.get("chromatogram") or {})
    return json.dumps([names, raw_sequences, aligned_sequences, chromatos])


def clustal_content(plasmid_name, reference_sequence, reads):
    records = [SeqRecord(Seq(reference_sequence), id=plasmid_name[:50])]
    for read in reads:
        alignment = read.get("alignment")
        if read.get("is_usable") and alignment:
            records.append(SeqRecord(Seq(alignment["reference_projection"]), id=read["name"][:50]))
    output = io.StringIO()
    alignment = MultipleSeqAlignment(records)
    from Bio import AlignIO
    AlignIO.write(alignment, output, "clustal")
    return output.getvalue()


def variants_csv(rows):
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["read", "coordinate_1based", "type", "expected", "observed", "quality", "low_quality"])
    for row in rows:
        writer.writerow([
            row.get("read", ""),
            row.get("coordinate", 0) + 1,
            row.get("type", ""),
            row.get("expected", ""),
            row.get("observed", ""),
            "" if row.get("quality") is None else row.get("quality"),
            row.get("low_quality", False),
        ])
    return output.getvalue()


def read_metrics_tsv(reads):
    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(["read", "formats", "source", "usable", "orientation", "start", "end", "identity", "raw_length", "trimmed_length", "warnings", "errors"])
    for read in reads:
        alignment = read.get("alignment") or {}
        quality = read.get("quality_metrics") or {}
        writer.writerow([
            read["name"],
            ",".join(read.get("formats", [])),
            read.get("selected_source", ""),
            read.get("is_usable", False),
            alignment.get("orientation", "unmapped"),
            alignment.get("start_display", ""),
            alignment.get("end_display", ""),
            alignment.get("identity", ""),
            quality.get("raw_length", ""),
            quality.get("trimmed_length", ""),
            "; ".join(read.get("warnings", [])),
            "; ".join(read.get("errors", [])),
        ])
    return output.getvalue()
