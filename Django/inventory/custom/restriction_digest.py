from dataclasses import dataclass
from itertools import combinations

from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.Restriction.Restriction_Dictionary import rest_dict


DEFAULT_MIN_FRAGMENTS = 2
DEFAULT_MAX_FRAGMENTS = 3
DEFAULT_MIN_BAND_DIFFERENCE_BP = 500
DEFAULT_MIN_FRAGMENT_SIZE_BP = 500
DEFAULT_MIN_BUFFER_ACTIVITY_PERCENT = 100
DEFAULT_MAX_ENZYMES = 2
DEFAULT_RESULT_LIMIT = 10

LEGACY_BUFFER_ACTIVITY_FIELDS = (
    ("activity_buffer_1_1", "NEB 1.1"),
    ("activity_buffer_2_1", "NEB 2.1"),
    ("activity_buffer_3_1", "NEB 3.1"),
    ("activity_buffer_CS", "NEB CutSmart"),
    ("activity_buffer_aari", "Thermo AarI"),
)


@dataclass(frozen=True)
class DigestRegion:
    start: int
    end: int


@dataclass(frozen=True)
class DigestConstraints:
    min_fragments: int = DEFAULT_MIN_FRAGMENTS
    max_fragments: int = DEFAULT_MAX_FRAGMENTS
    min_band_difference_bp: int = DEFAULT_MIN_BAND_DIFFERENCE_BP
    min_fragment_size_bp: int = DEFAULT_MIN_FRAGMENT_SIZE_BP
    min_buffer_activity_percent: int = DEFAULT_MIN_BUFFER_ACTIVITY_PERCENT
    max_enzymes: int = DEFAULT_MAX_ENZYMES
    limit: int = DEFAULT_RESULT_LIMIT
    required_regions: tuple = ()
    required_enzymes: tuple = ()


@dataclass(frozen=True)
class LabEnzyme:
    name: str
    display_name: str
    recognition_site: str
    fcut: int
    rcut: int
    temperature: int
    activities: dict


def normalize_regions(regions, sequence_length):
    normalized = []
    for region in regions or []:
        start = int(region["start"])
        end = int(region["end"])
        if start < 0 or end < 0 or start >= sequence_length or end >= sequence_length:
            raise ValueError("Digest regions must use coordinates inside the sequence.")
        normalized.append(DigestRegion(start=start, end=end))
    return tuple(normalized)


def lab_enzyme_from_model(enzyme):
    enzyme_info = rest_dict.get(enzyme.name) or {}
    return LabEnzyme(
        name=enzyme.name,
        display_name=str(enzyme),
        recognition_site=enzyme_info.get("site") or "",
        fcut=enzyme_info.get("fst5"),
        rcut=(enzyme_info.get("size") or 0) + enzyme_info.get("fst3") if enzyme_info.get("fst3") is not None else None,
        temperature=enzyme_info.get("opt_temp"),
        activities=enzyme_activity_map(enzyme),
    )


def enzyme_activity_map(enzyme):
    activity_map = getattr(enzyme, "buffer_activity_map", None)
    if callable(activity_map):
        return activity_map()
    if isinstance(activity_map, dict):
        return activity_map

    activities = {}
    for field_name, buffer_name in LEGACY_BUFFER_ACTIVITY_FIELDS:
        activity = getattr(enzyme, field_name, None)
        if activity is not None:
            activities[buffer_name] = activity
    return activities


def effective_cut_sites(sequence, enzyme, is_circular=True):
    sequence = str(sequence).upper()
    if not sequence:
        return []
    search_results = RestrictionBatch([enzyme.name]).search(Seq(sequence), linear=not is_circular)
    raw_positions = []
    for key, positions in search_results.items():
        if str(key) == enzyme.name:
            raw_positions = positions
            break
    cut_sites = []
    for position in raw_positions:
        cut_position = (int(position) - 1) % len(sequence) if is_circular else int(position) - 1
        if 0 <= cut_position < len(sequence):
            cut_sites.append({
                "position": cut_position,
                "position_label": cut_position + 1,
                "enzyme": enzyme.display_name,
                "enzyme_name": enzyme.name,
            })
    cut_sites.sort(key=lambda site: (site["position"], site["enzyme_name"]))
    return cut_sites


def enzymes_with_effective_cuts(sequence, enzymes, is_circular=True):
    sequence = str(sequence).upper()
    cutting_enzymes = []
    for enzyme in enzymes:
        if not rest_dict.get(enzyme.name):
            continue
        lab_enzyme = lab_enzyme_from_model(enzyme)
        if effective_cut_sites(sequence, lab_enzyme, is_circular=is_circular):
            cutting_enzymes.append(enzyme)
    return cutting_enzymes


def grouped_cut_sites(enzyme_cut_sites):
    grouped = {}
    for site in enzyme_cut_sites:
        grouped.setdefault(site["position"], {
            "position": site["position"],
            "position_label": site["position"] + 1,
            "enzymes": [],
        })
        grouped[site["position"]]["enzymes"].append(site["enzyme"])
    result = list(grouped.values())
    for site in result:
        site["enzymes"] = sorted(set(site["enzymes"]))
    return sorted(result, key=lambda site: site["position"])


def fragment_sizes(sequence_length, cut_positions, is_circular=True):
    positions = sorted(set(cut_positions))
    if is_circular:
        if not positions:
            return []
        if len(positions) == 1:
            return [sequence_length]
        return [
            (positions[(idx + 1) % len(positions)] - position) % sequence_length or sequence_length
            for idx, position in enumerate(positions)
        ]
    boundaries = [0] + positions + [sequence_length]
    return [
        boundaries[idx + 1] - boundaries[idx]
        for idx in range(len(boundaries) - 1)
    ]


def min_band_difference(fragment_lengths):
    if len(fragment_lengths) < 2:
        return None
    ordered = sorted(fragment_lengths)
    return min(ordered[idx + 1] - ordered[idx] for idx in range(len(ordered) - 1))


def region_contains_position(region, position, sequence_length):
    if region.start <= region.end:
        return region.start <= position <= region.end
    return position >= region.start or position <= region.end


def covered_regions(cut_positions, regions, sequence_length):
    covered = []
    for idx, region in enumerate(regions):
        is_covered = any(region_contains_position(region, position, sequence_length) for position in cut_positions)
        covered.append({
            "index": idx + 1,
            "start": region.start,
            "end": region.end,
            "start_label": region.start + 1,
            "end_label": region.end + 1,
            "covered": is_covered,
        })
    return covered


def compatible_temperature(enzyme_group):
    temperatures = [enzyme.temperature for enzyme in enzyme_group]
    if any(temperature is None for temperature in temperatures):
        return None, False
    return temperatures[0], all(temperature == temperatures[0] for temperature in temperatures)


def compatible_buffers(enzyme_group, min_activity):
    shared_buffer_names = None
    for enzyme in enzyme_group:
        enzyme_buffer_names = set(enzyme.activities.keys())
        if shared_buffer_names is None:
            shared_buffer_names = enzyme_buffer_names
        else:
            shared_buffer_names &= enzyme_buffer_names

    if not shared_buffer_names:
        return []

    buffers = []
    for buffer_name in sorted(shared_buffer_names):
        activities = []
        for enzyme in enzyme_group:
            activity = enzyme.activities.get(buffer_name)
            activities.append(activity)
        if any(activity is None for activity in activities):
            continue
        if all(activity >= min_activity for activity in activities):
            buffers.append({
                "key": buffer_name,
                "name": buffer_name,
                "activities": {
                    enzyme.display_name: enzyme.activities[buffer_name]
                    for enzyme in enzyme_group
                },
                "min_activity": min(activities),
                "avg_activity": round(sum(activities) / len(activities), 1),
            })
    buffers.sort(key=lambda buffer: (-buffer["min_activity"], -buffer["avg_activity"], buffer["name"]))
    return buffers


def violation(label, actual, target, missing, weight):
    return {
        "criterion": label,
        "actual": actual,
        "target": target,
        "missing": missing,
        "message": f"{label}: {actual}, missing {missing}",
        "score": weight * (missing / target if target else 1),
    }


def evaluate_digest(sequence, enzyme_group, constraints, is_circular=True):
    sequence_length = len(sequence)
    enzyme_sites = []
    per_enzyme_positions = {}
    for enzyme in enzyme_group:
        sites = effective_cut_sites(sequence, enzyme, is_circular=is_circular)
        enzyme_sites.extend(sites)
        per_enzyme_positions[enzyme.name] = tuple(site["position"] for site in sites)

    cut_sites = grouped_cut_sites(enzyme_sites)
    cut_positions = [site["position"] for site in cut_sites]
    enzyme_cut_counts = {
        enzyme.display_name: len(set(per_enzyme_positions.get(enzyme.name, ())))
        for enzyme in enzyme_group
    }
    fragments = fragment_sizes(sequence_length, cut_positions, is_circular=is_circular)
    sorted_fragments = sorted(fragments)
    min_fragment_size = min(fragments) if fragments else 0
    band_difference = min_band_difference(fragments)
    region_status = covered_regions(cut_positions, constraints.required_regions, sequence_length)
    temperature, temperature_ok = compatible_temperature(enzyme_group)
    buffers = compatible_buffers(enzyme_group, constraints.min_buffer_activity_percent)
    best_buffer = buffers[0] if buffers else None

    violations = []
    fragment_count = len(fragments)
    if fragment_count < constraints.min_fragments:
        missing = constraints.min_fragments - len(fragments)
        fragment_violation = violation("fragment count", fragment_count, constraints.min_fragments, missing, 0.8)
        fragment_violation["message"] = f"fragment count: {fragment_count}, below range {constraints.min_fragments}-{constraints.max_fragments}"
        violations.append(fragment_violation)
    if constraints.max_fragments and fragment_count > constraints.max_fragments:
        extra = fragment_count - constraints.max_fragments
        fragment_violation = violation("fragment count", fragment_count, constraints.max_fragments, extra, 0.8)
        fragment_violation["message"] = f"fragment count: {fragment_count}, above range {constraints.min_fragments}-{constraints.max_fragments}"
        violations.append(fragment_violation)
    if min_fragment_size < constraints.min_fragment_size_bp:
        missing = constraints.min_fragment_size_bp - min_fragment_size
        fragment_size_violation = violation("smallest fragment bp", min_fragment_size, constraints.min_fragment_size_bp, missing, 1.1)
        fragment_size_violation["message"] = f"smallest fragment: {min_fragment_size} bp"
        violations.append(fragment_size_violation)
    actual_difference = band_difference if band_difference is not None else 0
    if len(fragments) >= 2 and actual_difference < constraints.min_band_difference_bp:
        missing = constraints.min_band_difference_bp - actual_difference
        band_violation = violation("minimum band separation bp", actual_difference, constraints.min_band_difference_bp, missing, 0.9)
        band_violation["message"] = f"minimum band separation: {actual_difference} bp"
        violations.append(band_violation)
    uncovered = [region for region in region_status if not region["covered"]]
    if uncovered:
        violations.append({
            "criterion": "required regions",
            "actual": len(region_status) - len(uncovered),
            "target": len(region_status),
            "missing": len(uncovered),
            "message": f"required regions covered: {len(region_status) - len(uncovered)}/{len(region_status)}, missing {len(uncovered)}",
            "score": 2.0 * len(uncovered),
        })
    if not buffers:
        violations.append({
            "criterion": "shared buffer",
            "actual": "none",
            "target": f">= {constraints.min_buffer_activity_percent}%",
            "missing": "compatible buffer",
            "message": f"no shared buffer reaches {constraints.min_buffer_activity_percent}% for every enzyme",
            "score": 2.0,
        })
    if not temperature_ok:
        missing = "known shared temperature"
        actual = "unknown" if temperature is None else "incompatible"
        violations.append({
            "criterion": "reaction temperature",
            "actual": actual,
            "target": "same optimum",
            "missing": missing,
            "message": "reaction temperature is unknown or incompatible",
            "score": 1.6,
        })

    redundant = False
    if len(enzyme_group) > 1:
        union_positions = set(cut_positions)
        redundant = any(set(positions) == union_positions for positions in per_enzyme_positions.values())

    exact = not violations and not redundant
    if redundant:
        violations.append({
            "criterion": "redundant enzymes",
            "actual": "same cut pattern",
            "target": "added cut sites",
            "missing": "unique contribution",
            "message": "one enzyme does not add a new effective cut position",
            "score": 0.4,
        })

    closeness_score = round(sum(item["score"] for item in violations), 4)
    enzyme_names = [enzyme.display_name for enzyme in enzyme_group]
    return {
        "id": "weaver-digest-" + "-".join(enzyme.name for enzyme in enzyme_group),
        "enzymes": enzyme_names,
        "enzyme_names": [enzyme.name for enzyme in enzyme_group],
        "recognition_sites": {enzyme.display_name: enzyme.recognition_site for enzyme in enzyme_group},
        "enzyme_cut_counts": enzyme_cut_counts,
        "cut_sites": cut_sites,
        "cut_count": len(cut_sites),
        "fragment_count": len(fragments),
        "fragments_map_order": fragments,
        "fragments_by_size": sorted_fragments,
        "min_fragment_size": min_fragment_size,
        "min_band_difference": band_difference,
        "regions": region_status,
        "best_buffer": best_buffer,
        "compatible_buffers": buffers,
        "temperature": temperature if temperature_ok else None,
        "temperature_status": "compatible" if temperature_ok else "unknown" if temperature is None else "incompatible",
        "exact": exact,
        "status": "Exact match" if exact else "Closest match",
        "violations": violations,
        "closeness_score": closeness_score,
        "is_circular": is_circular,
    }


def result_sort_key(result, constraints):
    enzyme_count = len(result["enzymes"])
    best_buffer = result["best_buffer"] or {}
    extra_fragments = max(0, result["fragment_count"] - constraints.min_fragments)
    min_difference = result["min_band_difference"] if result["min_band_difference"] is not None else -1
    names = "+".join(result["enzymes"])
    if result["exact"]:
        return (
            0,
            enzyme_count,
            -(best_buffer.get("min_activity") or 0),
            -(best_buffer.get("avg_activity") or 0),
            -min_difference,
            extra_fragments,
            names,
        )
    return (
        1,
        result["closeness_score"],
        enzyme_count,
        names,
    )


def digest_candidates(sequence, enzymes, constraints=None, is_circular=True):
    constraints = constraints or DigestConstraints()
    sequence = str(sequence).upper()
    lab_enzymes = [lab_enzyme_from_model(enzyme) for enzyme in enzymes if rest_dict.get(enzyme.name)]
    max_enzymes = max(1, min(int(constraints.max_enzymes), len(lab_enzymes)))
    required_enzyme_names = {
        str(name).strip().lower()
        for name in constraints.required_enzymes
        if str(name).strip()
    }
    results = []
    seen_patterns = set()
    for enzyme_count in range(1, max_enzymes + 1):
        for enzyme_group in combinations(lab_enzymes, enzyme_count):
            if required_enzyme_names and not required_enzyme_names.issubset({
                enzyme.name.lower()
                for enzyme in enzyme_group
            }):
                continue
            result = evaluate_digest(sequence, enzyme_group, constraints, is_circular=is_circular)
            if any(count == 0 for count in result["enzyme_cut_counts"].values()):
                continue
            cut_pattern = tuple(site["position"] for site in result["cut_sites"])
            if enzyme_count > 1 and cut_pattern in seen_patterns and result["exact"]:
                continue
            seen_patterns.add(cut_pattern)
            results.append(result)
    results.sort(key=lambda result: result_sort_key(result, constraints))
    return results[:constraints.limit]


def serialize_digest_response(sequence, enzymes, constraints=None, is_circular=True):
    constraints = constraints or DigestConstraints()
    results = digest_candidates(sequence, enzymes, constraints, is_circular=is_circular)
    return {
        "results": results,
        "count": len(results),
        "exact_count": sum(1 for result in results if result["exact"]),
        "coordinate_system": "0-based internal, 1-based labels",
        "defaults": {
            "min_fragments": DEFAULT_MIN_FRAGMENTS,
            "max_fragments": DEFAULT_MAX_FRAGMENTS,
            "min_band_difference_bp": DEFAULT_MIN_BAND_DIFFERENCE_BP,
            "min_fragment_size_bp": DEFAULT_MIN_FRAGMENT_SIZE_BP,
            "min_buffer_activity_percent": DEFAULT_MIN_BUFFER_ACTIVITY_PERCENT,
            "max_enzymes": DEFAULT_MAX_ENZYMES,
        },
    }
