import re

from Bio.Seq import Seq

from inventory.custom.standards import assembly_standards
from inventory.custom.primer_dimers import PrimerDimerConditions
from inventory.custom.primer_dimers import PrimerDimerThresholds
from inventory.custom.primer_dimers import PrimerInput
from inventory.custom.primer_dimers import analyze_pair


TYPE_IIS_CLONING_SITES = ("CGTCTC", "GAGACG", "GGTCTC", "GAGACC", "GCTCTTC", "GAAGAGC")
MAX_CLONING_SITE_PREFIX = 8
TYPE_IIS_SPACER_LENGTHS = range(1, 9)
TYPE_IIS_OVERHANG_LENGTHS = (4, 3)
PREFERRED_TYPE_IIS_SPACER_LENGTH = 5
MIN_INFERRED_HYBRIDIZING_LENGTH = 15
MIN_PRIMER_3PRIME_BINDING_LENGTH = 15
YTK_OVERHANGS = {
    "8-1": {"name": "YTK Part 8 / Part 1", "overhang": "CCCT"},
    "1-2": {"name": "YTK Part 1 / Part 2", "overhang": "AACG"},
    "2-3": {"name": "YTK Part 2 / Part 3", "overhang": "TATG"},
    "3a-3b": {"name": "YTK Part 3a / Part 3b", "overhang": "TTCT"},
    "3-4": {"name": "YTK Part 3 / Part 4", "overhang": "ATCC"},
    "4a-4b": {"name": "YTK Part 4a / Part 4b", "overhang": "TGGC"},
    "4-5": {"name": "YTK Part 4 / Part 5", "overhang": "GCTG"},
    "5-6": {"name": "YTK Part 5 / Part 6", "overhang": "TACA"},
    "6-7": {"name": "YTK Part 6 / Part 7", "overhang": "GAGT"},
    "7-8": {"name": "YTK Part 7 / Part 8", "overhang": "CCGA"},
    "8a-8b": {"name": "YTK Part 8a / Part 8b", "overhang": "CAAT"},
}


def clean_dna(sequence):
    return re.sub(r"[^ACGTRYSWKMBDHVN]", "", str(sequence).upper())


def gc_content(sequence):
    sequence = clean_dna(sequence)
    if not sequence:
        return 0
    return (sequence.count("G") + sequence.count("C")) / len(sequence) * 100


def tm_value(sequence):
    sequence = clean_dna(sequence)
    if not sequence:
        return 0
    a = sequence.count("A")
    c = sequence.count("C")
    t = sequence.count("T")
    g = sequence.count("G")
    if len(sequence) < 14:
        return (a + t) * 2 + (c + g) * 4
    return 64.9 + 41 * (g + c - 16.4) / (a + t + c + g)


def recommended_annealing_tm(primer_tm, product_tm):
    return 0.3 * primer_tm + 0.7 * product_tm - 14.9


def known_cloning_overhangs():
    overhangs = set()
    for standard in assembly_standards.values():
        for level in standard.get("ohs", {}).values():
            for overhang in level.values():
                sequence = clean_dna(overhang.get("oh", ""))
                if sequence:
                    overhangs.add(sequence)
                    overhangs.add(str(Seq(sequence).reverse_complement()))
    return overhangs


KNOWN_CLONING_OVERHANGS = known_cloning_overhangs()


def ytk_overhang_lookup():
    lookup = {}
    for key, value in YTK_OVERHANGS.items():
        overhang = clean_dna(value["overhang"])
        reverse_complement = str(Seq(overhang).reverse_complement())
        lookup[overhang] = {
            "key": key,
            "name": value["name"],
            "canonical_overhang": overhang,
            "orientation": "forward",
        }
        lookup[reverse_complement] = {
            "key": key,
            "name": value["name"],
            "canonical_overhang": overhang,
            "orientation": "reverse_complement",
        }
    return lookup


YTK_OVERHANG_LOOKUP = ytk_overhang_lookup()


def classify_ytk_overhang(overhang):
    return YTK_OVERHANG_LOOKUP.get(clean_dna(overhang), {
        "key": "",
        "name": "",
        "canonical_overhang": "",
        "orientation": "",
    })


def infer_type_iis_overhang(sequence, min_hybridizing_length=MIN_INFERRED_HYBRIDIZING_LENGTH):
    sequence = clean_dna(sequence)
    if not sequence:
        return None

    candidates = []
    for site in TYPE_IIS_CLONING_SITES:
        site_index = sequence.find(site)
        if site_index < 0 or site_index > MAX_CLONING_SITE_PREFIX:
            continue

        for spacer_length in TYPE_IIS_SPACER_LENGTHS:
            for overhang_length in TYPE_IIS_OVERHANG_LENGTHS:
                split_index = site_index + len(site) + spacer_length + overhang_length
                if len(sequence) - split_index < min_hybridizing_length:
                    continue

                cloning_overhang = sequence[split_index - overhang_length:split_index]
                if cloning_overhang not in KNOWN_CLONING_OVERHANGS:
                    continue

                candidates.append({
                    "sequence_5": sequence[:split_index],
                    "sequence_3": sequence[split_index:],
                    "site": site,
                    "cloning_overhang": cloning_overhang,
                    "ytk": classify_ytk_overhang(cloning_overhang),
                    "spacer_length": spacer_length,
                    "overhang_length": overhang_length,
                })

    if not candidates:
        return None

    return sorted(candidates, key=lambda candidate: (
        abs(candidate["spacer_length"] - PREFERRED_TYPE_IIS_SPACER_LENGTH),
        -candidate["overhang_length"],
        -len(candidate["sequence_5"]),
        candidate["spacer_length"],
    ))[0]


def inferred_primer_parts(primer):
    sequence_5 = clean_dna(primer.sequence_5)
    sequence_3 = clean_dna(primer.sequence_3)
    if sequence_5 or not sequence_3:
        return {
            "sequence_5": sequence_5,
            "sequence_3": sequence_3,
            "inferred": False,
            "cloning_overhang": "",
            "site": "",
            "ytk": classify_ytk_overhang(""),
        }

    inferred = infer_type_iis_overhang(sequence_3)
    if not inferred:
        return {
            "sequence_5": sequence_5,
            "sequence_3": sequence_3,
            "inferred": False,
            "cloning_overhang": "",
            "site": "",
            "ytk": classify_ytk_overhang(""),
        }

    return {
        "sequence_5": inferred["sequence_5"],
        "sequence_3": inferred["sequence_3"],
        "inferred": True,
        "cloning_overhang": inferred["cloning_overhang"],
        "site": inferred["site"],
        "ytk": inferred["ytk"],
    }


def primer_full_sequence(primer):
    parts = inferred_primer_parts(primer)
    return clean_dna(f"{parts['sequence_5']}{parts['sequence_3']}")


def primer_dimer_analysis(primer_f, primer_r, annealing_temp_c=60.0):
    fwd_parts = inferred_primer_parts(primer_f)
    rev_parts = inferred_primer_parts(primer_r)
    primer_a = PrimerInput(
        name=display_primer_name(primer_f) or str(primer_f.name or "Forward primer"),
        sequence=primer_full_sequence(primer_f),
        hybridizing_sequence=fwd_parts["sequence_3"],
    )
    primer_b = PrimerInput(
        name=display_primer_name(primer_r) or str(primer_r.name or "Reverse primer"),
        sequence=primer_full_sequence(primer_r),
        hybridizing_sequence=rev_parts["sequence_3"],
    )
    return analyze_pair(
        primer_a,
        primer_b,
        "heterodimer",
        PrimerDimerConditions(annealing_temp_c=annealing_temp_c),
        PrimerDimerThresholds(),
    )


def display_primer_name(primer):
    return re.sub(r"^\d+[-_\s.]*", "", primer.name or "")


def display_primer_id(primer):
    match = re.match(r"^(\d+)", primer.name or "")
    if match:
        return int(match.group(1))
    return None


def find_literal_hits(sequence, query):
    query = clean_dna(query)
    if not query:
        return []

    hits = []
    start = 0
    sequence = sequence.upper()
    while True:
        index = sequence.find(query, start)
        if index == -1:
            break
        hits.append({
            "start": index,
            "end": index + len(query) - 1,
            "length": len(query),
        })
        start = index + 1
    return hits


def find_primer_binding_hits(sequence, primer_sequence, is_reverse=False, min_binding_length=MIN_PRIMER_3PRIME_BINDING_LENGTH):
    primer_sequence = clean_dna(primer_sequence)
    if not primer_sequence:
        return []

    minimum_length = min(int(min_binding_length), len(primer_sequence))
    for binding_length in range(len(primer_sequence), minimum_length - 1, -1):
        binding_sequence = primer_sequence[-binding_length:]
        query = str(Seq(binding_sequence).reverse_complement()) if is_reverse else binding_sequence
        hits = find_literal_hits(sequence, query)
        if not hits:
            continue

        unmatched_5 = primer_sequence[:-binding_length]
        for hit in hits:
            hit["binding_length"] = binding_length
            hit["binding_sequence"] = binding_sequence
            hit["query"] = query
            hit["unmatched_5"] = unmatched_5
        return hits

    return []


def primer_hit_partial_binding(hit):
    return bool(hit.get("unmatched_5"))


def partial_binding_warning(label, hit):
    unmatched = hit.get("unmatched_5", "")
    if not unmatched:
        return ""
    return f"{label} partial 3' binding: {hit.get('binding_length', 0)} bp aligned, {len(unmatched)} bp 5' unaligned"


def bases_are_complementary(base_a, base_b):
    return (
        (base_a == "A" and base_b == "T") or
        (base_a == "T" and base_b == "A") or
        (base_a == "C" and base_b == "G") or
        (base_a == "G" and base_b == "C")
    )


def matching_runs(sequence_a, sequence_b):
    runs = []
    for offset in range(-len(sequence_b) + 1, len(sequence_a)):
        current_start = None
        current_length = 0
        for a_index in range(max(0, offset), min(len(sequence_a), offset + len(sequence_b))):
            b_index = a_index - offset
            if bases_are_complementary(sequence_a[a_index], sequence_b[b_index]):
                if current_start is None:
                    current_start = (a_index, b_index)
                current_length += 1
            else:
                if current_length:
                    runs.append({
                        "length": current_length,
                        "a_start": current_start[0],
                        "a_end": a_index - 1,
                        "b_start": current_start[1],
                        "b_end": b_index - 1,
                    })
                current_start = None
                current_length = 0

        if current_length:
            a_end = min(len(sequence_a), offset + len(sequence_b)) - 1
            runs.append({
                "length": current_length,
                "a_start": current_start[0],
                "a_end": a_end,
                "b_start": current_start[1],
                "b_end": a_end - offset,
            })
    return runs


def run_alignment_preview(fwd_sequence, rev_antiparallel, run):
    offset = run["a_start"] - run["b_start"] if run else 0
    fwd_indent = max(0, -offset)
    rev_indent = max(0, offset)
    width = max(fwd_indent + len(fwd_sequence), rev_indent + len(rev_antiparallel))
    match_line = [" "] * width

    for fwd_index, fwd_base in enumerate(fwd_sequence):
        rev_index = fwd_index - offset
        if 0 <= rev_index < len(rev_antiparallel):
            if bases_are_complementary(fwd_base, rev_antiparallel[rev_index]):
                match_line[fwd_indent + fwd_index] = "|"

    fwd_aligned = " " * fwd_indent + fwd_sequence
    rev_aligned = " " * rev_indent + rev_antiparallel
    return {
        "fwd": f"FWD 5' {fwd_aligned.ljust(width)} 3'",
        "match": f"       {''.join(match_line)}",
        "rev": f"REV 3' {rev_aligned.ljust(width)} 5'",
    }


def best_complementarity_run(runs, fwd_length):
    if not runs:
        return None

    both_three_prime_runs = [
        run for run in runs
        if run["a_end"] == fwd_length - 1 and run["b_start"] == 0 and run["length"] >= 3
    ]
    if both_three_prime_runs:
        return sorted(both_three_prime_runs, key=lambda run: (-run["length"], run["a_start"], run["b_start"]))[0]

    three_prime_runs = [
        run for run in runs
        if (run["a_end"] == fwd_length - 1 or run["b_start"] == 0) and run["length"] >= 5
    ]
    if three_prime_runs:
        return sorted(three_prime_runs, key=lambda run: (-run["length"], run["a_start"], run["b_start"]))[0]

    return sorted(runs, key=lambda run: (-run["length"], run["a_start"], run["b_start"]))[0]


def primer_pair_complementarity(primer_f, primer_r):
    fwd_sequence = primer_full_sequence(primer_f)
    rev_sequence = primer_full_sequence(primer_r)
    if not fwd_sequence or not rev_sequence:
        return {
            "severity": "none",
            "max_contiguous": 0,
            "max_3prime_contiguous": 0,
            "max_both_3prime_contiguous": 0,
            "alignment": {},
            "warnings": [],
        }

    rev_antiparallel = rev_sequence[::-1]
    runs = matching_runs(fwd_sequence, rev_antiparallel)
    max_contiguous = max([run["length"] for run in runs] or [0])
    three_prime_runs = [
        run for run in runs
        if run["a_end"] == len(fwd_sequence) - 1 or run["b_start"] == 0
    ]
    both_three_prime_runs = [
        run for run in runs
        if run["a_end"] == len(fwd_sequence) - 1 and run["b_start"] == 0
    ]
    max_3prime_contiguous = max([run["length"] for run in three_prime_runs] or [0])
    max_both_3prime_contiguous = max([run["length"] for run in both_three_prime_runs] or [0])

    if max_both_3prime_contiguous >= 4 or max_3prime_contiguous >= 5:
        severity = "high"
    elif max_both_3prime_contiguous >= 3 or max_contiguous >= 6:
        severity = "medium"
    elif max_contiguous >= 4:
        severity = "low"
    else:
        severity = "none"

    warnings = []
    if max_both_3prime_contiguous >= 3:
        warnings.append(f"{max_both_3prime_contiguous} bp complementarity between both 3' ends")
    if max_3prime_contiguous >= 5 and max_3prime_contiguous > max_both_3prime_contiguous:
        warnings.append(f"{max_3prime_contiguous} bp complementarity involving a 3' end")
    if max_contiguous >= 4 and max_contiguous > max_3prime_contiguous:
        warnings.append(f"{max_contiguous} bp internal complementarity")
    alignment = run_alignment_preview(
        fwd_sequence,
        rev_antiparallel,
        best_complementarity_run(runs, len(fwd_sequence)),
    )

    return {
        "severity": severity,
        "max_contiguous": max_contiguous,
        "max_3prime_contiguous": max_3prime_contiguous,
        "max_both_3prime_contiguous": max_both_3prime_contiguous,
        "alignment": alignment,
        "warnings": warnings,
    }



def circular_distance(start, end, sequence_length):
    return (end - start) % sequence_length


def interval_intersection_length(a_start, a_end, b_start, b_end):
    start = max(a_start, b_start)
    end = min(a_end, b_end)
    if end < start:
        return 0
    return end - start + 1


def unwrap_hits(hits, sequence_length, window_start, window_end):
    unwrapped = []
    for hit in hits:
        for shift in (-sequence_length, 0, sequence_length):
            start = hit["start"] + shift
            end = hit["end"] + shift
            if end >= window_start and start <= window_end:
                unwrapped.append({
                    "start": start,
                    "end": end,
                    "original_start": hit["start"],
                    "original_end": hit["end"],
                    "length": hit["length"],
                    "binding_length": hit.get("binding_length", hit["length"]),
                    "binding_sequence": hit.get("binding_sequence", ""),
                    "unmatched_5": hit.get("unmatched_5", ""),
                })
    return unwrapped


def format_position(position, sequence_length):
    return (position % sequence_length) + 1


def format_range(start, end, sequence_length):
    start_display = format_position(start, sequence_length)
    end_display = format_position(end, sequence_length)
    if end >= sequence_length or start < 0 or start_display > end_display:
        return f"{start_display}..{end_display} (circular)"
    return f"{start_display}..{end_display}"


def circular_sequence_slice(sequence, start, end):
    if not sequence or end < start:
        return ""
    sequence_length = len(sequence)
    return "".join(sequence[position % sequence_length] for position in range(start, end + 1))


def classify_product(product_start, product_end, selection_start, selection_end, sequence_length):
    selection_length = selection_end - selection_start + 1
    product_length = product_end - product_start + 1
    coverage = interval_intersection_length(product_start, product_end, selection_start, selection_end)
    circular = product_start < 0 or product_end >= sequence_length or selection_start < 0 or selection_end >= sequence_length

    if product_start <= selection_start and product_end >= selection_end:
        category = "Full selection"
        rank = 0
    elif product_start >= selection_start and product_end <= selection_end:
        category = "Internal"
        rank = 3
    elif product_start < selection_start and product_end < selection_end and coverage > 0:
        category = "Left junction"
        rank = 1
    elif product_start > selection_start and product_end > selection_end and coverage > 0:
        category = "Right junction"
        rank = 2
    elif coverage > 0:
        category = "Partial"
        rank = 4
    else:
        category = "Nearby"
        rank = 5

    if circular:
        category = "Circular product" if category in ("Full selection", "Nearby") else f"{category} / circular"
        rank = min(rank, 1)

    return category, rank, coverage, selection_length, product_length


def primer_info(primer, hits):
    parts = inferred_primer_parts(primer)
    return {
        "primer": primer,
        "name": display_primer_name(primer),
        "display_id": display_primer_id(primer),
        "sequence_3": parts["sequence_3"],
        "sequence_5": parts["sequence_5"],
        "sequence_5_inferred": parts["inferred"],
        "cloning_overhang": parts["cloning_overhang"],
        "cloning_site": parts["site"],
        "hit_count": len(hits),
        "tm_3": tm_value(parts["sequence_3"]),
        "gc_3": gc_content(parts["sequence_3"]),
    }


def suggest_pcr_primers(
        sequence,
        primers,
        selection_start,
        selection_end,
        margin=300,
        max_results=80,
        min_product_size=100,
        max_tm_difference=5):
    sequence = clean_dna(sequence)
    sequence_length = len(sequence)
    if not sequence:
        return []

    selection_start = int(selection_start) % sequence_length
    selection_end = int(selection_end) % sequence_length
    if selection_end < selection_start:
        selection_end += sequence_length

    window_start = selection_start - int(margin)
    window_end = selection_end + int(margin)
    doubled_sequence = sequence + sequence

    fwd_candidates = []
    rev_candidates = []

    for primer in primers:
        sequence_3 = inferred_primer_parts(primer)["sequence_3"]
        if not sequence_3:
            continue

        if primer.fwd_or_rev == "r":
            hits = find_primer_binding_hits(doubled_sequence, sequence_3, is_reverse=True)
            hits = [hit for hit in hits if hit["start"] < sequence_length]
            unwrapped_hits = unwrap_hits(hits, sequence_length, window_start, window_end)
            if unwrapped_hits:
                rev_candidates.append({
                    "primer": primer,
                    "hits": unwrapped_hits,
                    "all_hits": hits,
                    "info": primer_info(primer, hits),
                })
        else:
            hits = find_primer_binding_hits(doubled_sequence, sequence_3)
            hits = [hit for hit in hits if hit["start"] < sequence_length]
            unwrapped_hits = unwrap_hits(hits, sequence_length, window_start, window_end)
            if unwrapped_hits:
                fwd_candidates.append({
                    "primer": primer,
                    "hits": unwrapped_hits,
                    "all_hits": hits,
                    "info": primer_info(primer, hits),
                })

    results = []
    seen = set()
    for fwd in fwd_candidates:
        for rev in rev_candidates:
            for f_hit in fwd["hits"]:
                for r_hit in rev["hits"]:
                    if r_hit["end"] < f_hit["start"]:
                        continue

                    product_start = f_hit["start"]
                    product_end = r_hit["end"]
                    product_length = product_end - product_start + 1
                    max_reasonable = (selection_end - selection_start + 1) + int(margin) * 2
                    if product_length <= 0 or product_length > max_reasonable:
                        continue

                    key = (
                        fwd["primer"].id,
                        rev["primer"].id,
                        f_hit["original_start"],
                        r_hit["original_start"],
                    )
                    if key in seen:
                        continue
                    seen.add(key)

                    category, rank, coverage, selection_length, template_product_length = classify_product(
                        product_start,
                        product_end,
                        selection_start,
                        selection_end,
                        sequence_length,
                    )

                    warnings = []
                    if fwd["info"]["hit_count"] > 1:
                        warnings.append(f"FWD has {fwd['info']['hit_count']} hits")
                    if rev["info"]["hit_count"] > 1:
                        warnings.append(f"REV has {rev['info']['hit_count']} hits")
                    for warning in (
                        partial_binding_warning("FWD", f_hit),
                        partial_binding_warning("REV", r_hit),
                    ):
                        if warning:
                            warnings.append(warning)
                    template_product_sequence = circular_sequence_slice(sequence, product_start, product_end)
                    less_stable_primer_tm = min(fwd["info"]["tm_3"], rev["info"]["tm_3"])
                    product_tm = tm_value(template_product_sequence)
                    annealing_tm = recommended_annealing_tm(less_stable_primer_tm, product_tm)
                    complementarity = primer_pair_complementarity(fwd["primer"], rev["primer"])
                    primer3_dimer = primer_dimer_analysis(fwd["primer"], rev["primer"], annealing_temp_c=annealing_tm)
                    warnings.extend(complementarity["warnings"])
                    if primer3_dimer["risk"] in ("MODERATE", "HIGH", "CALCULATION_ERROR"):
                        warnings.append(f"Primer3 dimer risk: {primer3_dimer['risk']}")

                    extra_left = max(0, selection_start - product_start)
                    extra_right = max(0, product_end - selection_end)
                    total_product_length = (
                        len(fwd["info"]["sequence_5"]) +
                        len(f_hit.get("unmatched_5", "")) +
                        template_product_length +
                        len(r_hit.get("unmatched_5", "")) +
                        len(rev["info"]["sequence_5"])
                    )
                    if total_product_length < min_product_size:
                        continue
                    tm_difference = abs(fwd["info"]["tm_3"] - rev["info"]["tm_3"])
                    if tm_difference > max_tm_difference:
                        continue

                    results.append({
                        "category": category,
                        "rank": rank,
                        "fwd": fwd["info"],
                        "rev": rev["info"],
                        "fwd_id": fwd["primer"].id,
                        "rev_id": rev["primer"].id,
                        "product_size": total_product_length,
                        "template_product_size": template_product_length,
                        "product_range": format_range(product_start, product_end, sequence_length),
                        "selection_coverage": coverage,
                        "selection_length": selection_length,
                        "coverage_percent": coverage / selection_length * 100 if selection_length else 0,
                        "extra_left": extra_left,
                        "extra_right": extra_right,
                        "tm_difference": tm_difference,
                        "fwd_partial_binding": primer_hit_partial_binding(f_hit),
                        "rev_partial_binding": primer_hit_partial_binding(r_hit),
                        "fwd_binding_length": f_hit.get("binding_length", f_hit["length"]),
                        "rev_binding_length": r_hit.get("binding_length", r_hit["length"]),
                        "fwd_unmatched_5": f_hit.get("unmatched_5", ""),
                        "rev_unmatched_5": r_hit.get("unmatched_5", ""),
                        "primer_complementarity": complementarity,
                        "primer3_dimer": primer3_dimer,
                        "warnings": warnings,
                    })

    results.sort(key=lambda result: (
        result["rank"],
        -result["coverage_percent"],
        result["extra_left"] + result["extra_right"],
        result["product_size"],
        result["fwd"]["name"],
        result["rev"]["name"],
    ))
    return results[:max_results]


def matching_primer_annotations(sequence, primers):
    sequence = clean_dna(sequence)
    sequence_length = len(sequence)
    if not sequence:
        return []

    doubled_sequence = sequence + sequence
    annotations = []
    for primer in primers:
        sequence_3 = inferred_primer_parts(primer)["sequence_3"]
        if not sequence_3:
            continue

        is_reverse = primer.fwd_or_rev == "r"
        hits = find_primer_binding_hits(doubled_sequence, sequence_3, is_reverse=is_reverse)
        hits = [hit for hit in hits if hit["start"] < sequence_length]

        for hit in hits:
            start = hit["start"]
            end = hit["end"]
            overlaps_self = end >= sequence_length
            primer_name = display_primer_name(primer)
            direction_label = "3'<-5'" if is_reverse else "5'->3'"
            annotations.append({
                "id": f"weaver-primer-{primer.id}-{start}-{'r' if is_reverse else 'f'}",
                "name": f"{primer_name} {direction_label}",
                "start": start,
                "end": end % sequence_length if overlaps_self else end,
                "forward": not is_reverse,
                "strand": -1 if is_reverse else 1,
                "arrowheadType": "BOTTOM" if is_reverse else "TOP",
                "primerBindsOn": "3prime",
                "type": "primer_bind",
                "annotationTypePlural": "primers",
                "color": "#6f9ceb" if is_reverse else "#54b86a",
                "notes": {
                    "weaver_primer_id": [str(display_primer_id(primer) or "")],
                    "weaver_primer_uuid": [str(primer.id)],
                    "weaver_hit_count": [str(len(hits))],
                    "direction": [direction_label],
                },
                "bases": sequence_3,
                "bindingBases": hit.get("binding_sequence", sequence_3),
                "overlapsSelf": overlaps_self,
            })

    annotations.sort(key=lambda annotation: (annotation["start"], annotation["end"], annotation["name"]))
    return annotations


def matching_amplicon_annotations(
        sequence,
        primers,
        min_product_size=100,
        max_product_size=None,
        max_results=250,
        max_tm_difference=5):
    sequence = clean_dna(sequence)
    sequence_length = len(sequence)
    if not sequence:
        return []

    max_product_size = max_product_size or sequence_length
    doubled_sequence = sequence + sequence
    fwd_hits = []
    rev_hits = []

    for primer in primers:
        parts = inferred_primer_parts(primer)
        sequence_3 = parts["sequence_3"]
        if not sequence_3:
            continue

        is_reverse = primer.fwd_or_rev == "r"
        hits = find_primer_binding_hits(doubled_sequence, sequence_3, is_reverse=is_reverse)
        hits = [hit for hit in hits if hit["start"] < sequence_length]
        primer_data = {
            "primer": primer,
            "name": display_primer_name(primer),
            "display_id": display_primer_id(primer),
            "sequence_5": parts["sequence_5"],
            "sequence_5_inferred": parts["inferred"],
            "cloning_overhang": parts["cloning_overhang"],
            "cloning_site": parts["site"],
            "tm_3": tm_value(sequence_3),
            "hit_count": len(hits),
        }

        for hit in hits:
            hit_data = {
                "primer": primer_data,
                "start": hit["start"],
                "end": hit["end"],
                "binding_sequence": hit.get("binding_sequence", sequence_3),
                "binding_length": hit.get("binding_length", len(sequence_3)),
                "unmatched_5": hit.get("unmatched_5", ""),
            }
            if is_reverse:
                rev_hits.append(hit_data)
            else:
                fwd_hits.append(hit_data)

    annotations = []
    seen = set()
    for f_hit in fwd_hits:
        for r_hit in rev_hits:
            r_start = r_hit["start"]
            r_end = r_hit["end"]
            if r_end < f_hit["start"]:
                r_start += sequence_length
                r_end += sequence_length

            template_size = r_end - f_hit["start"] + 1
            product_size = (
                len(f_hit["primer"]["sequence_5"]) +
                len(f_hit.get("unmatched_5", "")) +
                template_size +
                len(r_hit.get("unmatched_5", "")) +
                len(r_hit["primer"]["sequence_5"])
            )
            if product_size < min_product_size or product_size > max_product_size:
                continue
            tm_difference = abs(f_hit["primer"]["tm_3"] - r_hit["primer"]["tm_3"])
            if tm_difference > max_tm_difference:
                continue
            less_stable_primer_tm = min(f_hit["primer"]["tm_3"], r_hit["primer"]["tm_3"])
            template_product_sequence = doubled_sequence[f_hit["start"]:r_end + 1]
            amplicon_sequence = (
                f_hit["primer"]["sequence_5"] +
                f_hit.get("unmatched_5", "") +
                template_product_sequence +
                str(Seq(r_hit.get("unmatched_5", "")).reverse_complement()) +
                str(Seq(r_hit["primer"]["sequence_5"]).reverse_complement())
            )
            product_tm = tm_value(template_product_sequence)
            annealing_tm = recommended_annealing_tm(less_stable_primer_tm, product_tm)
            primer3_dimer = primer_dimer_analysis(
                f_hit["primer"]["primer"],
                r_hit["primer"]["primer"],
                annealing_temp_c=annealing_tm,
            )

            start = f_hit["start"]
            end = r_end
            display_end = end % sequence_length if end >= sequence_length else end
            overlaps_self = end >= sequence_length
            visual_start = start
            visual_end = display_end
            visual_size = template_size
            visual_is_alternative = False
            key = (f_hit["primer"]["display_id"], r_hit["primer"]["display_id"], start, display_end)
            if key in seen:
                continue
            seen.add(key)

            fwd_label = f_hit["primer"]["name"]
            rev_label = r_hit["primer"]["name"]
            warnings = []
            if f_hit["primer"]["hit_count"] > 1:
                warnings.append(f"FWD has {f_hit['primer']['hit_count']} hits")
            if r_hit["primer"]["hit_count"] > 1:
                warnings.append(f"REV has {r_hit['primer']['hit_count']} hits")
            for warning in (
                partial_binding_warning("FWD", f_hit),
                partial_binding_warning("REV", r_hit),
            ):
                if warning:
                    warnings.append(warning)
            complementarity = primer_pair_complementarity(f_hit["primer"]["primer"], r_hit["primer"]["primer"])
            warnings.extend(complementarity["warnings"])
            if primer3_dimer["risk"] in ("MODERATE", "HIGH", "CALCULATION_ERROR"):
                warnings.append(f"Primer3 dimer risk: {primer3_dimer['risk']}")

            primer3_tm = primer3_dimer.get("dimer_tm_c")
            primer3_dg = primer3_dimer.get("delta_g_kcal_mol")
            primer3_tm_note = f"{primer3_tm:.1f}" if primer3_tm is not None else ""
            primer3_dg_note = f"{primer3_dg:.2f}" if primer3_dg is not None else ""

            annotations.append({
                "id": f"weaver-amplicon-{f_hit['primer']['display_id']}-{r_hit['primer']['display_id']}-{start}-{display_end}",
                "name": f"{fwd_label} + {rev_label}",
                "start": start,
                "end": display_end,
                "forward": True,
                "strand": 1,
                "type": "misc_feature",
                "annotationTypePlural": "parts",
                "color": "#ac68cc",
                "overlapsSelf": overlaps_self,
                "notes": {
                    "weaver_amplicon": ["true"],
                    "amplicon_sequence": [amplicon_sequence],
                    "product_size": [str(product_size)],
                    "template_size": [str(template_size)],
                    "visual_start": [str(visual_start)],
                    "visual_end": [str(visual_end)],
                    "visual_size": [str(visual_size)],
                    "visual_is_alternative": ["true" if visual_is_alternative else "false"],
                    "fwd_tm": [f"{f_hit['primer']['tm_3']:.1f}"],
                    "rev_tm": [f"{r_hit['primer']['tm_3']:.1f}"],
                    "product_tm": [f"{product_tm:.1f}"],
                    "tm_difference": [f"{tm_difference:.1f}"],
                    "recommended_annealing_tm": [f"{annealing_tm:.1f}"],
                    "partial_primer_binding": [
                        "true" if primer_hit_partial_binding(f_hit) or primer_hit_partial_binding(r_hit) else "false"
                    ],
                    "fwd_partial_binding": ["true" if primer_hit_partial_binding(f_hit) else "false"],
                    "rev_partial_binding": ["true" if primer_hit_partial_binding(r_hit) else "false"],
                    "fwd_binding_length": [str(f_hit.get("binding_length", 0))],
                    "rev_binding_length": [str(r_hit.get("binding_length", 0))],
                    "fwd_unmatched_5": [f_hit.get("unmatched_5", "")],
                    "rev_unmatched_5": [r_hit.get("unmatched_5", "")],
                    # These are the positions at which the primer 3' ends
                    # anneal on the reference sequence (0-based, circular).
                    "fwd_3prime_position": [str(f_hit["end"] % sequence_length)],
                    "rev_3prime_position": [str(r_hit["start"] % sequence_length)],
                    "primer_complementarity_severity": [complementarity["severity"]],
                    "primer_complementarity_max": [str(complementarity["max_contiguous"])],
                    "primer_complementarity_3prime": [str(complementarity["max_3prime_contiguous"])],
                    "primer_complementarity_both_3prime": [str(complementarity["max_both_3prime_contiguous"])],
                    "primer_complementarity_alignment_fwd": [complementarity["alignment"].get("fwd", "")],
                    "primer_complementarity_alignment_match": [complementarity["alignment"].get("match", "")],
                    "primer_complementarity_alignment_rev": [complementarity["alignment"].get("rev", "")],
                    "primer_complementarity_warnings": complementarity["warnings"],
                    "primer3_dimer_risk": [primer3_dimer["risk"]],
                    "primer3_dimer_tm": [primer3_tm_note],
                    "primer3_dimer_dg": [primer3_dg_note],
                    "primer3_dimer_extendable_3prime": ["true" if primer3_dimer["extendable_3prime"] else "false"],
                    "primer3_dimer_3prime_run": [str(primer3_dimer["terminal_3prime_run"])],
                    "primer3_dimer_longest_run": [str(primer3_dimer["longest_complementary_run"])],
                    "primer3_dimer_last5": [str(primer3_dimer["paired_bases_last_5"])],
                    "primer3_dimer_last8": [str(primer3_dimer["paired_bases_last_8"])],
                    "primer3_dimer_reasons": [primer3_dimer["risk_reasons"]],
                    "primer3_dimer_recommendation": [primer3_dimer["recommendation"]],
                    "primer3_dimer_structure": [primer3_dimer["predicted_structure"]],
                    "primer3_dimer_error": [primer3_dimer["primer3_error"]],
                    "fwd_primer_id": [str(f_hit["primer"]["display_id"] or "")],
                    "rev_primer_id": [str(r_hit["primer"]["display_id"] or "")],
                    "fwd_primer": [fwd_label],
                    "rev_primer": [rev_label],
                    "warnings": warnings,
                },
            })

    annotations.sort(key=lambda annotation: (
        int(annotation["notes"]["product_size"][0]),
        annotation["start"],
        annotation["name"],
    ))
    return annotations[:max_results]


def primer_pair_amplicons(sequence, primer_f, primer_r, min_product_size=1, max_product_size=None):
    sequence = clean_dna(sequence)
    sequence_length = len(sequence)
    if not sequence:
        return []

    fwd_parts = inferred_primer_parts(primer_f)
    rev_parts = inferred_primer_parts(primer_r)
    fwd_sequence = fwd_parts["sequence_3"]
    rev_sequence = rev_parts["sequence_3"]
    if not fwd_sequence or not rev_sequence:
        return []

    max_product_size = max_product_size or sequence_length
    doubled_sequence = sequence + sequence
    fwd_hits = [
        hit for hit in find_primer_binding_hits(doubled_sequence, fwd_sequence)
        if hit["start"] < sequence_length
    ]
    rev_hits = [
        hit for hit in find_primer_binding_hits(doubled_sequence, rev_sequence, is_reverse=True)
        if hit["start"] < sequence_length
    ]

    amplicons = []
    seen = set()
    for f_hit in fwd_hits:
        for r_hit in rev_hits:
            r_start = r_hit["start"]
            r_end = r_hit["end"]
            if r_end < f_hit["start"]:
                r_start += sequence_length
                r_end += sequence_length

            template_size = r_end - f_hit["start"] + 1
            product_size = (
                len(fwd_parts["sequence_5"]) +
                len(f_hit.get("unmatched_5", "")) +
                template_size +
                len(r_hit.get("unmatched_5", "")) +
                len(rev_parts["sequence_5"])
            )
            if product_size < min_product_size or product_size > max_product_size:
                continue

            display_end = r_end % sequence_length if r_end >= sequence_length else r_end
            key = (f_hit["start"], display_end)
            if key in seen:
                continue
            seen.add(key)
            less_stable_primer_tm = min(tm_value(fwd_sequence), tm_value(rev_sequence))
            template_product_sequence = doubled_sequence[f_hit["start"]:r_end + 1]
            product_tm = tm_value(template_product_sequence)
            annealing_tm = recommended_annealing_tm(less_stable_primer_tm, product_tm)
            primer3_dimer = primer_dimer_analysis(primer_f, primer_r, annealing_temp_c=annealing_tm)

            amplicons.append({
                "start": f_hit["start"],
                "end": display_end,
                "product_range": format_range(f_hit["start"], r_end, sequence_length),
                "product_size": product_size,
                "template_product_size": template_size,
                "circular": r_end >= sequence_length,
                "fwd_hit_count": len(fwd_hits),
                "rev_hit_count": len(rev_hits),
                "tm_difference": abs(tm_value(fwd_sequence) - tm_value(rev_sequence)),
                "fwd_tm": tm_value(fwd_sequence),
                "rev_tm": tm_value(rev_sequence),
                "recommended_annealing_tm": annealing_tm,
                "partial_primer_binding": primer_hit_partial_binding(f_hit) or primer_hit_partial_binding(r_hit),
                "fwd_partial_binding": primer_hit_partial_binding(f_hit),
                "rev_partial_binding": primer_hit_partial_binding(r_hit),
                "fwd_binding_length": f_hit.get("binding_length", f_hit["length"]),
                "rev_binding_length": r_hit.get("binding_length", r_hit["length"]),
                "fwd_unmatched_5": f_hit.get("unmatched_5", ""),
                "rev_unmatched_5": r_hit.get("unmatched_5", ""),
                "primer3_dimer": primer3_dimer,
            })

    amplicons.sort(key=lambda amplicon: (
        amplicon["product_size"],
        amplicon["start"],
        amplicon["end"],
    ))
    return amplicons


def amplicon_segments(annotation, sequence_length):
    start = int(annotation["start"])
    end = int(annotation["end"])
    if annotation.get("overlapsSelf") or end < start:
        return [(start, sequence_length - 1), (0, end)]
    return [(start, end)]


def amplicons_overlap(left, right, sequence_length):
    for left_start, left_end in amplicon_segments(left, sequence_length):
        for right_start, right_end in amplicon_segments(right, sequence_length):
            if left_start <= right_end and right_start <= left_end:
                return True
    return False


def select_non_overlapping_amplicons(annotations, sequence_length):
    selected = []
    for annotation in annotations:
        if any(amplicons_overlap(annotation, chosen, sequence_length) for chosen in selected):
            continue
        selected.append(annotation)
    return selected
