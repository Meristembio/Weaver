import csv
import itertools
import re
from dataclasses import asdict
from dataclasses import dataclass
from io import StringIO

from Bio import SeqIO


VALID_DNA_RE = re.compile(r"^[ACGT]+$")
RISK_ORDER = {"LOW": 0, "MODERATE": 1, "HIGH": 2, "CALCULATION_ERROR": 3}


@dataclass
class PrimerDimerConditions:
    mv_conc_mM: float = 50.0
    dv_conc_mM: float = 1.5
    dntp_conc_mM: float = 0.6
    primer_conc_nM: float = 250.0
    annealing_temp_c: float | None = None
    calculation_temp_c: float = 37.0
    critical_3prime_bases: int = 8
    max_loop: int = 30


@dataclass
class PrimerDimerThresholds:
    high_3prime_run: int = 5
    high_both_3prime_run: int = 4
    high_tm_margin_c: float = 5.0
    high_delta_g_kcal_mol: float = -9.0
    moderate_3prime_run: int = 4
    moderate_paired_last5: int = 4
    moderate_paired_last8: int = 5
    moderate_tm_margin_c: float = 10.0
    moderate_delta_g_kcal_mol: float = -6.0
    low_internal_run: int = 4
    long_primer_warning_nt: int = 60


@dataclass
class PrimerInput:
    name: str
    sequence: str
    hybridizing_sequence: str = ""


class PrimerDimerInputError(ValueError):
    pass


def primer3_version():
    try:
        import primer3
        return getattr(primer3, "__version__", "unknown")
    except ImportError:
        return "not installed"


def normalize_sequence(sequence):
    return re.sub(r"\s+", "", str(sequence or "")).upper()


def validate_primers(primers, allow_iupac=False):
    if allow_iupac:
        raise PrimerDimerInputError("Ambiguous IUPAC bases are not supported by this analyzer yet.")

    seen_names = set()
    sequence_to_names = {}
    warnings = []
    validated = []
    for primer in primers:
        name = str(primer.name or "").strip()
        sequence = normalize_sequence(primer.sequence)
        if not name:
            raise PrimerDimerInputError("Primer with empty name.")
        if name in seen_names:
            raise PrimerDimerInputError(f"Duplicate primer name: {name}.")
        if not sequence:
            raise PrimerDimerInputError(f"Primer {name} has an empty sequence.")
        if not VALID_DNA_RE.match(sequence):
            raise PrimerDimerInputError(f"Primer {name} contains non-ACGT bases.")
        seen_names.add(name)
        sequence_to_names.setdefault(sequence, []).append(name)
        validated.append(PrimerInput(name=name, sequence=sequence, hybridizing_sequence=normalize_sequence(primer.hybridizing_sequence)))

    for sequence, names in sequence_to_names.items():
        if len(names) > 1:
            warnings.append(f"Identical sequence used by multiple primer names: {', '.join(names)}.")
    return validated, warnings


def primers_from_fasta(handle):
    records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        raise PrimerDimerInputError("No FASTA records found.")
    return [PrimerInput(name=record.id, sequence=str(record.seq)) for record in records]


def primers_from_table(handle):
    text = handle.read()
    if isinstance(text, bytes):
        text = text.decode("utf-8-sig")
    text = str(text)
    if not text.strip():
        raise PrimerDimerInputError("Empty primer table.")
    try:
        dialect = csv.Sniffer().sniff(text[:2048], delimiters=",\t;")
    except csv.Error:
        dialect = csv.excel_tab if "\t" in text.splitlines()[0] else csv.excel
    reader = csv.DictReader(text.splitlines(), dialect=dialect)
    fieldnames = [str(field or "").strip().lower() for field in (reader.fieldnames or [])]
    if "name" in fieldnames and "sequence" in fieldnames:
        primers = []
        for row in reader:
            normalized = {str(key or "").strip().lower(): str(value or "").strip() for key, value in row.items()}
            primers.append(PrimerInput(name=normalized.get("name", ""), sequence=normalized.get("sequence", "")))
        return primers
    if {"pair_id", "primer_role", "name", "sequence"}.issubset(set(fieldnames)):
        primers = []
        for row in reader:
            normalized = {str(key or "").strip().lower(): str(value or "").strip() for key, value in row.items()}
            primers.append(PrimerInput(name=normalized.get("name", ""), sequence=normalized.get("sequence", "")))
        return primers
    raise PrimerDimerInputError("Primer table must contain name and sequence columns.")


def read_primers(handle, input_format):
    if input_format == "fasta":
        return primers_from_fasta(handle)
    if input_format in ("csv", "tsv", "table"):
        return primers_from_table(handle)
    if input_format == "auto":
        content = handle.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8-sig")
        return primers_from_fasta(StringIO(content)) if str(content).lstrip().startswith(">") else primers_from_table(StringIO(content))
    raise PrimerDimerInputError(f"Unsupported input format: {input_format}.")


def primer_combinations(primers):
    for primer in primers:
        yield primer, primer, "homodimer"
    for primer_a, primer_b in itertools.combinations(primers, 2):
        yield primer_a, primer_b, "heterodimer"


def bases_are_complementary(base_a, base_b):
    return (
        (base_a == "A" and base_b == "T") or
        (base_a == "T" and base_b == "A") or
        (base_a == "C" and base_b == "G") or
        (base_a == "G" and base_b == "C")
    )


def scan_alignment(seq_a, seq_b, critical_3prime_bases):
    b_antiparallel = seq_b[::-1]
    alignments = []
    for offset in range(-len(b_antiparallel) + 1, len(seq_a)):
        pairs = []
        run = 0
        longest_run = 0
        terminal_a_run = 0
        terminal_b_run = 0
        for a_index in range(len(seq_a)):
            b_index = a_index - offset
            paired = 0 <= b_index < len(b_antiparallel) and bases_are_complementary(seq_a[a_index], b_antiparallel[b_index])
            if paired:
                pairs.append((a_index, b_index))
                run += 1
                longest_run = max(longest_run, run)
            else:
                run = 0

        pair_set = set(pairs)
        a_terminal = len(seq_a) - 1
        b_terminal = 0
        a_3prime_involved = any(a_index == a_terminal for a_index, _ in pairs)
        b_3prime_involved = any(b_index == b_terminal for _, b_index in pairs)
        while (a_terminal - terminal_a_run, (a_terminal - terminal_a_run) - offset) in pair_set:
            terminal_a_run += 1
        while (offset + terminal_b_run, terminal_b_run) in pair_set:
            terminal_b_run += 1

        a_last5_start = max(0, len(seq_a) - 5)
        a_last8_start = max(0, len(seq_a) - critical_3prime_bases)
        b_last5_end = min(len(b_antiparallel), 5)
        b_last8_end = min(len(b_antiparallel), critical_3prime_bases)
        paired_last5 = sum(1 for a_index, b_index in pairs if a_index >= a_last5_start or b_index < b_last5_end)
        paired_last8 = sum(1 for a_index, b_index in pairs if a_index >= a_last8_start or b_index < b_last8_end)
        alignments.append({
            "offset": offset,
            "total_pairs": len(pairs),
            "longest_run": longest_run,
            "a_3prime_involved": a_3prime_involved,
            "b_3prime_involved": b_3prime_involved,
            "both_3prime_involved": a_3prime_involved and b_3prime_involved,
            "terminal_a_run": terminal_a_run,
            "terminal_b_run": terminal_b_run,
            "terminal_3prime_run": max(terminal_a_run, terminal_b_run),
            "paired_last5": paired_last5,
            "paired_last8": paired_last8,
            "pairs": pairs,
        })
    return alignments


def best_alignment(alignments):
    if not alignments:
        return {}
    return sorted(alignments, key=lambda item: (
        -item["terminal_3prime_run"],
        -item["paired_last5"],
        -item["paired_last8"],
        -item["longest_run"],
        -item["total_pairs"],
        item["offset"],
    ))[0]


def risk_metrics_from_alignments(alignments):
    if not alignments:
        return {}
    metrics = dict(best_alignment(alignments))
    metrics["longest_run"] = max(alignment["longest_run"] for alignment in alignments)
    metrics["total_pairs"] = max(alignment["total_pairs"] for alignment in alignments)
    metrics["paired_last5"] = max(alignment["paired_last5"] for alignment in alignments)
    metrics["paired_last8"] = max(alignment["paired_last8"] for alignment in alignments)
    metrics["terminal_3prime_run"] = max(alignment["terminal_3prime_run"] for alignment in alignments)
    metrics["a_3prime_involved"] = any(alignment["a_3prime_involved"] for alignment in alignments)
    metrics["b_3prime_involved"] = any(alignment["b_3prime_involved"] for alignment in alignments)

    both_3prime_alignments = [
        alignment for alignment in alignments
        if alignment["both_3prime_involved"]
    ]
    if both_3prime_alignments:
        best_both = sorted(
            both_3prime_alignments,
            key=lambda alignment: (
                -min(alignment["terminal_a_run"], alignment["terminal_b_run"]),
                -alignment["terminal_3prime_run"],
                -alignment["paired_last5"],
                -alignment["paired_last8"],
            ),
        )[0]
        metrics["both_3prime_involved"] = True
        metrics["terminal_a_run"] = best_both["terminal_a_run"]
        metrics["terminal_b_run"] = best_both["terminal_b_run"]
    else:
        metrics["both_3prime_involved"] = False
    return metrics


def full_alignment_preview(seq_a, seq_b, alignment):
    b_antiparallel = seq_b[::-1]
    offset = alignment.get("offset", 0)
    a_indent = max(0, -offset)
    b_indent = max(0, offset)
    width = max(a_indent + len(seq_a), b_indent + len(b_antiparallel))
    match_line = [" "] * width
    for a_index, b_index in alignment.get("pairs", []):
        match_line[a_indent + a_index] = "|"
    a_aligned = (" " * a_indent + seq_a).ljust(width)
    b_aligned = (" " * b_indent + b_antiparallel).ljust(width)
    return "\n".join([
        f"A 5' {a_aligned} 3'",
        f"     {''.join(match_line)}",
        f"B 3' {b_aligned} 5'",
    ])


def thermo_value(result, attr):
    value = getattr(result, attr, None)
    return None if value is None else float(value)


def calc_primer3(primer_a, primer_b, interaction_type, conditions):
    try:
        import primer3
    except ImportError as error:
        return {"error": f"primer3-py is not installed: {error}"}

    kwargs = {
        "mv_conc": conditions.mv_conc_mM,
        "dv_conc": conditions.dv_conc_mM,
        "dntp_conc": conditions.dntp_conc_mM,
        "dna_conc": conditions.primer_conc_nM,
        "temp_c": conditions.calculation_temp_c,
        "max_loop": conditions.max_loop,
    }
    try:
        if interaction_type == "homodimer":
            result = primer3.bindings.calc_homodimer(primer_a.sequence, output_structure=True, **kwargs)
            end_a = primer3.bindings.calc_end_stability(primer_a.sequence, primer_b.sequence, **kwargs)
            end_b = end_a
        else:
            result = primer3.bindings.calc_heterodimer(primer_a.sequence, primer_b.sequence, output_structure=True, **kwargs)
            end_a = primer3.bindings.calc_end_stability(primer_a.sequence, primer_b.sequence, **kwargs)
            end_b = primer3.bindings.calc_end_stability(primer_b.sequence, primer_a.sequence, **kwargs)
    except Exception as error:
        return {"error": str(error)}

    return {
        "tm": thermo_value(result, "tm"),
        "dg_kcal_mol": thermo_value(result, "dg") / 1000 if thermo_value(result, "dg") is not None else None,
        "structure_found": bool(getattr(result, "structure_found", False)),
        "predicted_structure": getattr(result, "ascii_structure", "") or "",
        "end_stability_a_dg_kcal_mol": thermo_value(end_a, "dg") / 1000 if thermo_value(end_a, "dg") is not None else None,
        "end_stability_b_dg_kcal_mol": thermo_value(end_b, "dg") / 1000 if thermo_value(end_b, "dg") is not None else None,
    }


def classify_dimer(metrics, thermo, conditions, thresholds):
    if thermo.get("error"):
        return "CALCULATION_ERROR", [thermo["error"]]

    risk = "LOW"
    reasons = []
    terminal_run = metrics["terminal_3prime_run"]
    both_3prime_run = min(metrics["terminal_a_run"], metrics["terminal_b_run"]) if metrics["both_3prime_involved"] else 0
    extendable_3prime = metrics["a_3prime_involved"] or metrics["b_3prime_involved"]
    tm = thermo.get("tm")
    delta_g = thermo.get("dg_kcal_mol")
    tm_margin = conditions.annealing_temp_c - tm if conditions.annealing_temp_c is not None and tm is not None else None
    thermodynamic_high = tm_margin is not None and tm_margin <= thresholds.high_tm_margin_c
    thermodynamic_moderate = tm_margin is not None and tm_margin <= thresholds.moderate_tm_margin_c
    delta_g_high = delta_g is not None and delta_g <= thresholds.high_delta_g_kcal_mol
    delta_g_moderate = delta_g is not None and delta_g <= thresholds.moderate_delta_g_kcal_mol

    def set_risk(new_risk, reason):
        nonlocal risk
        if RISK_ORDER[new_risk] > RISK_ORDER[risk]:
            risk = new_risk
        reasons.append(reason)

    if both_3prime_run >= thresholds.high_both_3prime_run:
        set_risk("HIGH", f"Both 3' ends form an extendable {both_3prime_run} bp interaction.")
    if extendable_3prime and thermodynamic_high:
        set_risk("HIGH", f"Dimer Tm is within {tm_margin:.1f} C of annealing temperature.")
    if extendable_3prime and delta_g_high:
        set_risk("HIGH", f"3' interaction with delta G {delta_g:.2f} kcal/mol.")
    if terminal_run >= thresholds.high_3prime_run and extendable_3prime and (thermodynamic_moderate or delta_g_moderate):
        set_risk("HIGH", f"{terminal_run} bases consecutively complementary from a primer 3' end with thermodynamic support.")

    if risk != "HIGH":
        if both_3prime_run >= 3:
            set_risk("MODERATE", f"Both 3' ends form an extendable {both_3prime_run} bp interaction.")
        if terminal_run >= thresholds.moderate_3prime_run:
            set_risk("MODERATE", f"{terminal_run} bases consecutively complementary from a primer 3' end.")
        if metrics["paired_last5"] >= thresholds.moderate_paired_last5:
            set_risk("MODERATE", f"{metrics['paired_last5']} paired bases within the last 5 nt.")
        if metrics["paired_last8"] >= thresholds.moderate_paired_last8:
            set_risk("MODERATE", f"{metrics['paired_last8']} paired bases within the critical 3' region.")
        if thermodynamic_moderate:
            set_risk("MODERATE", f"Dimer Tm is {tm_margin:.1f} C below annealing temperature.")
        if delta_g_moderate:
            set_risk("MODERATE", f"Delta G {delta_g:.2f} kcal/mol.")
        if metrics["longest_run"] >= thresholds.low_internal_run and not extendable_3prime:
            set_risk("MODERATE", "Stable internal interaction without an extendable 3' end.")
        if extendable_3prime and tm_margin is not None and tm_margin > thresholds.moderate_tm_margin_c:
            reasons.append(f"Dimer Tm is {tm_margin:.1f} C below annealing temperature; thermodynamic stability is low at PCR annealing temperatures.")

    if not reasons:
        reasons.append("No relevant 3' complementarity under configured thresholds.")
    return risk, reasons


def recommendation_for_risk(risk):
    if risk == "HIGH":
        return "Redesign recommended, especially if 3' complementarity is extendable."
    if risk == "MODERATE":
        return "Review experimentally and consider redesign if small products or low efficiency appear."
    if risk == "CALCULATION_ERROR":
        return "Thermodynamic calculation failed; review input and Primer3 error."
    return "No important evidence of dimerization under the analyzed conditions."


def analyze_pair(primer_a, primer_b, interaction_type, conditions, thresholds):
    if len(primer_a.sequence) > thresholds.long_primer_warning_nt or len(primer_b.sequence) > thresholds.long_primer_warning_nt:
        length_warning = f"Primer length exceeds {thresholds.long_primer_warning_nt} nt; Primer3 thermodynamics may be less reliable."
    else:
        length_warning = ""
    alignments = scan_alignment(primer_a.sequence, primer_b.sequence, conditions.critical_3prime_bases)
    display_metrics = best_alignment(alignments)
    risk_metrics = risk_metrics_from_alignments(alignments)
    thermo = calc_primer3(primer_a, primer_b, interaction_type, conditions)
    risk, reasons = classify_dimer(risk_metrics, thermo, conditions, thresholds)
    if length_warning:
        reasons.append(length_warning)
    tm = thermo.get("tm")
    tm_margin = conditions.annealing_temp_c - tm if conditions.annealing_temp_c is not None and tm is not None else None
    return {
        "primer_a": primer_a.name,
        "primer_b": primer_b.name,
        "interaction_type": interaction_type,
        "risk": risk,
        "extendable_3prime": risk_metrics.get("a_3prime_involved", False) or risk_metrics.get("b_3prime_involved", False),
        "primer_a_3prime_involved": risk_metrics.get("a_3prime_involved", False),
        "primer_b_3prime_involved": risk_metrics.get("b_3prime_involved", False),
        "terminal_3prime_run": risk_metrics.get("terminal_3prime_run", 0),
        "longest_complementary_run": max([alignment["longest_run"] for alignment in alignments] or [0]),
        "paired_bases_last_5": risk_metrics.get("paired_last5", 0),
        "paired_bases_last_8": risk_metrics.get("paired_last8", 0),
        "total_paired_bases": max([alignment["total_pairs"] for alignment in alignments] or [0]),
        "dimer_tm_c": tm,
        "delta_g_kcal_mol": thermo.get("dg_kcal_mol"),
        "annealing_temp_c": conditions.annealing_temp_c,
        "tm_margin_c": tm_margin,
        "predicted_structure": thermo.get("predicted_structure", ""),
        "alignment": full_alignment_preview(primer_a.sequence, primer_b.sequence, display_metrics),
        "risk_reasons": "; ".join(reasons),
        "recommendation": recommendation_for_risk(risk),
        "primer3_error": thermo.get("error", ""),
    }


def analyze_primers(primers, conditions=None, thresholds=None):
    conditions = conditions or PrimerDimerConditions()
    thresholds = thresholds or PrimerDimerThresholds()
    validated, input_warnings = validate_primers(primers)
    results = [
        analyze_pair(primer_a, primer_b, interaction_type, conditions, thresholds)
        for primer_a, primer_b, interaction_type in primer_combinations(validated)
    ]
    results.sort(key=lambda row: (
        -RISK_ORDER.get(row["risk"], 0),
        not row["extendable_3prime"],
        -(row["dimer_tm_c"] if row["dimer_tm_c"] is not None else -999),
        row["delta_g_kcal_mol"] if row["delta_g_kcal_mol"] is not None else 999,
        -row["longest_complementary_run"],
        row["primer_a"],
        row["primer_b"],
    ))
    return {
        "conditions": asdict(conditions),
        "thresholds": asdict(thresholds),
        "primer3_version": primer3_version(),
        "input_warnings": input_warnings,
        "results": results,
        "summary_by_primer": summarize_by_primer(results),
    }


def summarize_by_primer(results):
    summary = {}
    for result in results:
        for primer_name in (result["primer_a"], result["primer_b"]):
            current = summary.get(primer_name)
            if current is None or RISK_ORDER[result["risk"]] > RISK_ORDER[current["risk"]]:
                summary[primer_name] = result
    return summary


RESULT_FIELDS = [
    "primer_a", "primer_b", "interaction_type", "risk", "extendable_3prime",
    "primer_a_3prime_involved", "primer_b_3prime_involved", "terminal_3prime_run",
    "longest_complementary_run", "paired_bases_last_5", "paired_bases_last_8",
    "total_paired_bases", "dimer_tm_c", "delta_g_kcal_mol", "annealing_temp_c",
    "tm_margin_c", "predicted_structure", "alignment", "risk_reasons",
    "recommendation", "primer3_error",
]


def write_results_csv(path, results):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        for row in results:
            writer.writerow({field: row.get(field, "") for field in RESULT_FIELDS})


def markdown_report(analysis):
    lines = [
        "# Primer Dimer Analysis",
        "",
        f"- Python/Primer3: primer3-py {analysis['primer3_version']}",
        "- Delta G from primer3-py is converted from cal/mol to kcal/mol.",
        "- Dimer predictions are computational risk estimates, not experimental validation.",
        "",
        "## Parameters",
        "",
    ]
    for key, value in analysis["conditions"].items():
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Thresholds", ""])
    for key, value in analysis["thresholds"].items():
        lines.append(f"- {key}: {value}")
    if analysis["input_warnings"]:
        lines.extend(["", "## Input Warnings", ""])
        lines.extend(f"- {warning}" for warning in analysis["input_warnings"])
    lines.extend(["", "## Results", "", "| primer_a | primer_b | type | risk | Tm C | dG kcal/mol | reasons |", "| --- | --- | --- | --- | --- | --- | --- |"])
    for result in analysis["results"]:
        lines.append(
            f"| {result['primer_a']} | {result['primer_b']} | {result['interaction_type']} | "
            f"{result['risk']} | {result['dimer_tm_c']} | {result['delta_g_kcal_mol']} | {result['risk_reasons']} |"
        )
    lines.extend(["", "## Summary By Primer", ""])
    for primer_name, result in analysis["summary_by_primer"].items():
        lines.append(f"- {primer_name}: worst interaction {result['risk']} with {result['primer_a']} / {result['primer_b']}.")
    return "\n".join(lines) + "\n"


def write_markdown_report(path, analysis):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(markdown_report(analysis))
