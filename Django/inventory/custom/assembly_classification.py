from dataclasses import dataclass
from dataclasses import field

from inventory.custom.type_iis import TypeIISEnzymeDefinition
from inventory.custom.type_iis import circular_slice
from inventory.custom.type_iis import interval_contains
from inventory.custom.type_iis import normalize_dna
from inventory.custom.type_iis import reverse_complement


CLASSIFIER_VERSION = "2026.08-phase1"


@dataclass(frozen=True)
class PartTypeDefinition:
    key: str
    name: str
    functional_category: str
    enzyme: str
    upstream_overhang: str
    downstream_overhang: str
    assembly_level: str
    plasmid_role: str = "insert"
    aliases: tuple = ()
    expected_feature_types: tuple = ()
    expected_feature_keywords: tuple = ()

    def __post_init__(self):
        object.__setattr__(self, "upstream_overhang", normalize_dna(self.upstream_overhang))
        object.__setattr__(self, "downstream_overhang", normalize_dna(self.downstream_overhang))
        object.__setattr__(
            self,
            "aliases",
            tuple(sorted({str(alias or "").strip().lower() for alias in self.aliases if str(alias or "").strip()})),
        )
        object.__setattr__(
            self,
            "expected_feature_types",
            tuple(sorted({str(value or "").strip().lower() for value in self.expected_feature_types if value})),
        )
        object.__setattr__(
            self,
            "expected_feature_keywords",
            tuple(sorted({str(value or "").strip().lower() for value in self.expected_feature_keywords if value})),
        )
        if not self.upstream_overhang or not self.downstream_overhang:
            raise ValueError(f"{self.key}: part definitions require upstream and downstream overhangs.")


@dataclass(frozen=True)
class AssemblyStandardDefinition:
    id: str
    name: str
    version: str
    enzymes: dict
    part_types: tuple
    circular_sequences_supported: bool = True


@dataclass(frozen=True)
class CandidateRegion:
    enzyme_name: str
    start: int
    end: int
    wraps_origin: bool
    upstream_overhang: str
    downstream_overhang: str
    sequence: str
    internal_site_count: int


@dataclass
class ClassificationResult:
    standard_id: str = ""
    standard_name: str = ""
    standard_version: str = ""
    classifier_version: str = CLASSIFIER_VERSION
    assembly_level: str = ""
    plasmid_role: str = ""
    part_type_key: str = ""
    part_name: str = ""
    functional_category: str = ""
    enzyme_name: str = ""
    upstream_overhang: str = ""
    downstream_overhang: str = ""
    start: int = 0
    end: int = 0
    wraps_origin: bool = False
    orientation: str = "forward"
    confidence: float = 0.0
    evidence: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    alternative_candidates: list = field(default_factory=list)
    source: str = "digest"

    @property
    def confidence_band(self):
        if self.confidence >= 0.90:
            return "HIGH"
        if self.confidence >= 0.70:
            return "MEDIUM"
        if self.confidence >= 0.50:
            return "LOW"
        return "UNKNOWN"

    @property
    def model_type_name(self):
        if self.plasmid_role == "receiver":
            return "Receiver"
        if self.plasmid_role == "insert":
            return "Insert"
        return None

    @property
    def model_level(self):
        if self.assembly_level in {"ENTRY_VECTOR", "LEVEL_0"}:
            return 0
        if self.assembly_level == "LEVEL_1":
            return 1
        if self.assembly_level == "LEVEL_2":
            return 2
        return None

    def as_dict(self):
        return {
            "standard_id": self.standard_id,
            "standard_name": self.standard_name,
            "standard_version": self.standard_version,
            "classifier_version": self.classifier_version,
            "assembly_level": self.assembly_level,
            "plasmid_role": self.plasmid_role,
            "part_type_key": self.part_type_key,
            "part_name": self.part_name,
            "functional_category": self.functional_category,
            "enzyme_name": self.enzyme_name,
            "upstream_overhang": self.upstream_overhang,
            "downstream_overhang": self.downstream_overhang,
            "start": self.start,
            "end": self.end,
            "wraps_origin": self.wraps_origin,
            "orientation": self.orientation,
            "confidence": self.confidence,
            "confidence_band": self.confidence_band,
            "evidence": list(self.evidence),
            "warnings": list(self.warnings),
            "alternative_candidates": list(self.alternative_candidates),
            "source": self.source,
            "model_type_name": self.model_type_name,
            "model_level": self.model_level,
        }


YTK_ENZYMES = {
    "BsaI": TypeIISEnzymeDefinition(
        name="BsaI",
        aliases=("BSAI",),
        recognition_site="GGTCTC",
        top_strand_cut_offset=7,
        bottom_strand_cut_offset=11,
        overhang_length=4,
        source="Lee et al. 2015 / Biopython Restriction",
        version="1.0",
    ),
    "BsmBI": TypeIISEnzymeDefinition(
        name="BsmBI",
        aliases=("BSMBI", "ESP3I"),
        recognition_site="CGTCTC",
        top_strand_cut_offset=7,
        bottom_strand_cut_offset=11,
        overhang_length=4,
        source="Lee et al. 2015 / Biopython Restriction",
        version="1.0",
    ),
}


YTK_PART_TYPES = (
    PartTypeDefinition(
        key="ytk_1",
        name="YTK Part 1",
        functional_category="connector",
        enzyme="BsaI",
        upstream_overhang="CCCT",
        downstream_overhang="AACG",
        assembly_level="LEVEL_0",
        aliases=("connector", "assembly connector"),
        expected_feature_keywords=("connector", "con", "assembly"),
    ),
    PartTypeDefinition(
        key="ytk_2",
        name="YTK Part 2",
        functional_category="promoter",
        enzyme="BsaI",
        upstream_overhang="AACG",
        downstream_overhang="TATG",
        assembly_level="LEVEL_0",
        expected_feature_types=("promoter", "regulatory"),
        expected_feature_keywords=("promoter",),
    ),
    PartTypeDefinition(
        key="ytk_3",
        name="YTK Part 3",
        functional_category="coding_sequence",
        enzyme="BsaI",
        upstream_overhang="TATG",
        downstream_overhang="ATCC",
        assembly_level="LEVEL_0",
        expected_feature_types=("cds", "gene"),
        expected_feature_keywords=("cds", "gene"),
    ),
    PartTypeDefinition(
        key="ytk_3a",
        name="YTK Part 3a",
        functional_category="coding_sequence",
        enzyme="BsaI",
        upstream_overhang="TATG",
        downstream_overhang="TTCT",
        assembly_level="LEVEL_0",
        expected_feature_types=("cds", "gene"),
        expected_feature_keywords=("cds", "gene", "fusion", "tag"),
    ),
    PartTypeDefinition(
        key="ytk_3b",
        name="YTK Part 3b",
        functional_category="coding_sequence",
        enzyme="BsaI",
        upstream_overhang="TTCT",
        downstream_overhang="ATCC",
        assembly_level="LEVEL_0",
        expected_feature_types=("cds", "gene"),
        expected_feature_keywords=("cds", "gene", "fusion", "tag"),
    ),
    PartTypeDefinition(
        key="ytk_4",
        name="YTK Part 4",
        functional_category="terminator",
        enzyme="BsaI",
        upstream_overhang="ATCC",
        downstream_overhang="GCTG",
        assembly_level="LEVEL_0",
        expected_feature_types=("terminator", "regulatory"),
        expected_feature_keywords=("terminator",),
    ),
    PartTypeDefinition(
        key="ytk_4a",
        name="YTK Part 4a",
        functional_category="tag",
        enzyme="BsaI",
        upstream_overhang="ATCC",
        downstream_overhang="TGGC",
        assembly_level="LEVEL_0",
        expected_feature_types=("cds", "misc_feature"),
        expected_feature_keywords=("tag", "fusion"),
    ),
    PartTypeDefinition(
        key="ytk_4b",
        name="YTK Part 4b",
        functional_category="terminator",
        enzyme="BsaI",
        upstream_overhang="TGGC",
        downstream_overhang="GCTG",
        assembly_level="LEVEL_0",
        expected_feature_types=("terminator", "regulatory"),
        expected_feature_keywords=("terminator",),
    ),
    PartTypeDefinition(
        key="ytk_5",
        name="YTK Part 5",
        functional_category="connector",
        enzyme="BsaI",
        upstream_overhang="GCTG",
        downstream_overhang="TACA",
        assembly_level="LEVEL_0",
        aliases=("connector", "assembly connector"),
        expected_feature_keywords=("connector", "con", "assembly"),
    ),
    PartTypeDefinition(
        key="ytk_6",
        name="YTK Part 6",
        functional_category="marker",
        enzyme="BsaI",
        upstream_overhang="TACA",
        downstream_overhang="GAGT",
        assembly_level="LEVEL_0",
        expected_feature_types=("cds", "gene"),
        expected_feature_keywords=("marker", "resistance", "camr", "ampr", "hygr", "kanr", "specr", "zeo"),
    ),
    PartTypeDefinition(
        key="ytk_7",
        name="YTK Part 7",
        functional_category="origin_or_homology",
        enzyme="BsaI",
        upstream_overhang="GAGT",
        downstream_overhang="CCGA",
        assembly_level="LEVEL_0",
        expected_feature_types=("rep_origin", "misc_feature"),
        expected_feature_keywords=("origin", "homology", "ars"),
    ),
    PartTypeDefinition(
        key="ytk_8",
        name="YTK Part 8",
        functional_category="bacterial_backbone",
        enzyme="BsaI",
        upstream_overhang="CCGA",
        downstream_overhang="CCCT",
        assembly_level="LEVEL_0",
        expected_feature_types=("rep_origin", "cds", "gene"),
        expected_feature_keywords=("origin", "resistance", "cole1", "camr", "ampr", "kanr", "zeo"),
    ),
    PartTypeDefinition(
        key="ytk_8a",
        name="YTK Part 8a",
        functional_category="bacterial_backbone",
        enzyme="BsaI",
        upstream_overhang="CCGA",
        downstream_overhang="CAAT",
        assembly_level="LEVEL_0",
        expected_feature_types=("rep_origin", "cds", "gene"),
        expected_feature_keywords=("origin", "resistance", "cole1", "camr", "ampr", "kanr", "zeo"),
    ),
    PartTypeDefinition(
        key="ytk_8b",
        name="YTK Part 8b",
        functional_category="homology",
        enzyme="BsaI",
        upstream_overhang="CAAT",
        downstream_overhang="CCCT",
        assembly_level="LEVEL_0",
        expected_feature_types=("misc_feature",),
        expected_feature_keywords=("homology",),
    ),
    PartTypeDefinition(
        key="ytk_234",
        name="YTK Part 234",
        functional_category="composite",
        enzyme="BsaI",
        upstream_overhang="AACG",
        downstream_overhang="GCTG",
        assembly_level="LEVEL_0",
        expected_feature_keywords=("promoter", "terminator", "cds"),
    ),
    PartTypeDefinition(
        key="ytk_234r",
        name="YTK Part 234r",
        functional_category="receiver_dropout",
        enzyme="BsaI",
        upstream_overhang="GCTG",
        downstream_overhang="AACG",
        assembly_level="LEVEL_1",
        plasmid_role="receiver",
        aliases=("dropout", "receiver", "gfp dropout", "rfp dropout"),
        expected_feature_types=("rep_origin", "cds", "gene", "misc_feature"),
        expected_feature_keywords=("origin", "resistance", "camr", "ampr", "kanr", "zeo", "spec", "dropout", "receiver"),
    ),
    PartTypeDefinition(
        key="ytk_678",
        name="YTK Part 678",
        functional_category="composite",
        enzyme="BsaI",
        upstream_overhang="TACA",
        downstream_overhang="CCCT",
        assembly_level="LEVEL_0",
        expected_feature_keywords=("marker", "origin"),
    ),
)


SUPPORTED_ASSEMBLY_STANDARDS = {
    "ytk": AssemblyStandardDefinition(
        id="ytk",
        name="Yeast Toolkit",
        version="1.0",
        enzymes=YTK_ENZYMES,
        part_types=YTK_PART_TYPES,
        circular_sequences_supported=True,
    )
}


def record_topology_is_circular(record):
    return str(record.annotations.get("topology") or "").strip().lower() == "circular"


def feature_metadata(record):
    metadata = []
    for feature in getattr(record, "features", []):
        labels = [str(label or "").strip() for label in feature.qualifiers.get("label", []) if str(label or "").strip()]
        label_text = " ".join(labels).lower()
        metadata.append({
            "type": str(getattr(feature, "type", "") or "").strip().lower(),
            "labels": labels,
            "label_text": label_text,
            "start": int(feature.location.start),
            "end": int(feature.location.end),
        })
    return metadata


def features_in_region(record, start, end, wraps_origin):
    sequence_length = len(record.seq)
    matches = []
    for feature in feature_metadata(record):
        probe_points = {feature["start"], max(feature["end"] - 1, feature["start"])}
        if wraps_origin:
            if any(interval_contains(point, start, end, sequence_length) for point in probe_points):
                matches.append(feature)
        elif feature["end"] > start and feature["start"] < end:
            matches.append(feature)
    return matches


def score_feature_evidence(part_definition, region_features):
    if not region_features:
        return 0.0, [], ["No supporting GenBank features found inside the detected Type IIS fragment."]

    evidence = []
    warnings = []
    compatible = False
    contradictory = False
    expected_types = set(part_definition.expected_feature_types)
    expected_keywords = set(part_definition.expected_feature_keywords)

    for feature in region_features:
        feature_type = feature["type"]
        label_text = feature["label_text"]
        if feature_type in expected_types or any(keyword in label_text for keyword in expected_keywords):
            compatible = True
        if expected_types:
            contradiction_map = {
                "promoter": {"terminator", "cds", "gene"},
                "coding_sequence": {"terminator", "promoter"},
                "terminator": {"promoter", "cds", "gene"},
            }
            for expected_type in expected_types:
                if feature_type in contradiction_map.get(expected_type, set()):
                    contradictory = True

    if compatible:
        evidence.append("Compatible feature annotation supports the overhang-based classification.")
        return 0.10, evidence, warnings

    if contradictory:
        warnings.append("Feature annotations contradict the overhang-derived YTK part type.")
        return -0.25, evidence, warnings

    warnings.append("GenBank features are present but do not support the detected YTK part type.")
    return -0.05, evidence, warnings


def score_alias_evidence(part_definition, region_features):
    labels = " ".join(
        label.lower()
        for feature in region_features
        for label in feature["labels"]
    )
    if part_definition.functional_category == "connector" and "con" in labels:
        return 0.05, ["Connector-like feature labels support this YTK connector call."]
    if any(alias in labels for alias in part_definition.aliases):
        return 0.05, ["Feature labels are compatible with this part definition."]
    return 0.0, []


def count_internal_sites(candidate_start, candidate_end, wraps_origin, sites, sequence_length):
    internal = 0
    for site in sites:
        if interval_contains(site.left_edge, candidate_start, candidate_end, sequence_length) and interval_contains(
                max(site.right_edge - 1, site.left_edge), candidate_start, candidate_end, sequence_length):
            internal += 1
    return internal


def candidate_regions_for_enzyme(sequence, enzyme_definition, circular=True):
    sequence = normalize_dna(sequence)
    sites = enzyme_definition.find_sites(sequence, circular=circular)
    sequence_length = len(sequence)
    candidates = []

    for first_site in sites:
        if first_site.orientation != "forward":
            continue
        for second_site in sites:
            if second_site.orientation != "reverse":
                continue
            candidate_start = first_site.right_edge
            candidate_end = second_site.left_edge
            wraps_origin = candidate_end < candidate_start
            fragment = circular_slice(sequence, candidate_start, candidate_end) if circular else sequence[candidate_start:candidate_end]
            if not fragment:
                continue
            internal_site_count = count_internal_sites(candidate_start, candidate_end, wraps_origin, [
                site for site in sites if site != first_site and site != second_site
            ], sequence_length)
            candidates.append(CandidateRegion(
                enzyme_name=enzyme_definition.name,
                start=candidate_start % sequence_length if circular else candidate_start,
                end=candidate_end % sequence_length if circular else candidate_end,
                wraps_origin=wraps_origin,
                upstream_overhang=first_site.overhang,
                downstream_overhang=second_site.overhang,
                sequence=fragment,
                internal_site_count=internal_site_count,
            ))

    return candidates


def standard_definition_for_id(standard_id):
    return SUPPORTED_ASSEMBLY_STANDARDS.get(str(standard_id or "").strip().lower())


def match_part_definition(candidate, part_definition):
    if candidate.upstream_overhang == part_definition.upstream_overhang and candidate.downstream_overhang == part_definition.downstream_overhang:
        return "forward"

    if reverse_complement(candidate.downstream_overhang) == part_definition.upstream_overhang and reverse_complement(
            candidate.upstream_overhang) == part_definition.downstream_overhang:
        return "reverse_complement"

    return ""


def digest_based_candidates(record, standard_definition):
    sequence = normalize_dna(record.seq)
    circular = record_topology_is_circular(record) and standard_definition.circular_sequences_supported
    features = feature_metadata(record)
    candidates = []

    for part_definition in standard_definition.part_types:
        enzyme_definition = standard_definition.enzymes.get(part_definition.enzyme)
        if not enzyme_definition:
            continue
        for candidate in candidate_regions_for_enzyme(sequence, enzyme_definition, circular=circular):
            orientation = match_part_definition(candidate, part_definition)
            if not orientation:
                continue

            region_features = features_in_region(record, candidate.start, candidate.end, candidate.wraps_origin)
            feature_score, feature_evidence, feature_warnings = score_feature_evidence(part_definition, region_features)
            alias_score, alias_evidence = score_alias_evidence(part_definition, region_features)
            evidence = [
                f"Overhang pair {candidate.upstream_overhang} -> {candidate.downstream_overhang} matches {part_definition.name}.",
                f"Detected from valid {candidate.enzyme_name} Type IIS sites flanking the candidate fragment.",
            ]
            warnings = list(feature_warnings)
            confidence = 0.50 + 0.15 + feature_score + alias_score

            if candidate.internal_site_count:
                confidence -= 0.40
                warnings.append("Additional internal Type IIS sites of the same enzyme make this digest ambiguous.")

            candidates.append(ClassificationResult(
                standard_id=standard_definition.id,
                standard_name=standard_definition.name,
                standard_version=standard_definition.version,
                assembly_level=part_definition.assembly_level,
                plasmid_role=part_definition.plasmid_role,
                part_type_key=part_definition.key,
                part_name=part_definition.name,
                functional_category=part_definition.functional_category,
                enzyme_name=candidate.enzyme_name,
                upstream_overhang=part_definition.upstream_overhang,
                downstream_overhang=part_definition.downstream_overhang,
                start=candidate.start,
                end=candidate.end,
                wraps_origin=candidate.wraps_origin,
                orientation=orientation,
                confidence=max(min(confidence, 0.99), 0.0),
                evidence=evidence + feature_evidence + alias_evidence,
                warnings=warnings,
                source="digest",
            ))

    return sorted(candidates, key=lambda result: result.confidence, reverse=True)


def legacy_label_fallback(record, preferred_standard_id="ytk"):
    labels = {
        str(label or "").strip().upper().replace("(1)", "")
        for feature in getattr(record, "features", [])
        for label in feature.qualifiers.get("label", [])
        if str(label or "").strip().upper().replace("(1)", "") in {"BSAI", "BSMBI", "SAPI"}
    }

    if labels == {"BSAI"}:
        return ClassificationResult(
            standard_id=preferred_standard_id,
            standard_name="Yeast Toolkit",
            standard_version="legacy",
            assembly_level="LEVEL_0",
            plasmid_role="insert",
            part_name="Legacy YTK Insert",
            enzyme_name="BsaI",
            confidence=0.45,
            evidence=["Fell back to legacy feature-label inference because no digest-based YTK match was found."],
            warnings=["Classification was inferred from GenBank labels only."],
            source="legacy_labels",
        )

    if labels == {"BSMBI"}:
        return ClassificationResult(
            standard_id=preferred_standard_id,
            standard_name="Yeast Toolkit",
            standard_version="legacy",
            assembly_level="ENTRY_VECTOR",
            plasmid_role="receiver",
            part_name="Legacy YTK Receiver",
            enzyme_name="BsmBI",
            confidence=0.45,
            evidence=["Fell back to legacy feature-label inference because no digest-based YTK match was found."],
            warnings=["Classification was inferred from GenBank labels only."],
            source="legacy_labels",
        )

    return None


class AssemblyStandardClassifier:
    def __init__(self, supported_standards=None):
        self.supported_standards = supported_standards or SUPPORTED_ASSEMBLY_STANDARDS

    def classify(self, record, standard_id=None, allow_legacy_fallback=True):
        standards = []
        if standard_id:
            standard_definition = standard_definition_for_id(standard_id)
            if standard_definition:
                standards.append(standard_definition)
        else:
            standards.extend(self.supported_standards.values())

        all_candidates = []
        for standard_definition in standards:
            all_candidates.extend(digest_based_candidates(record, standard_definition))

        all_candidates.sort(key=lambda candidate: candidate.confidence, reverse=True)
        if all_candidates:
            best = all_candidates[0]
            best.alternative_candidates = [
                {
                    "standard_id": candidate.standard_id,
                    "part_type_key": candidate.part_type_key,
                    "part_name": candidate.part_name,
                    "confidence": candidate.confidence,
                }
                for candidate in all_candidates[1:4]
            ]
            return best

        if allow_legacy_fallback:
            return legacy_label_fallback(record, preferred_standard_id=standard_id or "ytk")

        return None


def classify_assembly_record(record, standard_id=None, allow_legacy_fallback=True):
    return AssemblyStandardClassifier().classify(
        record,
        standard_id=standard_id,
        allow_legacy_fallback=allow_legacy_fallback,
    )


def assembly_metadata_from_classification(
        classification,
        confirmed_type_name=None,
        confirmed_level=None,
        confirmed_standard_id=None):
    if not classification:
        return {}

    return {
        "detected": classification.as_dict(),
        "confirmed": {
            "type_name": confirmed_type_name,
            "level": confirmed_level,
            "standard_id": confirmed_standard_id or classification.standard_id,
        },
    }
