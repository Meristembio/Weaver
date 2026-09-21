from dataclasses import dataclass

from Bio.Seq import Seq


def normalize_dna(sequence):
    return "".join(base for base in str(sequence or "").upper() if base in {"A", "C", "G", "T"})


def reverse_complement(sequence):
    return str(Seq(normalize_dna(sequence)).reverse_complement())


def circular_slice(sequence, start, end):
    sequence = normalize_dna(sequence)
    if not sequence:
        return ""
    start %= len(sequence)
    end %= len(sequence)
    if start == end:
        return ""
    if start < end:
        return sequence[start:end]
    return sequence[start:] + sequence[:end]


def interval_contains(point, start, end, sequence_length):
    if sequence_length <= 0:
        return False
    point %= sequence_length
    start %= sequence_length
    end %= sequence_length
    if start < end:
        return start <= point < end
    return point >= start or point < end


@dataclass(frozen=True)
class TypeIISSiteMatch:
    enzyme_name: str
    recognition_site: str
    site_start: int
    site_end: int
    orientation: str
    plus_cut: int
    minus_cut: int
    overhang: str

    @property
    def left_edge(self):
        return min(self.plus_cut, self.minus_cut)

    @property
    def right_edge(self):
        return max(self.plus_cut, self.minus_cut)


@dataclass(frozen=True)
class TypeIISEnzymeDefinition:
    name: str
    recognition_site: str
    top_strand_cut_offset: int
    bottom_strand_cut_offset: int
    overhang_length: int
    aliases: tuple = ()
    overhang_polarity: str = "five_prime"
    source: str = ""
    version: str = ""

    def __post_init__(self):
        recognition_site = normalize_dna(self.recognition_site)
        aliases = tuple(sorted({normalize_dna(alias) or str(alias).upper() for alias in self.aliases if alias}))
        object.__setattr__(self, "recognition_site", recognition_site)
        object.__setattr__(self, "aliases", aliases)
        if not recognition_site:
            raise ValueError("Type IIS enzyme definitions require a recognition site.")
        if self.overhang_length <= 0:
            raise ValueError("Type IIS enzyme definitions require a positive overhang length.")
        if abs(self.bottom_strand_cut_offset - self.top_strand_cut_offset) != self.overhang_length:
            raise ValueError(
                f"{self.name}: cut offsets must differ by exactly the overhang length."
            )

    @property
    def reverse_complement_site(self):
        return reverse_complement(self.recognition_site)

    def matches_name(self, enzyme_name):
        normalized = str(enzyme_name or "").strip().upper()
        return normalized == self.name.upper() or normalized in self.aliases

    def find_sites(self, sequence, circular=False):
        sequence = normalize_dna(sequence)
        if not sequence:
            return []

        recognition_site = self.recognition_site
        reverse_site = self.reverse_complement_site
        search_sequence = sequence if not circular else sequence + sequence[:len(recognition_site) - 1]
        matches = []

        max_start = len(sequence) if circular else max(len(sequence) - len(recognition_site) + 1, 0)
        for start in range(max_start):
            window = search_sequence[start:start + len(recognition_site)]
            if len(window) != len(recognition_site):
                continue

            if window == recognition_site:
                plus_cut = start + self.top_strand_cut_offset
                minus_cut = start + self.bottom_strand_cut_offset
                overhang = circular_slice(sequence, plus_cut, minus_cut) if circular else sequence[plus_cut:minus_cut]
                if len(overhang) == self.overhang_length:
                    matches.append(TypeIISSiteMatch(
                        enzyme_name=self.name,
                        recognition_site=recognition_site,
                        site_start=start % len(sequence),
                        site_end=(start + len(recognition_site)) % len(sequence),
                        orientation="forward",
                        plus_cut=plus_cut % len(sequence) if circular else plus_cut,
                        minus_cut=minus_cut % len(sequence) if circular else minus_cut,
                        overhang=overhang,
                    ))

            if reverse_site != recognition_site and window == reverse_site:
                plus_cut = start + len(recognition_site) - self.bottom_strand_cut_offset
                minus_cut = start + len(recognition_site) - self.top_strand_cut_offset
                overhang = circular_slice(sequence, plus_cut, minus_cut) if circular else sequence[plus_cut:minus_cut]
                if len(overhang) == self.overhang_length:
                    matches.append(TypeIISSiteMatch(
                        enzyme_name=self.name,
                        recognition_site=reverse_site,
                        site_start=start % len(sequence),
                        site_end=(start + len(recognition_site)) % len(sequence),
                        orientation="reverse",
                        plus_cut=plus_cut % len(sequence) if circular else plus_cut,
                        minus_cut=minus_cut % len(sequence) if circular else minus_cut,
                        overhang=overhang,
                    ))

        return sorted(matches, key=lambda match: (match.left_edge, match.right_edge, match.orientation))
