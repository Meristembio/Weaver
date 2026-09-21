# PCR and Primers

Weaver uses visible primer records and plasmid sequences to predict PCR products, design primers around a selected region, find amplicons on a map, search all visible plasmids, and assess primer-dimer risk.

## Primer Data Model

Create or import primers as described in [Inventory](inventory.md). The `Sequence (3' end)` field is the 3-prime hybridizing portion and is the sequence used to find template binding sites. The optional `Sequence (5' end / overhang)` field is appended to the full primer but is not required to match the template. This distinction matters for cloning primers: a 5-prime tail may remain unaligned while the 3-prime end binds.

Weaver can infer a Type IIS cloning overhang from a primer sequence when the explicit 5-prime field is empty and a recognized site, spacer, and known overhang are present. An explicit 5-prime field takes precedence. FASTA import can also populate these fields from header metadata.

## PCR Prediction for One Plasmid

Open a sequence-backed plasmid and choose its PCR action, or go directly to `/inventory/plasmid/pcr/<plasmid_id>`. Select a visible forward primer and reverse primer. If an inventory primer is not appropriate, leave the corresponding selector empty and enter its sequence in `Primer F sequence` or `Primer R sequence`; both primers must be supplied.

Submit the form. Weaver searches the plasmid sequence for the forward primer's 3-prime sequence and the reverse complement of the reverse primer's 3-prime sequence. It reports the product size, amplicon sequence, primer details, GC and approximate Tm-related values when available. It includes the 5-prime portions in the reported product size, while the binding search is based on the hybridizing portions. A missing template hit produces `FWD primer does not hit template` or `REV primer does not hit template`; an unavailable inventory primer and an omitted primer are reported as form errors.

## Selected-Region Suggestions

In the sequence viewer, select a target interval and open the PCR design action. Weaver opens the plasmid PCR page with the selected zero-based start and end coordinates, displays one-based coordinates to the user, and searches visible primers in a default 300 bp margin around the region. The request can also carry a non-negative margin and a maximum Tm difference; malformed coordinates produce `Bad PCR design coordinates`.

Review the ranked suggestions for product size, target coverage, coordinates, extra bases, Tm difference, and warnings. A suggestion is a candidate for review, not a replacement for experimental validation. Confirm that the primer 3-prime ends bind the intended sequence and that any 5-prime extension is intentional before ordering or using it.

## Amplicon Find in the Sequence Viewer

The viewer's primer-pair panel calls `/inventory/api/plasmid/<plasmid_id>/amplicon-matches/` for a readable plasmid. The default minimum product size is 100 bp; optionally set a maximum product size, one or more primer IDs, and required regions. A required region can be entered as one-based start and end coordinates or filled from the active map selection. Circular products and regions that cross the origin are supported. The default maximum primer Tm difference is 5 degrees C, and the default result set selects non-overlapping amplicons for clearer map display; disable non-overlapping selection when all candidates are needed.

Each candidate can include forward and reverse primer names and IDs, product and template sizes, circular coordinates, product sequence, recommended annealing temperature, Tm difference, primer3 dimer risk, a primer complementarity preview, and a warning when only the 3-prime portion binds and 5-prime bases remain unaligned. Select a candidate to add its Weaver-generated annotation to the map, clear those generated annotations, copy the amplicon sequence, or open the sequence in the local BLAST service.

If the plasmid has no readable sequence, the API returns an error instead of candidates. Bad numeric, JSON, or coordinate filter parameters return `Bad amplicon filter parameters`.

## Global PCR Search

Open `Services > PCR` at `/inventory/services/pcr/`. Select a forward and reverse primer and optionally enter minimum and maximum product sizes. Submit `Search plasmids` to scan plasmids you can read. The results table shows plasmid and project, product size, coordinates, forward and reverse hit counts, Tm difference, and primer3 dimer status. Circular products are labeled as such. Plasmids without a readable sequence are counted as skipped, and a search with no matching products reports that no visible plasmid generated an amplicon.

The dimer summary distinguishes `HIGH`, `MODERATE`, and lower-risk results and may include dimer Tm and Delta G. Primer complementarity is shown separately with the longest contiguous match, 3-prime match, both-3-prime match, warnings, and an alignment preview. Treat a high risk or a large Tm difference as a design warning requiring review, not as an automatic assay failure.

## Interpreting Tm and Dimer Data

The recommended annealing temperature is calculated from primer and product Tm values by the current PCR helper. The Tm difference is the absolute mismatch between the two primer Tm estimates. A small difference generally makes a pair easier to optimize, while a large difference flags a pair for review; the search defaults to a maximum difference of 5 degrees C. Primer3 dimer output evaluates heterodimer risk, including dimer Tm and Delta G when calculation succeeds. Complementarity at either 3-prime end, especially both 3-prime ends, is more concerning because it can support primer-dimer extension.

## Primer Dimer Command

The command-line tool `python manage.py analyze_primer_dimers ...` is separate from the web pages. It accepts FASTA, CSV, TSV, or auto-detected table input through the command's options and writes CSV or Markdown reports. Use the command help in the deployed checkout to confirm its required arguments. It does not create or update Inventory primer records.
