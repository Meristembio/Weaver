# Sanger Verification

Sanger verification compares one or more sequencing reads with a readable reference plasmid sequence and saves the parsed files, chromatograms, alignment evidence, metrics, and review decision as a verification run. Parsing and alignment are upload-time operations; saved-run listings and pagination read the persisted result and do not reopen or reanalyze AB1 files.

## Upload a Run

Open a plasmid and choose `Align > Sanger`, or use `/inventory/plasmid/align/sanger/<plasmid_id>`. The plasmid must have a readable reference sequence. Upload one or more `.ab1`, `.phd.1`, or `.seq` files. The form also accepts FASTA extensions for the separate FASTA alignment mode, but FASTA files cannot be mixed with Sanger trace files in one batch. At least one supported file is required.

Optionally enter a run label and notes, and enable saving a Clustal file. The default processing limit is 96 files, 8 MB per file, and 96 MB for the batch. Invalid extensions, an empty upload, mixed alignment modes, oversized files, or unreadable content are reported as validation or processing errors.

For batch uploads, the mapping CSV requires `ab1_file` and `plasmid_id`. The `primer_id` column is optional; an empty value or an omitted column stores the read without a primer assignment.

Batch uploads are processed and persisted one plasmid at a time inside a single transaction, so a later processing error rolls back the complete batch while keeping peak Python memory bounded by the current plasmid's result. The persistence boundary accepts either a list of uploaded-file records or a single record from the processing service.

Weaver groups files by base name when a read has more than one representation. For a group it prefers AB1, then PHD.1, then SEQ as the sequence source while retaining the uploaded files. SEQ supplies sequence text only; it cannot supply chromatogram traces or Phred quality values.

The saved-run `Reads` table shows the primer associated with the selected source file for each read. Users with write access to the plasmid can assign a primer, replace an incorrect assignment, or clear it when the metadata is missing.

## Read Orientation and Alignment

For each read, Weaver records the selected source, raw and trimmed sequence, parsing errors, warnings, quality metrics, and detected orientation. Orientation can be `Forward`, `Reverse complement`, `Ambiguous`, or `Unmapped`. The aligned read is compared with the plasmid reference, and the run aggregates coverage, alignment metrics, variants, gaps, and low-confidence regions. A read can be retained but marked unusable when its sequence or quality data do not support alignment.

The saved-run detail view is scoped to the selected run. When a run contains repeated reads with the same orientation and highly overlapping plasmid coverage, the detail view keeps the newest read using the sequencer date encoded in its filename, falling back to the date stored in the file. Same-orientation reads covering distinct regions remain visible. When a run contains complete Forward/Reverse evidence groups from different sequencing dates, the detail view uses the newest complete date group; same-orientation reads covering distinct regions within that group remain eligible. If usable evidence exists, failed or unmapped reads from the same upload are kept in a separate troubleshooting section with access to their chromatogram; they are excluded from the active sequencing evidence and metrics and remain hidden until their troubleshooting entry is selected. If no usable evidence exists, they remain visible as the run's review evidence. This is a display and review selection; it does not delete the persisted read records. The run selector displays the sequencing date encoded in the filename, falling back to AB1 metadata, and exposes failed reads as selectable `No reliable alignment` troubleshooting entries without exposing their filenames; same-date failed reads with consecutive sequencer IDs are grouped into one entry.

The Services table groups persisted runs by plasmid; it does not merge persisted runs. Within the selected run, it shows the newest complete Forward/Reverse sequencing-date group and omits unknown-orientation reads whenever reliable Forward or Reverse evidence is available, so older or failed evidence is not mixed into the active row. If a run has no reliable orientation, its unknown reads remain available in the summary for troubleshooting. The Date column shows the sequencing date (without time) encoded in the selected filename, falling back to the date stored in file metadata; the upload date (also without time) is available from the adjacent information icon. It shows the newest manually accepted run (`Verified / Approved`) by default when one exists. If none is manually accepted, it shows the newest automatically approved run (`PASS`), then the newest run needing review (`REVIEW`); only if neither exists does it fall back to the newest uploaded run. The run selector exposes every run for that plasmid so a different run can be reviewed explicitly.

Each inline AB1 chromatogram has its own `Auto-adjust` switch. Autoscaling is calculated only from that read's visible trace window, and changing one read's switch does not change the scale of any other inline chromatogram. Switching modes uses the same smooth scale transition as the detailed chromatogram viewer. The inline trace can also be dragged horizontally to move the visible alignment window.

When the alignment contains a deletion, the inline trace does not stretch the signal across the missing reference bases. It returns to the baseline through the deletion interval and resumes at the next aligned base; inline traces are drawn as lines without an artificial filled area.

Open a saved run at `/inventory/plasmid/align/sanger/<plasmid_id>/run/<run_id>`. The alignment browser provides a reference coordinate view, read coverage and orientation, feature context when the reference is annotated, variants, gaps, and quality regions. Use its read controls to move through the reference and inspect the evidence supporting a possible variant. Coordinates for stored variants are zero-based internally; use the displayed reference coordinate labels when communicating a finding.

The Sanger services header exposes batch upload as an upload icon with an `Upload AB1 batch` hover tooltip; counts and page metadata remain available through the table pagination controls rather than the header.

The one-base left and right controls advance once on click and repeat continuously while held, with keyboard Enter/Space support and automatic stop on release or cancellation.

## Quality, Coverage, Variants, and Gaps

Coverage shows which reference positions are supported by usable reads. A gap is a reference interval without sufficient aligned read support. A variant is an observed difference from the reference and may carry a quality value and a low-quality flag. Low-confidence and intermediate-confidence regions identify positions where the base call deserves manual review rather than automatic acceptance. Compare a variant with all covering reads and the chromatogram before deciding whether it is biological, a sequencing artifact, or unresolved.

## Chromatograms

For an AB1-backed read, choose its chromatogram action at `/inventory/plasmid/align/sanger/<plasmid_id>/run/<run_id>/read/<read_id>/chromatogram`. The viewer plots A, C, G, and T traces with base calls and quality coloring. Use `Standard` or `Autoscaled`, the left and right controls, the zoom slider, and the base number plus `Go` to inspect a region. Low and intermediate confidence regions are shaded. PHD.1 and SEQ reads do not provide the AB1 trace arrays needed for this viewer.

## Manual Decision and Validation State

The automated run classification can be `PASS`, `REVIEW`, `FAIL`, or `NO_DATA`. It is evidence for review, not the final laboratory decision. On the saved run, choose a manual decision of `Verified`, `Not verified`, `Inconclusive`, or `Pending`; enter an effective date and comment when appropriate, then save.

When a run is marked `Verified`, Weaver sets the associated plasmid's sequencing validation state to verified and stores the manual effective date as the plasmid sequencing date. Other decisions do not mark the plasmid as verified and clear the sequencing date when the saved decision is not verified. The run retains the reviewer, review timestamp, decision, effective date, comment, automated reasons, and any saved Clustal file.

## FASTA Alignment

For sequence text without trace data, choose `Align > Fasta` at `/inventory/plasmid/align/fasta/<plasmid_id>`. Paste FASTA text or upload one or more `.fa`, `.fas`, `.fasta`, or `.txt` files. Choose `Together` or `One at a time` in `View mode`, and optionally save a Clustal file. This mode provides sequence alignment evidence but no chromatogram or Phred interpretation.
