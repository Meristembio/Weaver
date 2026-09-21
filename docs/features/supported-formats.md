# Supported Formats

Use the format that matches the workflow. File extensions are validated by individual forms, so a format accepted by one page may not be accepted by another.

| Format | Used by | Important behavior |
| --- | --- | --- |
| GenBank `.gb`, `.gbk`, `.genbank` | Plasmid upload and batch import | Parsed as GenBank; feature annotations are retained when present. Batch import stores the original file on each plasmid. |
| FASTA `.fa`, `.fas`, `.fasta`, `.txt` | Primer import and FASTA alignment | FASTA alignment accepts text or files; primer import reads FASTA text and can derive metadata from headers. |
| AB1 `.ab1` | Sanger verification | Supplies sequence, quality values, peak positions, trace channels, and ABIF metadata when present. |
| PHD.1 `.phd.1` | Sanger verification | Supplies sequence and quality information when present, but not AB1 trace channels. |
| SEQ `.seq` | Sanger verification | Text sequence only; no chromatogram trace or Phred quality values. |
| Clustal `.clustal` | Saved alignment evidence | Clustal uploads are validated with a maximum size of 1 MB. Sanger and FASTA alignment can save Clustal output. |
| CSV and TSV | Primer dimer command | The command expects name and sequence columns or a compatible primer-pair table. |

## Plasmid Sequence Uploads

The manual plasmid form accepts a sequence file through the plasmid `sequence` field. Use GenBank when the viewer must show annotated features and FASTA when only the sequence is needed. A file that cannot be parsed is reported on the plasmid page or import log; Weaver does not infer missing feature annotations from FASTA.

## Sanger Upload Limits and Mixing Rules

The current Sanger processing limits are 96 files per batch, 8 MB per file, and 96 MB total. Upload at least one `.ab1`, `.phd.1`, or `.seq` file for Sanger mode. FASTA files belong to the separate alignment mode and cannot be mixed with Sanger trace files in a single batch. SEQ files can be aligned but cannot open the chromatogram viewer because they contain no trace data.

## Primer FASTA Metadata

Primer import accepts `.fa`, `.fas`, `.fasta`, and `.txt`. The import form lets the user choose whether the name comes from the FASTA ID or full description, select a default direction when the name does not indicate forward or reverse, and update existing primers with the same name. Headers can provide direction, overhang, and intended-use values. Invalid DNA/IUPAC characters or an empty sequence are skipped with a message.

## Command-Line Input

The `analyze_primer_dimers` management command supports FASTA, CSV, TSV, and auto-detected table input according to its command options. This command is distinct from web uploads and writes analysis reports; it does not change inventory records.
