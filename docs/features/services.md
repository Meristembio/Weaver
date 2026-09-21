# Services

The Services menu contains workflows that search across records or provide specialized tools. All service pages require login; inventory searches use projects that the signed-in user can read.

## BLAST

Open `Services > Blast` at `/inventory/services/blast/`. Enter a FASTA-formatted sequence in `Fasta Text Input (Preferred)` or upload a FASTA file, then choose a project or `All` in `Project to search in`. The optional `Use short input BLAST parameters?` flag changes the local search parameters for a short query. Submit `Do Blast`.

The service parses the query, searches readable plasmid sequences, and returns subject hits with alignment details. A malformed query produces `Input sequence not in FASTA format`; an empty request produces `No input sequence`. Plasmids without sequence files or with unreadable sequences are excluded and listed as not considered. The amplicon panel can open this service with a copied query sequence.

## Global PCR

Open `Services > PCR` at `/inventory/services/pcr/`. Choose a forward and reverse primer from visible primers, then set a minimum product size (default 100 bp) and optional maximum. The maximum must be greater than or equal to the minimum. `Search plasmids` scans readable plasmids in all projects where you are a member and reports products, coordinates, hit counts, Tm difference, dimer risk, and skipped records. See [PCR and Primers](pcr-and-primers.md) for interpretation.

## Batch Prints

Open `Services > Batch Prints` at `/inventory/services/batch-prints/`. Add one or more rows, choose `Plasmids` or `Glycerol stocks`, and enter one identifier per row. Plasmid lookup accepts an exact name, QR ID, UUID, or numeric plasmid index. Glycerol stock lookup accepts QR ID or UUID and can also use a numeric plasmid index for a stock associated with that plasmid. Enter the label date; plasmid rows also accept concentration in ng/ul, and rows can include a colony value.

Submit the rows to render the matching labels. The result keeps found labels and lists missing identifiers separately. Only records in projects readable by the signed-in user are considered.

## Experiments Map

Open `Plasmids > Experiments` at `/inventory/experiments`. Weaver lists experiments from projects you can read. Each experiment includes its description and a ReactFlow map built from its root plasmids, backbones, and inserts. Nodes expose plasmid name, index, type, level, working-colony status, dependencies, and a link to the plasmid detail page. The map statistics count total, validated, reference, pending, ready-to-build, and blocked nodes. Plasmids without an index are omitted from the map.

## Stats

Open `Services > Stats` at `/inventory/services/stats/`. Submit the `stats` refresh action to update the stored statistics, then read the plasmid count and summaries for sequence coverage, glycerol-stock coverage, plasmid type, assembly level, and monthly counts. The page is an inventory summary for visible data and is not a substitute for the individual record pages.

## L0 Designer and GTR

`Services > L0 Designer` opens `/inventory/services/l0d/`, the embedded Level 0 design tool. Its form requires sequence text, a ligation standard, L0 5-prime and L0 3-prime overhang choices, and a Type IIS enzyme. The tool returns generated output or an error; inspect the result before creating or validating a plasmid.

`/inventory/services/gtr/` opens the embedded GTR tool for Golden Gate or GoldenBraid-related work. It is a separate frontend service; the current Weaver routes expose the page but do not document its internal controls beyond the loaded tool. Use [Plasmid Design](plasmid-design.md) for the application-managed assembly workflows.
