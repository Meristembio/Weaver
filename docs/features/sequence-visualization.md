# Sequence Visualization

Weaver opens sequence-backed plasmids in the bundled Open Vector Editor (OVE) integration. GenBank files provide sequence plus feature annotations; FASTA provides sequence letters without GenBank feature metadata.

## Open a Sequence

Open `Inventory > Plasmids`, select a plasmid, and use its `View / Edit` sequence action, or open `/inventory/plasmid/view_edit/<plasmid_id>`. The plasmid must have a sequence file that Weaver can parse. A GenBank import keeps the uploaded file as the stored source, so its annotated features are available immediately. From the plasmid detail page you can also download the original file, a GenBank representation, or a FASTA representation.

If the record has no sequence, a user with write or admin access can use the sequence creation/edit action to attach or generate one. A reader can inspect a sequence but cannot create a missing sequence or edit a plasmid without project write access. Parsing or sequence-size failures are displayed on the plasmid detail page.

## Map and Sequence Work

Use OVE to inspect the circular map, sequence coordinates, and GenBank-derived features. Feature names, types, strands, and qualifiers come from the uploaded record; FASTA records do not contain equivalent annotations. Use the viewer's editing controls only after confirming the intended project and record, then save through Weaver's sequence workflow when the deployment exposes a save action. The plasmid detail page remains the source for the stored sequence and metadata.

The viewer-side panels can:

- list visible primer matches and run amplicon find with size, primer, and required-region filters;
- display selected amplicons on the map, clear Weaver-generated annotations, copy a product sequence, and hand it to local BLAST;
- request restriction digest candidates using enzyme availability, buffer, temperature, fragment, and required-region constraints.

These panels require a readable sequence and use the signed-in user's project visibility. Their detailed inputs and result interpretation are documented in [PCR and Primers](pcr-and-primers.md) and [Restriction Analysis](restriction-analysis.md).

## Sequence Downloads and Public Views

The plasmid detail page offers the stored original, GenBank, and FASTA download choices at `/inventory/plasmid/download/<plasmid_id>`. A separate public route, `/inventory/public/plasmid/<plasmid_id>/`, is available for records whose public visibility is enabled. Public visibility does not change the edit permissions of the project record.

## Relationship to Verification

Use the stored sequence as the reference for [Sanger Verification](sanger-verification.md), PCR prediction, and digest planning. Before interpreting an alignment or cut site, confirm that the reference is the intended construct and that any circular origin, annotated feature, or sequence edit is represented in the saved file.
