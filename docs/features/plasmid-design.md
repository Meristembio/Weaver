# Plasmid Design

Weaver supports manual plasmid records, a backbone-and-insert assembly wizard, and a Level 0 design form. All three create or update plasmid records in a project; they do not replace experimental validation.

## Choose a Workflow

- Use `Plasmids > Create > Form` at `/inventory/plasmid/create/` for a record whose metadata or sequence you are entering directly.
- Use `Plasmids > Create > Wizard` at `/inventory/plasmid/create/wizard` when a construct should be assembled from a backbone and inserts.
- Use `Plasmids > Create > L0 designer` or `Services > L0 Designer` when starting with a Level 0 sequence and assembly-standard overhangs.

Creating or editing a plasmid requires write or admin access to its project. Confirm that the destination project and any selected backbone or insert belong to records you are allowed to use; the form and server-side permission checks are authoritative for the deployment.

## Manual Plasmid Form

Enter a name and intended use. Add selectable markers, an optional sequence file, backbone, inserts, plasmid type, assembly level, description, destination project, ligation state, and creation date as appropriate. Save the form to create the record, then open the detail page to check the stored metadata, computed sequence size, and available sequence actions. Use GenBank when feature annotations must be retained; use FASTA when only sequence letters are required.

## Assembly Wizard

Open the wizard and follow its form controls to choose a backbone and compatible inserts, assembly standard, and construct metadata. The current frontend passes the selected data to Weaver's plasmid creation endpoint and then routes to the generated plasmid record. Review the resulting backbone, inserts, sequence, selectable markers, assembly level, and ligation state before using the construct in PCR, digest, or Sanger workflows.

Assembly standards are defined in the application and include Yeast Toolkit, Loop, and Golden Braid 2.0. The wizard's available choices and compatibility rules are the source of truth for a particular deployment. A missing component, incompatible overhang, or invalid construct is reported by the form or creation endpoint and must be corrected before a record can be completed.

## Level 0 Design

The Level 0 form accepts a sequence, a ligation standard, an L0 5-prime overhang, an L0 3-prime overhang, and the selected Type IIS enzyme. The supported standard-specific choices include the overhang sets defined for Yeast Toolkit, Loop, and Golden Braid 2.0. Enter sequence letters without relying on the viewer to repair them, select compatible end overhangs, and submit. Weaver returns the generated output or an error message; inspect the resulting plasmid sequence and metadata before continuing.

## After Creation

The detail page is the handoff point for downstream work. Attach or verify a readable sequence, inspect it in [Sequence Visualization](sequence-visualization.md), use [PCR and Primers](pcr-and-primers.md) or [Restriction Analysis](restriction-analysis.md) for in-silico checks, and record experimental results through the plasmid validation form and [Sanger Verification](sanger-verification.md). Derivative relationships can be used to rebuild dependent sequences from a changed parent when the user has write access.

## Current Limits

The application has a local `moclo` dependency, but the active assembly workflows use Weaver's own standard definitions and helpers. The web interface is therefore the authoritative source for supported standards and available overhangs in the current deployment.
