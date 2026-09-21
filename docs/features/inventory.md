# Inventory

The Inventory module stores the laboratory records used by the rest of Weaver: projects, plasmids, glycerol stocks, primers, restriction enzymes, strains, boxes, locations, and labels. Most inventory pages are available from the top navigation after sign-in.

## Sign In, Projects, and Permissions

Open `/accounts/login/` and submit the Django username and password. A successful login returns you to the requested page or to the profile page. Inventory pages require authentication.

Open `/organization/projects/` to see projects of which you are a member. Create a project at `/organization/project/create`. Project creation adds the creator as a member; the application assigns admin access when an existing membership is found and otherwise creates a write membership. A project has a name, optional description, public flag, and optional assembly standard.

Open a project at `/organization/project/<project_id>/`. The project page shows membership, record counts, visibility, and assembly standard. Select `Set current` to store the active project in the browser. Many create forms use this current project as their default. The navigation also provides a `show from all projects` toggle for list pages; this changes which records are displayed, but it does not grant access to projects where you are not a member.

Membership access has three levels:

- `read` allows a member to view records and use read-only analysis for that project.
- `write` allows record creation and edits in addition to reading.
- `admin` allows project administration, including membership changes and project deletion.

Only a member can read a project. Plasmid and glycerol-stock writes require write or admin access to the record's project. Primers are global inventory records; their writes require write or admin access in at least one project but are not tied to the active project. Project administration requires admin access, and the last project admin cannot be removed. A user who needs to work in a project must therefore be added to its memberships rather than relying on the active-project selector.

## Plasmids

Open `Inventory > Plasmids` at `/inventory/plasmids/`. The list can be scoped to the current project or to all projects visible to the signed-in user. Each record includes a name, numeric index when assigned, selectable markers, intended use, description, project, creation date, assembly type and level, sequence file, QR identifier, construction state, and validation fields.

The plasmid list is ordered with the most recent `created_on` date first and is
loaded with server-side pagination and search. Each response is limited to a
combined budget of 50 displayed rows: up to 20 recently viewed plasmids are
shown above the main list, and the remaining slots are filled by the main
list. Opening a plasmid records it in the browser's recent history; entries
remain subject to the selected project scope.

To create a record, use `Plasmids > Create > Form` or `/inventory/plasmid/create/`. Enter the name and intended use, then add any selectable markers, sequence file, backbone, inserts, type, assembly level, description, destination project, ligation state, and creation date offered by the form. Save only after confirming that the destination project is one where you have write or admin access. The resulting detail page links to the sequence viewer, downloads, PCR, digest, and FASTA or Sanger alignment tools when a sequence is present.

The list also exposes `Plasmids > Validation`, public plasmid views for records marked public, QR lookup, label rendering, and record editing or deletion subject to project permissions. A record without a readable sequence cannot be used for PCR, digest, or sequence alignment; the relevant page reports that the sequence could not be read.

## Batch GenBank Import

Use `Plasmids > Import GenBank` at `/inventory/plasmid/import/`. You need write or admin access to the destination project and one or more `.gb`, `.gbk`, or `.genbank` files. Choose one of these destination modes:

1. `Upload to an existing project`: select a writable project. The current project is selected by default when it is writable.
2. `Create a new project and upload there`: provide a project name, choose an assembly standard, and optionally mark the new project public. The new project is created with the importing user as a member and becomes the current project after a successful import.

Choose `Filename stem` or `GenBank LOCUS / record name` as the name source. Enable `Update existing plasmids with the same name` only when an import should replace the matching record in the destination project. Enable `Import as reference sequences` for records intended to be references, `Mark imported plasmids as public` when the public plasmid view is appropriate, and `Infer assembly type, level, and resistance` when YTK-like Type IIS annotations should be used for conservative metadata inference.

Each uploaded file becomes one plasmid. Weaver stores the original GenBank file as the plasmid sequence file, parses the sequence and annotated features, and reports created, updated, skipped, and error counts followed by an import log. Existing records are skipped unless updating was selected. Invalid or unreadable GenBank files appear as import errors rather than silently becoming empty records. Imported features remain available to the sequence viewer because the original annotated file is retained.

## Glycerol Stocks, Strains, Boxes, and Locations

Open `Stocks` at `/inventory/glycerolstocks/`. A glycerol stock records a strain, optional plasmid, optional parent stock, creation date, box, row, column, project, QR ID, and free-text details. To create one, use `Stocks > Create`, select the strain and plasmid when applicable, choose a box position, and save it to the appropriate project. The box-position controls use the configured row and column choices; the box itself links to a location.

Use `Stocks > Boxes` to inspect boxes and their locations, or `Stocks > From QR` at `/inventory/qr/g/` to look up a stock by QR identifier. A stock detail page links back to its plasmid when one exists and provides a label view. The stock list can be filtered to the current project or expanded to all projects visible to the user.

Strain, box, and location records do not have dedicated create links in the inventory navigation. When they are not already present, an administrator must create or maintain them through Django admin before a stock can reference them.

## Primers

Open `Components > Primers` at `/inventory/primers/`. A primer record contains a name, 3-prime hybridizing sequence, optional 5-prime extension or overhang, forward/reverse direction, intended use, and QR ID. Primers are global and the list does not use the project selector. The two sequence fields are written 5-prime to 3-prime: put only the region that must hybridize to the template in `Sequence (3' end)`, and put non-hybridizing added bases in `Sequence (5' end / overhang)`.

Create a primer at `/inventory/primer/create/` or import a FASTA file at `/inventory/primer/import/`. The importer accepts `.fa`, `.fas`, `.fasta`, or `.txt` text and offers the FASTA ID or full description as the name source, a default direction for names without `F` or `R`, and an option to update an existing primer with the same name globally. FASTA metadata can declare an overhang and intended use. Invalid DNA or an empty hybridizing sequence is skipped and reported; existing primers are skipped unless updating is enabled.

## Restriction Enzymes and Labels

Open `Components > Restriction Enzymes` at `/inventory/restrictionenzymes/`. Create an enzyme at `/inventory/restrictionenzyme/create/` by selecting a Biopython-known enzyme, optionally selecting its HF version, and adding a datasheet URL or lab description. The form prevents duplicate name and HF-version records. Add one or more buffer rows with an existing buffer or `+ New buffer...`, and enter an activity percentage from 0 to 100. The same buffer cannot be added twice to one enzyme. These records supply the enzyme and buffer data used by [Restriction Analysis](restriction-analysis.md).

For a single plasmid or stock, use its label action. For multiple labels, open `Services > Batch Prints` at `/inventory/services/batch-prints/`, choose `Plasmids` or `Glycerol stocks` per row, enter a name, QR ID, UUID, or supported numeric plasmid index, and provide the row-specific label details. Missing identifiers are reported in the result rather than silently printed.

## Validation Fields

The plasmid validation page records ligation state, working colony, an explicit `No colony` status for constructs used directly or obtained by synthesis, colony-PCR state/date/observations, digestion state/date/observations, sequencing state/date/observations, and an optional Clustal file. A blank colony with `No colony` unchecked remains `Not set`; it is not treated as a synthetic/direct-use construct. Use the validation form after experimental checks; it is separate from the automated PCR, digest, and Sanger analysis pages. Validation links accept either a numeric colony or `none`/`no-colony` in the colony position. See [Sanger Verification](sanger-verification.md) for how a saved sequencing decision updates the sequencing validation fields.
