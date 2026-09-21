# Restriction Analysis

Weaver uses lab restriction-enzyme records, Biopython recognition and cut-site metadata, and a plasmid's readable sequence to review cuts and plan digest combinations.

## Single-Plasmid Digest

Open a sequence-backed plasmid and choose its digest action, or use `/inventory/plasmid/digest/<plasmid_id>`. The page lists the restriction enzymes loaded in Weaver and their cut information. Select one or more enzyme names and submit `Digest`. Weaver treats the plasmid as circular, calculates cut positions and predicted fragments, and displays the selected enzymes, fragment sizes, and available buffer activity data. A plasmid without a readable sequence produces an error instead of a fragment table.

This page is a direct cut-site review. It uses the selected enzyme names and does not choose an optimal pair for you.

## OVE Digest Planning

The Open Vector Editor restriction panel calls `/inventory/api/plasmid/<plasmid_id>/restriction-digests/` for the current plasmid. Before planning, load the lab's enzymes under `Components > Restriction Enzymes` and add their buffer activity percentages. Enzymes absent from the database are not available to the planner.

The planner can evaluate one or two enzymes and accepts constraints for minimum and maximum fragment count, minimum fragment size, minimum band-size difference, required cut regions, required enzyme names, and a result limit. Required regions come from the active map selection or coordinate entries. The planner requires 100% shared activity in a buffer and a compatible known optimum temperature. It returns exact matches when all constraints are met and closest-match results when an exact plan is not available, so read the result classification before treating a combination as suitable.

The API validates numeric and JSON parameters, caps the enzyme count at two and the result limit at 50, rejects an unknown required enzyme, and rejects a maximum fragment count lower than the minimum. A bad request returns an error explaining which restriction-digest parameter needs correction. Circular coordinates may cross the origin; confirm the reported cut positions against the plasmid sequence and the laboratory protocol.

## Interpreting the Result

Use the fragment count and sizes to determine whether the digest distinguishes the expected construct. Use the required-region match to confirm that the cut sites cover the selected feature or interval. Use the shared-buffer and temperature fields to check that the proposed reaction conditions are compatible. The planner's approximation is a ranking aid, not a guarantee of enzyme performance, star activity, or complete digestion.
