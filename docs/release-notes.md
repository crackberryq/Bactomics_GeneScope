# GeneScope 0.1.0a15

- Replace whole-proteome category meters with clickable application-category bars and named protein comparisons for aliphatic index, pI, GRAVY and instability index.
- Move the whole-proteome scatter plot into collapsed background content at the end of the report.
- Add a complete chromosome density overview and a default 20 kb gene view. Dense intervals show all features as density rather than silently truncating the first 500 genes.
- Search saved gene names, protein descriptions and identifiers with suggestions, and navigate from category comparisons to gene locations.
- Prefer gene names, then protein names, on records and gene arrows; preserve identifiers in details. Preserve supplied GFF gene names in new result bundles.
- Hide empty operon and unknown-strand tracks. Supplied operon annotations remain distinct from gene proximity.
- Existing compatible saved results can be regenerated through Bactologix; missing labels use available protein annotations. No new reference search is required for these presentation changes.
