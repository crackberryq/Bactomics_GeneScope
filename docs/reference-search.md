# Reference protein search

Select a completed Database Builder `database/` folder as `reference_database` in the analysis configuration. Bactologix pins its checksum manifest when creating a draft. An empty setting retains descriptor-only analysis.

The native `reference_search` rule runs DIAMOND blastp 2.1.11 in `envs/reference-search.yml`, using very-sensitive mode and all target hits with one HSP per target. Defaults require at least 30% identity, 70% aligned coverage of both query and reference, and E-value at most 1e-5. These configurable discovery filters are not validated functional cutoffs. Near-best hits within 95% of the best bit score with differing gene/product labels are marked ambiguous. All qualifying hits remain available.

One terminal stop is removed. Unsupported residues cause an explicit skipped status, rather than a no-hit result. No qualifying hit does not establish biological absence. Supplied names and functions are retained separately from reference similarity evidence. Similarity does not establish expression, mineralization, complete operons or operating optima.

Outputs include `reference-search.json`, raw `reference-search.hits.tsv`, enriched `results.json`, `proteins.tsv`, and the collapsible HTML report. Database name/version, checksum and manifest hashes, DIAMOND version, command, filters and query hash are recorded. Database contents are checked before and after search. Linked review table rows and database bibliography remain inspectable as reference evidence rather than observations measured on the sample.
