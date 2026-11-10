# GeneScope 0.1.0a18

Fix the retention adapter's controller version check for native Conda installations.
The `snakemake-minimal` 9.6.2 package can report distribution metadata version
`0.0.0` while the executing engine correctly reports `snakemake.__version__` as
`9.6.2`. The adapter now validates that runtime identity. Missing distribution
metadata and the known `0.0.0` placeholder do not override a supported runtime;
unknown runtimes and conflicting concrete versions still stop execution.

Generated-data preservation from a17 remains unchanged. No search algorithms,
thresholds or scientific interpretation changed. CI now exercises native Conda
`snakemake-minimal=9.6.2` without replacing it with a pip Snakemake wheel, including
retention through failure/retry, reruns, report refresh and real-DIAMOND builds.

# GeneScope 0.1.0a17

Generated-data preservation is a core current-branch feature, applied to both
main and development. Search algorithms, thresholds and interpretation are unchanged.

- Preserve prior declared products, logs and benchmark records before Snakemake
  removes them on a rerun, retry or report refresh; retain incomplete outputs.
- Retain descriptor/report attempts and native DIAMOND query files, raw hit tables,
  commands, versions, standard output/error and temporary files. Use the supported
  `--no-unlink` flag and persistent temporary directories in DIAMOND 2.1.11.
- Preserve native Database Builder outputs and failed build attempts. Publishing
  replacements keeps the previous file or directory, including native sidecars.
- Refuse replacement when preservation fails or output ownership is unsafe.
  The adapter supports pinned Snakemake 9.6.2 and 9.27.1; it refuses other versions.

Current outputs remain at their existing paths. Previous workflow products are
under `output_dir/retained/` (or `build_dir/retained/` for databases); durable
sample-stage attempts are under `<sample>/micpscope/retained/`. Export the complete
output tree to include history. No existing lost files can be recovered by this
update. Cleanup requires a separate explicit user decision; completion is never
permission to delete generated evidence.

DIAMOND retention flags were checked against the pinned upstream source:
https://github.com/bbuchfink/diamond/blob/v2.1.11/src/basic/config.cpp

# GeneScope 0.1.0a16

This alpha release applies the pre-beta evidence-integrity fixes to the current
installable channel. The same corrections are retained in development for
parallel testing. It is not a declaration of analytical validation or a beta
promotion.

- Align all 164 MICP starter reference records with the annotation category vocabulary and reject unresolved category identifiers during database build/load. The corrected starter is **MICP 1.0.1**, catalog `2026.11-category.2`; reference sequences and search thresholds are unchanged.
- Preserve supplied labels separately from reference candidates. Conflicting labels/categories, ambiguous candidates and unknown categories remain unresolved, with visible counts at the top of the report and consistent JSON/TSV/HTML states. All reference categories are retained; the first declared category is the documented primary display group for an unambiguous, nonconflicting candidate.
- Check GFF features against declared region bounds and reject contradictory declarations. Missing bounds and identifier-only linkage are explicit warnings. Circular-origin features display as two intervals while counting as one source feature; multipart neighborhoods remain unresolved.
- Track renderer modules and entry scripts so report-only code changes regenerate HTML without rerunning protein analysis. Preserve explicit unassessed operon/promoter and engineering-performance states.
- Run regression, browser, real DIAMOND smoke and Snakemake dependency checks on both main and development. Browser CI uses the Playwright container's own Python to retain matching SSL libraries.

**Updating existing work:** install the new compatible source through Bactologix.
Rebuild an old MICP 1.0.0 database from the corrected bundled starter before new
searches, select the resulting 1.0.1 package for a new draft, and retain the old
run and database. Saved reports can be refreshed to show clearer unresolved
states, but refreshing does not rerun reference searches or replace the saved
reference database. No existing results should be deleted or overwritten.

See [integrity details](prebeta-integrity.md) for the exact interpretation
boundaries and regression scope. Software checks do not establish sensitivity,
specificity, functional activity, material performance or deployment safety.
