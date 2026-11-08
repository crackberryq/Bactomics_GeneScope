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
