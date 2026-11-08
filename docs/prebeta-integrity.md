# Pre-beta evidence integrity fixes

This release carries the development fixes addressing the November 2026 review. Passing these
regressions does not qualify analytical sensitivity, specificity or biological
performance.

The bundled MICP database is now **1.0.1**, with category vocabulary
`2026.11-category.2`. Its knowledge catalog uses the same declared module IDs as
the reference records. Database build and load reject unknown or duplicate
category identifiers. Existing 1.0.0 packages with inconsistent category metadata
must be rebuilt from the corrected specification; their checksums must not be
changed in place.

Supplied annotations and reference similarity evidence are recorded separately.
Ambiguous candidates and disagreements between annotation and reference category
or gene labels remain unresolved, and are excluded from resolved category
summaries. A label disagreement requests review; it does not prove which label
is correct. No-hit, skipped, not-requested and unresolved states are distinct in
JSON, TSV, collapsed protein headings and visual summaries. Refreshing older
saved reports recomputes display grouping without modifying the saved evidence.

All qualifying reference categories are retained in `reference_categories`.
For an unambiguous candidate without annotation disagreement, the first category
declared by the reference is the primary display group. This ordering is explicit
package metadata, not a biological priority or a completeness inference.
Near-tied references with different category sets remain unresolved even if
their gene labels agree. This preserves one protein record per report while
retaining secondary category associations in exports and record evidence.

GFF features must agree with declared sequence-region bounds. Contradictory
declarations reject. A contig without a region declaration remains usable with
an explicit unverified bounds warning, including partially declared files.
A circular region marked `Is_circular=true` may have a feature
crossing the origin; original coordinates are retained and the viewer displays
both intervals. Multipart CDS records remain separate. Circular/multipart
neighborhood inference remains unassessed. Imported coordinates are linked only
by identifiers: no sequence pairing verification or operon inference is claimed.

The workflow tracks rendering modules and entry scripts explicitly. Changing
only the explorer, report assets or report script regenerates HTML without
recomputing protein data. Reference-search code dependencies include the FASTA
parser. Development pushes run the same CI checks without publishing releases.

Validation includes all 164 starter references with neutral synthetic query
identifiers, conflicting and near-tied evidence, JSON/TSV/HTML agreement, invalid
and circular coordinates, browser behavior, and an actual Snakemake rendering
dependency check. These fixtures test computational contracts, not empirical
functional assignment accuracy.

On hosts where Snakemake's benchmark process inspection raises
`psutil.NoSuchProcess`, `GENESCOPE_TEST_NO_BENCHMARK=1` disables benchmarking only
in the dependency test's temporary workflow copy. It does not modify production
rules or validate production resource benchmarking. CI uses the normal rules.
