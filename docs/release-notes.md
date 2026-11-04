## GeneScope 0.1.0a11

- Clearly identify descriptor-only runs where no application database search was requested.
- Add a compact coverage panel and in-app next steps to select a completed application database.
- Distinguish missing searches from negative findings; retain collapsed protein records and scientific interpretation limits.

## GeneScope 0.1.0a10

- Publish from the renamed canonical repository, crackberryq/Bactomics_GeneScope.
- Update package repository/issue metadata and migration documentation.
- Preserve historical release tags, workflow aliases, Python imports and output paths. Analysis behavior and application database schemas are unchanged.

## GeneScope 0.1.0a9

- Separate the reusable GeneScope engine from the bundled MICP application package.
- Add application schema 2.0 with id, display name, version and description; derive package filenames, report titles and database provenance from this metadata.
- Read interpretation and category definitions from the selected package. Generic runs no longer receive implicit MICP annotations.
- Preserve schema 1.0 databases, historical paths and MICPscope.smk through compatibility interfaces; add GeneScope.smk and genescope CLI.
- Validate custom applications, named archives, legacy builds and escaped report labels.

## MICP Scope 0.1.0a8

- Search a completed reference database with a dedicated Snakemake DIAMOND rule and pinned environment.
- Verify database checksums before and after search; record database identity, search settings, tool version and query hash.
- Display candidate matches, ambiguity, alignment coverage and linked review evidence in collapsed report sections; retain supplied annotations separately.
- Export reference scores in TSV. Distinguish no qualifying hits from unsupported inputs.
- Test the full bundled-database search and report with real DIAMOND in CI.

## MICP Scope 0.1.0a7

- Bundle the default literature starter with 164 reviewed protein sequences and 207 table entries from the supplied processes and regulation reviews.
- Default DatabaseBuilder.smk to this portable specification when reference_spec is omitted.
- Preserve citations, relative regulatory positions, uncertainties and unresolved selections; no absolute promoter coordinates are invented.
- Verify an offline starter build with real DIAMOND in CI.

# MICP Scope 0.1.0a6

Adds sample-free DatabaseBuilder.smk for user-specified references, local protein and regulatory sequence imports, pinned accession downloads, DIAMOND indexes, packaging and checksum verification. Includes dedicated Conda environment, stage benchmarks and synthetic build tests. See docs/database-builder.md. This release builds reference packages; the protein report does not yet search them.
