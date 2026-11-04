# GeneScope application packages

GeneScope is the analysis engine. The selected database defines the application. The bundled MICP starter is one application, not a different pipeline.

```yaml
schema_version: "2.0"
application:
  id: micp
  name: MICP
  version: "1.0.0"
  description: Microbially induced calcium carbonate precipitation
proteins: [] # Supply at least one documented reference record
regulatory: []
modules: []
annotations: {}
references: []
```

The display name and description accept plain text; id and version use filesystem-safe identifiers. A build creates `micp-1.0.0.zip`, and its completion receipt records all four application fields. GeneScope reports use **GeneScope — MICP** and record application and software versions independently. Package names are never inferred from a local directory name.

Supply sequences through the existing local FASTA or pinned accession fields. Optional modules define functional categories; protein `categories` link references to them. An optional `annotations.knowledge_catalog` supplies the conservative header-annotation catalog (categories, entries, sources and evidence labels). The MICP starter includes this catalog. No arbitrary Python or rule code is loaded from packages. Gene groups and regulatory observations remain descriptive metadata until their corresponding analysis capabilities are implemented.

Schema 1.0 databases remain readable using their existing name and version. Unknown schemas and conflicting identity fields are rejected. Existing frozen runs remain unchanged. The repository is now [crackberryq/Bactomics_GeneScope](https://github.com/crackberryq/Bactomics_GeneScope). The Python package, persisted Bactologix pipeline ID and output folder retain their historical `micpscope` identifiers for compatibility. `GeneScope.smk` and the `genescope` CLI are the new public entry points; `MICPscope.smk` remains an alias and the old CLI is an explicit legacy annotation interface.

The builder validates files, metadata and integrity, not the biological truth of supplied claims. Similarity evidence does not establish activity, gene expression or operating optima.
