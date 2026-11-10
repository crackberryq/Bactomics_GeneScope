"""Keep each descriptor attempt, including partial results on failure."""
from pathlib import Path
import sys

root, sample, proteins, data, table, *extra = sys.argv[1:]
evidence = extra[:1] if extra and extra[0] != "--reference-results" else []
reference_args = extra[len(evidence):]
sys.path.insert(0, str(Path(root) / "src"))
from micpscope.cli import main
from micpscope.retention import retained_attempt, publish

parent = Path(data).parent
parent.mkdir(parents=True, exist_ok=True)
with retained_attempt(parent, "properties") as staging:
    output = Path(staging) / "analysis"
    args = ["--proteins", proteins, "--sample", sample, "--out", str(output), "--data-only"]
    if evidence:
        args += ["--evidence", evidence[0]]
    main(args + reference_args)
    publish(output / "results.json", data)
    publish(output / "proteins.tsv", table)

