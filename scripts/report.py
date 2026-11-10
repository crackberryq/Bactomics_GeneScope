import json
from pathlib import Path
import sys

root, source, destination = sys.argv[1:]
sys.path.insert(0, str(Path(root) / "src"))
from micpscope.cli import render

from micpscope.retention import retained_attempt, publish

with retained_attempt(Path(destination).parent, "report") as attempt:
    report = render(json.loads(Path(source).read_text()))
    generated = attempt / "report.html"
    generated.write_text(report, encoding="utf-8")
    publish(generated, destination)
print("Report generated:", destination)
