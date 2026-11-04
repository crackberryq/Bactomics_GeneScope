import argparse
from pathlib import Path
import sys

parser = argparse.ArgumentParser()
parser.add_argument('root')
parser.add_argument('proteins')
parser.add_argument('database')
parser.add_argument('output')
parser.add_argument('--threads', type=int, default=1)
parser.add_argument('--expected', default='')
parser.add_argument('--identity', type=float, default=30)
parser.add_argument('--coverage', type=float, default=70)
parser.add_argument('--evalue', type=float, default=1e-5)
args = parser.parse_args()
sys.path.insert(0, str(Path(args.root) / 'src'))
from micpscope.reference import search
search(args.proteins, args.database, args.output, args.threads, args.expected,
       args.identity, args.coverage, args.evalue)
