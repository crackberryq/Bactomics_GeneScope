"""Conservative interpretation of supplied annotation, not sequence annotation."""
import hashlib
import json
import re
from pathlib import Path


def load_catalog():
    raw = (Path(__file__).parent / 'data' / 'micp_catalog.json').read_bytes()
    catalog = json.loads(raw)
    return catalog, hashlib.sha256(raw).hexdigest()


def annotation_context(identifier, description, catalog):
    # FASTA identifier is retained separately. Never substring-match locus tags.
    product = description.split(maxsplit=1)[1] if len(description.split(maxsplit=1)) == 2 else ''
    explicit = re.findall(r'(?:\bgene=|\bGN=)([A-Za-z0-9_-]+)', description, re.I)
    tokens = set(re.findall(r'(?<![\w-])[A-Za-z][A-Za-z0-9]*(?![\w-])', product))
    tokens.add(identifier)
    symbols = {x.casefold() for x in (explicit or tokens)}
    candidates = []
    for entry in catalog['entries']:
        gene_match = any(x.casefold() in symbols for x in entry['gene_symbols'])
        phrase_match = not explicit and any(
            re.search(r'(?<!\w)' + re.escape(x) + r'(?!\w)', product, re.I)
            for x in entry['product_aliases'])
        if gene_match or phrase_match:
            candidates.append(entry)
    result = dict(gene_name=', '.join(explicit) or 'unassigned',
                  protein_name=product or 'unannotated protein',
                  function=product or 'Function not supplied', category='unassigned',
                  annotation_status='supplied annotation only' if product else 'unannotated',
                  catalog_ids=[], references=[], evidence_level='unassessed', literature_scope='Not assessed',
                  annotation_basis='Input FASTA header; not independently validated')
    if len(candidates) == 1:
        entry = candidates[0]
        result.update(gene_name=entry['gene_symbols'][0] if entry['gene_symbols'] else 'family only',
                      protein_name=entry['protein_name'], function=entry['function'],
                      category=entry['category'], annotation_status='catalog-linked annotation',
                      catalog_ids=[entry['id']], references=entry['references'], evidence_level=entry['evidence_level'],
                      literature_scope=entry['scope'], annotation_basis='Header label matched to literature catalog; not sequence-confirmed')
    elif candidates:
        result.update(annotation_status='ambiguous annotation',
                      catalog_ids=[x['id'] for x in candidates],
                      literature_scope='Multiple catalog labels matched; no category assigned')
    return result


def interpret_descriptors(row):
    result = dict(temperature_indication='Undetermined: no validated aliphatic-index-to-temperature conversion',
                  ph_indication='Undetermined: theoretical pI is not optimal activity pH',
                  instability_class='not calculated', hydropathy_class='not calculated')
    if row['status'] == 'calculated':
        index = row['instability_index']
        result['instability_class'] = ('below 40: stable-class heuristic' if index < 40 else
                                       'above 40: unstable-class heuristic' if index > 40 else
                                       '40: boundary; not classified')
        gravy = row['gravy']
        result['hydropathy_class'] = ('negative mean hydropathy' if gravy < 0 else
                                      'positive mean hydropathy' if gravy > 0 else 'zero mean hydropathy')
    return result


def enrich(proteins, catalog):
    for row in proteins:
        row.update(annotation_context(row['protein_id'], row['description'], catalog))
        row['supplied_annotation'] = {k: row[k] for k in (
            'gene_name', 'protein_name', 'function', 'category', 'annotation_status', 'annotation_basis')}
        row['annotation_category'] = row['category']
        row['reference_search_status'] = 'not_requested'
        row['evidence_status'] = 'not_requested'
        row['category_basis'] = 'Supplied annotation only; no sequence reference search requested'
        row.update(interpret_descriptors(row))
    return proteins
