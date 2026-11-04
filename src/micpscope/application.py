"""Application package identity and compatibility, independent of analysis code."""
import re

SAFE = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}')


def application(spec):
    schema = spec.get('schema_version')
    if schema not in {'1.0', '2.0'}:
        raise ValueError('Unsupported application specification schema; expected 1.0 or 2.0')
    app = spec.get('application')
    if app is None and schema == '1.0':
        app = dict(id=spec.get('name'), name=spec.get('name'), version=spec.get('version'),
                   description=spec.get('description', 'Legacy reference package'))
    if not isinstance(app, dict) or set(app) != {'id', 'name', 'version', 'description'}:
        raise ValueError('application requires id, name, version and description')
    for key in ('id', 'version'):
        if not isinstance(app[key], str) or not SAFE.fullmatch(app[key]) or '..' in app[key]:
            raise ValueError('Application ' + key + ' must be a safe identifier (1–80 characters)')
    for key, maximum in (('name', 160), ('description', 4000)):
        if not isinstance(app[key], str) or not app[key].strip() or len(app[key]) > maximum or any(ord(c)<32 for c in app[key]):
            raise ValueError('Application ' + key + ' must be nonempty plain text')
    if schema == '2.0' and any(key in spec for key in ('name', 'version', 'description')):
        raise ValueError('Put identity fields inside application, not at the top level')
    return dict(app)


def archive_name(spec):
    app = application(spec)
    return app['id'] + '-' + app['version'] + '.zip'


def catalog(spec):
    """Optional application-supplied annotation catalog; no implicit MICP labels."""
    annotations = spec.get('annotations', {})
    supplied = annotations.get('knowledge_catalog') if isinstance(annotations, dict) else None
    if supplied:
        return validate_catalog(supplied)
    modules = spec.get('modules', [])
    categories = [dict(id=m['id'], name=m.get('name', m['id']), description=m.get('description', m.get('interpretation', 'Reference-defined functional group')))
                  for m in modules if isinstance(m, dict) and m.get('id') != 'unassigned'] if isinstance(modules, list) else []
    categories.append(dict(id='unassigned', name='Unassigned proteins', description='No application category assigned; not an absence claim'))
    return dict(catalog_version=application(spec)['version'], scope=application(spec)['description'],
                categories=categories, entries=[], sources=[], coverage_gaps=[])


def validate_catalog(value):
    """Reject malformed catalogs and executable link schemes before packaging."""
    if not isinstance(value, dict) or any(not isinstance(value.get(k, []), list) for k in ('categories','entries','sources')):
        raise ValueError('Application catalog must contain category, entry and source lists')
    categories = value.get('categories', [])
    ids = [c.get('id') for c in categories if isinstance(c, dict)]
    if len(ids) != len(categories) or any(not isinstance(x,str) for x in ids) or len(set(ids)) != len(ids) or 'unassigned' not in ids:
        raise ValueError('Catalog categories require unique IDs including unassigned')
    for category in categories:
        if any(not isinstance(category.get(k), str) for k in ('id','name','description')):
            raise ValueError('Catalog categories require id, name and description')
    for source in value.get('sources', []):
        if not isinstance(source, dict) or any(not isinstance(source.get(k),str) for k in ('id','url','title','type')) or not source['url'].startswith(('https://','http://')):
            raise ValueError('Catalog sources require IDs, titles, types and HTTP(S) URLs')
    for entry in value.get('entries', []):
        if not isinstance(entry, dict) or entry.get('category') not in ids:
            raise ValueError('Catalog entry must select a declared category')
        if any(not isinstance(entry.get(k),str) for k in ('id','protein_name','function','evidence_level','scope')) or any(not isinstance(entry.get(k),list) or any(not isinstance(x,str) for x in entry[k]) for k in ('gene_symbols','product_aliases','references')):
            raise ValueError('Malformed catalog entry')
    return value
