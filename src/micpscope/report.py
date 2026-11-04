"""Standalone, offline report with category summaries and inspectable evidence."""
import html
import json
from .context import neighborhood


def esc(value):
    return html.escape(str(value if value is not None else '—'), quote=True)


def table(rows, columns):
    return '<div class="scroll"><table><thead><tr>' + ''.join(
        '<th scope="col">' + esc(c.replace('_', ' ')) + '</th>' for c in columns
    ) + '</tr></thead><tbody>' + ''.join(
        '<tr>' + ''.join('<td>' + esc(row.get(c, '—')) + '</td>' for c in columns) + '</tr>'
        for row in rows) + '</tbody></table></div>'


def protein_context(protein, context):
    parts = []
    domains = protein.get('domains', [])
    if domains:
        width = 700
        diagram = ['<svg viewBox="0 0 760 ' + str(36 + len(domains)*30) + '" role="img" aria-label="Recorded protein signature positions">']
        for i, d in enumerate(domains):
            x = 30 + (d['start']-1) / max(1, protein['length']) * width
            w = max(1, (d['end']-d['start']+1) / max(1, protein['length']) * width)
            y = 12 + i*30
            diagram.append(f'<line x1="30" x2="730" y1="{y+7}" y2="{y+7}" stroke="#d3e2e6"/>')
            diagram.append(f'<rect x="{x:.2f}" y="{y}" width="{w:.2f}" height="14" rx="3" fill="#087c78"><title>' + esc(d['accession']) + ': ' + esc(d['description']) + '</title></rect>')
        diagram.append('</svg>')
        parts.append('<details class="subpanel"><summary>Protein signatures · ' + str(len(domains)) + '</summary>' + ''.join(diagram) + '<p>Imported signatures can describe families, sites or domains; they are not all structural domains. Coordinates are 1-based and inclusive.</p>' + table(domains, ['analysis','accession','description','start','end','score','interpro']) + '</details>')
    nearby = neighborhood(protein, context)
    if nearby:
        parts.append('<details class="subpanel"><summary>Genome neighborhood · ' + esc(protein['coordinates'][0]['contig']) + '</summary><p>Focal protein: ' + esc(protein['protein_id']) + '. Up to two adjacent CDS records on either side; recorded coordinates do not establish co-transcription or a complete operon.</p>' + table(nearby, ['feature_id','protein_id','start','end','strand','product']) + '</details>')
    counts = protein.get('amino_acid_counts', {})
    if counts:
        composition = [dict(residue=a, count=n, percent=round(100*n/protein['length'],2)) for a,n in sorted(counts.items())]
        parts.append('<details class="subpanel"><summary>Amino-acid composition</summary>' + table(composition,['residue','count','percent']) + '</details>')
    copies = protein.get('identical_sequence_records', [])
    if len(copies) > 1:
        parts.append('<p><b>Identical sequence records:</b> ' + esc(', '.join(copies)) + '. Sequence equality alone does not establish separate genomic copies.</p>')
    parts.append('<details class="subpanel"><summary>Assessment coverage</summary><p>Signatures: ' + esc(protein.get('domain_status','not assessed')) + '<br>Coordinates: ' + esc(protein.get('coordinate_status','not assessed')) + '<br>Localization, disorder, structure and operon prediction: not assessed.</p></details>')
    return ''.join(parts)


def reference_details(protein, search):
    if not search:
        return ''
    hits = protein.get('reference_hits', [])
    rows, evidence = [], []
    records = {r['id']: r for r in search.get('literature_records', [])}
    for hit in hits:
        ref = search['references'][hit['reference_id']]
        rows.append(dict(hit, gene=ref.get('gene', ''), protein=ref.get('protein_name', ''),
                         organism=ref.get('organism', ''), categories=', '.join(ref.get('categories', []))))
        linked = [records[r] for r in ref.get('references', []) if r in records]
        for record in linked:
            evidence.append('<details class="subpanel"><summary>' + esc(ref['id'] + ' · ' + record['id']) +
                '</summary>' + table(record.get('fields', []), ['heading', 'text']) +
                '<p>' + esc(record.get('curation_status', '')) + '</p></details>')
    return ('<details class="subpanel"><summary>Reference similarity · ' +
            esc(protein.get('reference_search_status', 'unassessed')) + '</summary><p>' +
            esc(search['interpretation']) + '</p>' + table(rows, ['reference_id', 'gene', 'protein', 'organism',
            'categories', 'identity', 'query_coverage', 'reference_coverage', 'evalue', 'bitscore']) +
            ''.join(evidence) + '</details>')


def reference_summary(search):
    if not search:
        return '<details class="card"><summary>Reference search · not requested</summary><p>No reference database was selected for this run.</p></details>'
    matched = sum(bool(q['hits']) for q in search['queries'].values())
    return ('<details class="card" id="reference-search"><summary><strong>Reference search · ' +
            str(matched) + ' proteins with candidates</strong></summary><p>' + esc(search['interpretation']) +
            '</p>' + table([search['database']], ['name', 'version', 'description', 'checksums_sha256', 'manifest_sha256']) +
            '<p>Tool: ' + esc(search['tool_version']) + '</p>' + table([search['thresholds']],
            ['identity', 'query_coverage', 'reference_coverage', 'evalue', 'ambiguity_bitscore_fraction']) +
            '<details><summary>Database bibliography and provenance</summary><pre style="white-space:pre-wrap;overflow-wrap:anywhere">' +
            esc(json.dumps(search.get('bibliography', []), indent=2, ensure_ascii=False)) + '</pre></details></details>')


def render_report(result):
    app = result.get('application') or {}
    application_title = 'GeneScope' + (' — ' + str(app['name']) if app.get('name') else '')
    proteins = result['proteins']
    search_results = result.get('reference_search')
    catalog = result.get('knowledge_catalog', {})
    categories = catalog.get('categories', [{'id':'unassigned','name':'Protein properties','description':''}])
    sources = {s['id']:s for s in catalog.get('sources', [])}
    linked = sum(p.get('annotation_status') == 'catalog-linked annotation' for p in proteins)
    calculated = sum(p['status'] == 'calculated' for p in proteins)
    summary = []
    sections = []
    empty_sections = []
    context = result.get('context', {})
    flags = [dict(protein_id=p['protein_id'], flags=', '.join(p.get('quality_flags', [])) or p.get('notes','')) for p in proteins if p.get('quality_flags') or p['status'] != 'calculated']
    export_columns = ['protein_id','gene_name','protein_name','category','function','length','aliphatic_index','theoretical_pi','gravy','instability_index','annotation_status']
    export_data = json.dumps({'columns':export_columns,'rows':[{k:p.get(k,'') for k in export_columns} for p in proteins]},ensure_ascii=False).replace('<', r'\u003c').replace('>', r'\u003e').replace('&', r'\u0026')
    index_by_id = {p['protein_id']:i for i,p in enumerate(proteins)}
    for cat in categories:
        rows = [p for p in proteins if p.get('category', 'unassigned') == cat['id']]
        summary.append(dict(category=cat['name'], annotation_records=len(rows),
                            interpretation='Supplied or reference-associated groups; inspect record evidence' if rows else 'No assigned annotation; not an absence claim'))
        details = []
        for p in rows:
            refs = ' · '.join('<a href="' + esc(sources[r]['url']) + '">' + esc(r) + '</a>'
                              for r in p.get('references', []) if r in sources)
            metrics = table([p], ['length','status','aliphatic_index','theoretical_pi','gravy',
                                  'instability_index','molecular_weight_da','aromaticity'])
            details.append('<details class="protein" data-index="' + str(index_by_id[p['protein_id']]) + '"><summary><strong>' + esc(p.get('gene_name', 'unassigned')) +
                           '</strong> · ' + esc(p.get('protein_name', p.get('description'))) +
                           ' <span class="muted">[' + esc(p['protein_id']) + ']</span></summary>' +
                           '<p><b>Function:</b> ' + esc(p.get('function', 'Not supplied')) + '</p>' +
                           '<p><b>Annotation:</b> ' + esc(p.get('annotation_status', 'unassessed')) + ' · ' +
                           esc(p.get('annotation_basis', 'Not assessed')) + '</p>' + metrics +
                           '<p><b>Evidence level:</b> ' + esc(p.get('evidence_level', 'unassessed')) +
                           '<br><b>Literature scope:</b> ' + esc(p.get('literature_scope', 'Not assessed')) + '</p>' +
                           '<p><b>Descriptor interpretation:</b> ' + esc(p.get('instability_class', 'Not assessed')) +
                           '; ' + esc(p.get('hydropathy_class', 'Not assessed')) + '.</p>' +
                           '<p><b>Temperature:</b> ' + esc(p.get('temperature_indication', 'Undetermined')) +
                           '<br><b>Activity pH:</b> ' + esc(p.get('ph_indication', 'Undetermined')) + '</p>' +
                           '<p><b>Input annotation:</b> ' + esc(p.get('description', '')) + '<br><b>Quality:</b> ' +
                           esc(p.get('notes') or 'No residue flags; completeness remains unassessed') + '</p>' +
                           '<p><b>References:</b> ' + (refs or 'No catalog reference assigned') + '</p>' + reference_details(p, search_results) + protein_context(p, context) + '</details>')
        section = ('<details class="card category" data-category="' + esc(cat['id']) + '" data-empty="' + ('false' if rows else 'true') + '"><summary><strong>' + esc(cat['name']) + '</strong><span class="count">' + str(len(rows)) + ' records</span></summary><p>' + esc(cat['description']) + '</p>' + (''.join(details) or '<p>No assigned records; this does not establish biological absence.</p>') + '</details>')
        (sections if rows else empty_sections).append(section)
    options = ''.join('<option value="' + esc(c['id']) + '">' + esc(c['name']) + '</option>' for c in categories)
    quality = '<details class="card" id="quality"><summary><strong>Input quality</strong><span class="count">' + str(len(flags)) + ' flagged records</span></summary><p>Residue checks only; no inferred truncation, frameshift or biological completeness.</p>' + (table(flags,['protein_id','flags']) if flags else '<p>No residue flags recorded.</p>') + ''.join('<p>' + esc(w) + '</p>' for w in context.get('warnings',[])) + '</details>'
    evidence_columns = ['sample_id','subject_id','evidence_type','property','value','unit','source','notes']
    evidence = table(result['evidence'], evidence_columns) if result['evidence'] else '<p>No measured or literature observations supplied.</p>'
    bibliography = ''.join('<li id="ref-' + esc(s['id']) + '"><a href="' + esc(s['url']) + '">' +
                           esc(s['title']) + '</a> · ' + esc(s['type']) + '</li>' for s in sources.values())
    prov = result['provenance']
    version = result.get('software_version', 'unknown')
    return '''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>''' + esc(application_title) + ''' · ''' + esc(result['sample_id']) + '''</title>
<style>
:root{font-family:system-ui,sans-serif;color:#173443;background:#eef3f5}*{box-sizing:border-box}
body{max-width:1320px;margin:auto;padding:28px}header{background:#123c49;color:white;padding:30px;border-radius:16px}
h1{font-size:36px;margin:8px 0}h2{margin-top:32px}.eyebrow{letter-spacing:3px;color:#88e3ca}
.card{background:white;border:1px solid #ccdce2;border-radius:12px;padding:20px;margin-top:18px}
.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:12px;border-bottom:1px solid #dbe5e9;vertical-align:top}th{background:#e4f1ee}
.muted{color:#526974}a{color:#066e70}nav{display:flex;gap:18px;flex-wrap:wrap;margin:22px 0}
.stats{display:flex;gap:16px;flex-wrap:wrap}.stats span{flex:1;min-width:160px;padding:16px;background:#e4f1ee;border-radius:8px}
.protein{border-top:1px solid #ccdce2;padding:16px 0;overflow-wrap:anywhere}summary{cursor:pointer;line-height:1.6}
select{font:inherit;padding:12px;max-width:100%;border-radius:8px} .count{float:right;color:#526974;font-size:13px;margin-left:12px}.category>summary{padding:3px}.category p{font-size:14px}.subpanel{margin:12px 0;padding:12px;background:#f5f8f9;border-radius:8px}svg{width:100%;height:auto}.toolbar{display:flex;flex-wrap:wrap;gap:8px;align-items:center}.card>summary{font-size:17px}#overview h2{margin-top:0} .empty-categories{margin-top:18px}
input{font:inherit;padding:12px;width:100%;border:1px solid #75949d;border-radius:8px}button{font:inherit;padding:9px 16px;margin:10px 8px 0 0;cursor:pointer}
[hidden]{display:none!important}li{margin:12px 0}p{line-height:1.6}code{overflow-wrap:anywhere}
@media(max-width:600px){body{padding:12px}.card,header{padding:16px}h1{font-size:28px}}
@media print{body{padding:0;background:white}input,button{display:none}.scroll{overflow:visible}details{break-inside:avoid}}
</style></head><body><header><div class="eyebrow">BACTOMICS / GENESCOPE</div><p>''' + esc(application_title) + ''' · ''' + esc(app.get("version", "No application selected")) + '''</p><h1>''' + esc(result['sample_id']) + '''</h1>
<p>Gene context, protein properties & application evidence · ''' + esc(version) + '''</p></header>
<nav aria-label="Report sections"><a href="#overview">Overview</a><a href="#proteins">Genes & proteins</a><a href="#conditions">Temperature & pH</a><a href="#evidence">Evidence</a><a href="#methods">Methods & sources</a></nav>
<section class="card" id="overview"><h2>Overview</h2><div class="stats"><span><b>''' + str(len(proteins)) + '''</b> protein records</span><span><b>''' + str(linked) + '''</b> catalog-linked annotations</span><span><b>''' + str(calculated) + '''</b> calculated descriptors</span></div>
<p>Explore functional categories below. Names are supplied annotations; operating optima remain undetermined.</p><details><summary>Scope and interpretation</summary><p>''' + esc(result.get('limitations', '')) + '</p>' + table(summary, ['category','annotation_records','interpretation']) + '''</details></section>''' + quality + reference_summary(search_results) + '''
<section id="proteins"><h2>Genes & proteins</h2><p>Names and functions come from supplied annotation and a cited catalog. Database similarity candidates, when searched, are shown separately with their scores and reference evidence. Expand a record for descriptors, context and references.</p>
<label for="protein-search">Find a gene, protein, function or identifier</label><input id="protein-search" type="search" placeholder="Search annotations…">
<div class="toolbar"><label for="category-filter">Category</label><select id="category-filter"><option value="">All categories</option>''' + options + '''</select><button type="button" id="expand">Expand visible records</button><button type="button" id="collapse">Collapse all</button><button type="button" id="export">Export filtered TSV</button><button type="button" id="print">Print / Save PDF</button></div>
<p id="search-status" role="status" aria-live="polite"></p>''' + ''.join(sections) + '<details class="empty-categories" id="empty-categories"><summary>Show unassessed categories (' + str(len(empty_sections)) + ')</summary>' + ''.join(empty_sections) + '''</details></section>
<details class="card" id="conditions"><summary><strong>Temperature & pH indications</strong></summary>
<p><b>Best operating temperature: undetermined. Best operating pH: undetermined.</b></p>
<p>Aliphatic index describes side-chain composition and has a reported association with thermostability of globular proteins. The reviewed method supplies no universal cutoff that converts this index into a temperature or temperature range. Thermostability and activity optimum are different measurements.</p>
<p>Theoretical pI estimates the charge-neutral point; it is not an enzyme's activity optimum. GRAVY describes mean hydropathy. The instability-index threshold is a sequence heuristic: below 40 is the stable class, above 40 the unstable class, and exactly 40 is left unclassified here. None establishes formulation stability or shelf life.</p>
<p>Source-linked measured temperature and pH observations can be recorded below with protein/strain identity, assay endpoint and preparation context in notes. Conflicting observations remain visible; values are not averaged or transferred to a different sample.</p></details>
<details class="card" id="evidence"><summary><strong>Application evidence</strong></summary><p>Supplied observations are displayed as recorded and have not been independently verified. A literature value is not automatically applicable to this sample.</p>''' + evidence + '''</details>
<details class="card" id="methods"><summary><strong>Methods & sources</strong></summary>
<p>Biopython ''' + esc(prov.get('biopython_version', 'unknown')) + ''' ProteinAnalysis calculates mass, theoretical pI, GRAVY, instability and aromaticity. Aliphatic index = 100 × [A + 2.9V + 3.9(I + L)] / sequence length, using residue counts. One terminal stop is removed; unsupported residues are flagged without descriptor calculation.</p>
<p>Annotation labels are matched conservatively to catalog gene symbols or product phrases. Ambiguous matches stay unresolved. A selected reference database enables a separate DIAMOND similarity search. No pathway-completeness or gene-expression inference is performed. Generic function annotations are retained even without a catalog match.</p>
<p>Catalog ''' + esc(catalog.get('catalog_version', 'not recorded')) + ' · ' + esc(catalog.get('scope', '')) + '''</p>
<p>Coverage gaps: ''' + esc('; '.join(catalog.get('coverage_gaps', []))) + '''.</p><ul>''' + bibliography + '''</ul>
<p>Input SHA-256: <code>''' + esc(prov['input_sha256']) + '''</code><br>Catalog SHA-256: <code>''' + esc(prov.get('catalog_sha256', 'not recorded')) + '''</code><br>Generated: ''' + esc(prov['created_utc']) + '''</p></details>
<script id="export-data" type="application/json">''' + export_data + '''</script>
<script type="text/javascript">
const search=document.getElementById('protein-search'),select=document.getElementById('category-filter');
const records=[...document.querySelectorAll('.protein')],categories=[...document.querySelectorAll('.category')];
const empties=document.getElementById('empty-categories');
function filter(){
 const q=search.value.toLocaleLowerCase().trim(), cat=select.value;let n=0;
 records.forEach(r=>{const match=(!cat||r.closest('.category').dataset.category===cat)&&(!q||r.textContent.toLocaleLowerCase().includes(q));r.hidden=!match;r.open=!!q&&match;if(match)n++});
 categories.forEach(c=>{const items=[...c.querySelectorAll('.protein')];c.hidden=items.length?!items.some(r=>!r.hidden):(!!q||!!cat&&c.dataset.category!==cat);c.open=!!q&&!c.hidden;});
 empties.hidden=!!q; if(cat&&!q&&categories.some(c=>c.dataset.category===cat&&c.dataset.empty==='true'))empties.open=true;
 document.getElementById('search-status').textContent=n+' of '+records.length+' records match';
}
search.addEventListener('input',filter);select.addEventListener('change',filter);filter();
document.getElementById('expand').addEventListener('click',()=>{categories.filter(c=>!c.hidden&&c.dataset.empty==='false').forEach(c=>c.open=true);records.filter(r=>!r.hidden).forEach(r=>r.open=true)});
document.getElementById('collapse').addEventListener('click',()=>document.querySelectorAll('details').forEach(d=>d.open=false));
document.querySelectorAll('nav a').forEach(a=>a.addEventListener('click',()=>{const d=document.querySelector(a.getAttribute('href'));if(d&&d.tagName==='DETAILS')d.open=true}));
document.getElementById('export').addEventListener('click',()=>{
 const data=JSON.parse(document.getElementById('export-data').textContent);
 const clean=v=>{let s=String(v??'').replace(/[\\t\\r\\n]/g,' ');return /^[=+@-]/.test(s)?"'"+s:s};
 const rows=records.filter(r=>!r.hidden).map(r=>data.rows[Number(r.dataset.index)]);
 const text=[data.columns.join('\\t'),...rows.map(r=>data.columns.map(k=>clean(r[k])).join('\\t'))].join('\\n');
 const url=URL.createObjectURL(new Blob([text],{type:'text/tab-separated-values;charset=utf-8'}));const a=document.createElement('a');a.href=url;a.download='micp-proteins-filtered.tsv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
});
let printState=[];
window.addEventListener('beforeprint',()=>{printState=[...document.querySelectorAll('details')].map(d=>[d,d.open]);printState.forEach(([d])=>{if(!d.hidden)d.open=true})});
window.addEventListener('afterprint',()=>printState.forEach(([d,open])=>d.open=open));
document.getElementById('print').addEventListener('click',()=>window.print());
</script></body></html>'''

