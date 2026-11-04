"""Standalone, offline report with category summaries and inspectable evidence."""
import html
import json
from .context import neighborhood
from .explorer import explorer
from .report_assets import CSS


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


def assessment_summary(result):
    search = result.get('reference_search')
    app = result.get('application') or {}
    proteins = result['proteins']
    if not search:
        title = 'Descriptor-only run — application assessment not performed'
        explanation = ('No reference database search was requested. Unassigned proteins are not negative findings. '
                       'This report contains supplied annotations, calculated protein properties and available coordinates.')
        next_step = ('In Bactologix, build or register an application database under Databases & models → GeneScope. '
                     'Select the completed database for a new or edited draft, review the configuration, and run it. '
                     'Building a database does not attach it to an existing frozen run.')
    else:
        matched = sum(bool(q.get('hits')) for q in search.get('queries', {}).values())
        title = 'Application reference search — ' + str(app.get('name') or search.get('database', {}).get('name') or 'selected database')
        explanation = (str(matched) + ' proteins have qualifying reference candidates. '
                       'Similarity candidates require review; proteins without qualifying hits are not proof of biological absence.')
        next_step = 'Review the reference scores, coverage, competing matches and source evidence in the expandable protein records.'
    rows = [
        {'evidence':'Protein descriptors', 'coverage':str(sum(p['status'] == 'calculated' for p in proteins)) + ' calculated'},
        {'evidence':'Application reference search', 'coverage':'Recorded' if search else 'Not requested'},
        {'evidence':'Gene coordinates', 'coverage':str(sum(bool(p.get('coordinates')) for p in proteins)) + ' proteins with coordinates'},
        {'evidence':'Imported signatures', 'coverage':str(sum(bool(p.get('domains')) for p in proteins)) + ' proteins with recorded matches'},
        {'evidence':'Operon / promoter prediction', 'coverage':'Not assessed by this workflow'},
        {'evidence':'Measured or literature observations', 'coverage':str(len(result.get('evidence', []))) + ' supplied records'},
    ]
    return ('<section class="card" id="assessment-status" role="status"><h2>' + esc(title) +
            '</h2><p>' + esc(explanation) + '</p><p><strong>Next step:</strong> ' + esc(next_step) +
            '</p><details><summary>Analysis coverage and missing evidence</summary>' +
            table(rows, ['evidence','coverage']) + '</details></section>')


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
    unassigned_sections = []
    def unmatched(p):
        return (bool(search_results) and not p.get("reference_hits")) or p.get("category", "unassigned") == "unassigned"
    groups = [dict(c, _unassigned=False) for c in categories]
    groups.append(dict(id="unassigned", name="Unassigned proteins", description="No qualifying match to the selected database or no assigned application category. Existing annotations are retained; this is not evidence of biological absence.", _unassigned=True))
    empty_sections = []
    context = result.get('context', {})
    flags = [dict(protein_id=p['protein_id'], flags=', '.join(p.get('quality_flags', [])) or p.get('notes','')) for p in proteins if p.get('quality_flags') or p['status'] != 'calculated']
    export_columns = ['protein_id','gene_name','protein_name','category','function','length','aliphatic_index','theoretical_pi','gravy','instability_index','annotation_status']
    export_data = json.dumps({'columns':export_columns,'rows':[{k:p.get(k,'') for k in export_columns} for p in proteins]},ensure_ascii=False).replace('<', r'\u003c').replace('>', r'\u003e').replace('&', r'\u0026')
    index_by_id = {p['protein_id']:i for i,p in enumerate(proteins)}
    for cat in groups:
        rows = [p for p in proteins if unmatched(p)] if cat["_unassigned"] else [p for p in proteins if not unmatched(p) and p.get("category", "unassigned") == cat["id"]]
        if cat["id"] == "unassigned" and not cat["_unassigned"]:
            continue
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
                           '<button type="button" class="locate-protein" data-protein="' + esc(p['protein_id']) + '">View in genome explorer</button>' +
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
        if cat['_unassigned']:
            unassigned_sections.append(section)
        else:
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
<style>''' + CSS + r'''
main{overflow-wrap:anywhere}.card{background:#fff;border:1px solid var(--line);border-radius:14px;padding:24px;margin-top:20px;box-shadow:var(--shadow);min-width:0}.card>summary{font-size:16px;font-weight:650;cursor:pointer}.card[open]>summary{margin-bottom:20px}#proteins{margin-top:36px}.scroll{overflow:auto}.protein{border-top:1px solid var(--line);padding:16px 0}.protein>summary{cursor:pointer;font-size:13px}.subpanel{margin:12px 0;padding:15px;background:#f6f8fa;border:1px solid var(--line);border-radius:9px}.count{float:right;color:var(--muted);font-size:11px}.toolbar{display:flex;flex-wrap:wrap;align-items:end;gap:10px;margin:16px 0}.toolbar label{display:flex;flex-direction:column;font-size:11px;font-weight:650;gap:5px}input,select{max-width:100%;min-height:40px;padding:9px 12px;border:1px solid #cbd8df;border-radius:7px;background:#fff;color:var(--ink)}input[type=search]{width:100%}button{border:1px solid #cbd8df;border-radius:7px;padding:9px 14px;background:#fff;color:#294a62;font-size:12px}button:hover{background:var(--teal-soft)}button:disabled{opacity:.5;cursor:default}.locate-protein{margin:14px 0}.stats{grid-template-columns:repeat(3,minmax(0,1fr))}.stats>span{padding:22px;border:1px solid var(--line);border-radius:12px;font-size:11px;color:var(--muted);background:#fff}.stats b{display:block;color:var(--navy);font-size:32px;line-height:1.3;margin-bottom:8px}svg{max-width:100%;height:auto}#gx-list{max-height:300px}#gx-map{width:100%;min-width:780px}#gx-properties{width:100%}#gx-start,#gx-end{width:125px}#gx-find{flex:1;min-width:180px}.gx-scroll{overflow:auto;border:1px solid var(--line);border-radius:9px;background:#fff;margin:16px 0}#gx-selection{background:#eef7f5;border:1px solid #cee0e5;border-radius:10px;padding:18px}#gx-results{max-height:240px;overflow:auto}#gx-results button{display:block;width:100%;text-align:left;margin:5px 0}#gx-categories p{display:flex;justify-content:space-between;gap:12px;font-size:12px}meter{width:45%}#unassigned-records{margin-top:40px}.empty-categories{margin:20px 0;color:var(--muted)}summary{cursor:pointer}svg [role=button]{cursor:pointer}svg [role=button]:focus{stroke:#d89925;stroke-width:3}.hero{display:block}.main{padding-top:28px}.sidebar nav a{display:block;color:#bdcddd;padding:10px 12px;text-decoration:none;border-radius:9px;font-size:13px}.sidebar nav a:hover,.sidebar nav a[aria-current]{background:#223b53;color:#fff}.sidebar .brand{color:#fff;font-size:25px;font-weight:800;margin-bottom:30px}.sidebar-foot{font-size:11px;color:#a8bbc9}.hero p{color:#dbe8ef}#overview h2{margin-bottom:20px}#search-status,#gx-status{font-size:12px;color:var(--muted)}
@media(max-width:900px){.app{display:block}.sidebar{position:static;height:auto;padding:20px}.sidebar nav{flex-direction:row;flex-wrap:wrap;gap:4px}.sidebar .brand{margin-bottom:12px}.sidebar-foot{display:none}.main{padding:20px 16px}.hero{padding:26px}.stats{grid-template-columns:1fr}.card{padding:18px}.count{float:none;display:block}.toolbar{align-items:stretch}}
@media print{.app{display:block}.sidebar,.toolbar,.locate-protein{display:none}.main{padding:0}.card{break-inside:avoid}.scroll,.gx-scroll{overflow:visible}#gx-map{min-width:0}}
</style></head><body><a class="skip-link" href="#overview">Skip to report</a><div class="app"><aside class="sidebar"><a class="brand" href="#overview">GeneScope</a><p class="nav-label">Application report</p><nav class="report-nav" aria-label="Report sections"><a href="#overview">Overview</a><a href="#genome-explorer">Genome explorer</a><a href="#visual-summary">Visual summary</a><a href="#proteins">Genes & proteins</a><a href="#conditions">Temperature & pH</a><a href="#evidence">Evidence</a><a href="#methods">Methods & sources</a><a href="#unassigned-records">Unassigned proteins</a></nav><div class="sidebar-foot">Offline report · Bactomics</div></aside><main class="main"><header class="hero"><div class="eyebrow">BACTOMICS / GENESCOPE</div><p>''' + esc(application_title) + ''' · ''' + esc(app.get("version", "No application selected")) + '''</p><h1>''' + esc(result['sample_id']) + '''</h1>
<p>Gene context, protein properties & application evidence · ''' + esc(version) + '''</p></header>

''' + assessment_summary(result) + '''
<section class="card" id="overview"><h2>Overview</h2><div class="stats"><span><b>''' + str(len(proteins)) + '''</b> protein records</span><span><b>''' + str(linked) + '''</b> catalog-linked annotations</span><span><b>''' + str(calculated) + '''</b> calculated descriptors</span></div>
<p>Explore functional categories below. Names are supplied annotations; operating optima remain undetermined.</p><details><summary>Scope and interpretation</summary><p>''' + esc(result.get('limitations', '')) + '</p>' + table(summary, ['category','annotation_records','interpretation']) + '''</details></section>''' + explorer(result) + quality + reference_summary(search_results) + '''
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
<section id="unassigned-records" aria-label="Unassigned proteins">''' + ''.join(unassigned_sections) + '''</section>
<script id="export-data" type="application/json">''' + export_data + '''</script>
<script type="text/javascript">
const search=document.getElementById('protein-search'),select=document.getElementById('category-filter');
const records=[...document.querySelectorAll('.protein')],categories=[...document.querySelectorAll('.category')];
const empties=document.getElementById('empty-categories');
function filter(){
 const q=search.value.toLocaleLowerCase().trim(), cat=select.value;let n=0;
 records.forEach(r=>{const match=(!cat||r.closest('.category').dataset.category===cat)&&(!q||r.textContent.toLocaleLowerCase().includes(q));r.hidden=!match;r.open=!!q&&match;if(match)n++});
 categories.forEach(c=>{const items=[...c.querySelectorAll('.protein')];c.hidden=items.length?!items.some(r=>!r.hidden):(!!q||!!cat&&c.dataset.category!==cat);c.open=(!!q||!!cat)&&!c.hidden;});
 empties.hidden=!!q; if(cat&&!q&&categories.some(c=>c.dataset.category===cat&&c.dataset.empty==='true'))empties.open=true;
 document.getElementById('search-status').textContent=n+' of '+records.length+' records match';
}
search.addEventListener('input',filter);select.addEventListener('change',filter);filter();
document.getElementById('expand').addEventListener('click',()=>{categories.filter(c=>!c.hidden&&c.dataset.empty==='false').forEach(c=>c.open=true);records.filter(r=>!r.hidden).forEach(r=>r.open=true)});
document.getElementById('collapse').addEventListener('click',()=>document.querySelectorAll('details').forEach(d=>d.open=false));
document.querySelectorAll('nav a').forEach(a=>a.addEventListener('click',()=>{document.querySelectorAll('nav a').forEach(n=>n.removeAttribute('aria-current'));a.setAttribute('aria-current','location');const d=document.querySelector(a.getAttribute('href'));if(d&&d.tagName==='DETAILS')d.open=true}));
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
</script></main></div></body></html>'''

