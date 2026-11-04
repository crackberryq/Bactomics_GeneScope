"""Offline coordinate explorer; supplied operons are never inferred from proximity."""
import json


def explorer(result):
    context = result.get('context', {})
    data = dict(features=context.get('features', []), operons=context.get('operons', []),
                regions=context.get('sequence_regions', {}),
                proteins=[{k: p.get(k) for k in ('protein_id', 'gene_name', 'protein_name', 'category', 'aliphatic_index', 'theoretical_pi')}
                          for p in result['proteins']])
    encoded = json.dumps(data, ensure_ascii=False).replace('<', r'\u003c').replace('>', r'\u003e').replace('&', r'\u0026')
    return '''<details class="card" id="genome-explorer"><summary><strong>Genome explorer · coordinates & operon annotations</strong></summary>
<p>Coordinates are 1-based, inclusive. Arrows show CDS strand; the ? track indicates unknown strand. Grey CDS marks have no matching protein record. Teal marks have a matching protein record. Purple tracks are supplied operon annotations, not independently validated predictions. Proximity alone does not establish co-transcription or operon completeness.</p>
<div class="toolbar"><label>Contig <select id="gx-contig"></select></label><label>Start <input id="gx-start" type="number" min="1"></label><label>End <input id="gx-end" type="number" min="1"></label><button id="gx-go" type="button">Show region</button><button id="gx-reset" type="button">Whole contig</button><button id="gx-left" type="button">← Pan left</button><button id="gx-in" type="button">Zoom in</button><button id="gx-out" type="button">Zoom out</button><button id="gx-right" type="button">Pan right →</button></div>
<div class="toolbar"><label for="gx-find">Find gene name, protein ID or feature ID</label><input id="gx-find" type="search" placeholder="Gene name or identifier…"><button id="gx-locate" type="button">Locate gene</button></div><div id="gx-results" aria-label="Matching genes"></div><p id="gx-status" role="status"></p><div class="gx-scroll"><svg id="gx-map" viewBox="0 0 1000 210" role="img" aria-label="Genome coordinates and strands"></svg></div>
<p id="gx-selection" role="status">Select a gene arrow or an annotation below.</p><details class="subpanel"><summary>Visible feature table</summary><div id="gx-list" class="scroll"></div></details>
</details><details class="card" id="visual-summary"><summary><strong>Visual summary · annotation groups & protein properties</strong></summary><p>Counts reflect assigned annotation groups, not pathway completeness. Theoretical pI and aliphatic index are sequence descriptors, not operating pH or temperature predictions.</p><div id="gx-categories"></div><svg id="gx-properties" viewBox="0 0 1000 330" role="img" aria-label="Theoretical pI versus aliphatic index"></svg><p id="gx-properties-status"></p></details>
<script id="gx-data" type="application/json">''' + encoded + '''</script><script type="text/javascript">
(()=>{
const data=JSON.parse(document.getElementById('gx-data').textContent),byId=new Map(data.proteins.map((p,i)=>[p.protein_id,{...p,index:i}]));
let focused=null;
const el=id=>document.getElementById(id),ns='http://www.w3.org/2000/svg';
function svg(tag,attrs,text){const n=document.createElementNS(ns,tag);Object.entries(attrs).forEach(([k,v])=>n.setAttribute(k,v));if(text!==undefined)n.textContent=text;return n}
function jump(p){if(!p)return;el('protein-search').value='';el('category-filter').value='';el('protein-search').dispatchEvent(new Event('input'));const r=document.querySelector('.protein[data-index="'+p.index+'"]');if(r){r.hidden=false;r.closest('.category').hidden=false;r.closest('.category').open=true;r.open=true;r.scrollIntoView({block:'start'});}}
const names=[...new Set([...data.features,...data.operons].map(f=>f.contig))].sort();names.forEach(name=>{const o=document.createElement('option');o.value=name;o.textContent=name;el('gx-contig').append(o)});
function bounds(){const c=el('gx-contig').value,fs=[...data.features,...data.operons].filter(f=>f.contig===c),decl=data.regions[c];return decl||[1,Math.max(1,...fs.map(f=>f.end))]}
function selectFeature(f){focused=f;draw();const p=byId.get(f.protein_id);el('gx-selection').textContent=(f.feature_id||'Feature')+' · '+f.start+'–'+f.end+' · strand '+f.strand+' · '+(f.product||'')+(f.source?' · source '+f.source:'');if(p)jump(p)}
function draw(){const c=el('gx-contig').value,[lo,hi]=bounds(),start=Number(el('gx-start').value),end=Number(el('gx-end').value);if(!Number.isInteger(start)||!Number.isInteger(end)||start<lo||end>hi||start>end){el('gx-status').textContent='Enter an integer interval within '+lo+'–'+hi+'.';return}
const map=el('gx-map');map.replaceChildren();const x=v=>50+900*(v-start)/(end-start+1),visible=data.features.filter(f=>f.contig===c&&f.end>=start&&f.start<=end),ops=data.operons.filter(f=>f.contig===c&&f.end>=start&&f.start<=end);
for(const [y,label] of [[65,'+'],[115,'−'],[155,'?'],[190,'Operon']]){map.append(svg('line',{x1:50,x2:950,y1:y,y2:y,stroke:'#cedbe0'}),svg('text',{x:2,y:y+4,'font-size':12},label))}
map.append(svg('text',{x:50,y:22,'font-size':14},String(start)),svg('text',{x:950,y:22,'text-anchor':'end','font-size':14},String(end)));
const limit=500,shown=visible.slice(0,limit);for(const f of [...shown,...ops.slice(0,100)]){const operon=ops.includes(f),a=x(Math.max(start,f.start)),b=Math.max(a+2,x(Math.min(end+1,f.end+1))),y=operon?190:f.strand==='+'?65:f.strand==='-'?115:155,h=operon?5:10,t=Math.min(8,(b-a)/2);let points=`${a},${y-h} ${b},${y-h} ${b},${y+h} ${a},${y+h}`;if(!operon&&f.strand==='+')points=`${a},${y-h} ${b-t},${y-h} ${b},${y} ${b-t},${y+h} ${a},${y+h}`;if(!operon&&f.strand==='-')points=`${a},${y} ${a+t},${y-h} ${b},${y-h} ${b},${y+h} ${a+t},${y+h}`;const n=svg('polygon',{points,stroke:f===focused?'#d89925':'none','stroke-width':3,fill:operon?'#7858a6':f.protein_id?'#087c78':'#84949b',tabindex:0,role:'button','aria-label':`${f.feature_id} ${f.start}–${f.end} strand ${f.strand}`});n.append(svg('title',{},`${f.feature_id}: ${f.product||''} (${f.start}–${f.end})`));n.addEventListener('click',()=>selectFeature(f));n.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();selectFeature(f)}});map.append(n);if(!operon&&b-a>35){const label=byId.get(f.protein_id)?.gene_name;map.append(svg('text',{x:(a+b)/2,y:y-15,'text-anchor':'middle','font-size':11,fill:'#3d586b'},label&&label!=='unassigned'?label:f.feature_id))}}
el('gx-status').textContent=`${c}: ${start}–${end}. ${data.regions[c]?'Declared sequence-region bounds.':'Extent ends at the last annotated feature; full contig length unknown.'} ${visible.length} CDS; ${ops.length} supplied operon annotations. ${visible.length>limit?'Showing first 500 CDS; narrow the interval to inspect the rest.':''} ${ops.length>100?'Showing first 100 operons; narrow the interval.':''}`;
const list=el('gx-list');list.replaceChildren();const table=document.createElement('table'),head=document.createElement('tr');for(const label of ['Feature','Coordinates','Strand','Annotation']){const th=document.createElement('th');th.textContent=label;head.append(th)}table.append(head);for(const f of [...ops.slice(0,100),...shown]){const row=document.createElement('tr'),cell=document.createElement('td'),button=document.createElement('button');button.textContent=f.feature_id;button.onclick=()=>selectFeature(f);cell.append(button);row.append(cell);for(const value of [f.start+'–'+f.end,f.strand,f.product||'']){const td=document.createElement('td');td.textContent=value;row.append(td)}table.append(row)}list.append(table);}
function reset(){const [a,b]=bounds();el('gx-start').value=a;el('gx-end').value=b;draw()}
function locate(f){
 focused=f;el('genome-explorer').open=true;el('gx-contig').value=f.contig;
 const [a,b]=bounds();el('gx-start').value=Math.max(a,f.start-5000);el('gx-end').value=Math.min(b,f.end+5000);
 draw();el('gx-selection').textContent='Located '+(byId.get(f.protein_id)?.gene_name||f.feature_id)+' · '+f.feature_id+' · '+f.contig+':'+f.start+'–'+f.end+' · strand '+f.strand;
 el('genome-explorer').scrollIntoView({block:'start'});
}
function showMatches(matches,q){
 const holder=el('gx-results');holder.replaceChildren();
 if(matches.length===1){locate(matches[0]);return}
 el('gx-status').textContent=matches.length?matches.length+' matches. Choose a coordinate below.':'No matching gene coordinates for '+q+'.';
 for(const f of matches.slice(0,100)){const b=document.createElement('button');b.type='button';b.textContent=(byId.get(f.protein_id)?.gene_name||f.feature_id)+' · '+f.feature_id+' · '+f.contig+':'+f.start+'–'+f.end+' ('+f.strand+')';b.onclick=()=>locate(f);holder.append(b)}
 if(matches.length>100){const note=document.createElement('p');note.textContent='First 100 matches shown; refine your search.';holder.append(note)}
}
function findGene(){
 const q=el('gx-find').value.trim().toLocaleLowerCase();if(!q){el('gx-results').replaceChildren();el('gx-status').textContent='Enter a gene name or identifier.';return}
 const labels=f=>[f.protein_id,f.feature_id,f.gene,f.gene_name,f.name,byId.get(f.protein_id)?.gene_name].filter(Boolean).map(v=>String(v).toLocaleLowerCase());
 const exact=data.features.filter(f=>labels(f).includes(q));
 showMatches(exact.length?exact:data.features.filter(f=>labels(f).some(v=>v.includes(q))),q);
}
el('gx-locate').onclick=findGene;el('gx-find').onkeydown=e=>{if(e.key==='Enter'){e.preventDefault();findGene()}};
document.addEventListener('click',e=>{const b=e.target.closest('.locate-protein');if(!b)return;el('genome-explorer').open=true;el('gx-find').value=b.dataset.protein;showMatches(data.features.filter(f=>f.protein_id===b.dataset.protein||f.feature_id===b.dataset.protein),b.dataset.protein);el('genome-explorer').scrollIntoView({block:'start'})});
function move(scale,shift){const [lo,hi]=bounds(),a=Number(el('gx-start').value),b=Number(el('gx-end').value),span=Math.min(hi-lo+1,Math.max(20,Math.round((b-a+1)*scale))),center=(a+b)/2+(b-a+1)*shift;const start=Math.max(lo,Math.min(hi-span+1,Math.round(center-span/2)));el('gx-start').value=start;el('gx-end').value=start+span-1;draw()}
el('gx-left').onclick=()=>move(1,-.5);el('gx-right').onclick=()=>move(1,.5);el('gx-in').onclick=()=>move(.5,0);el('gx-out').onclick=()=>move(2,0);
el('gx-contig').onchange=reset;el('gx-go').onclick=draw;el('gx-reset').onclick=reset;if(names.length)reset();else{el('gx-status').textContent='No GFF coordinates supplied. Add the matching annotation GFF to view genome locations.';for(const id of ['gx-contig','gx-start','gx-end','gx-go','gx-reset','gx-left','gx-right','gx-in','gx-out'])el(id).disabled=true}
const counts=new Map();for(const p of data.proteins)counts.set(p.category||'unassigned',(counts.get(p.category||'unassigned')||0)+1);for(const [name,count] of [...counts].sort((a,b)=>b[1]-a[1])){const line=document.createElement('p'),label=document.createElement('span'),bar=document.createElement('meter');label.textContent=name+' · '+count+' ';bar.min=0;bar.max=Math.max(1,data.proteins.length);bar.value=count;bar.setAttribute('aria-label',name+' record count');line.append(label,bar);el('gx-categories').append(line)}
const points=data.proteins.filter(p=>typeof p.theoretical_pi==='number'&&Number.isFinite(p.theoretical_pi)&&typeof p.aliphatic_index==='number'&&Number.isFinite(p.aliphatic_index)),plot=el('gx-properties');const maxAI=Math.max(100,...points.map(p=>p.aliphatic_index)),maxPI=Math.max(14,...points.map(p=>p.theoretical_pi));plot.append(svg('path',{d:'M60 20V280H960',fill:'none',stroke:'#526974'}),svg('text',{x:460,y:320,'font-size':16},'Theoretical pI (0–'+maxPI+')'),svg('text',{x:65,y:16,'font-size':14},'Aliphatic index (0–'+maxAI.toFixed(1)+')'));
for(let i=0;i<=4;i++){const x=60+i*225,y=280-i*62.5;plot.append(svg('text',{x,y:300,'text-anchor':'middle','font-size':12},(maxPI*i/4).toFixed(1)),svg('text',{x:50,y:y+4,'text-anchor':'end','font-size':12},(maxAI*i/4).toFixed(1)))}
for(const p of points){const n=svg('circle',{cx:60+p.theoretical_pi/maxPI*900,cy:280-p.aliphatic_index/maxAI*250,r:3,fill:'#087c78',opacity:.45});n.append(svg('title',{},`${p.protein_id}: pI ${p.theoretical_pi}, AI ${p.aliphatic_index}`));n.onclick=()=>jump(byId.get(p.protein_id));plot.append(n)}el('gx-properties-status').textContent=points.length+' proteins plotted; '+(data.proteins.length-points.length)+' without finite descriptors omitted. Hover for values; click a point for its protein record. All values are also available in the searchable protein tables.';
})();</script>'''
