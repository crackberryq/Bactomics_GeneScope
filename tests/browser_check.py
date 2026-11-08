"""Exercise the generated offline report using synthetic, nonfunctional sequences."""
from pathlib import Path
import tempfile
import hashlib
import json
from micpscope.report import render_report

from playwright.sync_api import sync_playwright
from micpscope.cli import legacy_main as main

with tempfile.TemporaryDirectory() as directory:
    base = Path(directory)
    fasta = base / 'test.faa'
    fasta.write_text('>p [gene=ureC] synthetic label fixture\nAVIL\n>q hypothetical protein\nMX\n')
    (base/'test.gff').write_text('##gff-version 3\n##sequence-region ctg 1 1000\nctg\tfixture\toperon\t5\t50\t.\t+\t.\tID=op;Name=Example\nctg\tfixture\tCDS\t10\t21\t.\t+\t0\tID=p;product=synthetic\n')
    digest=hashlib.md5(b'AVIL').hexdigest()
    (base/'test.interpro.tsv').write_text(f'p\t{digest}\t4\tFixture\tTEST001\tSynthetic signature\t2\t3\t-\tT\t29-10-2026\n')
    main(['--proteins', str(fasta), '--sample', 'synthetic-demo', '--out', str(base/'out')])
    with sync_playwright() as runtime:
        browser = runtime.chromium.launch()
        page = browser.new_page(viewport={'width':1440, 'height':1000})
        errors = []
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.goto((base/'out/report.html').as_uri())
        assert page.locator('.protein').count() == 2
        assert page.locator('details[open]').count() == 0
        assert page.locator('#unassigned-records .protein').count() == 1
        assert page.locator('#unassigned-records details[open]').count() == 0
        assert page.evaluate("document.querySelector('#unassigned-records').compareDocumentPosition(document.querySelector('#methods')) & Node.DOCUMENT_POSITION_PRECEDING")
        page.get_by_label('Find a gene, protein, function or identifier').fill('ureC')
        assert page.locator('.protein:visible').count() == 1
        page.get_by_role('button', name='Expand visible records').click()
        assert page.locator('.protein[open]').count() == 1
        assert page.get_by_text('Catalytic component of the multimeric urea-hydrolysing enzyme.', exact=False).is_visible()
        page.get_by_role('button', name='Collapse all').click()
        assert page.locator('details[open]').count() == 0
        page.get_by_label('Find a gene, protein, function or identifier').fill('no-such-label')
        assert page.locator('.protein:visible').count() == 0
        page.get_by_label('Find a gene, protein, function or identifier').fill('')
        page.get_by_role('button', name='Expand visible records').click()
        page.locator('#category-filter').select_option('unassigned')
        page.get_by_role('button', name='Expand visible records').click()
        assert page.locator('.protein:visible').count() == 1
        with page.expect_download() as download:
            page.get_by_role('button', name='Export filtered TSV').click()
        content=Path(download.value.path()).read_text()
        assert '\nq\t' in content and '\np\t' not in content
        page.locator('#category-filter').select_option('')
        page.get_by_role('button', name='Collapse all').click()
        page.get_by_label('Find a gene, protein, function or identifier').fill('ureC')
        page.get_by_text('Protein signatures · 1', exact=True).click()
        assert page.get_by_role('img', name='Recorded protein signature positions').is_visible()
        page.get_by_text('Genome neighborhood · ctg', exact=True).click()
        assert page.get_by_text('Focal protein: p.', exact=False).is_visible()
        page.get_by_label('Find a gene, protein, function or identifier').fill('')
        page.get_by_role('button', name='Collapse all').click()
        page.get_by_role('link', name='Genome explorer', exact=True).click()
        assert page.locator('#genome-explorer').get_attribute('open') is not None
        assert '1 supplied operon annotations' in page.locator('#gx-status').inner_text()
        assert page.locator('#gx-map polygon').count() == 2
        page.locator('#gx-find').fill('ureC')
        page.locator('#gx-find').press('Enter')
        assert 'Located ureC' in page.locator('#gx-selection').inner_text()
        page.locator('#gx-start').fill('100')
        page.locator('#gx-go').click()
        assert page.locator('#gx-map polygon').count() == 0
        page.locator('#gx-reset').click()
        page.locator('#gx-find').fill('p')
        page.locator('#gx-locate').click()
        assert 'Located ureC' in page.locator('#gx-selection').inner_text()
        page.locator('#gx-map polygon').first.press('Enter')
        assert page.locator('.protein[data-index="0"]').get_attribute('open') is not None
        page.locator('.protein[data-index="0"] .locate-protein').click()
        assert page.locator('#genome-explorer').get_attribute('open') is not None
        assert 'Located ureC' in page.locator('#gx-selection').inner_text()
        page.get_by_role('link', name='Visual summary', exact=True).click()
        assert page.locator('#gx-properties circle').count() == 1
        Path('ui-previews').mkdir(exist_ok=True)
        page.get_by_role('link', name='Genome explorer', exact=True).click()
        page.evaluate('window.scrollTo({top:0,behavior:"instant"})')
        page.screenshot(path='ui-previews/genescope-explorer.png', full_page=True)
        page.get_by_role('button', name='Collapse all').click()
        Path('ui-previews').mkdir(exist_ok=True)
        page.evaluate('window.scrollTo({top:0,behavior:"instant"})')
        page.screenshot(path='ui-previews/micp-report-desktop.png', full_page=True)
        page.set_viewport_size({'width':390, 'height':844})
        assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
        page.evaluate('window.scrollTo({top:0,behavior:"instant"})')
        page.screenshot(path='ui-previews/micp-report-mobile.png', full_page=True)
        assert not errors, errors
        # Duplicate symbols must offer every locus, never silently select one.
        result=json.loads((base/'out/results.json').read_text())
        result['proteins'][1]['gene_name']='ureC'
        result['context']['features'].append(dict(result['context']['features'][0], feature_id='q', protein_id='q', start=100, end=110))
        (base/'duplicates.html').write_text(render_report(result))
        page.goto((base/'duplicates.html').as_uri())
        page.get_by_role('link', name='Genome explorer', exact=True).click()
        page.locator('#gx-find').fill('ureC')
        page.locator('#gx-find').press('Enter')
        assert page.locator('#gx-results button').count()==2
        page.locator('#gx-results button').nth(1).click()
        assert 'q · ctg:100–110' in page.locator('#gx-selection').inner_text()
        page.locator('#gx-find').fill('no-such-gene')
        page.locator('#gx-find').press('Enter')
        assert 'No matching gene coordinates' in page.locator('#gx-status').inner_text()
        assert page.locator('#gx-results button').count()==0
        # Database-unmatched records retain annotations, but appear only at the end.
        result['reference_search']={'queries':{'p':{'hits':[]},'q':{'hits':[]}},'database':{},'interpretation':'fixture','tool_version':'fixture','thresholds':{}}
        (base/'unmatched.html').write_text(render_report(result))
        page.goto((base/'unmatched.html').as_uri())
        assert page.locator('#proteins .protein').count()==0
        assert page.locator('#unassigned-records .protein').count()==2
        assert page.locator('#unassigned-records details[open]').count()==0
        assert not errors, errors
        # Conflicting and ambiguous candidates are visible before opening protein details.
        hit = dict(reference_id='ref', identity=99, query_coverage=100,
                   reference_coverage=100, evalue=1e-10, bitscore=100)
        result['reference_search'] = dict(
            queries={'p': {'status': 'reference_candidate', 'hits': [hit]},
                     'q': {'status': 'ambiguous_candidates', 'hits': [hit]}},
            references={'ref': {'id': 'ref', 'gene': 'tasA', 'protein_name': 'Synthetic reference', 'categories': ['matrix']}},
            knowledge_catalog=result['knowledge_catalog'], database={},
            interpretation='Synthetic state fixture', tool_version='fixture', thresholds={})
        (base/'conflict.html').write_text(render_report(result))
        page.goto((base/'conflict.html').as_uri())
        assert page.locator('#proteins .protein').count() == 0
        assert page.locator('#unassigned-records .protein').count() == 2
        assert 'Annotation / reference conflict' in page.locator('.protein[data-index="0"] > summary').inner_text()
        assert 'Ambiguous reference candidates' in page.locator('.protein[data-index="1"] > summary').inner_text()
        assert page.locator('#gx-categories button').count() == 0
        with page.expect_download() as download:
            page.get_by_role('button', name='Export filtered TSV').click()
        content = Path(download.value.path()).read_text()
        assert 'annotation_reference_conflict' in content and 'ambiguous_candidates' in content
        assert not errors, errors
        # A circular-origin CDS has two drawable intervals but remains one source CDS.
        circular = json.loads((base/'out/results.json').read_text())
        feature = dict(circular['context']['features'][0], start=995, end=1006,
                       circular_wrap=True, display_segments=[[995, 1000], [1, 6]])
        circular['context']['features'] = [feature]
        circular['context']['operons'] = []
        circular['proteins'][0]['coordinates'] = [feature]
        (base/'circular.html').write_text(render_report(circular))
        page.goto((base/'circular.html').as_uri())
        page.get_by_role('link', name='Genome explorer', exact=True).click()
        assert '1 CDS;' in page.locator('#gx-status').inner_text()
        assert '2 drawable CDS intervals' in page.locator('#gx-status').inner_text()
        assert page.locator('#gx-map polygon').count() == 2
        assert '1 genes' in page.locator('#gx-overview').text_content()
        assert not errors, errors
        # Large chromosome: late genes are included, names are readable and category bars exclude background.
        result['reference_search']=None
        result['proteins'][0]['gene_name']='unassigned'
        result['proteins'][0]['protein_name']='Named example protein'
        result['proteins'][1]['gene_name']='unknown'
        result['context']['sequence_regions']={'ctg':[1,4529685]}
        template=result['context']['features'][0]
        result['context']['features']=[dict(template,feature_id='feature-'+str(i),protein_id=None,start=100+i*6000,end=400+i*6000) for i in range(700)]
        result['context']['features'].append(dict(template,feature_id='p',protein_id='p',start=4200100,end=4201100))
        result['context']['operons']=[]
        result['proteins'][0]['coordinates']=[result['context']['features'][-1]]
        result['proteins'][1]['coordinates']=[]
        (base/'large.html').write_text(render_report(result))
        page.goto((base/'large.html').as_uri())
        page.get_by_role('link',name='Genome explorer',exact=True).click()
        assert int(page.locator('#gx-end').input_value())==20000
        assert '701 genes' in page.locator('#gx-overview').text_content()
        assert 'Operon' not in page.locator('#gx-map').text_content()
        page.locator('#gx-find').fill('Named example protein')
        page.locator('#gx-find').press('Enter')
        assert int(page.locator('#gx-start').input_value())>4000000
        assert 'Named example protein' in page.locator('#gx-selection').inner_text()
        page.locator('#gx-reset').click()
        assert 'Density view includes all CDS' in page.locator('#gx-status').inner_text()
        assert page.locator('#gx-map rect').count()>100
        page.get_by_role('link',name='Visual summary',exact=True).click()
        assert 'unassigned' not in page.locator('#gx-categories').inner_text()
        assert page.locator('#gx-background').get_attribute('open') is None
        assert 'Named example protein' in page.locator('#gx-comparison').inner_text()
        page.locator('#gx-metric').select_option('gravy')
        page.locator('#gx-categories button').first.click()
        assert page.locator('.protein[data-index="0"]').is_visible()
        assert 'unassigned' not in page.locator('.protein[data-index="0"] > summary').inner_text()
        page.get_by_role('link',name='Visual summary',exact=True).click()
        page.screenshot(path='ui-previews/genescope-focused-summary.png',full_page=True)
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
        assert not errors, errors
        browser.close()

