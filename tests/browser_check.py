"""Exercise the generated offline report using synthetic, nonfunctional sequences."""
from pathlib import Path
import tempfile
import hashlib

from playwright.sync_api import sync_playwright
from micpscope.cli import legacy_main as main

with tempfile.TemporaryDirectory() as directory:
    base = Path(directory)
    fasta = base / 'test.faa'
    fasta.write_text('>p [gene=ureC] synthetic label fixture\nAVIL\n>q hypothetical protein\nMX\n')
    (base/'test.gff').write_text('##gff-version 3\nctg\tfixture\tCDS\t10\t21\t.\t+\t0\tID=p;product=synthetic\n')
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
        Path('ui-previews').mkdir(exist_ok=True)
        page.screenshot(path='ui-previews/micp-report-desktop.png', full_page=True)
        page.set_viewport_size({'width':390, 'height':844})
        assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth')
        page.screenshot(path='ui-previews/micp-report-mobile.png', full_page=True)
        assert not errors, errors
        browser.close()

