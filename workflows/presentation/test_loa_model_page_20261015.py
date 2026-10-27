"""Bounded compute-node browser checks for the exact two Page fragments."""
import functools
import http.server
import json
import base64
import gzip
import math
import os
from pathlib import Path
import sys
import struct
import threading

BASE = Path('/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1')
ROOT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_model_pilot_20261015_v1/page_v3')
assert os.environ.get('SLURM_JOB_ID')
sys.path.insert(0, str(BASE/'playwright_runtime'))
os.environ['PLAYWRIGHT_BROWSERS_PATH'] = str(BASE/'browser_runtime')
from playwright.sync_api import sync_playwright

class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass

server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(Quiet, directory=str(ROOT)))
threading.Thread(target=server.serve_forever, daemon=True).start()
url = f'http://127.0.0.1:{server.server_port}/'
errors, external, checks = [], [], []
try:
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={'width':1200,'height':1000}, accept_downloads=True)
        def onerror(e):
            errors.append(str(e))
            print('BROWSER_ERROR: '+str(e), flush=True)
        page.on('pageerror', onerror)
        def guard(route):
            if route.request.url.startswith(url):
                route.continue_()
            else:
                external.append(route.request.url)
                route.abort()
        page.route('**/*', guard)
        for name in ['field_NGC', 'field_SGC', 'properties']:
            print('CHECKING '+name, flush=True)
            page.goto(url+name+'.html')
            page.wait_for_function('window.__ready || window.__error', timeout=30000)
            assert page.evaluate('window.__ready'), page.evaluate('window.__error')
            assert not page.evaluate('document.documentElement.scrollWidth > innerWidth')
            if name.startswith('field_'):
                assert page.evaluate('window.__fieldState.meshes.every(m=>m.valid>0 && m.triangles>0)')
                for cap in ['0']:
                    page.select_option('#cap', cap)
                    for model in map(str, range(6)):
                        page.select_option('#a', model)
                        assert page.evaluate('window.__fieldState.meshes[0].valid>0')
                    for mode in ['environment','entropy','density']:
                        page.select_option('#mode', mode)
                        assert page.evaluate('window.__fieldState.meshes.every(m=>Number.isFinite(m.min))')
                page.fill('#zmin','0.4');page.locator('#zmin').dispatch_event('change')
                assert page.evaluate('window.__fieldState.meshes.every(m=>m.valid===0)')
                page.fill('#zmin','0.15');page.locator('#zmin').dispatch_event('change')
                page.select_option('#a','2');page.select_option('#cap','0')
                checks.append(name+': six fits, three fields and empty redshift selection')
            else:
                assert page.evaluate('window.__propertyState.count>0')
                # Independently decode the binary rows and reproduce displayed
                # statistics, rather than merely checking that marks exist.
                html = (ROOT/'properties.html').read_text()
                data = json.loads(html.split('const DATA=',1)[1].split(',$=id=>',1)[0])
                raw = gzip.decompress(base64.b64decode(data['data']))
                source_rows = list(struct.iter_unpack('<7B3f',raw))
                actual = page.evaluate('window.__propertyState.shown')
                for mi in [2,4]:
                    for env in range(4):
                        subset=[r for r in source_rows if r[1]==mi and r[5]==env and r[6]==6]
                        n=sum(r[7] for r in subset); n2=sum(r[8] for r in subset); sy=sum(r[9] for r in subset)
                        row=next(r for r in actual if r['view']=='Environment' and r['model']==data['names'][mi] and r['group']==['Void','Wall','Filament','Cluster'][env])
                        assert math.isclose(row['n'],n,rel_tol=1e-12,abs_tol=1e-12)
                        if n:
                            assert math.isclose(row['mean'],sy/n,rel_tol=1e-12,abs_tol=1e-12)
                            assert math.isclose(row['neff'],n*n/n2,rel_tol=1e-12,abs_tol=1e-12)
                checks.append('independent Python/JavaScript mean, weight and effective-count parity for eight displayed groups')
                for prop in map(str,range(7)):
                    page.select_option('#prop',prop)
                    assert page.evaluate('window.__propertyState.count>0')
                page.select_option('#mass','2');page.select_option('#cap','0')
                assert page.evaluate('window.__propertyState.shown.filter(r=>r.mean!==null).every(r=>r.n>0 && r.neff>0)')
                with page.expect_download() as download:
                    page.click('#download')
                download.value.save_as(str(ROOT/'QA_filter_export.csv'))
                page.select_option('#z','3')
                assert page.evaluate('window.__propertyState.count===0')
                page.click('#reset');page.select_option('#prop','6')
                checks.append('properties: seven properties, intersected filters, CSV export and empty redshift selection')
            page.screenshot(path=str(ROOT/(name+'_wide.png')), full_page=True)
            page.set_viewport_size({'width':390,'height':844})
            assert not page.evaluate('document.documentElement.scrollWidth > innerWidth')
            page.screenshot(path=str(ROOT/(name+'_narrow.png')), full_page=True)
            page.set_viewport_size({'width':1200,'height':1000})
        assert not errors, errors
        assert not external, external
        browser.close()
    result={'passed':True,'job':os.environ['SLURM_JOB_ID'],'checks':checks,'errors':errors,'external_requests':external}
    (ROOT/'BROWSER_QA.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
finally:
    server.shutdown()
