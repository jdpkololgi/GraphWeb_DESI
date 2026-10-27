"""Offline browser and payload reconciliation for the full-cap Page exports."""
import base64
import functools
import gzip
import http.server
import json
import os
from pathlib import Path
import sys
import threading

ROOT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1/page_fullcap_20261026_v2')
RUNTIME = Path('/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1')
sys.path.insert(0, str(RUNTIME/'playwright_runtime'))
os.environ['PLAYWRIGHT_BROWSERS_PATH'] = str(RUNTIME/'browser_runtime')
from playwright.sync_api import sync_playwright

class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass

assert os.environ.get('SLURM_JOB_ID')
server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(Quiet, directory=str(ROOT)))
threading.Thread(target=server.serve_forever, daemon=True).start()
url = f'http://127.0.0.1:{server.server_port}/'
errors, external, checks = [], [], []
try:
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={'width':1200,'height':1000})
        page.on('pageerror', lambda e: errors.append(str(e)))
        def guard(route):
            if route.request.url.startswith(url):
                route.continue_()
            else:
                external.append(route.request.url)
                route.abort()
        page.route('**/*', guard)
        for cap in ('NGC', 'SGC'):
            fragment = (ROOT/(cap+'.html')).read_text()
            payload = fragment.split("const payload='",1)[1].split("';",1)[0]
            data = json.loads(gzip.decompress(base64.b64decode(payload)))
            assert data['tileCount'] == len(data['tiles'])
            assert len(base64.b64decode(data['density'])) == len(data['tiles'])*96
            page.goto(url+cap+'.preview.html')
            frame = page.frames[1]
            frame.wait_for_function('window.__ready || window.__error', timeout=30000)
            assert frame.evaluate('window.__ready'), frame.evaluate('window.__error')
            snapshot = frame.evaluate('window.__loaState')
            assert snapshot['retained'] == data['retained'] and snapshot['tiles'] == data['tileCount']
            assert snapshot['pixels'] > 0
            assert frame.locator('#fractions tr').count() == 6
            for m in range(6):
                frame.select_option('#fit-a',str(m))
                assert frame.locator('#label-a').inner_text() == data['names'][m]
                counts = data['summary'][data['models'][m]]['class_counts']
                actual = frame.locator('#fractions tr').nth(m).locator('td').all_text_contents()
                assert actual == [f'{100*x/sum(counts):.2f}' for x in counts]
            for shell in range(4):
                frame.select_option('#shell',str(shell))
                assert frame.evaluate('window.__loaState.shown') > 0
            frame.select_option('#shell','-1')
            layer = str(data['tiles'][0][0][2])
            frame.select_option('#layer',layer)
            assert frame.evaluate('window.__loaState.shown') == sum(t[0][2]==int(layer) for t in data['tiles'])
            frame.select_option('#layer','-1')
            frame.select_option('#tile','0')
            assert f"{data['tiles'][0][3]:,} retained" in frame.locator('#selected').inner_text()
            frame.evaluate("document.getElementById('yaw').value=30;document.getElementById('yaw').dispatchEvent(new Event('change'))")
            assert frame.locator('#yaw-value').inner_text() == '30°'
            # External collaborator state must update controls without a save loop.
            frame.evaluate("window.dispatchEvent(new CustomEvent('openai:set_globals',{detail:{globals:{widgetState:{version:1,cap:window.__loaState.cap,privateContent:{a:2,b:4,yaw:0,tilt:25,shell:-1,layer:-1,tile:-1}}}}}))")
            assert frame.locator('#fit-a').input_value() == '2'
            for width in (1200, 390, 320):
                page.set_viewport_size({'width':width,'height':1000})
                page.wait_for_timeout(400)
                assert not frame.evaluate('document.documentElement.scrollWidth > innerWidth'), width
                assert frame.locator('canvas').first.bounding_box()['width'] > 100
                page.screenshot(path=str(ROOT/f'{cap}_{width}.png'), full_page=True)
            checks.append({'cap':cap,'tiles':data['tileCount'],'retained':data['retained'],'six_fits':True,'four_shells':True,'layers':True,'shared_state':True,'widths':[1200,390,320]})
            page.set_viewport_size({'width':1200,'height':1000})
        assert not errors, errors
        # The bundled preview shell requests optional tooltip/icon libraries.
        # They are blocked throughout this offline test; the fragment needs none.
        allowed_shell = ('https://unpkg.com/@floating-ui/', 'https://unpkg.com/lucide@')
        assert all(u.startswith(allowed_shell) for u in external), external
        browser.close()
    result={'passed':True,'job':os.environ['SLURM_JOB_ID'],'checks':checks,'errors':errors,'blocked_optional_preview_requests':external,'external_requests_completed':0}
    (ROOT/'BROWSER_QA.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
finally:
    server.shutdown()
