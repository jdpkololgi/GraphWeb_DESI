"""Offline headless interaction/render checks for the generated local viewer."""
import argparse
import json
import os
import sys
from pathlib import Path


def main(a):
    if 'SLURM_JOB_ID' not in os.environ:raise RuntimeError('Browser QA requires compute allocation')
    if a.playwright_path:sys.path.insert(0,str(a.playwright_path))
    from playwright.sync_api import sync_playwright
    a.out.mkdir(parents=True,exist_ok=True)
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True,args=['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
        page=browser.new_page(viewport={'width':1440,'height':1080},device_scale_factor=1)
        errors=[];requests=[];checks=[]
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.on('request',lambda r:requests.append(r.url))
        page.goto(a.html.resolve().as_uri(),wait_until='load',timeout=120000)
        page.wait_for_function('window.viewerReady === true',timeout=120000)
        def ready():
            page.wait_for_function('!busy',timeout=120000)
            assert 'Could not render' not in page.locator('#status').inner_text()
            assert page.locator('#plot canvas').count()>0
        ready();page.screenshot(path=str(a.out/'NGC_detail_desktop.png'),full_page=True)
        checks.append({'view':'NGC detail','traces':page.evaluate('document.getElementById("plot").data.map(t=>t.type)')})
        page.select_option('#draw','16');ready()
        page.uncheck('#galaxies');ready();assert page.evaluate('document.getElementById("plot").data.length')==1
        assert page.evaluate('document.getElementById("plot").data[0].i.length')>0
        assert page.evaluate('document.getElementById("plot")._fullLayout.scene.xaxis.range[1]-document.getElementById("plot")._fullLayout.scene.xaxis.range[0]')>100
        page.screenshot(path=str(a.out/'NGC_density_only.png'),full_page=True)
        page.check('#galaxies');ready()
        page.select_option('#view','SGC detail');ready()
        page.screenshot(path=str(a.out/'SGC_detail_desktop.png'),full_page=True)
        checks.append({'view':'SGC detail','traces':page.evaluate('document.getElementById("plot").data.map(t=>t.type)')})
        page.select_option('#view','NGC overview');ready()
        assert page.locator('#draw').is_disabled()
        page.screenshot(path=str(a.out/'NGC_overview_desktop.png'),full_page=True)
        for name in ['NGC atlas','SGC atlas']:
            if page.locator('#view option').filter(has_text=name).count():
                page.select_option('#view',name);ready();page.select_option('#draw','8');ready()
                page.screenshot(path=str(a.out/(name.replace(' ','_')+'.png')),full_page=True)
                checks.append({'view':name,'traces':page.evaluate('document.getElementById("plot").data.map(t=>t.type)')})
        page.set_viewport_size({'width':430,'height':932});page.select_option('#view','SGC detail');ready()
        page.screenshot(path=str(a.out/'SGC_detail_mobile.png'),full_page=True)
        overflow=page.evaluate('document.documentElement.scrollWidth > innerWidth')
        remote=[r for r in requests if r.startswith(('http:','https:'))]
        report={'passed':not errors and not remote and not overflow,'page_errors':errors,'remote_requests':remote,'mobile_horizontal_overflow':overflow,'checks':checks}
        (a.out/'BROWSER_CHECK.json').write_text(json.dumps(report,indent=2)+'\n')
        browser.close();print(json.dumps(report),flush=True)
        if not report['passed']:raise RuntimeError('Viewer browser check failed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--html',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--playwright-path',type=Path);main(p.parse_args())
