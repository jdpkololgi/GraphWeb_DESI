"""Decode, inspect and browser-play the completed scientific movie on CPU."""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
from PIL import Image

def main():
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);args=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    root=args.directory;video=root/'loa-ngc-density-tour.webm'
    report=json.loads((root/'COMPLETE.json').read_text())
    assert report['passed'] and not report['smoke']
    assert hashlib.sha256(video.read_bytes()).hexdigest()==report['video_sha256']
    probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-count_frames','-select_streams','v:0',
        '-show_entries','stream=codec_name,width,height,r_frame_rate,nb_read_frames:format=duration,size','-of','json',str(video)],text=True))
    stream=probe['streams'][0]
    assert (stream['width'],stream['height'])==(1280,720)
    assert stream['r_frame_rate']=='24/1' and int(stream['nb_read_frames'])==1920
    assert abs(float(probe['format']['duration'])-80)<.1
    assert len(report['timeline'])==1920
    comparisons=[r for r in report['timeline'] if 56<=r['time']<72]
    assert {r['draw'] for r in comparisons}=={0,1,2,3}
    assert len({(r['cut'],r['span']) for r in comparisons})==1
    assert {r['draw'] for r in report['timeline'] if r['time']>=72}=={4}
    # Decode sampled video frames to ensure nonblank pictures and actual draw
    # changes, excluding labels/sidebar from the comparison.
    sampled=[];frames={}
    for t in [5,25,37,50,58,62,66,70,76]:
        raw=subprocess.check_output(['ffmpeg','-v','error','-ss',str(t),'-i',str(video),'-frames:v','1',
            '-f','image2pipe','-c:v','png','pipe:1'])
        a=np.asarray(Image.open(io.BytesIO(raw)).convert('RGB'))
        assert a.shape==(720,1280,3)
        volume=a[145:600,200:850].astype(float)
        if t==50:
            # 600-Mpc box at span 620 occupies x=243..765 in the movie.
            # Interpolation-padding pixels farther left must remain masked.
            assert np.max(a[150:600,210:230])<5
        # The maximum-cutoff frame is deliberately sparse, not normalized to
        # match the brightness of lower thresholds.
        assert np.mean(volume)>(.02 if t==37 else 1.),(t,float(np.mean(volume)))
        frames[t]=volume
        sampled.append(dict(time=t,mean_rgb=float(volume.mean()),nonblack=int(np.any(volume>8,axis=2).sum())))
    differences=[float(np.mean(np.abs(frames[a]-frames[b]))) for a,b in [(58,62),(62,66),(66,70),(70,76)]]
    assert min(differences)>.01,differences
    base=Path('/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1')
    sys.path.insert(0,str(base/'playwright_runtime'));os.environ['PLAYWRIGHT_BROWSERS_PATH']=str(base/'browser_runtime')
    from playwright.sync_api import sync_playwright
    with sync_playwright() as pw:
        browser=pw.chromium.launch(headless=True,args=['--no-sandbox'])
        page=browser.new_page(viewport={'width':1280,'height':760})
        page.goto(video.as_uri(),timeout=60000)
        page.wait_for_function('document.querySelector("video")?.readyState>=2',timeout=60000)
        duration=page.evaluate('document.querySelector("video").duration')
        assert abs(duration-80)<.1
        page.evaluate('document.querySelector("video").currentTime=35')
        page.wait_for_function('!document.querySelector("video").seeking')
        page.evaluate('document.querySelector("video").play()')
        page.wait_for_timeout(1200)
        assert page.evaluate('document.querySelector("video").currentTime')>35.3
        page.evaluate('document.querySelector("video").pause()')
        page.screenshot(path=str(root/'playback-check.png'))
        browser.close()
    result=dict(passed=True,probe=probe,sampled_frames=sampled,draw_image_mean_absolute_differences=differences,
        browser_playback=True,sha256=report['video_sha256'],job=os.environ['SLURM_JOB_ID'])
    (root/'VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
