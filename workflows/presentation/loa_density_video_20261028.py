"""CPU-only, mask-preserving volume-raycast tour of existing shared field draws.

No inference or density smoothing. Full-cap views pool 2^3 native cells;
close-ups and slices use native cells. Interpolation requires all eight
supported corners. All rendering settings and source hashes are recorded.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import subprocess
import time

import h5py
import numpy as np
from numba import njit, prange
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import uniform_filter

ROOT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1/production_v1/NGC_seed17')
SOURCE = ROOT/'unet_x_17_FIELD.h5'
CELL = 6.766
WIDTH, HEIGHT = 960, 540
FPS, SECONDS = 24, 80

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

@njit
def colour(v):
    t=min(1.,max(0.,(math.log2(max(v,.01))+1.)/4.))
    # Continuous fixed log-density colour mapping: blue / violet / pink / ivory.
    knots=((.12,.20,.48),(.43,.18,.55),(.93,.23,.54),(1.,.57,.55),(1.,.94,.80))
    a=min(3,int(t*4)); f=t*4-a
    return (knots[a][0]*(1-f)+knots[a+1][0]*f,
            knots[a][1]*(1-f)+knots[a+1][1]*f,
            knots[a][2]*(1-f)+knots[a+1][2]*f)

@njit
def sample(a,x,y,z):
    i,j,k=int(math.floor(x)),int(math.floor(y)),int(math.floor(z))
    if i<0 or j<0 or k<0 or i+1>=a.shape[0] or j+1>=a.shape[1] or k+1>=a.shape[2]:return -1.
    value=0.
    for dx in range(2):
        for dy in range(2):
            for dz in range(2):
                v=a[i+dx,j+dy,k+dz]
                if v<0:return -1.
                value+=v*((x-i) if dx else (1-x+i))*((y-j) if dy else (1-y+j))*((z-k) if dz else (1-z+k))
    return value

@njit(parallel=True)
def raycast(a,origin,cell,center,right,up,forward,span,cut,opacity,cliplo,cliphi,width,height):
    out=np.zeros((height,width,3),np.uint8)
    step=cell*.35
    for py in prange(height):
        for px in range(width):
            uv=(px+.5-width/2)*span/height
            vv=(height/2-py-.5)*span/height
            p=center+right*uv+up*vv
            lo=-1e9;hi=1e9
            for ax in range(3):
                if abs(forward[ax])<1e-9:
                    if p[ax]<cliplo[ax] or p[ax]>cliphi[ax]:hi=-1e10
                else:
                    t0=(cliplo[ax]-p[ax])/forward[ax];t1=(cliphi[ax]-p[ax])/forward[ax]
                    lo=max(lo,min(t0,t1));hi=min(hi,max(t0,t1))
            rgb=np.zeros(3);trans=1.
            t=lo
            while t<hi and trans>.005:
                dt=min(step,hi-t)
                q=(p+forward*(t+dt*.5)-origin)/cell
                v=sample(a,q[0],q[1],q[2])
                if v>cut:
                    strength=min(6.,max(0.,math.log2(v))**2)*min(1.,(v-cut)/.25)
                    alpha=1.-math.exp(-opacity*strength*dt/CELL)
                    c=colour(v)
                    for ch in range(3):rgb[ch]+=trans*alpha*c[ch]
                    trans*=1-alpha
                t+=dt
            for ch in range(3):out[py,px,ch]=min(255,int(rgb[ch]*255+.5))
    return out

@njit(parallel=True)
def slice_image(a,origin,cell,center,span,z,width,height):
    out=np.zeros((height,width,3),np.uint8)
    for py in prange(height):
        for px in range(width):
            x=center[0]+(px+.5-width/2)*span/height
            y=center[1]+(height/2-py-.5)*span/height
            # The loaded array has interpolation padding beyond Region A;
            # do not display that padding as part of the labelled 600-Mpc box.
            if abs(x-center[0])>300 or abs(y-center[1])>300:continue
            v=sample(a,(x-origin[0])/cell,(y-origin[1])/cell,(z-origin[2])/cell)
            if v>=0:
                c=colour(v)
                for ch in range(3):out[py,px,ch]=int(c[ch]*255)
    return out

def basis(angle):
    eye=np.array([math.cos(angle),math.sin(angle),.6]);eye/=np.linalg.norm(eye)
    forward=-eye;right=np.cross(forward,[0.,0.,1.]);right/=np.linalg.norm(right)
    return right,np.cross(right,forward),forward

def pooled(a):
    x,y,z=a.shape
    return a.reshape(x//2,2,y//2,2,z//2,2).mean((1,3,5))

def load():
    receipt=json.loads((ROOT/'COMPLETE.json').read_text())
    record=next(r for r in receipt['results'] if r['model']=='unet_x_17')
    digest=sha(SOURCE)
    assert receipt['passed'] and digest==record['field']['sha256']
    with h5py.File(SOURCE) as f:
        origin=np.asarray(f.attrs['origin_mpc_h'])+CELL/2
        assert float(f.attrs['cell_mpc_h'])==CELL and f['rho'].shape[0]==4
        valid=(f['support'][:]>=.5)&(f['coverage'][:]>0)
        # Pick a support-rich region nearest the supported geometry centroid,
        # independently of density or choice of draw.
        coords=np.array(np.where(valid)).T
        mid=np.median(coords,axis=0)
        support_fraction=uniform_filter(valid.astype('f4'),size=65,mode='constant')
        candidates=np.array(np.where(support_fraction>.98)).T
        if not len(candidates):raise ValueError('No well-supported native detail region')
        chosen=candidates[np.argmin(np.sum((candidates-mid)**2,axis=1))]
        radius=48;low=np.maximum(chosen-radius,0);high=np.minimum(chosen+radius+1,valid.shape)
        crop=tuple(slice(int(l),int(h)) for l,h in zip(low,high))
        detail=[];full=None
        pooled_valid=(pooled(f['support'][:])>=.5)&(pooled((f['coverage'][:]>0).astype('f4'))==1)
        for d in range(4):
            rho=f['rho'][d]
            detail.append(np.where(valid[crop],rho[crop],-1).astype('f4'))
            if d==0:full=np.where(pooled_valid,pooled(rho),-1).astype('f4')
        # Arithmetic mean in linear density, never a mean rendered image.
        detail.append(np.where(valid[crop],np.mean(np.stack(detail),axis=0),-1).astype('f4'))
    support_indices=np.array(np.where(pooled_valid))
    bounds=np.array([support_indices.min(1),support_indices.max(1)]).astype(float)*CELL*2+origin+CELL/2
    return dict(full=full,detail=detail,origin=origin,full_origin=origin+CELL/2,
                detail_origin=origin+low*CELL,roi=origin+chosen*CELL,bounds=bounds,
                sha256=digest,roi_support_fraction=float(support_fraction[tuple(chosen)]),
                chosen=chosen.tolist(),crop=[low.tolist(),high.tolist()])

def font(size):
    import matplotlib
    path=Path(matplotlib.get_data_path())/'fonts/ttf/DejaVuSans.ttf'
    return ImageFont.truetype(str(path),size)

def test():
    a=np.ones((8,8,8),dtype='f4')*2
    assert abs(sample(a,2.25,3.4,4.1)-2)<1e-6
    a[2,3,4]=-1;assert sample(a,2.25,3.4,4.1)==-1
    a[:]=2
    origin=np.zeros(3);center=np.array([3.5,3.5,3.5]);right=np.array([1.,0.,0.]);up=np.array([0.,1.,0.]);forward=np.array([0.,0.,1.])
    out=raycast(a,origin,1.,center,right,up,forward,2.,1.,.1,np.zeros(3),np.ones(3)*7,1,1)
    # Exact optical path length, including any final fractional ray step.
    expected=np.asarray(colour(2.))*(1-np.exp(-.1*7/CELL))*255
    assert np.max(np.abs(out[0,0]-expected))<1.1,(out,expected)
    assert not raycast(a,origin,1.,center,right,up,forward,2.,3.,.1,np.zeros(3),np.ones(3)*7,1,1).any()
    return dict(trilinear_constant=True,mask_rejects_invalid_corner=True,uniform_slab_alpha=True,threshold=True)

def state(t,data):
    start=math.radians(35);angle=start
    global_center=data['bounds'].mean(0);global_span=float(np.linalg.norm(np.diff(data['bounds'],axis=0))) *.8
    roi=data['roi'];roi_span=900.;cut=1.25;draw=0;opacity=.10
    center=global_center;span=global_span;cliplo=data['bounds'][0];cliphi=data['bounds'][1];native=False;z=None
    if t<18:
        chapter='01 / ORBIT';title='One draw. One fixed density-to-opacity mapping.'
        angle=start+2*np.pi*t/18
    elif t<30:
        chapter='02 / ZOOM';title='Zoom and crop to Region A: a support-selected 600 Mpc/h box.'
        q=(t-18)/12;q=q*q*(3-2*q);native=True
        center=global_center*(1-q)+roi*q;span=global_span*(1-q)+roi_span*q
        # Continuous clipping transition explicitly reported on screen.
        cliplo=data['bounds'][0]*(1-q)+(roi-300)*q
        cliphi=data['bounds'][1]*(1-q)+(roi+300)*q
    elif t<44:
        chapter='03 / DENSITY SWEEP';title='Fixed camera. Only the opacity cutoff changes.'
        native=True;center=roi;span=roi_span;cliplo=roi-300;cliphi=roi+300
        cut=.75*(3.2/.75)**(math.sin(np.pi*(t-30)/14)**2)
    elif t<56:
        chapter='04 / SLICES';title='Native xy slices reveal structures hidden in projection.'
        native=True;center=roi;span=620.;cliplo=roi-300;cliphi=roi+300
        z=roi[2]-240+480*(t-44)/12
    elif t<72:
        chapter='05 / DRAW COMPARISON';title='Same region, camera and transfer function. No morphing.'
        native=True;center=roi;span=roi_span;cliplo=roi-300;cliphi=roi+300;draw=min(3,int((t-56)/4))
    else:
        chapter='06 / ARITHMETIC MEAN';title='Mean of four density fields — not another posterior draw.'
        native=True;center=roi;span=roi_span;cliplo=roi-300;cliphi=roi+300;draw=4
    # Zoom uses the complete native first draw until the ROI fills the view.
    return dict(chapter=chapter,title=title,center=center,span=span,angle=angle,cut=cut,
                opacity=opacity,draw=draw,native=native,cliplo=cliplo,cliphi=cliphi,z=z)

def annotate(pixels,s,t,data,fonts):
    im=Image.new('RGB',(1280,720),(7,9,17));d=ImageDraw.Draw(im)
    im.paste(Image.fromarray(pixels),(24,115))
    pink=(255,102,173);white=(235,239,249);muted=(151,164,185)
    d.text((24,19),'LOA  /  MATTER FIELD',font=fonts[28],fill=white)
    d.text((24,57),'NGC  ·  Expanded U-Net  ·  seed 17  ·  saved shared-state inference',font=fonts[17],fill=muted)
    d.text((24,89),s['title'],font=fonts[16],fill=pink)
    d.text((1005,24),s['chapter'],font=fonts[15],fill=pink)
    drawlabel='MEAN OF 4' if s['draw']==4 else 'DRAW '+str(s['draw']+1)+' / 4'
    d.text((1005,120),drawlabel,font=fonts[25],fill=white)
    d.text((1005,159),'NATIVE SLICE' if s['z'] is not None else 'TRANSLUCENT VOLUME',font=fonts[13],fill=muted)
    d.text((1005,192),'Density  ρ / ρ̄',font=fonts[19],fill=white)
    for yy in range(180):
        v=2**(3-4*yy/179);c=tuple(int(x*255) for x in colour(v))
        d.line((1008,226+yy,1028,226+yy),fill=c)
    for v in [.5,1,2,4,8]:
        y=226+int((3-math.log2(v))/4*179)
        d.text((1040,y-8),str(v)+' × mean',font=fonts[15],fill=muted)
    if s['z'] is None:
        d.text((1005,430),f'Opacity cutoff: {s["cut"]:.2f}',font=fonts[16],fill=white)
        d.text((1005,457),'Colour scale fixed',font=fonts[14],fill=muted)
        d.text((1005,481),'Opacity strength: 0.10',font=fonts[14],fill=muted)
    else:
        d.text((1005,430),f'z = {s["z"]:.1f} Mpc/h',font=fonts[16],fill=white)
        d.text((1005,457),'x →    y ↑',font=fonts[16],fill=muted)
        d.text((1005,481),'Colour = local density',font=fonts[14],fill=muted)
    spacing=CELL if s['native'] else 2*CELL
    for y,text in [(526,f'Cells: {spacing:.3f} Mpc/h'),(553,'Black: transparent / masked'),(580,'Region A: 600 Mpc/h' if s['native'] else 'Full-cap overview'),(607,'No galaxy markers')]:
        d.text((1005,y),text,font=fonts[13],fill=muted)
    if s['z'] is None:
        right,up,_=basis(s['angle']);cx,cy=910,610
        for j,axis in enumerate('xyz'):
            end=(cx+right[j]*38,cy-up[j]*38)
            d.line((cx,cy,*end),fill=muted,width=2)
            d.text((end[0]+3,end[1]-8),axis,font=fonts[13],fill=white)
    # Orthographic projection: exact physical scale independent of depth.
    scale=100 if s['span']<1000 else 500
    length=scale/s['span']*HEIGHT
    d.line((52,622,52+length,622),fill=white,width=3)
    d.line((52,617,52,627),fill=white,width=2);d.line((52+length,617,52+length,627),fill=white,width=2)
    d.text((52,632),f'{scale} Mpc/h',font=fonts[14],fill=white)
    if s['native']:
        x,y,z=data['roi']
        d.text((250,637),f'Region A centre: ({x:.1f}, {y:.1f}, {z:.1f}) Mpc/h',font=fonts[13],fill=muted)
    else:
        right,up,_=basis(s['angle']);corners=[]
        for a in (-1,1):
            for b in (-1,1):
                for c in (-1,1):
                    q=data['roi']+np.array([a,b,c])*300-s['center']
                    corners.append((24+WIDTH/2+q@right/s['span']*HEIGHT,115+HEIGHT/2-q@up/s['span']*HEIGHT))
        for i in range(8):
            for bit in (1,2,4):
                if i<(i^bit):d.line((*corners[i],*corners[i^bit]),fill=(145,153,174),width=1)
        d.text((corners[7][0]+8,corners[7][1]),'Region A',font=fonts[14],fill=white)
    d.text((24,670),'Four draws illustrate variation, not convergence or calibrated real-DESI uncertainty.',font=fonts[15],fill=muted)
    d.text((1165,670),f'{int(t):02d} / 80 s',font=fonts[14],fill=muted)
    d.rectangle((24,704,1256,707),fill=(28,33,46));d.rectangle((24,704,24+1232*t/80,707),fill=pink)
    return im

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--smoke',action='store_true');args=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID') and socket.gethostname().startswith('nid')
    args.out.mkdir(exist_ok=False,parents=True)
    start=time.monotonic();tests=test();data=load()
    # Native full draw is needed during the continuous zoom, then discardable.
    with h5py.File(SOURCE) as f:
        v=(f['support'][:]>=.5)&(f['coverage'][:]>0)
        native_full=np.where(v,f['rho'][0],-1).astype('f4')
    fonts={n:font(n) for n in [13,14,15,16,17,19,25,28]}
    times=[0.,18.,29.9,35.,49.,56.,60.,64.,68.,75.] if args.smoke else np.arange(FPS*SECONDS)/FPS
    encoder=None;lastkey=None;lastpixels=None;timeline=[]
    if not args.smoke:
        encoder=subprocess.Popen(['ffmpeg','-hide_banner','-loglevel','error','-f','rawvideo','-pixel_format','rgb24','-video_size','1280x720','-framerate',str(FPS),'-i','pipe:0','-an','-c:v','libvpx-vp9','-threads','8','-cpu-used','4','-crf','24','-b:v','0','-pix_fmt','yuv420p',str(args.out/'loa-ngc-density-tour.webm')],stdin=subprocess.PIPE)
    for number,t in enumerate(times):
        s=state(float(t),data);right,up,forward=basis(s['angle'])
        if s['native']:
            a=native_full if t<30 else data['detail'][s['draw']]
            origin=data['origin'] if t<30 else data['detail_origin'];cell=CELL
        else:a=data['full'];origin=data['full_origin'];cell=2*CELL
        lo=np.maximum(s['cliplo'],origin);hi=np.minimum(s['cliphi'],origin+(np.array(a.shape)-1)*cell)
        key=(s['chapter'],s['draw'],s['angle'],tuple(s['center']),s['span'],s['cut'],s['z'],tuple(lo),tuple(hi))
        tick=time.monotonic()
        if key==lastkey:pixels=lastpixels
        elif s['z'] is not None:pixels=slice_image(a,origin,cell,s['center'],s['span'],s['z'],WIDTH,HEIGHT)
        else:pixels=raycast(a,origin,cell,s['center'],right,up,forward,s['span'],s['cut'],s['opacity'],lo,hi,WIDTH,HEIGHT)
        lastkey=key;lastpixels=pixels
        im=annotate(pixels,s,float(t),data,fonts)
        if args.smoke or number%120==0:
            im.save(args.out/f'frame-{number:04d}.png')
            print(json.dumps(dict(frame=number,time=float(t),render_seconds=time.monotonic()-tick,elapsed=time.monotonic()-start,nonblack=int(np.any(pixels>0,axis=2).sum()))),flush=True)
        if encoder:encoder.stdin.write(im.tobytes())
        timeline.append(dict(frame=number,time=float(t),chapter=s['chapter'],draw=s['draw'],cut=s['cut'],span=s['span']))
    if encoder:
        encoder.stdin.close();assert encoder.wait()==0
    report=dict(passed=True,smoke=args.smoke,job=os.environ['SLURM_JOB_ID'],node=socket.gethostname(),
        source_field=str(SOURCE),field_sha256=data['sha256'],renderer_sha256=sha(__file__),
        tests=tests,roi_center_mpc_h=data['roi'].tolist(),roi_selection='support-rich region nearest supported-mask centroid; independent of density',
        roi_support_fraction=data['roi_support_fraction'],native_cell_mpc_h=CELL,overview_cell_mpc_h=2*CELL,
        mask='support >= 0.5 and coverage > 0; trilinear interpolation rejects any missing corner',
        opacity='alpha = 1-exp(-0.10 * min(6,max(0,log2(rho))^2) * min(1,(rho-cut)/0.25) * ray_step/6.766), rho>cut; front-to-back compositing',
        duration_seconds=SECONDS,fps=FPS,dimensions=[1280,720],timeline=timeline,elapsed_seconds=time.monotonic()-start,
        posterior_qualified=False,video_sha256=sha(args.out/'loa-ngc-density-tour.webm') if encoder else None)
    (args.out/'COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='timeline'}),flush=True)

if __name__=='__main__':main()
