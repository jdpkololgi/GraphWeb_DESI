#!/usr/bin/env python3
"""Build an animated line-of-sight skewer through the DESI wedge driven by REAL NPE
posteriors (iteration 1). Mirrors the Desktop mockup layout (galaxy strip + width
ribbon + class strip + 3 eigenvalue posterior densities + class-probability bar +
scrubber), but every curve is computed from the model's posterior samples.

Skewer: galaxies are projected onto the wedge's long axis (PC1 of comoving XYZ);
we pan a position window along it. At each step the galaxies inside the window
(within a transverse tube) pool their posterior λ samples into KDEs for λ1/λ2/λ3,
and the class bar/width come from those galaxies. A galaxy-distribution panel
(PC1 vs PC2, coloured by inferred class) shows the structure being traversed with a
moving cursor band.

Reads desi_wedge_flowjax_preds.npz (needs lambda_samples_subset + sample_subset_idx).
Writes skewer_posterior_animation_real.html (self-contained, dark theme).
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.stats import gaussian_kde
from astropy.cosmology import Planck18 as cosmo

CC = {"void": "#A1FCDD", "wall": "#4E84F7", "filament": "#EB336F", "cluster": "#F5C144"}
ORDER = ["void", "wall", "filament", "cluster"]
EC = ["#3A86FF", "#FF006E", "#D62828"]  # λ1, λ2, λ3 curve colours
TH = 0.2


def _mixture(means, stds, xg):
    """Per-bin eigenvalue density = mixture over galaxies of N(mean_g, std_g).
    Uses every galaxy in the bin (mean+std available for all), so it is smooth and
    correctly widens where the posteriors widen."""
    stds = np.maximum(np.asarray(stds), 1e-3)
    z = (xg[None, :] - means[:, None]) / stds[:, None]
    return (np.exp(-0.5 * z * z) / (stds[:, None] * np.sqrt(2 * np.pi))).mean(0)


def main(args):
    d = np.load(args.preds_npz)
    ra, dec, z = d["ra"], d["dec"], d["z"]
    classprob, hard = d["classprob"], d["hard_class"]
    lam_mean, lam_std = d["lambda_mean"], d["lambda_std"]
    width = lam_std.mean(1)

    dist = cosmo.comoving_distance(z).value
    rar, decr = np.deg2rad(ra), np.deg2rad(dec)

    # ---- radial line-of-sight pencil beam at the richest inferred CLUSTER's RA/Dec ----
    # (a thin 3D line through a shell-like wedge would cut the empty interior; the
    # physical "skewer" is a fixed-RA/Dec sightline binned in comoving distance, so the
    # cluster — localized in distance — appears as a spike, with void/wall around it.)
    clu = classprob[:, 3]
    P = np.vstack([dist * np.cos(decr) * np.cos(rar), dist * np.cos(decr) * np.sin(rar), dist * np.sin(decr)]).T
    near = np.linalg.norm(P - P[np.argmax(clu)], axis=1) < args.anchor_r
    ra_c = float(np.average(ra[near], weights=clu[near])); dec_c = float(np.average(dec[near], weights=clu[near]))
    rac, dcc = np.deg2rad(ra_c), np.deg2rad(dec_c)
    ang = np.degrees(np.arccos(np.clip(np.sin(decr) * np.sin(dcc) + np.cos(decr) * np.cos(dcc) * np.cos(rar - rac), -1, 1)))
    tube = ang < args.theta_deg                       # pencil beam
    t_mpc = np.sign(dec - dec_c) * np.deg2rad(ang) * dist   # signed transverse offset (Mpc) for the panel
    dmin, dmax = np.percentile(dist[tube], [1, 99])
    print(f"LOS at RA={ra_c:.2f} Dec={dec_c:.2f}; beam θ<{args.theta_deg}° -> {int(tube.sum())} galaxies; "
          f"distance {dmin:.0f}-{dmax:.0f} Mpc", flush=True)

    centers = np.linspace(dmin, dmax, args.n_frames)
    binw = (dmax - dmin) / args.n_frames * 1.6
    xg = np.linspace(args.lmin, args.lmax, args.grid)
    frames = []
    for c in centers:
        m = tube & (np.abs(dist - c) < binw)
        n = int(m.sum())
        if n >= args.min_gal:
            kdes = [_mixture(lam_mean[m, k], lam_std[m, k], xg).tolist() for k in range(3)]
            probs = {c2: float(classprob[m, ORDER.index(c2)].mean()) for c2 in ORDER}
            wid = float(width[m].mean())
        else:
            kdes = [[0.0] * args.grid] * 3; probs = {c2: 0.0 for c2 in ORDER}; wid = 0.0
        frames.append({"s": float(c), "n": n, "kde": kdes, "probs": probs, "width": wid})

    rng = np.random.default_rng(0)
    gi = np.where(tube)[0]
    gi = rng.choice(gi, min(len(gi), args.scatter_points), replace=False)
    gal = {"s": dist[gi].round(2).tolist(), "t": t_mpc[gi].round(2).tolist(), "c": [ORDER[h] for h in hard[gi]]}

    payload = {"frames": frames, "xg": xg.tolist(), "gal": gal, "th": TH, "cc": CC, "ec": EC,
               "order": ORDER, "smin": float(centers[0]), "smax": float(centers[-1]),
               "tmin": float(np.percentile(t_mpc[gi], 2)), "tmax": float(np.percentile(t_mpc[gi], 98))}

    html = _HTML.replace("__DATA__", json.dumps(payload))
    out = Path(args.out)
    out.write_text(html, encoding="utf-8")
    print(f"frames={len(frames)} median galaxies/bin={int(np.median([f['n'] for f in frames]))}")
    print("Saved:", out)


_HTML = r"""<!doctype html><html><head><meta charset="utf-8"><style>
body{background:#000;color:#F2F2F2;font-family:'IBM Plex Sans',system-ui,sans-serif;margin:0;padding:16px}
canvas{display:block;width:100%}.row{display:flex;align-items:center;gap:12px;margin-top:12px}
button{background:#1a1a1a;color:#F2F2F2;border:1px solid #555;border-radius:6px;padding:4px 12px;cursor:pointer}
input[type=range]{flex:1}.ro{font-size:13px;color:#a6a69f;min-width:120px}h3{font-weight:600;margin:0 0 6px}
.leg{font-size:12px;color:#a6a69f;margin-top:6px}.sw{display:inline-block;width:11px;height:11px;border-radius:2px;margin:0 5px 0 12px;vertical-align:-1px}
</style></head><body>
<h3>Line-of-sight skewer through the DESI wedge — real NPE posteriors</h3>
<canvas id="gal" height="150"></canvas>
<canvas id="strip" height="86"></canvas>
<canvas id="post" height="280"></canvas>
<div id="bar" style="display:flex;height:22px;border-radius:6px;overflow:hidden;border:0.5px solid #555;margin-top:12px"></div>
<div class="leg" id="leg"></div>
<div class="row"><button id="play">pause</button><input type="range" id="scrub" min="0" max="100" value="0"><span class="ro" id="pos"></span></div>
<script>
var D=__DATA__;var F=D.frames,XG=D.xg,TH=D.th,CC=D.cc,EC=D.ec,ORD=D.order;
var dpr=Math.max(1,Math.min(2,window.devicePixelRatio||1));
function ctx(id,h){var c=document.getElementById(id);c.width=c.clientWidth*dpr;c.height=h*dpr;var x=c.getContext('2d');x.setTransform(dpr,0,0,dpr,0,0);x._w=c.clientWidth;x._h=h;return x;}
var cg=ctx('gal',150),cs=ctx('strip',86),cp=ctx('post',280);
var L=46,R=16;
function sx(W,s){return L+(s-D.smin)/(D.smax-D.smin)*(W-L-R);}
function blend(p){var r=0,g=0,b=0;ORD.forEach(function(k){var h=CC[k],pp=Math.max(0,p[k]);r+=pp*parseInt(h.substr(1,2),16);g+=pp*parseInt(h.substr(3,2),16);b+=pp*parseInt(h.substr(5,2),16);});return 'rgb('+(r|0)+','+(g|0)+','+(b|0)+')';}
function drawGal(i){var W=cg._w,H=cg._h;cg.clearRect(0,0,W,H);
  var g=D.gal;for(var j=0;j<g.s.length;j++){var X=sx(W,g.s[j]);var Y=10+(g.t[j]-D.tmin)/(D.tmax-D.tmin)*(H-20);cg.fillStyle=CC[g.c[j]];cg.globalAlpha=0.5;cg.fillRect(X,Y,1.4,1.4);}cg.globalAlpha=1;
  var cx=sx(W,F[i].s);cg.strokeStyle='#f2f2f2';cg.lineWidth=1.5;cg.beginPath();cg.moveTo(cx,2);cg.lineTo(cx,H-2);cg.stroke();
  cg.fillStyle='#a6a69f';cg.font='11px sans-serif';cg.fillText('galaxies along skewer (colour = inferred class)',L,12);}
function drawStrip(i){var W=cs._w,H=cs._h;cs.clearRect(0,0,W,H);
  var sigMin=0.03,sigMax=0.12,baseR=30;cs.beginPath();cs.moveTo(L,baseR);
  for(var k=0;k<F.length;k++){var X=sx(W,F[k].s);var h=Math.max(0,Math.min(20,(F[k].width-sigMin)/(sigMax-sigMin)*20));cs.lineTo(X,baseR-h);}
  cs.lineTo(sx(W,F[F.length-1].s),baseR);cs.closePath();cs.fillStyle='rgba(127,119,221,.30)';cs.fill();
  cs.fillStyle='#a6a69f';cs.font='11px sans-serif';cs.fillText('posterior width',L,12);
  var top=40,bh=26;for(var m=0;m<F.length;m++){var X2=sx(W,F[m].s);var X3=(m<F.length-1)?sx(W,F[m+1].s):X2+2;cs.fillStyle=blend(F[m].probs);cs.fillRect(X2,top,Math.max(1.5,X3-X2+0.6),bh);}
  cs.strokeStyle='rgba(255,255,255,.22)';cs.strokeRect(L,top,sx(W,F[F.length-1].s)-L,bh);
  var cx=sx(W,F[i].s);cs.strokeStyle='#f2f2f2';cs.lineWidth=1.5;cs.beginPath();cs.moveTo(cx,top-4);cs.lineTo(cx,top+bh+4);cs.stroke();}
function drawPost(i){var W=cp._w,H=cp._h,PT=14,PB=30,pH=H-PT-PB;cp.clearRect(0,0,W,H);
  var lmin=XG[0],lmax=XG[XG.length-1];function px(l){return L+(l-lmin)/(lmax-lmin)*(W-L-R);}
  var dmax=0;F[i].kde.forEach(function(c){c.forEach(function(v){if(v>dmax)dmax=v;});});dmax=Math.max(dmax,0.5);
  function py(v){return PT+pH-v/dmax*pH;}
  cp.fillStyle='rgba(255,255,255,.05)';cp.fillRect(px(TH),PT,px(lmax)-px(TH),pH);
  cp.strokeStyle='rgba(255,255,255,.22)';cp.beginPath();cp.moveTo(L,PT+pH);cp.lineTo(W-R,PT+pH);cp.stroke();
  cp.fillStyle='#a6a69f';cp.font='12px sans-serif';cp.textAlign='center';
  [-1,0,1,2,3,4].forEach(function(t){if(t<lmin||t>lmax)return;var X=px(t);cp.fillText(t,X,PT+pH+16);});
  cp.fillText('eigenvalue λ',(L+W-R)/2,H-2);
  cp.setLineDash([5,4]);cp.strokeStyle='#9a9a93';cp.beginPath();cp.moveTo(px(TH),PT);cp.lineTo(px(TH),PT+pH);cp.stroke();cp.setLineDash([]);
  cp.fillStyle='#F2F2F2';cp.textAlign='left';cp.fillText('λ_th=0.2',px(TH)+5,PT+12);
  var dash=[[],[7,4],[2,4]],nm=['λ₁','λ₂','λ₃'];
  for(var k=0;k<3;k++){var c=F[i].kde[k];cp.beginPath();for(var j=0;j<XG.length;j++){var X=px(XG[j]),Yv=py(c[j]);if(j===0)cp.moveTo(X,Yv);else cp.lineTo(X,Yv);}
    cp.lineTo(px(lmax),PT+pH);cp.lineTo(px(lmin),PT+pH);cp.closePath();cp.globalAlpha=0.13;cp.fillStyle=EC[k];cp.fill();cp.globalAlpha=1;
    cp.beginPath();cp.setLineDash(dash[k]);for(var j2=0;j2<XG.length;j2++){var X2=px(XG[j2]),Y2=py(c[j2]);if(j2===0)cp.moveTo(X2,Y2);else cp.lineTo(X2,Y2);}cp.strokeStyle=EC[k];cp.lineWidth=1.8;cp.stroke();cp.setLineDash([]);
    cp.fillStyle='#F2F2F2';cp.fillText(nm[k],px(lmax)-40,PT+12+k*15);}}
function bar(i){var p=F[i].probs,el=document.getElementById('bar');el.innerHTML='';
  ORD.forEach(function(k){var s=document.createElement('div');s.style.background=CC[k];s.style.width=(Math.max(0,p[k])*100)+'%';el.appendChild(s);});
  document.getElementById('leg').innerHTML=ORD.map(function(k){return '<span class="sw" style="background:'+CC[k]+'"></span>'+k+' '+Math.round(Math.max(0,p[k])*100)+'%';}).join('');}
var i=0,playing=true,last=0;
function render(){drawGal(i);drawStrip(i);drawPost(i);bar(i);document.getElementById('scrub').value=Math.round(i/(F.length-1)*100);document.getElementById('pos').textContent='s = '+Math.round(F[i].s)+' Mpc  ('+F[i].n+' gal)';}
function loop(ts){if(playing){if(ts-last>120){i=(i+1)%F.length;last=ts;}}render();requestAnimationFrame(loop);}
document.getElementById('play').onclick=function(){playing=!playing;this.textContent=playing?'pause':'play';};
document.getElementById('scrub').oninput=function(){playing=false;document.getElementById('play').textContent='play';i=Math.round(this.value/100*(F.length-1));render();};
window.addEventListener('resize',function(){cg=ctx('gal',150);cs=ctx('strip',86);cp=ctx('post',280);render();});
requestAnimationFrame(loop);
</script></body></html>"""


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--preds-npz", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--n-frames", type=int, default=60)
    ap.add_argument("--theta-deg", type=float, default=1.5, help="pencil-beam angular radius (deg) around the sightline")
    ap.add_argument("--anchor-r", type=float, default=15.0, help="radius (Mpc) to locate the cluster centroid")
    ap.add_argument("--min-gal", type=int, default=8)
    ap.add_argument("--grid", type=int, default=120)
    ap.add_argument("--lmin", type=float, default=-1.0); ap.add_argument("--lmax", type=float, default=2.0)
    ap.add_argument("--scatter-points", type=int, default=8000)
    a = ap.parse_args()
    if a.out is None:
        a.out = str(Path(a.preds_npz).parent / "skewer_posterior_animation_real.html")
    main(a)
