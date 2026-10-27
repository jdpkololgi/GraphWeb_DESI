#!/usr/bin/env python3
"""Talk slides from saved mock examples and matched cached DESI summaries."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,hashlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle,Patch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import PercentFormatter
ROOT=Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010');OUT=ROOT/'talk_wall_overlay_extended_20261017';OUT.mkdir(exist_ok=True)
SRC=Path('/global/u2/d/dkololgi/GraphWeb_DESI/docs/evidence/filament_mock_examples_20261017')
C=['#9CF6D5','#5B8FF9','#E85D9E','#F2C14E'];N=['Void','Wall','Filament','Cluster'];G='#B9C2D0';W='#F5F5F5';R=7/.6766
plt.rcParams.update({'font.family':'DejaVu Sans','text.color':W,'axes.labelcolor':W,'xtick.color':G,'ytick.color':G,'axes.edgecolor':'#596270','pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'black'})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(f,name,pdf):
 for ext in ('png','pdf','svg'):f.savefig(OUT/(name+'.'+ext),dpi=240 if ext=='png' else 120,facecolor='black')
 pdf.savefig(f,facecolor='black');plt.close(f)
def axes(f,rect):
 ax=f.add_axes(rect,facecolor='black')
 for side in ('top','right'):ax.spines[side].set_visible(False)
 return ax
m=json.loads((SRC/'EXAMPLES.json').read_text());(OUT/'data').mkdir(exist_ok=True)
with PdfPages(OUT/'main_wall.pdf') as mainpdf, PdfPages(OUT/'backup_slides.pdf') as backpdf:
 for row in m['examples']:
  key=row['key'];p=SRC/(key+'.npz');assert sha(p)==row['array_sha256'];d=dict(np.load(p));xyz=d['delta_Mpc'];central=d['central'];probs=d['probability'];truth=d['truth_class']
  np.testing.assert_array_equal(probs[central][0],row['probability']);np.testing.assert_allclose(d['truth'][central][0],row['truth_eigenvalues'],atol=1e-6,rtol=0)
  use=np.linalg.norm(xyz,axis=1)<=R;assert use.sum()==row['neighbors_R']+1
  pred=probs.argmax(1);f=plt.figure(figsize=(16,9),dpi=120,facecolor='black')
  main=key=='true_wall_line';err=key=='filament_error'
  title='A galaxy string can live in a tidal wall' if main else ('Backup · a wall with sheet-like galaxy geometry' if not err else 'Backup · a genuine wall–filament error')
  f.text(.055,.935,title,fontsize=27)
  f.text(.055,.872,'Galaxy geometry: where the tracers lie.    Tidal environment: the smoothed gravitational field.',fontsize=17,color=G)
  vertical=2 if key=='true_wall_sheet' else 1
  for j,label in enumerate(('Galaxy positions','Tidal truth at each galaxy','P12-A most probable class')):
   ax=axes(f,[.065+j*.31,.36,.25,.43]);q=xyz[use]
   if j==0:ax.scatter(q[:,0],q[:,vertical],s=44,c=G,lw=0)
   else:
    vals=truth[use] if j==1 else pred[use]
    for k in range(4):
     z=vals==k;ax.scatter(q[z,0],q[z,vertical],s=45,c=C[k],marker=['o','s','^','D'][k],lw=0)
   center_color=W if j==0 else C[int(truth[central][0]) if j==1 else int(pred[central][0])]
   ax.scatter([0],[0],marker='*',s=230,c=center_color,edgecolors='white',linewidths=1,zorder=5)
   ax.add_patch(Circle((0,0),R,fill=False,edgecolor='#657080',lw=1,ls=':'))
   ax.set(xlim=(-11.6,11.6),ylim=(-11.6,11.6),aspect='equal',xticks=[-10,0,10],yticks=[-10,0,10]);ax.tick_params(labelsize=12)
   ax.set_title(label,fontsize=19,pad=12);ax.set_xlabel('Δx [Mpc]',fontsize=13)
   if j==0:ax.set_ylabel('Δheight [Mpc]' if vertical==2 else 'Δy [Mpc]',fontsize=13)
   else:ax.tick_params(labelleft=False)
  handles=[Line2D([],[],marker=['o','s','^','D'][k],color=C[k],ls='',label=N[k],ms=7) for k in range(4)]
  handles.append(Line2D([],[],marker='*',color=W,ls='',label='Central galaxy',ms=10))
  f.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.292),ncol=5,frameon=False,fontsize=14)
  lam=np.array(row['truth_eigenvalues']);pw=row['probability'][1]
  if main:
   f.text(.055,.218,'True wall: one tidal eigenvalue exceeds 0.2',fontsize=20,color=C[1])
   f.text(.055,.165,'Its galaxy neighbourhood remains line-dominated in 3D.',fontsize=19)
   f.text(.945,.218,f'P(wall) = {pw:.1%}',ha='right',fontsize=20,color=C[1])
   f.text(.945,.165,'A string of galaxies need not be a tidal filament.',ha='right',fontsize=16,color=G)
  elif not err:
   f.text(.055,.218,'True wall · sheet-dominated full 3D neighbourhood',fontsize=20,color=C[1])
   f.text(.055,.165,'The original projected slab also passed the same string-like selection.',fontsize=18,color=G)
   f.text(.945,.218,f'P(wall) = {pw:.1%}',ha='right',fontsize=20,color=C[1])
  else:
   f.text(.055,.218,f'True filament: λ₂ = {lam[1]:.3f} > 0.2',fontsize=20,color=C[2])
   f.text(.055,.165,'Near the class boundary, P12-A ranks wall highest.',fontsize=18,color=G)
   f.text(.945,.218,f'P(wall) = {pw:.1%} · P(filament) = {row["probability"][2]:.1%}',ha='right',fontsize=18)
  f.text(.5,.095,f'Identical galaxies across panels · full sphere, R = 7 h⁻¹ Mpc · {row["neighbors_R"]} neighbours + centre',ha='center',fontsize=14,color=G)
  f.text(.5,.05,'ph006 development mock · native host truth at observed positions · per-galaxy labels, not a reconstructed field',ha='center',fontsize=13,color=G)
  save(f,'main_string_like_wall' if main else 'backup_'+key,mainpdf if main else backpdf)
  shutil.copy2(p,OUT/'data'/p.name)
shutil.copy2(SRC/'EXAMPLES.json',OUT/'data/EXAMPLES.json')
# Matched continuous and discrete estimates from one verified wedge sample.
cache=ROOT/'density_checks_20261014/matched_galaxies.npz';cm=json.loads((cache.parent/'MANIFEST.json').read_text());assert sha(cache)==cm['files'][cache.name]
a=dict(np.load(cache));x=a['trace'];y=a['log_ssfr'];_,cell=np.unique(a['mass_z_cells'],return_inverse=True);_,reg=np.unique(a['region'],return_inverse=True);nc=cell.max()+1;nr=reg.max()+1
edges=np.sort(np.r_[np.quantile(x,np.linspace(0,1,9)),1.75]);groups={'continuous':np.searchsorted(edges[1:-1],x,side='right'),'classes':a['P_WEB'].argmax(1)};sizes={'continuous':9,'classes':4}
counts={k:np.bincount(g*nc+cell,minlength=sizes[k]*nc).reshape(sizes[k],nc) for k,g in groups.items()};keep=np.logical_and.reduce([(n>=20).all(0) for n in counts.values()]);mask=keep[cell];assert keep.any()
ref=np.bincount(cell[mask],minlength=nc)[keep].astype(float);ref/=ref.sum();mult=np.random.default_rng(202610171).multinomial(nr,np.ones(nr)/nr,size=128)
stats={};boots={};valids={};checks=[]
for key,g in groups.items():
 ng=sizes[key];ix=(g*nc+cell)*nr+reg
 def agg(v):return np.bincount(ix,weights=v,minlength=ng*nc*nr).reshape(ng,nc,nr)[:,keep,:]
 n=agg(np.ones(len(x)));sy=agg(y);sx=agg(x);total=n.sum(2);mx=(sx.sum(2)/total)@ref;my=(sy.sum(2)/total)@ref
 bn=np.einsum('gcr,kr->kgc',n,mult);by=np.einsum('gcr,kr->kgc',sy,mult);valids[key]=(bn>0).all((1,2));boots[key]=(by,bn)
 for j in range(ng):
  selected=mask&(g==j);lookup=np.zeros(nc);lookup[keep]=ref/total[j];w=lookup[cell[selected]]
  np.testing.assert_allclose(np.bincount(cell[selected],weights=w,minlength=nc)[keep],ref,atol=1e-13)
  np.testing.assert_allclose([w@x[selected],w@y[selected]],[mx[j],my[j]],atol=1e-12)
  np.testing.assert_allclose(w@cell[selected],ref@np.flatnonzero(keep),atol=1e-12)
 stats[key]={'x':mx.tolist(),'mean':my.tolist(),'counts':np.bincount(g[mask],minlength=ng).tolist()}
valid=valids['continuous']&valids['classes'];assert valid.sum()>=100 # Require at least 100 usable replicates; report actual count without replacement.
for k,(by,bn) in boots.items():stats[k]['interval']=np.quantile((by[valid]/bn[valid])@ref,[.16,.84],axis=0).tolist()
orig=json.loads((ROOT/'a07_combined/plotted_statistics.json').read_text());theme=json.loads((ROOT/'theme.json').read_text());colors=[theme['class_palette'][k] for k in ['void','sheet','filament','knot']]
f=plt.figure(figsize=(1920/144,1080/144),dpi=144,facecolor='black');left=axes(f,[.08,.31,.385,.43]);right=axes(f,[.605,.31,.35,.43])
for ax in (left,right):
 ax.tick_params(axis='x',labelsize=17,pad=11,length=0,colors='#C8CDD4');ax.tick_params(axis='y',labelsize=16,pad=7,colors='#C8CDD4')
 for side in ('left','bottom'):ax.spines[side].set_color('#77818E')
left.set_xticks(range(4),orig['class_order']);f.text(.2725,.858,'Environment fractions',ha='center',fontsize=22);f.text(.78,.858,'sSFR and environment',ha='center',fontsize=22)
series=list(orig['class_fractions']);alphas=[.42,.70,1.];width=.23
for j,name in enumerate(series):
 values=orig['class_fractions'][name];xs=np.arange(4)+(j-1)*width;left.bar(xs,values,width,color=colors,alpha=alphas[j],edgecolor='#E2E7ED',linewidth=.6)
 for xx,v in zip(xs,values):left.text(xx,v+.01,f'{v*100:.0f}',ha='center',fontsize=12,color='#E2E7ED')
left.set(ylim=(0,.53),xlim=(-.6,3.6),yticks=[0,.1,.2,.3,.4,.5]);left.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));left.set_ylabel('Fraction of galaxies',fontsize=18,labelpad=9)
left.legend(handles=[Patch(facecolor='#E2E7ED',alpha=al,label=s) for s,al in zip(series,alphas)],loc='lower center',bbox_to_anchor=(.5,1.015),ncol=3,fontsize=12,frameon=False,handlelength=1,columnspacing=.8)
r=stats['continuous'];xx=np.array(r['x']);yy=np.array(r['mean']);ci=np.array(r['interval']);right.plot(xx,yy,'o-',color='#DEE4EC',ms=4,lw=1.5);right.fill_between(xx,*ci,color='#DEE4EC',alpha=.17,lw=0)
r=stats['classes'];ci=np.array(r['interval'])
for k in range(4):
 yy=r['mean'][k];right.errorbar(r['x'][k],yy,yerr=[[max(0,yy-ci[0,k])],[max(0,ci[1,k]-yy)]],fmt='D',color=colors[k],mec='black',mew=.8,ms=9,capsize=4,zorder=5)
 right.annotate(N[k],(r['x'][k],yy),xytext=(0,13 if k%2==0 else -21),textcoords='offset points',ha='center',fontsize=12,color=colors[k])
right.set(ylim=(-11.25,-10.65),yticks=[-11.2,-11.0,-10.8]);right.set_ylabel('Mean log₁₀(sSFR / yr⁻¹)',fontsize=18,labelpad=11);right.set_xlabel('Inferred tidal trace',fontsize=16,labelpad=9)
f.text(.78,.785,'Same galaxies · matched stellar mass and redshift',ha='center',fontsize=12,color=G)
f.text(.2725,.209,'0.20 ≤ z < 0.30',ha='center',fontsize=18)
f.text(.78,.183,'Line: trace bins · diamonds: class averages',ha='center',fontsize=12,color=G)
f.text(.52,.148,'Classes overlap in trace; diamonds sit at each class’s mean trace.',ha='center',fontsize=15,color=G)
f.text(.52,.105,f'Right: {mask.sum():,}-galaxy wedge · 16–84% spatial bootstrap · DESI calibration unverified',ha='center',fontsize=14,color=G)
# Match original output resolution despite original A07 figure dimensions.
for ext in ('png','pdf','svg'):f.savefig(OUT/('A07_continuous_and_classes.'+ext),dpi=288 if ext=='png' else 144,facecolor='black')
plt.close(f)
record=dict(statistics=stats,input_n=len(x),retained_n=int(mask.sum()),common_cells=int(keep.sum()),valid_paired_bootstraps=int(valid.sum()),trace_edges=edges.tolist(),source_cache_sha256=sha(cache),method='Original eight trace quantile bins with upper bin split at trace=1.75, and four P_WEB argmax classes; >=20 rows in each group of both partitions per .1dex mass/.01z cell. Common rows, fixed pooled reference, same spatial resamples; intervals retain only jointly valid draws. x and y are independently standardized to the identical reference mix.',limitations='Descriptive cached-wedge results, not full-survey sSFR. Hard classes discard posterior uncertainty; no SED/selection/model systematics. Trace does not uniquely determine class. Class summaries overlap and are not samples of the bin curve.',checks='Exact joint-cell balance, independent direct weighted x/y parity, cell-only outcome equality; source hashes and original central probability/truth/count parity.')
(OUT/'statistics.json').write_text(json.dumps(record,indent=2)+'\n');np.savez_compressed(OUT/'data/overlay_rows.npz',TARGETID=a['TARGETID'],trace=x,log_ssfr=y,mass=a['mass'],z=a['z'],class_index=groups['classes'],trace_bin=groups['continuous'],cell=cell,region=reg,retained=mask)
shutil.copy2(cache.parent/'provenance.json',OUT/'data/wedge_provenance.json');shutil.copy2(ROOT/'a07_combined/plotted_statistics.json',OUT/'data/original_fullsurvey_statistics.json')
(OUT/'scripts').mkdir(exist_ok=True);shutil.copy2(__file__,OUT/'scripts'/Path(__file__).name)
print(json.dumps(record,indent=2))
