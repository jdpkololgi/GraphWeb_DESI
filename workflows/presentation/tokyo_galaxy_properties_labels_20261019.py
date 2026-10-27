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
ROOT=Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010');OUT=ROOT/'notation_corrections_20261019/F_galaxy_properties';OUT.mkdir(exist_ok=True)
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
record=json.loads((ROOT/'talk_wall_overlay_extended_20261017/statistics.json').read_text());stats=record['statistics']
retained_n=record['retained_n']
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
right.set(ylim=(-11.25,-10.65),yticks=[-11.2,-11.0,-10.8]);right.set_ylabel('Mean log₁₀(sSFR / yr⁻¹)',fontsize=18,labelpad=11);right.set_xlabel('Inferred tidal trace\n'+r'$\mathrm{Tr}(T)=\lambda_1+\lambda_2+\lambda_3\propto\delta$',fontsize=14,labelpad=7)
f.text(.78,.785,'Same galaxies · matched stellar mass and redshift',ha='center',fontsize=12,color=G)
f.text(.2725,.209,'0.20 ≤ z < 0.30',ha='center',fontsize=18)
f.text(.78,.155,'Line: trace bins · diamonds: class averages',ha='center',fontsize=12,color=G)
f.text(.52,.113,'Classes overlap in trace; diamonds sit at each class’s mean trace.',ha='center',fontsize=15,color=G)
f.text(.52,.065,f'Right: {retained_n:,}-galaxy wedge · 16–84% spatial bootstrap · DESI calibration unverified',ha='center',fontsize=14,color=G)
# Match original output resolution despite original A07 figure dimensions.
for ext in ('png','pdf','svg'):f.savefig(OUT/('A07_continuous_and_classes.'+ext),dpi=288 if ext=='png' else 144,facecolor='black')
plt.close(f)
shutil.copy2(__file__,OUT/Path(__file__).name)
shutil.copy2(ROOT/'talk_wall_overlay_extended_20261017/statistics.json',OUT/'statistics.json')
