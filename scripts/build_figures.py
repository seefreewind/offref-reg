"""Draw only frozen data. No correlations, intervals or inferential tests are recomputed."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter, PercentFormatter
from PIL import Image
from theme import *
apply()
plt.rcParams.update({'svg.fonttype':'none','pdf.fonttype':42})
R=Path(__file__).resolve().parents[1];S=R/'source_data';O=R/'figures';QA=R/'qa'
B=json.loads((QA/'frozen_bundle.json').read_text())
h=B['harmonization'];ref=B['reference'];pc=B['placement_context'];gap=B['gap'];gn=B['gnrs']
metrics=pd.read_csv(S/'frozen_metrics.tsv',sep='\t',float_precision='round_trip')
panels=[];exports=[];fig_no=0
def data(name):return pd.read_csv(S/name,sep='\t',float_precision='round_trip')
def canvas(n,height):
    global fig_no
    fig_no=n;return plt.figure(figsize=(WIDTH_MM/25.4,height/25.4))
def panel(fig,rect,letter,title,off=False,source=None,claim=None):
    ax=fig.add_axes(rect)
    fig.text(rect[0]-.025,rect[1]+rect[3]+.018,letter,fontsize=9.5,weight='bold',va='bottom')
    ax.set_title(title,loc='left',pad=9,fontsize=8.5)
    if off:ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    panels.append({'figure':str(fig_no),'panel':letter,'title':title,'source':source or f'Figure{fig_no}{letter}_source.tsv','claim':claim or title,'rect':rect,'font_min_pt':7.2})
    return ax
def text(ax,x,y,s,size=7.5,**kw):return ax.text(x,y,s,fontsize=size,transform=ax.transAxes,**kw)
def arrow(ax,a,b,ls='-'):
    ax.annotate('',xy=b,xytext=a,xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.8,'color':BORZOI,'linestyle':ls})
def box(ax,x,y,w,ht,label,color=REFERENCE,ls='-',size=7.5):
    ax.add_patch(Rectangle((x,y),w,ht,transform=ax.transAxes,fc='white',ec=color,lw=.7,ls=ls))
    text(ax,x+w/2,y+ht/2,label,size,ha='center',va='center',color=color)
def allele(ax,y,alt,label):
    text(ax,0,y+.105,label,weight='bold')
    for x in [.17,.68]:ax.add_patch(Rectangle((x,y),.30,.055,fc=LIGHT,ec='none',transform=ax.transAxes))
    ax.plot([.47,.68],[y+.028,y+.028],lw=.7,c=BORZOI,transform=ax.transAxes)
    if alt:ax.add_patch(Rectangle((.47,y-.013),.21,.082,fc=INSERTION,ec='none',transform=ax.transAxes))
def zero(ax):ax.axhline(0,c=LIGHT,lw=.6,zorder=0);ax.axvline(0,c=LIGHT,lw=.6,zorder=0)
def scatter(ax,d,x,y,n,rho,color=INSERTION,xlab='AlphaGenome ΔRNA',ylab='Borzoi ΔRNA',annotate=True):
    valid=np.isfinite(d[x])&np.isfinite(d[y]);assert int(valid.sum())==n
    ax.scatter(d.loc[valid,x],d.loc[valid,y],s=8,alpha=.42,c=color,edgecolors='none',zorder=2)
    zero(ax);ax.set(xlabel=xlab,ylabel=ylab)
    for axis in [ax.xaxis,ax.yaxis]:
        sf=ScalarFormatter(useMathText=False);sf.set_powerlimits((-2,2));axis.set_major_formatter(sf);axis.get_offset_text().set_size(7.2)
    ax.tick_params(length=3,width=.7)
    if annotate:text(ax,.04,.94,f'n = {n}\nρ = {rho:.3f}',va='top')
def point_ci(ax,x,lo,hi,y,color=INSERTION,marker='o',horizontal=True):
    if horizontal:ax.plot([lo,hi],[y,y],c=color,lw=.9);ax.plot(x,y,marker,c=color,ms=4,clip_on=False)
    else:ax.plot([y,y],[lo,hi],c=color,lw=.9);ax.plot(y,x,marker,c=color,ms=4)
def forest(ax,rows,labels,sign=False,colors=None,markers=None,xlim=None,values=True):
    for i,(_,z) in enumerate(rows.iterrows()):
        v=z.sign_fraction*100 if sign else z.rho
        lo=z.sign_low*100 if sign else z.rho_low;hi=z.sign_high*100 if sign else z.rho_high
        c=colors[i] if colors else INSERTION;mk=markers[i] if markers else 'o'
        point_ci(ax,v,lo,hi,i,c,mk)
        if values:text(ax,1.04,1-(i+.5)/len(rows),f'{v:.1f}%' if sign else f'{v:.3f}',va='center',color=c)
    ax.set(yticks=range(len(rows)),yticklabels=labels,ylim=(len(rows)-.5,-.5))
    ax.axvline(50 if sign else 0,c=LIGHT,lw=.6,ls=':');ax.set_xlabel('Sign concordance (%)' if sign else 'Spearman ρ')
    ax.tick_params(axis='y',length=0,pad=5)
    if xlim:ax.set_xlim(xlim)
    if sign:ax.xaxis.set_major_formatter(PercentFormatter(100,decimals=0))
def finish(fig,n,height,supp=False):
    dest=R/'supplementary' if supp else O;stem=f'Figure_{n}' if supp else f'Figure{n}'
    fig.savefig(dest/f'{stem}.pdf');fig.savefig(dest/f'{stem}.svg')
    for dpi in [300,600]:fig.savefig(dest/f'{stem}_{dpi}dpi.png',dpi=dpi)
    fig.savefig(QA/f'{stem}_preview.png',dpi=300)
    # Grayscale and a deterministic deuteranopia simulation, performed in Python.
    im=Image.open(QA/f'{stem}_preview.png').convert('RGB')
    im.convert('L').save(QA/f'{stem}_grayscale.png')
    a=np.asarray(im)/255;lin=np.where(a<=.04045,a/12.92,((a+.055)/1.055)**2.4)
    mat=np.array([[.367322,.860646,-.227968],[.280085,.672501,.047413],[-.011820,.042940,.968881]])
    sim=np.clip(lin@mat.T,0,1);sim=np.where(sim<=.0031308,12.92*sim,1.055*sim**(1/2.4)-.055)
    Image.fromarray((np.clip(sim,0,1)*255).astype('uint8')).save(QA/f'{stem}_deuteranopia.png')
    if not supp:exports.append({'figure':n,'width_mm':170,'height_mm':height,'PDF':'Embedded Arial TrueType','SVG':'Editable text','raster_dpi':[300,600]})
    plt.close(fig)

# Figure 1: overview, with the comparison as the dominant lower-right panel.
fig=canvas(1,127)
a=panel(fig,[.065,.57,.27,.32],'a','Research question',True,'Figure1_source.tsv')
allele(a,.66,False,'REF');allele(a,.38,True,'ALT');text(a,0,.22,'Shared flanks + complete insertion',7.2)
text(a,0,.075,'Does reference-space S2F\nperformance transfer to complete\nnon-reference alleles?',7.5,va='top')
b=panel(fig,[.39,.57,.23,.32],'b','Frozen benchmark',True,'Figure1_source.tsv')
for y,label in [(.78,'32,969 Class B loci'),(.48,'400 frozen core loci'),(.18,'318 jointly scorable')]:box(b,.05,y,.9,.18,label)
arrow(b,(.5,.77),(.5,.67));arrow(b,(.5,.47),(.5,.37))
b.plot([.95,1.035,1.035],[.57,.57,-.055],c=LIGHT,lw=.6,transform=b.transAxes,clip_on=False)
text(b,.05,-.07,'82 no common exon\nRetained in denominator',7.2,color=BORZOI)
c=panel(fig,[.685,.57,.26,.32],'c','Paired scoring',True,'Figure1_source.tsv')
text(c,.5,.84,'REF + ALT',8,ha='center');box(c,0,.47,.46,.20,'AlphaGenome',AG,size=7.2);box(c,.54,.47,.46,.20,'Borzoi',BORZOI,size=7.2)
arrow(c,(.45,.80),(.23,.69));arrow(c,(.55,.80),(.77,.69))
text(c,.5,.23,'ΔRNA = ALT − REF\nSame homologous exon bases',7.5,ha='center',va='top')
d=panel(fig,[.065,.08,.31,.32],'d','Evidence hierarchy',True,'Figure1_source.tsv')
for i,label in enumerate(['Technical reproducibility','Cross-model reliability','Biological calibration']):
    y=.72-i*.28;text(d,0,y,str(i+1),9,weight='bold');text(d,.12,y,label,8)
    if i<2:arrow(d,(.04,y-.07),(.04,y-.19))
e=panel(fig,[.445,.075,.50,.325],'e','Core finding',True,'Figure1_source.tsv')
for x,label,value,col in [(0,'REFERENCE','ρ = 0.626',REFERENCE),(.36,'OFF-REFERENCE','ρ = −0.111',INSERTION),(.77,'GAP','Δρ = 0.737',REFERENCE)]:
    text(e,x,.77,label,7.2,color=col);text(e,x,.58,value,10,color=col)
text(e,0,.34,'Nonmatched benchmarks · GNRS signed truth: 9/33 usable',7.2)
text(e,0,.06,'Reference-space coherence did not transfer\nto the evaluated complete insertions.',8.5)
finish(fig,1,127)

# Figure 2: retain every jointly finite raw observation; an explicitly labelled inset
# shows the crowded centre without removing or changing any observations.
fig=canvas(2,147);d2=data('Figure2a_source.tsv')
a=panel(fig,[.125,.51,.78,.375],'a','Complete-insertion effects')
scatter(a,d2,'AlphaGenome_original','Borzoi_original',318,h['original']['rho'])
inset=a.inset_axes([.69,.10,.27,.30]);scatter(inset,d2,'AlphaGenome_original','Borzoi_original',318,h['original']['rho'],annotate=False,xlab='',ylab='')
inset.set(xlim=(-.001,.001),ylim=(-.003,.003));inset.set_title('Centre zoom',fontsize=7.2,pad=2);inset.tick_params(labelsize=7)
b=panel(fig,[.215,.095,.21,.245],'b','Primary metrics')
b.axis('off')
rank=fig.add_axes([.215,.265,.21,.065]);forest(rank,metrics.iloc[[0]],['Rank'],xlim=(-.35,.12))
sgn=fig.add_axes([.215,.115,.21,.065]);forest(sgn,metrics.iloc[[0]],['Sign'],sign=True,xlim=(30,65))
c=panel(fig,[.65,.095,.23,.245],'c','Harmonization sensitivity')
forest(c,metrics.iloc[:5],['Original','Cached exact','New strict','Near family','Full RNA'],xlim=(-.35,.12),values=False)
for i,(_,z) in enumerate(metrics.iloc[:5].iterrows()):text(c,1.05,1-(i+.5)/5,f'{z.sign_fraction*100:.1f}%',va='center',color=BORZOI)
text(c,1.05,1.03,'Sign',7.2,color=BORZOI)
finish(fig,2,147)

# Figure 3: within-model robustness is distinct from cross-model reliability.
fig=canvas(3,131)
a=panel(fig,[.07,.68,.25,.18],'a','Placement',True)
for i,(x,label) in enumerate([(0,'−32.8 kb'),(.37,'Baseline'),(.74,'+32.8 kb')]):
    text(a,x,.74,label,7.2,ha='left');a.plot([x,x+.22],[.48,.48],c=LIGHT,lw=3,transform=a.transAxes)
    a.plot(x+.11,.48,'D',c=INSERTION,ms=4,transform=a.transAxes)
text(a,0,.10,'Same complete allele\nPrespecified ±32,768 bp shifts',7.2)
b=panel(fig,[.48,.64,.46,.22],'b','Within-model placement robustness',True)
for j,m in enumerate(['AlphaGenome','Borzoi']):
    ax=fig.add_axes([.47+j*.255,.68,.205,.13]);z=metrics[(metrics.analysis_label=='PLACEMENT')&(metrics.model==m)]
    forest(ax,z,['−32.8 kb','+32.8 kb'] if j==0 else ['',''],colors=[COLORS[m]]*2,markers=[MARKERS[m]]*2,xlim=(.5,1.03),values=False)
    ax.set_title(m,fontsize=7.8,pad=8)
    for k,(_,row) in enumerate(z.iterrows()):ax.annotate(f'ρ {row.rho:.3f} · sign {row.sign_fraction*100:.1f}%',(.5,k),xytext=(0,8),textcoords='offset points',fontsize=7.2)
text(b,0,-.36,'n = 72/100 per model and shift',7.2)
c=panel(fig,[.125,.11,.31,.34],'c','Official-member stability')
mm=data('Figure3c_source.tsv')
for j,m in enumerate(['AlphaGenome','Borzoi']):
    d=mm[(mm.model==m)&mm.valid_members_N.eq(4)];v=d.dominant_sign_fraction.to_numpy();assert len(v)==77
    c.scatter(j+np.linspace(-.14,.14,len(v)),v,s=12,c=COLORS[m],marker=MARKERS[m],alpha=.45,edgecolors='none')
    c.plot([j-.20,j+.20],[.75,.75],c=COLORS[m],lw=1.2)
    text(c,.03+j*.51,.99,'87.0%' if j==0 else '90.9%',7.2,va='top')
c.set(xticks=[0,1],xticklabels=['AlphaGenome','Borzoi'],ylim=(.44,1.12),xlim=(-.4,1.4),ylabel='Dominant sign fraction',yticks=[.5,.75,1])
text(c,0,-.22,'Median = 0.75 · n = 77/100 each\nTop labels: non-unanimous fraction',7.2)
d=panel(fig,[.59,.12,.35,.33],'d','Robustness synthesis',True)
text(d,0,.81,'WITHIN-MODEL PLACEMENT',7.5,weight='bold');text(d,0,.63,'High stability',9,color=REFERENCE)
text(d,0,.37,'CROSS-MODEL INSERTIONS',7.5,weight='bold');text(d,0,.19,'Limited agreement',9,color=INSERTION)
text(d,0,-.08,'Official members are fitted models,\nnot biological replicates.',7.2,color=BORZOI)
finish(fig,3,131)

# Figure 4: categorical context steps; native points are identities, not estimates.
fig=canvas(4,151)
a=panel(fig,[.07,.64,.38,.24],'a','Distinct context interventions',True)
for y,m,width,desc in [(.75,'AlphaGenome','1,048.6 kb','Native shorter inputs'),(.30,'Borzoi','524.3 kb','N-mask within native input')]:
    text(a,0,y,m,8,weight='bold',color=COLORS[m]);text(a,0,y-.13,'16.4 → 131.1 → '+width,7.5)
    text(a,0,y-.26,desc,7.2,color=BORZOI)
b=panel(fig,[.635,.64,.285,.24],'b','Rank stability vs native')
c=panel(fig,[.125,.115,.31,.33],'c','Sign stability vs native')
for ax,sign in [(b,False),(c,True)]:
    for j,m in enumerate(['AlphaGenome','Borzoi']):
        z=metrics[(metrics.analysis_label=='CONTEXT')&(metrics.model==m)];xs=np.array([0,1,2])+(-.035 if j==0 else .035)
        vals=list(z.sign_fraction*100 if sign else z.rho)+[100 if sign else 1]
        ax.plot(xs,vals,c=COLORS[m],ls='-' if j==0 else '--',lw=.8)
        for k,(_,row) in enumerate(z.iterrows()):
            lo=row.sign_low*100 if sign else row.rho_low;hi=row.sign_high*100 if sign else row.rho_high
            point_ci(ax,vals[k],lo,hi,xs[k],COLORS[m],MARKERS[m],horizontal=False)
        ax.plot(xs[2],vals[2],MARKERS[m],c=COLORS[m],mfc='white',ms=4)
    ax.set(xticks=[0,1,2],xticklabels=['16 kb','131 kb','Native'],xlim=(-.2,2.2),xlabel='Context condition',ylabel='Sign concordance (%)' if sign else 'Spearman ρ')
    ax.axhline(50 if sign else 0,c=LIGHT,lw=.6,ls=':')
    ax.set_ylim((15,110) if sign else (-.55,1.14));text(ax,.02,.97,'n = 29/100',va='top',size=7.2)
    if sign:ax.yaxis.set_major_formatter(PercentFormatter(100,decimals=0))
for m in ['AlphaGenome','Borzoi']:b.plot([],[],MARKERS[m]+('-' if m=='AlphaGenome' else '--'),c=COLORS[m],label=m,ms=4)
b.legend(loc='lower right',fontsize=7.2)
d=panel(fig,[.635,.115,.285,.33],'d','Raw effects · AlphaGenome 16 kb')
cx=data('Figure4d_source.tsv');cx=cx[cx.plotted]
scatter(d,cx,'native_1048576','restricted_16384',29,.05566502463054186,AG,xlab='Native ΔRNA',ylab='16 kb ΔRNA')
finish(fig,4,151)

# Figure 5: centerpiece, two raw scatters and two metric contrasts.
fig=canvas(5,157)
a=panel(fig,[.085,.83,.85,.08],'a','Nonmatched benchmark comparator',True)
text(a,0,.55,'REFERENCE',8,weight='bold',color=REFERENCE);text(a,0,-.10,'500 SNVs · 473 genes',8)
text(a,.52,.55,'OFF-REFERENCE',8,weight='bold',color=INSERTION);text(a,.52,-.10,'400 insertions · 318 jointly scorable',8)
text(a,.34,.53,'≠',12,ha='center',color=BORZOI)
b=panel(fig,[.13,.485,.31,.23],'b','Reference effects')
scatter(b,data('Figure5b_source.tsv'),'AlphaGenome_TIER1','Borzoi_TIER1',500,ref['strict_cross_model']['rho'],REFERENCE)
c=panel(fig,[.625,.485,.31,.23],'c','Off-reference effects')
scatter(c,d2,'AlphaGenome_original','Borzoi_original',318,h['original']['rho'],INSERTION)
d=panel(fig,[.22,.19,.22,.16],'d','Rank agreement')
rr=metrics.iloc[[5,0]];forest(d,rr,['Reference','Off-reference'],colors=[REFERENCE,INSERTION],markers=['o','D'],xlim=(-.35,.85),values=False)
for i,(_,z) in enumerate(rr.iterrows()):text(d,.02,.88-i*.5,f'{z.rho:.3f} [{z.rho_low:.3f}, {z.rho_high:.3f}]',7.2)
text(d,-.4,-.52,f"Δρ = {gap['rho_difference']:.3f} [{gap['rho_difference_CI95'][0]:.3f}, {gap['rho_difference_CI95'][1]:.3f}]",7.5)
e=panel(fig,[.715,.19,.22,.16],'e','Direction agreement')
forest(e,rr,['Reference','Off-reference'],sign=True,colors=[REFERENCE,INSERTION],markers=['o','D'],xlim=(35,80),values=False)
for i,(_,z) in enumerate(rr.iterrows()):text(e,.02,.88-i*.5,f'{z.sign_fraction*100:.1f}% [{z.sign_low*100:.1f}, {z.sign_high*100:.1f}]',7.2)
text(e,-.4,-.52,f"Δ = {gap['sign_difference']*100:.1f} pp [{gap['sign_difference_CI95'][0]*100:.1f}, {gap['sign_difference_CI95'][1]*100:.1f}]",7.5)
f=panel(fig,[.085,.015,.85,.035],'f','Benchmark generalization gap',True)
text(f,0,.45,'Nonmatched: variant class and ascertainment differ',7.5)
text(f,0,-.05,'The contrast does not isolate a causal effect of reference absence.',7.5,color=BORZOI)
finish(fig,5,157)

# Figure 6: calibration remains incomplete; no visual pass badge for model agreement.
fig=canvas(6,165)
a=panel(fig,[.07,.73,.40,.14],'a','GNRS recovery boundary',True)
box(a,0,.58,.40,.32,'33 frozen loci');box(a,.60,.58,.40,.32,'9 jointly usable');arrow(a,(.41,.74),(.59,.74))
box(a,.13,.04,.74,.28,'13 theoretical Borzoi maximum',BORZOI,ls='--',size=7.2)
arrow(a,(.20,.57),(.40,.36),ls='--')
b=panel(fig,[.655,.75,.265,.12],'b','GNRS direction · descriptive only')
z=metrics[metrics.analysis_label=='GNRS descriptive'];forest(b,z,['AlphaGenome','Borzoi'],sign=True,colors=[AG,BORZOI],markers=['o','s'],xlim=(-2,102),values=False)
for i,(_,v) in enumerate(z.iterrows()):text(b,1.03,1-(i+.5)/2,f'{int(v.concordant_N)}/9\n{v.sign_fraction*100:.1f}%',7.2,va='center')
text(b,-.38,-.70,'Below prespecified N = 20 gate',7.2,color=BORZOI)
c=panel(fig,[.07,.39,.33,.21],'c','External truth attrition',True)
for y,label in [(.76,'13,272 association rows'),(.47,'1,994 anchor/length candidate rows'),(.18,'558 Class B catalogue loci')]:
    box(c,0,y,1,.20,label,size=7.2)
for y in [.75,.46]:arrow(c,(.5,y),(.5,y-.08))
text(c,0,-.06,'Exact sequence / orientation /\nindependent signed truth: unresolved',7.2,color=BORZOI)
d=panel(fig,[.48,.37,.46,.23],'d','Evidence established and unresolved',True)
text(d,0,.92,'ESTABLISHED',7.2,weight='bold');text(d,.52,.92,'UNRESOLVED',7.2,weight='bold',color=BORZOI)
left=['+ Complete paired alleles','+ Homology-aware scoring','+ Limited model agreement','+ Stronger reference agreement','+ Context sensitivity']
right=['○ Biological model ranking','○ Population-level accuracy','○ Matched reference-SV test','○ Independent signed truth']
for i,label in enumerate(left):text(d,0,.72-i*.17,label,7.2)
for i,label in enumerate(right):text(d,.52,.72-i*.17,label,7.2,color=BORZOI)
e=panel(fig,[.07,.055,.87,.16],'e','Five levels of annotation evidence',True)
labels=['Sequence\nvalidity','Coordinate\nvalidity','Technical\nrobustness','Cross-model\nreliability','Biological\ncalibration']
states=['DOCUMENTED','DOCUMENTED','ASSESSED','ASSESSED\nLimited agreement','INCOMPLETE']
for i,(label,state) in enumerate(zip(labels,states)):
    x=i*.205;col=REFERENCE if i<3 else (INSERTION if i==3 else BORZOI)
    text(e,x,.94,f'LEVEL {i+1}',7.2,weight='bold',color=col);text(e,x,.65,label,7.5,va='top',color=col)
    e.plot([x,x+.17],[.27,.27],c=col if i<4 else LIGHT,lw=.8,ls='-' if i<4 else '--',transform=e.transAxes)
    text(e,x,.18,state,7,va='top',color=col)
finish(fig,6,165)

# Material actually relocated, avoiding placeholder destinations.
fig=canvas('S12',92);ex=pd.read_csv(R/'supplementary/Figure_S12_source.tsv',sep='\t')
ax=panel(fig,[.15,.20,.72,.62],'a','Frozen post-hoc illustrative loci',source='supplementary/Figure_S12_source.tsv')
xs=np.arange(len(ex))
for m in ['AlphaGenome','Borzoi']:
    dx=-.07 if m=='AlphaGenome' else .07;ax.plot(xs+dx,ex[m+'_z'],MARKERS[m],c=COLORS[m],ms=5,label=m)
ax.axhline(0,c=LIGHT,lw=.6);ax.set(xticks=xs,xticklabels=['L1','L2','L3','L4'],ylabel='Standardized ΔRNA');ax.legend()
finish(fig,'S12',92,True)
fig=canvas('S13',141)
for i,(m,setting) in enumerate([('AlphaGenome',16384),('AlphaGenome',131072),('Borzoi',16384),('Borzoi',131072)]):
    cx=data('Figure4_raw_complete.tsv');cx=cx[cx.model==m].pivot(index='OFFREF_ID',columns='setting',values='delta');cx=cx[cx.notna().all(axis=1)]
    base=1048576 if m=='AlphaGenome' else 524288
    ax=panel(fig,[.135+(i%2)*.49,.58-(i//2)*.46,.30,.30],chr(97+i),m+f' {setting/1000:.1f} kb',source='Figure4_raw_complete.tsv')
    zz=next(z for z in pc if z['model']==m and z['analysis']=='CONTEXT' and z['alternative']==setting)
    scatter(ax,cx,base,setting,29,zz['rho'],COLORS[m],xlab='Native ΔRNA',ylab='Restricted ΔRNA')
finish(fig,'S13',141,True)
(QA/'PANEL_MANIFEST.json').write_text(json.dumps(panels,indent=2))
(QA/'EXPORT_MANIFEST.json').write_text(json.dumps(exports,indent=2))
print('SIX_MAIN_AND_TWO_RELOCATED_FIGURES_EXPORTED')
