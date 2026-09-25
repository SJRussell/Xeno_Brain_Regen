"""Presentation-only revisions; canonical numeric inputs and tests are unchanged."""
from pathlib import Path
import sys, os, re, json, io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import pdist
from PIL import Image, ImageDraw, ImageFont
from pptx import Presentation  # read-only geometry and original assets

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
SRC=BASE/'figure_rebuild_2026-07'
OUT=ROOT/'figures'; OUT.mkdir(exist_ok=True)
PAN=ROOT/'panels'; PAN.mkdir(exist_ok=True)
sys.path.insert(0,str(SRC))
os.environ['XENO_FS_SCALE']='1.8'
import _style as S
original_save=Figure.savefig
def save(self,path,*args,**kwargs):
    dest=PAN/Path(path).name
    original_save(self,dest,*args,**kwargs)
Figure.savefig=save

def load(name,scale=1.8):
    S.FS_TITLE=10*scale; S.FS_AXIS=9*scale; S.FS_TICK=9*scale
    S.FS_LEGEND=8*scale; S.FS_ANNOT=8*scale; S.apply()
    source=(SRC/(name+'.py')).read_text(encoding='utf-8')
    source=source.split('if __name__ == "__main__":')[0]
    if name=='render_fig2':
        source=source.replace('float(d["nlp"].max()) * 1.12','float(d["nlp"].max()) * 1.30')
        source=source.replace('expand=(1.2, 1.5), ensure_inside_axes=True','expand=(1.2, 1.5), ensure_inside_axes=True, iter_lim=600, force_text=(0.5,0.8), force_static=(0.2,0.4)')
        source=source.replace('auto_label_n=6)','auto_label_n=4)')
        source=source.replace('left=0.13, right=0.90, top=0.97','left=0.13, right=0.90, top=0.85')
        source=source.replace('ax.legend(loc=best, **lkw)','fig.legend(*ax.get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.52,1.0), ncol=2, fontsize=S.FS_LEGEND - 1, frameon=False, columnspacing=0.8, handletextpad=0.2)')
    if name=='render_fig6':
        source=source.replace('fontsize=S.FS_AXIS - 2, fontstyle="italic"','fontsize=S.FS_AXIS, fontstyle="italic"')
        source=source.replace('textwrap.fill(d[:1].upper() + d[1:], 26)','textwrap.fill(d[:1].upper() + d[1:], 24)')
        source=source.replace('fig.text(0.55, 0.965, "Synaptic and axon GO BP: 24 HPT vs 7 DPT (brain)"','fig.text(0.55, 0.965, "Synaptic and axon-associated GO terms"')
        source=source.replace('lm, bm, gx, gy = 0.115, 0.125, 0.075, 0.125','lm, bm, gx, gy = 0.12, 0.15, 0.08, 0.15')
        source=source.replace('ax.text(0.03, 0.04,','ax.text(0.5, 1.08,').replace('ha="left", va="bottom", fontsize=S.FS_ANNOT - 2.5','ha="center", va="bottom", fontsize=S.FS_ANNOT - 0.5')
        source=source.replace('fontstyle="italic", color="#333333", pad=2','fontstyle="italic", color="#333333", pad=22')
        source=source.replace('rm, tm = 0.02, 0.135','rm, tm = 0.035, 0.21')
        source=source.replace('fontsize=S.FS_ANNOT - 0.5','fontsize=S.FS_ANNOT - 2.5')
        source=source.replace('Brain synaptic and structural gene trajectories','Synaptic and structural gene trajectories')
        source=source.replace('Liao zbtb18+ excitatory-neuron signature in brain samples','zbtb18+ excitatory-neuron signature')
        source=source.replace('ax.set_xlim(-2.5, 2.5); ax.set_ylim(-0.6, len(top) - 0.4)','ax.set_xlim(-2.5, 2.5); ax.set_xticks([-2, 0, 2]); ax.set_ylim(-0.6, len(top) - 0.4)')
        source=source.replace('NES (MDT vs ND)','NES')
        source=source.replace('fig.text(0.022, cy','fig.text(0.003, cy')
    if name=='render_fig7':
        source=source.replace('uninjured baseline (','baseline (')
        source=source.replace('etv1+ signature vs baseline','etv1+ signature vs baseline')
    if name=='render_fig9':
        source=source.replace('fontsize=S.FS_TICK - 2)','fontsize=S.FS_TICK, rotation=90)')
    ns={'__file__':str(SRC/(name+'.py')),'__name__':name+'_revision'}
    exec(compile(source,str(SRC/(name+'.py')),'exec'),ns)
    ns['OUT']=PAN;ns['FPNG']=PAN
    return ns

GROUPS={
 'Metabolism': 'OXIDATIVE_PHOSPHORYLATION GLYCOLYSIS FATTY_ACID_METABOLISM CHOLESTEROL_HOMEOSTASIS BILE_ACID_METABOLISM XENOBIOTIC_METABOLISM HEME_METABOLISM ADIPOGENESIS PEROXISOME',
 'Immune and inflammatory': 'ALLOGRAFT_REJECTION COMPLEMENT INFLAMMATORY_RESPONSE INTERFERON_ALPHA_RESPONSE INTERFERON_GAMMA_RESPONSE IL2_STAT5_SIGNALING IL6_JAK_STAT3_SIGNALING TNFA_SIGNALING_VIA_NFKB',
 'Stress and survival': 'HYPOXIA REACTIVE_OXYGEN_SPECIES_PATHWAY UNFOLDED_PROTEIN_RESPONSE DNA_REPAIR P53_PATHWAY APOPTOSIS UV_RESPONSE_UP UV_RESPONSE_DN',
 'Growth and cell cycle': 'MTORC1_SIGNALING PI3K_AKT_MTOR_SIGNALING MYC_TARGETS_V1 MYC_TARGETS_V2 E2F_TARGETS G2M_CHECKPOINT MITOTIC_SPINDLE PROTEIN_SECRETION SPERMATOGENESIS',
 'Development and signaling': 'EPITHELIAL_MESENCHYMAL_TRANSITION ANGIOGENESIS MYOGENESIS APICAL_JUNCTION APICAL_SURFACE COAGULATION HEDGEHOG_SIGNALING NOTCH_SIGNALING WNT_BETA_CATENIN_SIGNALING TGF_BETA_SIGNALING KRAS_SIGNALING_UP KRAS_SIGNALING_DN ANDROGEN_RESPONSE ESTROGEN_RESPONSE_EARLY ESTROGEN_RESPONSE_LATE PANCREAS_BETA_CELLS'
}
GCOL=['#E6F0F4','#EDF1E4','#F8ECE2','#ECE8F3','#F3E8EA']
CAT={f'HALLMARK_{g}':cat for cat,genes in GROUPS.items() for g in genes.split()}
assert len(CAT)==50
def short(t):
    s=t.replace('HALLMARK_','').replace('_',' ').capitalize()
    return {'Epithelial mesenchymal transition':'Epithelial–mesenchymal (EMT)', 'Reactive oxygen species pathway':'Reactive oxygen species', 'Tnfa signaling via nfkb':'TNFα signaling via NFκB', 'Il6 jak stat3 signaling':'IL6–JAK–STAT3 signaling', 'Pi3k akt mtor signaling':'PI3K–AKT–mTOR signaling', 'Wnt beta catenin signaling':'Wnt / β-catenin signaling', 'Mtorc1 signaling':'mTORC1 signaling', 'Myc targets v1':'MYC targets V1', 'Myc targets v2':'MYC targets V2', 'E2f targets':'E2F targets', 'G2m checkpoint':'G2M checkpoint', 'Il2 stat5 signaling':'IL2–STAT5 signaling', 'Dna repair':'DNA repair', 'P53 pathway':'p53 pathway', 'Tgf beta signaling':'TGFβ signaling'}.get(s,s)
def groupsort(terms): return sorted(terms,key=lambda t:(list(GROUPS).index(CAT[t]),t))

def grouped7b(ns):
    m=ns['_hallmark_vs_baseline'](); d=m[m.sig_either].reindex(groupsort(m[m.sig_either].index))
    fig=plt.figure(figsize=(7.2,7.3)); ax=fig.add_axes([.61,.19,.36,.72])
    yy=[]; pos=0; bands=[]
    for cat in GROUPS:
        members=[x for x in d.index if CAT[x]==cat]
        if not members: continue
        start=pos
        for t in members: yy.append((t,pos));pos+=1
        bands.append((cat,start,pos));pos+=.65
    for cat,a,b in bands:
        ax.axhspan(a-.5,b-.5,color=GCOL[list(GROUPS).index(cat)],zorder=0)
        ax.text(-.04,(a+b-1)/2,'',transform=ax.get_yaxis_transform())
    ys=np.array([y for _,y in yy]); d=d.loc[[t for t,_ in yy]]
    ax.axvline(0,color='#555555',lw=.6)
    for y,r in zip(ys,d.itertuples()): ax.plot([r.NES_ND,r.NES_MDT],[y,y],color='#888888',lw=1.4)
    for c,p,label,sz,z in [('NES_ND','padj_ND','ND',95,3),('NES_MDT','padj_MDT','MDT',35,4)]:
        color=S.COL_ND if label=='ND' else S.COL_MDT
        for sig in [True,False]:
            sel=(d[p]<.05)==sig
            ax.scatter(d.loc[sel,c],ys[sel.values],s=sz,facecolor=color if sig else 'white',edgecolor='black' if sig else color,linewidths=.7,zorder=z)
    ax.set_yticks(ys);ax.set_yticklabels([short(t) for t in d.index],fontsize=18)
    ax.invert_yaxis(); ax.set_ylim(pos-.2,-.7);ax.tick_params(axis='y',length=0)
    ax.set_xlabel('NES vs baseline',fontsize=17)
    ax.tick_params(axis='x',labelsize=14); ax.grid(axis='x',alpha=.15)
    fig.text(.52,.97,'Hallmark shifts from uninjured baseline',ha='center',fontsize=17)
    handles=[Line2D([],[],marker='o',ls='',mfc=c,mec='black',ms=s,label=l) for c,s,l in [(S.COL_ND,9,'ND'),(S.COL_MDT,6,'MDT')]]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.75,.035),ncol=2,frameon=False,fontsize=16)
    fig.text(.05,.045,'Groups as in Fig. 9D',fontsize=13,color='#555555')
    fig.savefig(PAN/'fig7B_hallmark_vs_baseline.png',dpi=300);plt.close(fig)
    (ROOT/'hallmark_display_groups.json').write_text(json.dumps(CAT,indent=2))

def grouped9d(ns):
    data=[]
    for path,label in ns['COMP']:
        d=pd.read_csv(path,sep='\t')[['term','NES','padj']];d['column']=label;data.append(d)
    long=pd.concat(data); rows=groupsort(long.term.unique());cols=[label for _,label in ns['COMP']]
    nes=long.pivot(index='term',columns='column',values='NES').reindex(rows)[cols]
    q=long.pivot(index='term',columns='column',values='padj').reindex(rows)[cols]
    Z=linkage(pdist(np.nan_to_num(nes.to_numpy().T),metric='euclidean'),method='ward')
    order=leaves_list(Z); values=nes.iloc[:,order].to_numpy();qs=q.iloc[:,order].to_numpy()
    fig=plt.figure(figsize=(8.0,15.3));L=.46;W=.47
    axd=fig.add_axes([L,.923,W,.042]); dendrogram(Z,ax=axd,no_labels=True,color_threshold=0,above_threshold_color='#444444');axd.axis('off')
    ax=fig.add_axes([L,.155,W,.76]);im=ax.imshow(values,cmap='RdBu_r',vmin=-3,vmax=3,aspect='auto')
    for i in range(50):
        for j in range(8):
            if np.isfinite(qs[i,j]) and qs[i,j]<.05: ax.text(j,i,f'{values[i,j]:+.1f}',ha='center',va='center',fontsize=12)
    ax.set_yticks(range(50));ax.set_yticklabels([short(t) for t in rows],fontsize=15)
    disp={ns['XEN24']:'X. laevis 24 HPT',ns['XEN7']:'X. laevis 7 DPT','Mouse cortex burn vs nt':'Cortex burn vs untreated','Mouse cortex formalin vs nt':'Cortex formalin vs untreated','Mouse cortex burn vs formalin':'Cortex burn vs formalin'}
    ax.set_xticks(range(8)); ax.set_xticklabels([disp.get(cols[i],cols[i]) for i in order],rotation=55,ha='right',fontsize=14)
    ax.tick_params(length=0);ax.set_xticks(np.arange(-.5,8,1),minor=True);ax.set_yticks(np.arange(-.5,50,1),minor=True);ax.grid(which='minor',color='#EEEEEE',lw=.4);ax.tick_params(which='minor',length=0)
    bounds=[]
    for k,cat in enumerate(GROUPS):
        ix=[i for i,t in enumerate(rows) if CAT[t]==cat]; a=min(ix);b=max(ix)+1
        if a:ax.axhline(a-.5,color='white',lw=3)
        # Labeled bands outside the labels; biological grouping changes row order only.
        y=.155+.76*(1-b/50);h=.76*(b-a)/50
        fig.patches.append(Rectangle((.013,y),.040,h,transform=fig.transFigure,color=GCOL[k],ec='none'))
        fig.text(.033,y+h/2,cat.replace(' and ',' / '),rotation=90,ha='center',va='center',fontsize=13)
        bounds.append((cat,a,b))
    cax=fig.add_axes([.12,.950,.24,.012]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_ticks([-3,0,3]);cb.ax.tick_params(labelsize=14);cb.set_label('NES (limits ±3)',fontsize=14)
    fig.text(.60,.986,'Hallmark profiles across comparator contrasts',ha='center',fontsize=17)
    fig.savefig(PAN/'fig9D_hallmark_comparator_heatmap.png',dpi=300);plt.close(fig)
    (ROOT/'figure9_numeric_check.json').write_text(json.dumps({'rows':rows,'columns':[cols[i] for i in order],'significant_cells':int((qs<.05).sum()),'NES':values.tolist()},indent=2))

FONT='C:/Windows/Fonts/arial.ttf';BOLD='C:/Windows/Fonts/arialbd.ttf'
def compose(name,placements,size,letters=None):
    dpi=300;canvas=Image.new('RGB',tuple(round(x*dpi) for x in size),'white')
    for path,x,y,w,h in placements:
        im=Image.open(path).convert('RGB');scale=min(w*dpi/im.width,h*dpi/im.height)
        im=im.resize((round(im.width*scale),round(im.height*scale)),Image.Resampling.LANCZOS)
        canvas.paste(im,(round(x*dpi),round(y*dpi)))
    draw=ImageDraw.Draw(canvas)
    for letter,x,y in letters or []:draw.text((x*dpi,y*dpi),letter,font=ImageFont.truetype(BOLD,round(21*dpi/72)),fill='black')
    canvas.save(OUT/f'{name}.png',dpi=(dpi,dpi))

def main():
    # The four editor-targeted figures are reconstructed from unchanged data.
    n6=load('render_fig6',1.95)
    n6['panel_6A_expanded'](figsize=(7.1,10.0));n6['panel_6B'](figsize=(6.8,7.0));n6['panel_6C'](figsize=(6.8,4.3))
    compose('figure06',[(PAN/'fig6A_synaptic_axon_24v7.png',.15,.3,7.1,10),(PAN/'fig6B_synaptic_gene_trajectories.png',7.4,.3,6.8,7),(PAN/'fig6C_liao_zbtb18_ucell.png',7.4,7.5,6.8,4.3)],(14.3,11.9),[('A',.1,.08),('B',7.35,.08),('C',7.35,7.25)])
    n7=load('render_fig7',1.8)
    n7['panel_7A']();grouped7b(n7);n7['panel_7C'](figsize=(4.5,4.5));n7['panel_7D'](figsize=(5.7,4.5));n7['panel_7E'](figsize=(4.0,4.5))
    compose('figure07',[(PAN/'fig7a_process_state.png',.1,.3,7.4,5.2),(PAN/'fig7B_hallmark_vs_baseline.png',7.5,.3,7.2,7.3),(PAN/'fig7C_baseline_concordance.png',.1,7.95,4.5,4.5),(PAN/'fig7D_gdf15_baseline_anchor.png',4.7,7.95,5.7,4.5),(PAN/'fig7E_walker_etv1_baseline_anchor.png',10.5,7.95,4,4.5)],(14.8,12.55),[('A',.08,.05),('B',7.45,.05),('C',.08,7.7),('D',4.65,7.7),('E',10.45,7.7)])
    fp=OUT/'figure07.png';im=Image.open(fp);dr=ImageDraw.Draw(im)
    for i,(cat,col) in enumerate(zip(GROUPS,GCOL)):
        c,r=divmod(i,3);x=.65+c*3.8;y=6.0+r*.48
        dr.rectangle((x*300,y*300,(x+.25)*300,(y+.25)*300),fill=col)
        dr.text(((x+.34)*300,(y+.12)*300),cat,anchor='lm',font=ImageFont.truetype(FONT,round(14*300/72)),fill='#333333')
    im.save(fp,dpi=(300,300))
    n9=load('render_fig9',1.65)
    n9['panel_9A'](figsize=(5.9,4.0));n9['panel_9B'](figsize=(5.9,4.0));n9['panel_9C'](figsize=(5.9,4.0));grouped9d(n9)
    compose('figure09',[(PAN/'fig9A_walker_umap.png',.15,.3,5.9,4),(PAN/'fig9B_walker_ucell_ranked_bar.png',.15,4.65,5.9,4),(PAN/'fig9C_xlaevis_walker_ucell.png',.15,9,5.9,4),(PAN/'fig9D_hallmark_comparator_heatmap.png',6.25,.3,8,15.3)],(14.35,15.75),[('A',.08,.05),('B',.08,4.4),('C',.08,8.75),('D',6.18,.05)])
    print('Figures 6, 7, 9 rebuilt')
if __name__=='__main__': main()
