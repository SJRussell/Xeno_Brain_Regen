from rebuild_figures import *
import importlib.util

n=load('render_fig2',1.60)
# Re-render the heatmap from its locked row/column order; enlarge text, not bitmap pixels.
S.ZMAX_HEATMAP=4; S.TP_COL={'24HPA':'#F0A860','7DPA':'#4C9E7C'}
S.FS_SCALE=1.5;S.DECK_FS_TICK=12.5;S.DECK_FS_LEGEND=13.5;S.DECK_FS_TITLE=16
spec=importlib.util.spec_from_file_location('heatmap_for_revision',BASE/'figures_src/_heatmap.py')
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
expr,meta=h.load_brain();genes=pd.read_csv(BASE/'tables/brain/brain_heatmap_CTtemporal_DEGs_gene_order.tsv',sep='\t').gene_id.tolist()
h.render_heatmap(expr,meta,genes,'MDT temporal DEGs (924 displayed)',PAN/'fig2F_temporal_heatmap.png',lock_cols=h.CANON_2F_COLS,cluster_rows=False,gene_labels=False,figsize=(6.2,7.2),title_size=16,xtick_size=12.5,legend_size=13.5)

compose('figure02',[
 (PAN/'fig2A_brain_volcano_24hpt_deck.png',.65,.75,3.15,2.73),
 (PAN/'fig2B_brain_volcano_7dpt_deck.png',.65,3.85,3.15,3.01),
 (PAN/'fig2C_brain_volcano_MDT_7v24_deck.png',.65,7.18,3.15,2.73),
 (PAN/'fig2D_brain_volcano_ND_7v24_deck.png',.65,10.27,3.15,3.01),
 (PAN/'fig2F_temporal_heatmap.png',3.90,.45,6.2,7.2),
 (PAN/'fig2E_temporal_venn.png',4.0,8.0,6.1,4.95)
],(10.2,13.45),[('A',.58,.08),('B',.58,6.65),('C',3.85,.08),('D',3.95,7.65)])
fp=OUT/'figure02.png';im=Image.open(fp);draw=ImageDraw.Draw(im);DPI=300
for text_,y,col in [('MDT vs ND, 24 HPT',.50,'#F4D3CC'),('MDT vs ND, 7 DPT',3.60,'#F4D3CC'),('MDT: 7 DPT vs 24 HPT',6.94,'#D8E7F5'),('ND: 7 DPT vs 24 HPT',10.03,'#D8E7F5')]:
    box=(.88*DPI,y*DPI,3.7*DPI,(y+.24)*DPI);draw.rectangle(box,fill=col)
    draw.text(((.88+3.7)/2*DPI,(y+.12)*DPI),text_,anchor='mm',font=ImageFont.truetype(FONT,round(14*DPI/72)),fill='black')
for label,y,hh in [('Treatment contrasts',.5,6.25),('Timepoint contrasts',6.94,6.3)]:
    tag=Image.new('RGB',(round(hh*DPI),round(.46*DPI)),'#EAEAEA');d=ImageDraw.Draw(tag)
    d.text((tag.width/2,tag.height/2),label+'  (−log10 adjusted P)',anchor='mm',font=ImageFont.truetype(FONT,round(14*DPI/72)),fill='black')
    im.paste(tag.rotate(90,expand=True),(round(.04*DPI),round(y*DPI)))
im.save(fp,dpi=(300,300))
print('Figure 2 rebuilt with neutral contrast labels, retained colored headers, larger chart fonts')
