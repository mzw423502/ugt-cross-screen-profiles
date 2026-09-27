"""Publication figures from frozen reaction units, at an explicit 180-mm width.

Sources use neutral symbols; chemical series use heading/line colors; paired
reaction states use a separate palette plus patterns. No raw data are altered.
"""
from pathlib import Path
import os
for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import json, io, re, xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D
from PIL import Image
import pymupdf

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'05_REPRODUCIBILITY/input'; GEN=ROOT/'05_REPRODUCIBILITY/generated'
OUT=ROOT/'02_FIGURES'
INK='#24282d'; MUTED='#636b72'; GRID='#d8dde1'
YCOL=MCOL=INK
COLORS={'HCA':'#76578b','Coumarins':'#247d76','D':INK}
STATE_COLORS=['#41526b','#d5ab5f','#8fc4dc','#eef0f2','#ffffff']
STATE_HATCH=['','///','..','','xx']
STATE_NAMES=['Positive in both','Y18 only','M25 only','Neither positive','Y18 untested']
W=180

def apply_figure_style():
    # Role-mapped style adapted from figure-style; fixed canvas size.
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':9,'axes.labelsize':9,'axes.titlesize':9,'legend.fontsize':8,
        'xtick.labelsize':8,'ytick.labelsize':8,'axes.linewidth':.6,
        'xtick.direction':'out','ytick.direction':'out','xtick.major.size':2.5,
        'ytick.major.size':2.5,'xtick.major.width':.6,'ytick.major.width':.6,
        'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,
        'svg.fonttype':'path','savefig.dpi':300,'savefig.bbox':None,
        'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,
        'hatch.linewidth':.35})
apply_figure_style()

def canvas(height):
    f=plt.figure(figsize=(W/25.4,height/25.4),dpi=150)
    bg=f.add_axes([0,0,1,1]); bg.set(xlim=(0,W),ylim=(height,0));bg.axis('off')
    f.mm_height=height
    return f,bg

def axes(f,x,y,w,h):
    return f.add_axes([x/W,1-(y+h)/f.mm_height,w/W,h/f.mm_height])

def tx(ax,x,y,s,size=9,**kw):
    return ax.text(x,y,s,fontsize=size,va='center',**kw)

def head(bg,x,y,letter,title):
    tx(bg,x,y,letter,11,weight='bold');tx(bg,x+6,y,title)

def clean_axes(ax):
    ax.spines[['top','right']].set_visible(False);ax.tick_params(length=2.5,width=.6)

def source_legend(bg,x,y):
    bg.plot(x,y,marker='o',ms=4,color=INK);tx(bg,x+3,y,'Y18',8)
    bg.plot(x+15,y,marker='s',ms=4,mfc='white',mec=INK,mew=.9);tx(bg,x+18,y,'M25',8)

def arrow(bg,x1,y1,x2,y2):
    bg.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=7,lw=.7,color=MUTED,shrinkA=0,shrinkB=0))

def swatch(bg,x,y,state,w=3.5,h=2.7):
    bg.add_patch(Rectangle((x,y-h/2),w,h,facecolor=STATE_COLORS[state],edgecolor=MUTED,lw=.35,hatch=STATE_HATCH[state]))

def load_data():
    s3=pd.read_csv(INPUT/'S3_ReactionUnits600.csv')
    s3['y_pos']=s3.y_status.eq('TESTED_ACTIVE');s3['m_pos']=s3.m_status.eq('TESTED_ACTIVE')
    s3['tested']=s3.y_status.isin(['TESTED_ACTIVE','TESTED_INACTIVE']) & s3.m_status.isin(['TESTED_ACTIVE','TESTED_INACTIVE'])
    s3['state']=np.select([~s3.tested,s3.y_pos&s3.m_pos,s3.y_pos&~s3.m_pos,~s3.y_pos&s3.m_pos],[4,0,1,2],default=3)
    focal=['Caffeic acid','Ferulic acid','Sinapic acid','Esculetin','Scopoletin','Umbelliferone']
    assert len(s3)==600 and s3.tested.sum()==578
    return s3,focal

def save(fig,stem,overlays=None):
    OUT.mkdir(exist_ok=True);fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
    for t in fig.findobj(plt.Text):
        if t.get_visible() and t.get_text().strip():
            b=t.get_window_extent(renderer)
            if b.x0 < -1 or b.y0 < -1 or b.x1 > fig.bbox.width+1 or b.y1 > fig.bbox.height+1:outside.append(t.get_text())
    if outside:raise ValueError(f'{stem}: labels outside canvas: {outside}')
    fig.savefig(OUT/f'{stem}.pdf',facecolor='white');fig.savefig(OUT/f'{stem}.svg',facecolor='white');plt.close(fig)
    if overlays:
        doc=pymupdf.open(OUT/f'{stem}.pdf');svgroot=ET.parse(OUT/f'{stem}.svg').getroot()
        for svg,box in overlays:
            x,y,w,h=box;pmm=72/25.4
            src=pymupdf.open(stream=svg.encode(),filetype='svg');vector=pymupdf.open(stream=src.convert_to_pdf(),filetype='pdf')
            doc[0].show_pdf_page(pymupdf.Rect(x*pmm,y*pmm,(x+w)*pmm,(y+h)*pmm),vector,0)
            sr=ET.fromstring(svg)
            sr.attrib.update({'x':str(x*pmm),'y':str(y*pmm),'width':str(w*pmm),'height':str(h*pmm),
                              'viewBox':'0 0 360 280','preserveAspectRatio':'xMidYMid meet'})
            svgroot.append(sr)
        tmp=OUT/f'{stem}_tmp.pdf';doc.save(tmp,garbage=4,deflate=True);doc.close();tmp.replace(OUT/f'{stem}.pdf')
        ET.ElementTree(svgroot).write(OUT/f'{stem}.svg',encoding='utf-8',xml_declaration=True)
    doc=pymupdf.open(OUT/f'{stem}.pdf')
    for dpi,ext in [(300,'png'),(1000,'tif')]:
        pix=doc[0].get_pixmap(matrix=pymupdf.Matrix(dpi/72,dpi/72),alpha=False)
        im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
        im.save(OUT/f'{stem}.{ext}',dpi=(dpi,dpi),**({'compression':'tiff_lzw'} if ext=='tif' else {}));del pix,im
    doc.close()

def protein(bg,x,y,lysate=False):
    for dx,dy,r in [(-1.6,0,1.7),(0,-1,1.8),(1.7,.5,1.5),(0,1.2,1.6)]:
        bg.add_patch(Circle((x+dx,y+dy),r,facecolor='#aeb8c3',edgecolor=INK,lw=.4))
    if lysate:
        for dx,dy,r in [(-4,-2,.5),(-4,2,.65),(4,-1.5,.55),(4,2.3,.6),(.5,-4,.45)]:
            bg.add_patch(Circle((x+dx,y+dy),r,facecolor='white',edgecolor=MUTED,lw=.45))

def vessel(bg,x,y,pool=False):
    bg.plot([x-4,x-4,x+4,x+4],[y-3,y+3,y+3,y-3],color=MUTED,lw=.7)
    if pool:
        for i,(dx,dy) in enumerate([(-2,-1),(0,1),(2,-1.3),(-2,1.5),(1.8,1.5)]):
            bg.plot(x+dx,y+dy,marker=['o','s','^'][i%3],ms=2.3,color=INK,mfc='white')
    else:bg.plot(x,y,marker='o',ms=3,color=INK,mfc='white')

def spectrum(bg,x,y,msms=False):
    bg.plot([x-5,x+5],[y+3,y+3],lw=.55,color=MUTED)
    for dx,hh in [(-3,2),(-1,5), (1,3),(3,4 if msms else 1.5)]:bg.plot([x+dx,x+dx],[y+3,y+3-hh],lw=.8,color=INK)

def fig1(s3,focal):
    acc=focal+['Baicalein','Fisetin','Kaempferol','Quercetin','Genistein','β-Sitosterol','Dihydrojasmonic acid','Gibberellin A3','Kinetin']
    enz=sorted(s3.enzyme_id.unique());f,bg=canvas(214)
    head(bg,3,5,'A','Distinct assays mapped through enzyme and acceptor identity')
    for y,source in [(19,'Y18'),(41,'M25')]:
        tx(bg,3,y,source,9,weight='bold');protein(bg,29,y,source=='M25');vessel(bg,64,y,source=='M25');spectrum(bg,99,y,source=='M25')
        arrow(bg,38,y,53,y);arrow(bg,73,y,88,y)
        tx(bg,29,y+7,'Purified enzyme' if source=='Y18' else r'$\it{E.\ coli}$ lysate',8,ha='center')
        tx(bg,64,y+7,'Single acceptor' if source=='Y18' else '40-acceptor pool',8,ha='center')
        tx(bg,99,y+7,'MS / GAR score' if source=='Y18' else 'LC–MS/MS',8,ha='center')
    bg.plot([109,115,115,109],[19,19,41,41],color=MUTED,lw=.7);arrow(bg,115,30,124,30)
    tx(bg,151,15,'Enzyme sequence',8,ha='center');tx(bg,151,22,'Acceptor connectivity',8,ha='center')
    bg.plot([128,175],[27,27],color=GRID,lw=.7)
    tx(bg,151,34,'40 UGTs × 15 acceptors',9,ha='center',weight='bold');tx(bg,151,41,'578 tested pairs',8,ha='center')
    tx(bg,151,48,'22 Y18 untested',8,ha='center',color=MUTED)
    head(bg,3,59,'B','Positive fractions across all 15 acceptors');source_legend(bg,145,59)
    dat=s3[s3.tested].groupby('acceptor_label').agg(y=('y_pos','mean'),m=('m_pos','mean')).reindex(acc)
    for start,stop,xlabel,xplot in [(0,6,3,38),(6,11,64,99),(11,15,125,160)]:
        aa=acc[start:stop];yy=np.linspace(0,5,len(aa));ax=axes(f,xplot,68,18,41);ax.set(xlim=(-.07,1.07),ylim=(5.6,-.6))
        for y,a in zip(yy,aa):
            r=dat.loc[a];i=acc.index(a);ax.plot([r.y,r.m],[y-.09,y+.09],color='#aeb7bf',lw=.8)
            ax.plot(r.y,y-.09,marker='o',ms=3.4,color=INK);ax.plot(r.m,y+.09,marker='s',ms=3.4,mfc='white',mec=INK,mew=.85)
            y_mm=68+(y+.6)/6.2*41;color=COLORS['HCA' if i<3 else 'Coumarins'] if i<6 else INK
            tx(bg,xlabel,y_mm,f'{i+1:02d}',8,color=color)
            tx(bg,xlabel+5,y_mm,a.replace('Dihydrojasmonic acid','Dihydrojasmonic\nacid'),8,color=color,linespacing=1.05)
        ax.set_xticks([0,.5,1],['0','0.5','1']);ax.set_yticks([]);ax.spines[['left','top','right']].set_visible(False);ax.tick_params(axis='x',pad=2)
    head(bg,3,123,'C','Paired reaction states across the 40 UGTs')
    matrix=s3.pivot(index='enzyme_id',columns='acceptor_label',values='state').reindex(index=enz,columns=acc)
    for block,(labelx,mx) in enumerate([(3,27),(94,118)]):
        cw=3.95;rh=3.1;top=139
        tx(bg,mx+1.5*cw,131,'HCA',8,ha='center',color=COLORS['HCA']);tx(bg,mx+4.5*cw,131,'Coumarins',8,ha='center',color=COLORS['Coumarins'])
        tx(bg,mx+10.5*cw,131,'Other acceptors',8,ha='center',color=MUTED)
        for j in range(15):tx(bg,mx+(j+.5)*cw,136,str(j+1),8,ha='center',color=COLORS['HCA' if j<3 else 'Coumarins'] if j<6 else MUTED)
        for i,e in enumerate(enz[block*20:(block+1)*20]):
            tx(bg,labelx,top+(i+.5)*rh,e,8,weight='bold' if e=='UGT84A2' else 'normal')
            for j,state in enumerate(matrix.loc[e]):
                bg.add_patch(Rectangle((mx+j*cw,top+i*rh),cw,rh,facecolor=STATE_COLORS[int(state)],edgecolor='white',lw=.35))
                if state in (1,2,4):
                    cx=mx+(j+.5)*cw;cy=top+(i+.5)*rh
                    if state==1:bg.plot([cx-.55,cx+.55],[cy+.65,cy-.65],lw=.5,color='#71552c')
                    elif state==2:bg.plot(cx,cy,marker='.',ms=2.3,color='#3b7189')
                    else:
                        bg.plot([cx-.55,cx+.55],[cy-.55,cy+.55],color='#7c8389',lw=.45);bg.plot([cx-.55,cx+.55],[cy+.55,cy-.55],color='#7c8389',lw=.45)
        for j in (0,3,6,15):bg.plot([mx+j*cw,mx+j*cw],[top,top+20*rh],color='white',lw=1)
    for x,i in zip([3,41,72,103,142],range(5)):swatch(bg,x,208,i);tx(bg,x+5,208,STATE_NAMES[i],8)
    save(f,'Figure_1')

def fig2(s3,focal):
    f,bg=canvas(177);head(bg,3,5,'A','Changes within each chemical series');head(bg,94,5,'B','Positive counts and reaction identity');source_legend(bg,45,14)
    dat=s3[s3.tested].groupby('acceptor_label').agg(y=('y_pos','mean'),m=('m_pos','mean')).reindex(focal)
    ax=axes(f,35,23,47,58);ys=[0,1,2,4,5,6]
    for i,(a,y) in enumerate(zip(focal,ys)):
        r=dat.loc[a];ax.plot([r.y,r.m],[y-.1,y+.1],lw=1,color='#aeb7bf');ax.plot(r.y,y-.1,marker='o',ms=4.3,color=INK)
        ax.plot(r.m,y+.1,marker='s',ms=4.3,mfc='white',mec=INK,mew=1)
    ax.set(xlim=(-.06,1.06),ylim=(6.6,-.6),xlabel='Positive fraction');ax.set_yticks(ys,focal);ax.set_xticks([0,.5,1],['0','0.5','1'])
    for i,t in enumerate(ax.get_yticklabels()):t.set_color(COLORS['HCA' if i<3 else 'Coumarins'])
    ax.axhline(3,color=GRID,lw=.7);clean_axes(ax)
    for x,state,label in [(98,0,'Both'),(122,1,'Y18 only'),(150,2,'M25 only')]:swatch(bg,x,14,state);tx(bg,x+5,14,label,8)
    ax=axes(f,98,22,78,59);ax.set(xlim=(0,120),ylim=(2.10,-.85))
    for i,s in enumerate(['HCA','Coumarins']):
        q=s3[s3.tested&s3.focal_series.eq(s)];left=0; row_y=i*1.5
        for state in (0,1,2):
            n=int(q.state.eq(state).sum())
            if n:
                ax.barh(row_y,n,left=left,height=.36,color=STATE_COLORS[state],edgecolor=MUTED,hatch=STATE_HATCH[state],lw=.4)
                if n>=10:ax.text(left+n/2,row_y,str(n),ha='center',va='center',fontsize=9,color='white' if state==0 else INK,bbox={'facecolor':STATE_COLORS[state],'edgecolor':'none','pad':.7})
                else:ax.annotate(str(n),(left+n/2,row_y-.18),xytext=(left+n/2+4,row_y-.4),fontsize=8,ha='center',arrowprops={'arrowstyle':'-','lw':.5,'color':MUTED})
            left+=n
        y=int(q.y_pos.sum());m=int(q.m_pos.sum());neither=int(q.state.eq(3).sum())
        ax.text(0,row_y-.60,f'{s}: {y} → {m} positive calls',color=COLORS[s],fontsize=9)
        ax.text(0,row_y+.45,f'Neither positive: {neither}',fontsize=8,color=MUTED)
    ax.set_yticks([]);ax.spines[['left','right','top']].set_visible(False);ax.set_xticks([0,40,80,120]);ax.set_xlabel('Pairs positive in either screen')
    head(bg,3,103,'C','Enzyme-weighted series contrasts');head(bg,94,103,'D','Positive-enzyme overlap')
    boot=pd.read_csv(GEN/'Fig2_effects.csv');rows=pd.concat([boot.iloc[:3],boot.iloc[5:6]]);ax=axes(f,35,115,47,45)
    for i,(_,r) in enumerate(rows.iterrows()):
        color=COLORS.get(r.metric_id,COLORS['Coumarins']);yy=[0,1,2.3,3.7][i]
        ax.plot([r.lower,r.upper],[yy,yy],color=color,lw=1.5 if i<3 else 1);ax.plot(r.estimate,yy,marker='D',ms=4.2,mfc=color if i<3 else 'white',mec=color,mew=1)
    ax.set_yticks([0,1,2.3,3.7],['ΔHCA','Δcoumarin','D · 39 UGTs','D · 30 exact']);ax.set(ylim=(4.3,-.6),xlim=(-.68,.25),xlabel='Change in positive fraction')
    ax.set_xticks([-.6,-.3,0],['−0.6','−0.3','0']);ax.axvline(0,color=GRID,lw=.8);ax.axhline(3,color=GRID,lw=.7);clean_axes(ax)
    rows=[]
    for a,q in s3[s3.tested].groupby('acceptor_label'):
        ys1=set(q.loc[q.y_pos,'enzyme_id']);ms=set(q.loc[q.m_pos,'enzyme_id']);rows.append((a,len(ys1&ms)/len(ys1|ms) if ys1|ms else np.nan,q.focal_series.iloc[0]))
    jd=pd.DataFrame(rows,columns=['acceptor','jaccard','series']);jd.to_csv(GEN/'acceptor_jaccard.csv',index=False);ax=axes(f,133,115,43,45)
    for i,a in enumerate(focal):
        r=jd.set_index('acceptor').loc[a];color=COLORS[r.series];ax.plot([0,r.jaccard],[ys[i],ys[i]],color=color,lw=1);ax.plot(r.jaccard,ys[i],marker='D',ms=3.7,color=color)
    ax.set_yticks(ys,focal);ax.set(xlim=(-.04,1.04),ylim=(6.6,-.6),xlabel='Jaccard index');ax.set_xticks([0,.5,1],['0','0.5','1'])
    ax.axhline(3,color=GRID,lw=.7);clean_axes(ax);save(f,'Figure_2')

def molecular_svg(smiles,expected_key):
    # Explicit, reviewed scaffold/substitution maps for the six recorded graphs.
    # Guard against drawing a different input; no stereochemical or product-site assignment.
    specifications={
        'QAIPRVGONGVQAS':('O=C(O)C=CC1=CC=C(O)C(O)=C1','HCA',{2:'OH',3:'OH'}),
        'KSEBMYQBYZTDHS':('O=C(O)C=CC1=CC=C(O)C(OC)=C1','HCA',{2:'MeO',3:'OH'}),
        'PCMORTLOPMLEFB':('COC1=CC(C=CC(O)=O)=CC(OC)=C1O','HCA',{2:'MeO',3:'OH',4:'MeO'}),
        'ILEDWLMCKZNDJK':('O=C1C=CC2=C(O1)C=C(O)C(O)=C2','Coumarins',{2:'OH',3:'OH'}),
        'RODXRVNMMDRFIK':('O=C1OC2=C(C=C(C(O)=C2)OC)C=C1','Coumarins',{2:'MeO',3:'OH'}),
        'ORHBXUUXSCNDEV':('O=C1C=CC2=C(O1)C=C(O)C=C2','Coumarins',{3:'OH'}),
    }
    recorded,kind,groups=specifications[expected_key];assert smiles==recorded
    elements=[]
    def line(a,b):elements.append(f'<path d="M {a[0]:.2f},{a[1]:.2f} L {b[0]:.2f},{b[1]:.2f}" fill="none" stroke="#20252b" stroke-width="2.6"/>')
    def atom(x,y,label,anchor='middle'):
        # Vector glyphs prevent unembedded fallback fonts in molecule overlays.
        glyph=TextPath((0,0),label,size=30,prop=FontProperties(family='Arial'))
        b=glyph.get_extents()
        left=x-b.width/2 if anchor=='middle' else (x-b.width if anchor=='end' else x)
        shape=glyph.transformed(Affine2D().scale(1,-1).translate(left-b.x0,y+(b.y0+b.y1)/2))
        bounds=shape.get_extents()
        elements.append(f'<rect x="{bounds.x0-2}" y="{bounds.y0-2}" width="{bounds.width+4}" height="{bounds.height+4}" fill="white"/>')
        parts=[]
        for coords,code in shape.iter_segments(curves=True,simplify=False):
            if code==1:parts.append('M '+' '.join(f'{v:.3f}' for v in coords))
            elif code==2:parts.append('L '+' '.join(f'{v:.3f}' for v in coords))
            elif code==3:parts.append('Q '+' '.join(f'{v:.3f}' for v in coords))
            elif code==4:parts.append('C '+' '.join(f'{v:.3f}' for v in coords))
            elif code==79:parts.append('Z')
        elements.append('<path d="'+' '.join(parts)+'" fill="#20252b"/>')
    if kind=='HCA':
        center=np.array([132.,140.]);ang=np.arange(6)*np.pi/3
        ring=[center+47*np.array([np.cos(a),np.sin(a)]) for a in ang]
        for i in range(6):line(ring[i],ring[(i+1)%6])
        for i in (0,2,4):line(center+.81*(ring[i]-center),center+.81*(ring[(i+1)%6]-center))
        chain=[ring[0],np.array([205.,124.]),np.array([231.,140.]),np.array([257.,124.])]
        for a,b in zip(chain,chain[1:]):line(a,b)
        line(chain[1]+[0,6],chain[2]+[0,6]);line([254,121],[254,91]);line([260,121],[260,91]);atom(257,81,'O')
        line(chain[-1],[281,138]);atom(285,142,'OH','start')
    else:
        center=np.array([140.,140.]);ang=np.pi/6+np.arange(6)*np.pi/3
        ring=[center+53*np.array([np.cos(a),np.sin(a)]) for a in ang]
        for i in range(6):line(ring[i],ring[(i+1)%6])
        for i in (1,3,5):line(center+.82*(ring[i]-center),center+.82*(ring[(i+1)%6]-center))
        oxygen=np.array([231.8,87.]);carbonyl=np.array([277.7,113.5]);c3=np.array([277.7,166.5]);c4=np.array([231.8,193.])
        chain=[ring[5],oxygen,carbonyl,c3,c4,ring[0]]
        for a,b in zip(chain,chain[1:]):line(a,b)
        line(c3+[-5,-3],c4+[-5,-3]);line(carbonyl+[0,-3],[307,93]);line(carbonyl+[3,3],[310,99]);atom(325,86,'O')
        atom(*oxygen,'O')
    for i,label in groups.items():
        a=ring[i];v=(a-center)/np.linalg.norm(a-center);b=a+v*26;line(a,b)
        atom(b[0]-4,b[1]+3,'HO' if label=='OH' else label,'end')
    return '<svg xmlns="http://www.w3.org/2000/svg" width="360" height="280" viewBox="0 0 360 280">'+''.join(elements)+'</svg>'

def fig3(s3,focal):
    f,bg=canvas(212);overlays=[];head(bg,3,5,'A','UGT84A2 retains coumarin reactions as hydroxycinnamate calls decline')
    tx(bg,45,19,'Hydroxycinnamates',9,weight='bold',ha='center',color=COLORS['HCA']);tx(bg,136,19,'Coumarins',9,weight='bold',ha='center',color=COLORS['Coumarins'])
    bg.plot([4,86],[24,24],color=COLORS['HCA'],lw=1.4);bg.plot([95,177],[24,24],color=COLORS['Coumarins'],lw=1.4)
    q=s3[s3.enzyme_id.eq('UGT84A2')].set_index('acceptor_label').reindex(focal);identity=pd.read_csv(INPUT/'S2_AcceptorIdentity.csv').set_index('key14');assert q.tested.all()
    for j,a in enumerate(focal):
        xx=4 if j<3 else 95; yy=27+(j%3)*28; key=q.loc[a,'key14']
        svg=molecular_svg(identity.loc[key,'y_smiles'],key);overlays.append((svg,(xx,yy,39,28)))
        tx(bg,xx+46,yy+8,a,9);state=int(q.loc[a,'state'])
        swatch(bg,xx+46,yy+17,state,w=4,h=3)
        tx(bg,xx+52,yy+17,'Both screens' if state==0 else 'Y18 only',8)
    tx(bg,45,116,'All three positive only in Y18',8,ha='center');tx(bg,136,116,'Two remain positive in both screens',8,ha='center')
    head(bg,3,130,'B','Series profiles distinguish selective and broad contraction');source_legend(bg,141,139)
    tx(bg,77,148,'Hydroxycinnamates',9,ha='center',color=COLORS['HCA']);tx(bg,151,148,'Coumarins',9,ha='center',color=COLORS['Coumarins'])
    enzymes=['UGT84A2','UGT84A1','UGT76E12','UGT73C5','UGT72D1'];y0=156;hh=40;ylimits=(4.6,-.6)
    highlight_y=y0+.6/5.2*hh;bg.add_patch(Rectangle((3,highlight_y-3.3),174,6.6,facecolor='#f0f2f4',edgecolor='none'))
    for i,e in enumerate(enzymes):tx(bg,6,y0+(i+.6)/5.2*hh,e,9,weight='bold' if e=='UGT84A2' else 'normal')
    for x,series in [(53,'HCA'),(127,'Coumarins')]:
        ax=axes(f,x,y0,48,hh);ax.set_facecolor('none')
        for i,e in enumerate(enzymes):
            z=s3[s3.enzyme_id.eq(e)&s3.tested&s3.focal_series.eq(series)];assert len(z)==3
            y=int(z.y_pos.sum());m=int(z.m_pos.sum());ax.plot([y,m],[i-.11,i+.11],lw=1.2 if i==0 else .9,color='#aab3bc')
            ax.plot(y,i-.11,'o',ms=4.6 if i==0 else 4.1,color=INK);ax.plot(m,i+.11,'s',ms=4.6 if i==0 else 4.1,mfc='white',mec=INK,mew=1)
        ax.set(xlim=(-.16,3.16),ylim=ylimits);ax.set_xticks([0,1,2,3]);ax.set_yticks([]);ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='x',pad=3)
    tx(bg,111,207,'Positive acceptors (of three tested per series)',9,ha='center');save(f,'Figure_3',overlays)

def main():
    s3,focal=load_data();fig1(s3,focal);fig2(s3,focal);fig3(s3,focal)
    (GEN/'figure_visual_encoding.json').write_text(json.dumps({'width_mm':180,'font':'Arial with DejaVu Sans fallback','font_roles_pt':[9,8],
        'sources':{'Y18':'black filled circle','M25':'black open square'},'series':COLORS,'reaction_states':dict(zip(STATE_NAMES,STATE_COLORS)),
        'missing':'crossed white cell, excluded from binary denominators','molecules':'Explicit scaffold/substitution depictions, guarded by recorded SMILES and connectivity keys',
        'minimum_intended_width_mm':169,'single_column':'not intended; use full-width figures'},indent=2),encoding='utf-8')
    print('Built three 180-mm vector figures, with 300-dpi PNG and 1000-dpi TIFF.')

if __name__=='__main__':main()
