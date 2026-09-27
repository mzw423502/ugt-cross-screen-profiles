"""Build readable Supplementary Information tables and figures."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.lines import Line2D
from build_figures import apply_figure_style, COLORS, INK

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "05_REPRODUCIBILITY" / "input"
OUT = ROOT / "03_SUPPLEMENT"
GEN = ROOT / "05_REPRODUCIBILITY" / "generated"
TABLES = []

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.linewidth": 0.7, "xtick.direction": "out", "ytick.direction": "out",
    "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
})

apply_figure_style()

def clean_axes(ax):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.tick_params(length=3)

def savefig(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", facecolor="white", bbox_inches="tight")
    Image.open(OUT / f"{stem}.png").convert("RGB").save(OUT / f"{stem}.tif", dpi=(300,300), compression="tiff_lzw")
    plt.close(fig)

def md_table(headers, rows):
    TABLES.append({'headers':headers,'rows':rows})
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    out.extend("| " + " | ".join(str(x).replace('|', '\\|') for x in row) + " |" for row in rows)
    return "\n".join(out)

def fmt(x, digits=3):
    if pd.isna(x): return "—"
    if isinstance(x, (float, np.floating)): return f"{x:.{digits}f}"
    return str(x)

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s1 = pd.read_csv(INPUT / "S1_EnzymeIdentity.csv")
    s2 = pd.read_csv(INPUT / "S2_AcceptorIdentity.csv")
    s3 = pd.read_csv(INPUT / "S3_ReactionUnits600.csv")
    acc = pd.read_csv(INPUT / "acceptor_paired_summary_085.csv")
    cls = pd.read_csv(INPUT / "chemical_class_paired_summary_085.csv")
    enz = pd.read_csv(GEN / "per_enzyme_series.csv")
    loo = pd.read_csv(INPUT / "leave_one_compound_out_series_sensitivity.csv")
    boot = pd.read_csv(GEN / "bootstrap_series_final.csv")
    exact = pd.read_csv(INPUT / "exact31_series_summary.csv")
    eff = pd.read_csv(INPUT / "effect_sizes_family_cluster_bootstrap.csv")
    th = pd.read_csv(INPUT / "threshold_sensitivity_identity_fixed.csv")
    evidence = pd.read_csv(INPUT / "biochemical_evidence_table.csv")
    prov = pd.read_csv(INPUT / "source_manifest.csv")

    jacc = []
    for _, row in acc.iterrows():
        r = s3[(s3.acceptor_label == row.acceptor_label) &
               s3.y_status.isin(["TESTED_ACTIVE","TESTED_INACTIVE"]) &
               s3.m_status.isin(["TESTED_ACTIVE","TESTED_INACTIVE"])]
        ys = set(r.loc[r.y_status.eq("TESTED_ACTIVE"), "enzyme_id"])
        ms = set(r.loc[r.m_status.eq("TESTED_ACTIVE"), "enzyme_id"])
        union = ys | ms
        jacc.append(len(ys & ms) / len(union) if union else np.nan)
    acc["jaccard"] = jacc

    lines = [
        "# Supplementary information",
        "",
        "## Supplementary methods",
        "",
        "This supplement presents sequence matching, connectivity-level acceptor matching, reaction summaries, sensitivity analyses and published biochemical context for the Arabidopsis comparison. Y18 denotes Yang et al. (2018), and M25 denotes Sirirungruang et al. (2025). Tables S1–S16 and Figures S1–S3 accompany Supplementary_Data.xlsx, which also contains the complete reaction matrix and full-precision analysis results.",
        "",
        "The comparison contains 40 common UGT constructs, 15 connectivity-matched acceptors and 600 enzyme–acceptor units. Twenty-two Y18 units are NOT_TESTED; the binary comparison therefore contains 578 tested pairs. Positive-in-both, Y18-only, M25-only and no-qualifying-product states are mutually exclusive; NOT_TESTED remains separate.",
        "",
        "For enzyme e and series s, Δe,s = (M25 positives − Y18 positives)/number of acceptors tested in both screens. The primary estimates average these contrasts over the 39 enzymes with at least two tested members of both focal series; D = mean ΔHCA − mean Δcoumarin. The exact-sequence analysis applies the same coverage rule to 31 exact names, retaining 30 enzymes.",
        "",
        "Uncertainty uses 500 shared family-prefix cluster-bootstrap draws per analysis set (seed 20260919; NumPy default_rng). The prefix includes the UGT number and following letter, such as UGT84A. The primary and exact sets contain 17 and 13 groups. Each draw samples the original number of groups with replacement, retaining all enzymes with cluster multiplicity. Both series and D use the same draw. Intervals are the 2.5th and 97.5th percentiles with linear interpolation. These intervals describe variation among enzyme groups; the screening calls do not provide replicate-level measurement uncertainty.",
        "",
        "For each acceptor, Y and M are the positive-enzyme sets restricted to enzymes tested in both sources. Jaccard overlap is |Y ∩ M| / |Y ∪ M|; an empty union is undefined, not zero. Agreement is the number of concordant positive or negative calls divided by tested pairs. Compound leave-one-out summaries pool the remaining tested pairs within a series. Threshold comparisons hold enzyme identity, acceptor identity and Y18 tested status fixed.",
        "",
        "## Supplementary Table S1. Enzyme construct matching",
        "",
        "The table gives the sequence evidence used to link the two screen records. Exact matches have identical sequence strings; high-identity matches meet at least 98.5% aligned identity and 95% coverage.",
        "",
    ]
    rows=[]
    for _,r in s1.iterrows():
        match = "exact" if r.y18_sequence_sha256 == r.m25_sequence_sha256 else "high identity"
        rows.append([r.enzyme_id, match, fmt(r.aligned_residue_identity*100,2)+"%", fmt(r.aligned_core_coverage*100,2)+"%", f'{r.y18_length} / {r.m25_length}'])
    lines += [md_table(["Enzyme","Match class","Identity","Coverage","Y18 / M25 length"],rows),"",
              "Sequence hashes are retained in S1_EnzymeIdentity.csv as provenance fields.",
              "",
              "## Supplementary Table S2. Acceptor identity and connectivity matching",
              "",
              "The 14-character InChIKey block defines the matched molecular connectivity. Stereochemical equivalence remains unresolved at this level; the excluded collision is listed explicitly.",
              ""]
    rows=[]
    for _,r in s2.iterrows():
        label=s3.loc[s3.key14.eq(r.key14),'acceptor_label']
        rows.append([label.iloc[0] if len(label) else 'Catechin / epicatechin / cianidanol', str(r.key14), str(r.m_names), "included" if bool(r.identity_include) else "excluded: stereoisomer collision"])
    lines += [md_table(["Acceptor","Connectivity key","M25 name","Decision"],rows),"",
              "## Supplementary Table S3. Reaction-state summary",
              "",
              "The 600-unit matrix is supplied as S3_ReactionUnits600.csv. The summary below reports the 0.85 M25 layer after enzyme and acceptor matching.",
              ""]
    tested=s3[s3.comparable_status=="comparable_tested"]
    state_counts={
        "positive in both": int(((tested.y_status=="TESTED_ACTIVE")&(tested.m_status=="TESTED_ACTIVE")).sum()),
        "Y18-only": int(((tested.y_status=="TESTED_ACTIVE")&(tested.m_status=="TESTED_INACTIVE")).sum()),
        "M25-only": int(((tested.y_status=="TESTED_INACTIVE")&(tested.m_status=="TESTED_ACTIVE")).sum()),
        "no qualifying product in either": int(((tested.y_status=="TESTED_INACTIVE")&(tested.m_status=="TESTED_INACTIVE")).sum()),
    }
    rows=[[k,v,f"{v/len(tested):.1%}"] for k,v in state_counts.items()]
    lines += [md_table(["Reaction state","Units","Fraction of tested pairs"],rows),"",
              "Y18 NOT_TESTED units (n=22) remain in the machine-readable matrix and are excluded from the tested-pair denominator.",
              "",
              "## Supplementary Table S4. Complete 15-acceptor landscape",
              "",
              "This table is the full acceptor context used to identify the biologically interpretable HCA–coumarin contrast. Positive fractions use each acceptor's tested denominator; Jaccard overlap is the intersection divided by the union of positive enzyme sets.",
              ""]
    rows=[]
    for _,r in acc.sort_values(["series","acceptor_label"]).iterrows():
        rows.append([r.acceptor_label,r.series,int(r.n_comparable),int(r.y18_active_n),int(r.m25_active_n),fmt(r.agreement_fraction),fmt(r.jaccard)])
    lines += [md_table(["Acceptor","Chemical series","Tested pairs","Y18+","M25+","Agreement","Jaccard"],rows),"",
              "## Supplementary Table S5. Chemical-series summaries","",
              "Class-level summaries retain the tested denominator for each chemical class.",""]
    rows=[]
    for _,r in cls.sort_values("series").iterrows():
        rows.append([r.series,int(r.acceptor_n),int(r.n_comparable),int(r.y18_active_n),int(r.m25_active_n),fmt(r.y18_fraction),fmt(r.m25_fraction),fmt(r.agreement_fraction)])
    lines += [md_table(["Series","Acceptors","Tested pairs","Y18+","M25+","Y18 fraction","M25 fraction","Agreement"],rows),"",
              "## Supplementary Table S6. Enzyme-level series contrasts","",
              "All enzymes except UGT71D1 meet the two-series coverage rule. Fractions are positive calls divided by tested acceptors in each series; D = delta HCA minus delta coumarin.",""]
    rows=[]
    for _,r in enz.iterrows():
        rows.append([r.enzyme_id,f"{int(r.y_HCA)}/{int(r.n_HCA)} → {int(r.m_HCA)}/{int(r.n_HCA)}" if r.n_HCA else 'Not tested',fmt(r.delta_HCA),f"{int(r.y_Coumarins)}/{int(r.n_Coumarins)} → {int(r.m_Coumarins)}/{int(r.n_Coumarins)}",fmt(r.delta_Coumarins),fmt(r.D_HCA_minus_Coumarins)])
    lines += [md_table(["Enzyme","HCA Y18 → M25","ΔHCA","Coumarin Y18 → M25","Δcoumarin","D"],rows),"",
              "Arrows join counts over the same tested denominator. UGT71D1 is the coverage-limited row and has no D estimate.","",
              "## Supplementary Table S7. Compound leave-one-out context","",
              "Each row removes one acceptor from its focal series and recomputes the pooled series fractions.",""]
    rows=[]
    for _,r in loo.iterrows():
        label=s3.loc[s3.key14.eq(r.left_out_key),'acceptor_label'].iloc[0]
        rows.append([r.series,label,int(r.n),fmt(r.y18_fraction),fmt(r.m25_fraction),fmt(r.delta)])
    lines += [md_table(["Series","Omitted acceptor","Pairs","Y18 fraction","M25 fraction","Δ"],rows),"",
              "## Supplementary Table S8. Exact-sequence sensitivity","",
              "The exact-sequence set contains 31 names and 30 enzymes meeting the series coverage rule. This table pools all available exact-sequence reactions, including the coverage-limited construct in the coumarin denominator. Equal-enzyme contrasts in Table S9 use only the 30 enzymes meeting the coverage rule.",""]
    rows=[]
    for _,r in exact.iterrows():
        rows.append([r.focal_series,int(r.n_comparable),int(r.enzyme_n),int(r.y18_active_n),int(r.m25_active_n),fmt(r.y18_active_fraction),fmt(r.m25_active_fraction),fmt(r.delta_m25_minus_y18)])
    lines += [md_table(["Series","Tested pairs","Enzymes","Y18+","M25+","Y18 fraction","M25 fraction","delta"],rows),"",
              "## Supplementary Table S9. Shared cluster bootstrap","",
              "The same family-prefix resample supplies delta HCA, delta coumarin and D. Complete resample rows are available in generated/bootstrap_full.csv.",""]
    rows=[]
    for _,r in boot.iterrows():
        rows.append([r.metric,int(r.n_enzymes),int(r.n_clusters),fmt(r.raw_point_estimate),fmt(r.bootstrap_low),fmt(r.bootstrap_high),fmt(r.bootstrap_mean)])
    lines += [md_table(["Metric","Enzymes","Clusters","Estimate","Lower","Upper","Bootstrap mean"],rows),"",
              "## Supplementary Figure S1. Threshold sensitivity","",
              "The four M25 product-calling thresholds are evaluated with fixed enzyme matching and fixed tested denominators.",
              "![Supplementary Figure S1](Supplementary_Figure_S1.png)",
              "Supplementary Figure S1. Agreement and M25 positive calls across the four cosine thresholds. Each threshold uses the same 578 tested pairs. Visualization code assisted by OpenAI Codex.","",
              "## Supplementary Table S10. M25 threshold sensitivity",""]
    rows=[]
    for _,r in th.iterrows():
        rows.append([f"{r.threshold:.2f}",int(r.n_identity_included_tested),int(r.y18_active_n),int(r.m25_active_n),int(r.agreement_n),fmt(r.agreement_fraction)])
    lines += [md_table(["M25 threshold","Tested pairs","Y18+","M25+","Agreement calls","Agreement fraction"],rows),"",
              "## Supplementary Figure S2. Exact-sequence sensitivity","",
              "The exact-sequence set preserves the direction of the two-series contrast.",
              "![Supplementary Figure S2](Supplementary_Figure_S2.png)",
              "Supplementary Figure S2. Estimates and 95% family-prefix cluster intervals for HCA, coumarins and D. Filled diamonds denote the primary 39-enzyme set; open diamonds denote the 30 exact-sequence enzymes meeting the coverage rule. Each analysis uses 500 shared resamples. Visualization code assisted by OpenAI Codex.","",
              "## Supplementary Table S11. Product-row aggregation","",
              md_table(["Step","Definition"],[
                  ["Unit","One enzyme × connectivity-matched acceptor"],
                  ["Threshold","Use the product table corresponding to the selected cosine-similarity threshold"],
                  ["Aggregation","Group qualifying single- and double-glycosylated entries by enzyme, acceptor name and canonical SMILES; take the maximum AUC"],
                  ["M25 call","Positive if the maximum AUC exceeds zero; otherwise the screened unit has no qualifying positive-area product"],
              ]),"",
              "## Supplementary Table S12. Source provenance","",
              "The source manifest preserves DOI, file names and SHA-256 values for the author-released datasets.",""]
    rows=[['Yang et al. (2018)','10.1038/s41589-018-0154-9','Article, supplementary acceptor GAR matrix and sequence records'],['Sirirungruang et al. (2025)','10.1038/s41467-025-61530-6','Article, Source Data and author-released threshold product tables']]
    lines += [md_table(["Source","DOI","Records"],rows),"",
              "## Supplementary Figure S3. Enzyme coverage across the focal series","",
              "Coverage is shown for all 40 common enzyme constructs.",
              "![Supplementary Figure S3](Supplementary_Figure_S3.png)",
              "Supplementary Figure S3. Number of tested HCA and coumarin acceptors among the 40 common enzymes. Enzymes sharing a coverage combination are aggregated; marker area is proportional to the labelled count. Filled markers represent the 39 enzymes meeting the coverage rule; the open marker represents UGT71D1. Visualization code assisted by OpenAI Codex.","",
              "## Supplementary Table S13. Focal acceptor identity",""]
    focal = json.loads((INPUT/'recompute_config.json').read_text())['focal_acceptors']
    rows=[]
    for _,r in s2.iterrows():
        if r.key14 in focal:
            rows.append([focal[r.key14]['label'],str(r.key14),str(r.y_smiles)])
    lines += [md_table(["Acceptor","Connectivity key","Canonical SMILES"],rows),"",
              "## Supplementary Table S14. Published biochemistry and cross-screen substrate profiles","",
              "Published enzyme assays and plant studies provide biochemical context for the observed profiles. Their experimental conditions differ from those of the two screens.",""]
    evidence_headers=['Enzyme','Published biochemical evidence','Cross-screen profile','Interpretation']
    rows=[
      ['UGT84A1','Hydroxycinnamate glucose-ester activity, including caffeic acid; UV-B-responsive UGT84A biology (Lim et al., 2001; Meißner et al., 2008).','12/14 positives in Y18; 2/14 in M25. All three HCA and all three coumarins are Y18-only.','Contraction spans both series.'],
      ['UGT84A2','Sinapate glucose-ester activity; UGT84A2-derived sinapoylglucose participates in anthocyanin modification (Lim et al., 2001; Yonekura-Sakakibara et al., 2012).','12/14 → 4/14 positives. All three HCA are Y18-only; esculetin and umbelliferone are positive in both.','Different series responses within an identical recorded enzyme sequence.'],
      ['UGT84A3','Hydroxycinnamate glucose-ester activity and phenylpropanoid redundancy in the UGT84A family (Meißner et al., 2008).','All three HCA are Y18-only; all three coumarins are positive in both.','Shared HCA contraction with persistent coumarin calls.'],
      ['UGT84A4','Hydroxycinnamate glucosyltransferase assays and UV-B response experiments (Meißner et al., 2008).','All three HCA are Y18-only; two Y18 coumarin-positive calls are absent in M25.','Coumarin behavior differs among UGT84A enzymes.'],
      ['UGT76E12','The two large screens provide the evidence for these six acceptors.','9/15 → 1/15 positives. Caffeic acid and all three coumarins are Y18-only; ferulic and sinapic acids are negative in both.','Contraction spans both series.'],
      ['UGT73C5','Brassinosteroid glucosylation established in planta (Poppenberger et al., 2005).','All three coumarins are positive in both screens; all three HCA are negative in both.','Stable biochemical screen profile; physiological evidence concerns brassinosteroids.'],
      ['UGT72D1','Purified His-tagged Q9ZU72 (At2g18570; 470 aa) forms sinapic-acid and coniferyl-aldehyde 4-O-glucosides; sinapic-acid conversion depends on pH (Li et al., 2024).','13/15 pairwise calls agree. Esculetin and scopoletin are positive in both; all three HCA are negative in both.','Published sinapic-acid activity coexists with negative calls under both screen conditions.'],
    ]
    lines += [md_table(evidence_headers,rows),"",
              "## Supplementary Table S15. Experimental contexts","",
              "Y18 and M25 are independent studies; the purified validation is part of M25. pH and concentration values follow the cross-study comparison in Sirirungruang et al. (2025); purified validation details follow its Source Data. Conditions co-vary and are not isolated causal effects.","",
              md_table(["Feature","Y18","M25","M25 purified validation"],[
                  ["Protein","Recombinant enzymes; GST-affinity purification","Recombinant E. coli lysate","Purified GST-fusion enzymes"],
                  ["Acceptor format","Single acceptor","Pools of 40 acceptors","Single acceptor"],
                  ["Reported pH","7.8","6.8","7.6"],
                  ["Reported acceptor concentration","0.1 mg/mL","10 µM each","50 µM"],
                  ["Reported UDP-glucose","177 µM","83 µM","UDP-glucose donor"],
                  ["Detection","Mass-spectrometric GAR score","LC–MS/MS product calls","Purified-enzyme product signal"],
                  ["Replication","Single high-throughput measurements","Product rows aggregated into units","Three technical replicates in recovered Source Data"],
                  ["Use here","Main matched source","Main matched source","Same-study context; no rows added to main matrix"],
              ]),"",
              "## Supplementary Table S16. Threshold effects with fixed denominators","",
              "Counts and equal-enzyme series contrasts are recomputed at each threshold below.","",
              "The corresponding machine-readable threshold tables, reaction matrix and identity tables are supplied in 05_REPRODUCIBILITY/input/ and Supplementary_Data.xlsx.",
    ]
    threshold_rows=[]
    for threshold in [.75,.80,.85,.90]:
        q=pd.read_csv(INPUT/f'q2_strict_overlap_{threshold:.2f}.csv')
        q=q[q.identity_include & q.comparable_status.eq('comparable_tested')].copy()
        q['y']=q.y_status.eq('TESTED_ACTIVE');q['m']=q.m_status.eq('TESTED_ACTIVE')
        pools=[];deltas=[]
        for series in ['HCA','Coumarins']:
            z=q[q.focal_series.eq(series)];pools.append(f'{int(z.y.sum())}/{len(z)} → {int(z.m.sum())}/{len(z)}')
        for _,z in q.groupby('enzyme_id'):
            groups=[z[z.focal_series.eq(s)] for s in ['HCA','Coumarins']]
            if all(len(g)>=2 for g in groups):deltas.append([g.m.mean()-g.y.mean() for g in groups])
        dh,dc=np.array(deltas).mean(axis=0)
        threshold_rows.append([f'{threshold:.2f}',*pools,fmt(dh),fmt(dc),fmt(dh-dc)])
    lines += [md_table(['Threshold','HCA Y18 → M25','Coumarin Y18 → M25','ΔHCA','Δcoumarin','D'],threshold_rows),'',
              'Each contrast averages over the same 39 enzymes meeting the coverage rule. Arrows show pooled counts over matched tested pairs.','',
              '## Supplementary references','',
              'Yang et al. (2018), Nature Chemical Biology 14, 1109–1117. https://doi.org/10.1038/s41589-018-0154-9','',
              'Sirirungruang et al. (2025), Nature Communications 16, 6366. https://doi.org/10.1038/s41467-025-61530-6','',
              'Biochemical sources for Table S14: Lim et al. (2001), https://doi.org/10.1074/jbc.M007263200; Meißner et al. (2008), https://doi.org/10.1007/s00425-008-0768-3; Yonekura-Sakakibara et al. (2012), https://doi.org/10.1111/j.1365-313X.2011.04779.x; Poppenberger et al. (2005), https://doi.org/10.1073/pnas.0504279102; Li et al. (2024), https://doi.org/10.1038/s42004-024-01231-1.','']
    # Standalone source references use the same corrected records as the main paper.
    ref_start=lines.index('## Supplementary references')
    references=(ROOT/'06_QC/REFERENCES_CURRENT.md').read_text(encoding='utf-8')
    required=['10.1038/s41589-018-0154-9','10.1038/s41467-025-61530-6','10.1074/jbc.M007263200','10.1007/s00425-008-0768-3','10.1111/j.1365-313X.2011.04779.x','10.1073/pnas.0504279102','10.1038/s42004-024-01231-1']
    lines=lines[:ref_start]+['## Supplementary references','']
    for record in references.splitlines():
        if any(doi.lower() in record.lower() for doi in required):lines.extend([record,''])
    (OUT / "Supplementary_Information.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    (GEN/'supplement_readable_tables.json').write_text(json.dumps(TABLES,ensure_ascii=False,indent=2),encoding='utf-8')

    fig, axs = plt.subplots(1,2,figsize=(6.8,2.7),constrained_layout=True)
    axs[0].plot(th.threshold,th.agreement_fraction*100,marker="o",color=INK,lw=1.5)
    axs[0].set(xlabel="M25 cosine threshold",ylabel="Agreement (%)",title="Agreement remains stable")
    axs[0].set_xticks(th.threshold); axs[0].set_ylim(68,73); clean_axes(axs[0])
    axs[1].plot(th.threshold,th.m25_active_n,marker="o",color=INK,lw=1.5)
    axs[1].set(xlabel="M25 cosine threshold",ylabel="M25 positive calls",title="Calls decline with threshold")
    axs[1].set_xticks(th.threshold); axs[1].set_ylim(160,260); clean_axes(axs[1])
    fig.suptitle("Threshold sensitivity with fixed identity and tested denominators",x=0.02,ha="left",fontsize=9)
    savefig(fig,"Supplementary_Figure_S1")

    fig, ax = plt.subplots(figsize=(6.3,2.8))
    fig.subplots_adjust(left=.18,right=.97,bottom=.25,top=.77)
    metrics=["HCA","coumarins","D"]; color={"HCA":COLORS["HCA"],"coumarins":COLORS["Coumarins"],"D":COLORS["D"]}
    full=boot[boot.scope.eq("all_39_series_adequate")].set_index("metric")
    exact_data=boot[boot.scope.eq('exact31_30_series_adequate')].set_index('metric')
    for i,metric_set in enumerate(["39 enzymes","30 exact enzymes"]):
        for j,m in enumerate(metrics):
            if i==0:
                r=full.loc[m]; val,lo,hi=r.raw_point_estimate,r.bootstrap_low,r.bootstrap_high
            else:
                r=exact_data.loc[m]; val,lo,hi=r.raw_point_estimate,r.bootstrap_low,r.bootstrap_high
            yy=j+(i-.5)*.23
            ax.plot([lo,hi],[yy,yy],color=color[m],lw=1.5)
            ax.scatter(val,yy,marker='D',facecolor=color[m] if i==0 else 'white',edgecolor=color[m],s=28,zorder=3)
    ax.axvline(0,color="#999",lw=.7)
    ax.set_yticks([0,1,2]); ax.set_yticklabels(["ΔHCA","Δcoumarin","D"],fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Change in positive fraction; 95% cluster interval")
    fig.legend([Line2D([],[],marker='D',color=INK,ls='none'),Line2D([],[],marker='D',mfc='white',color=INK,ls='none')],['39 enzymes','30 exact-sequence enzymes'],loc='upper center',ncol=2,bbox_to_anchor=(.57,.96),frameon=False)
    ax.set_xlim(-.72,.3); clean_axes(ax); savefig(fig,"Supplementary_Figure_S2")

    fig, ax = plt.subplots(figsize=(4.6,3.7),constrained_layout=True)
    inc=enz.included_in_series_analysis.astype(bool)
    grouped = enz.groupby(["n_HCA","n_Coumarins","included_in_series_analysis"], dropna=False).size().reset_index(name="n")
    for _,r in grouped.iterrows():
        adequate = bool(r.included_in_series_analysis)
        ax.scatter(r.n_HCA,r.n_Coumarins,s=18*r.n,color=INK if adequate else "white",edgecolors=INK,label="Included in series contrast" if adequate and r.n==grouped[grouped.included_in_series_analysis].n.max() else ("Insufficient series coverage" if not adequate else None))
        ax.annotate(f"n={int(r.n)}",(r.n_HCA,r.n_Coumarins),xytext=(12,12),textcoords='offset points',fontsize=7)
    ax.set(xlabel="Tested HCA acceptors",ylabel="Tested coumarin acceptors",title="Coverage of the 40 common enzymes")
    ax.set_xlim(-.2,3.55); ax.set_ylim(-.2,3.55); ax.set_xticks([0,1,2,3]); ax.set_yticks([0,1,2,3])
    handles=[Line2D([0],[0],marker="o",color="none",markerfacecolor=INK,markeredgecolor=INK,markersize=5,label="Included in series contrast"), Line2D([0],[0],marker="o",color="none",markerfacecolor="white",markeredgecolor=INK,markersize=5,label="Insufficient series coverage")]
    ax.legend(handles=handles,frameon=False,fontsize=6,loc="upper left"); clean_axes(ax); savefig(fig,"Supplementary_Figure_S3")

if __name__ == "__main__":
    main()
