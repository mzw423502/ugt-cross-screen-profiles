"""Assemble current readable tables and full-precision data for XLSX export."""
from pathlib import Path
import csv, json, re

ROOT=Path(__file__).resolve().parent
def typed(x):
    if x is None or not isinstance(x,str):return x
    if re.fullmatch(r'-?\d+(?:\.\d+)?%',x):return float(x[:-1])/100
    if re.fullmatch(r'-?\d+',x):return int(x)
    if re.fullmatch(r'-?\d+\.\d+(?:e[+-]?\d+)?',x,re.I):return float(x)
    if x in ('True','False'):return x=='True'
    return x

def main():
    labels=['Enzyme_matches','Acceptor_matches','Reaction_states','Acceptor_landscape','Chemical_series','Enzyme_contrasts','Compound_leaveout','Exact_sequence','Cluster_bootstrap','Threshold_counts','Product_aggregation','Sources','Focal_structures','Biochemical_context','Assay_conditions','Threshold_effects']
    tables=json.loads((ROOT/'generated/supplement_readable_tables.json').read_text(encoding='utf-8'))
    sheets=[]
    for i,(label,t) in enumerate(zip(labels,tables),1):
        sheets.append({'name':f'S{i}_{label}','rows':[t['headers']]+[[typed(v) for v in row] for row in t['rows']]})
    files={
      'Reaction_units':'input/S3_ReactionUnits600.csv',
      'Sequence_records':'input/S1_EnzymeIdentity.csv',
      'Acceptor_records':'input/S2_AcceptorIdentity.csv',
      'Enzyme_full_precision':'generated/per_enzyme_series.csv',
      'Bootstrap_draws':'generated/bootstrap_full.csv',
      'Bootstrap_full_precision':'generated/bootstrap_series_final.csv',
      'Source_manifest':'input/source_manifest.csv',
    }
    for name,file in files.items():
        with (ROOT/file).open(encoding='utf-8-sig',newline='') as f:
            rows=list(csv.reader(f))
        sheets.append({'name':name,'rows':[rows[0]]+[[typed(v) if v!='' else None for v in r] for r in rows[1:]]})
    notes=[['Item','Definition'],
      ['Article','Cross-screen substrate profiles of Arabidopsis UDP-glycosyltransferases diverge by chemical series'],
      ['Summary tables','S1–S16 correspond exactly to the numbered readable tables in Supplementary Information. Summary decimal values are rounded for display.'],
      ['Experimental unit','One enzyme construct and one connectivity-matched acceptor; multiple products do not create independent reactions.'],
      ['Coverage','600 units; 578 tested pairs; 22 Y18 NOT_TESTED units. NOT_TESTED is missing, not zero.'],
      ['Status semantics','TESTED_ACTIVE = qualifying source call; TESTED_INACTIVE = tested without a qualifying call; NOT_TESTED and AMBIGUOUS are excluded from binary denominators.'],
      ['Source labels','y = Yang et al. 2018; m = Sirirungruang et al. 2025. Main M25 threshold is 0.85.'],
      ['Bootstrap_draws','500 draws per set; primary 39 enzymes/17 clusters and exact 30 enzymes/13 clusters; both contrasts use the same resample. Seed 20260919.'],
      ['Full precision','The three full-precision/draw sheets preserve numerical outputs. The executable calculations and CSVs are in 05_REPRODUCIBILITY.'],
      ['Blank / em dash','Missing or undefined, never silently coded as zero.'],
      ['Y18 source','https://doi.org/10.1038/s41589-018-0154-9'],
      ['M25 source','https://doi.org/10.1038/s41467-025-61530-6'],
      ['Rebuild','python prepare_workbook_data.py, then node build_workbook.mjs (requires @oai/artifact-tool).'],
    ]
    sheets.append({'name':'README','rows':notes})
    (ROOT/'generated/workbook_content.json').write_text(json.dumps(sheets,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    print(f'{len(sheets)} synchronized workbook sheets assembled')
if __name__=='__main__':main()
