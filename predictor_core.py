from __future__ import annotations
import math, json
from pathlib import Path
from collections import Counter
import numpy as np
import xgboost as xgb
try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except Exception:
    ProteinAnalysis = None

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"; HYDRO=set("AVILMFWY")
AA_MASS={'A':89.0935,'C':121.159,'D':133.1032,'E':147.1299,'F':165.19,'G':75.0669,'H':155.1552,'I':131.1736,'K':146.1882,'L':131.1736,'M':149.2124,'N':132.1184,'P':115.131,'Q':146.1451,'R':174.2017,'S':105.093,'T':119.1197,'V':117.1469,'W':204.2262,'Y':181.1894}
KD={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3}
CODONS=[a+b+c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT']; DINUCS=[a+b for a in 'ACGT' for b in 'ACGT']
ROOT=Path(__file__).resolve().parent; MODEL_DIR=ROOT/'models'
PROTEIN_FEATURE_NAMES=['length','log_length','mw_approx','gravy','charge_proxy','aromatic_frac','hydrophobic_frac','positive_frac','negative_frac','charged_frac','polar_frac','small_frac','cys_count','cys_frac','odd_cys','gly_frac','pro_frac','glypro_frac','unknown_frac','entropy','max_aa_frac','longest_hydrophobic_run']+[f'frac_{a}' for a in AA_ORDER]
DNA_FEATURE_NAMES=['dna_length','dna_log_length','gc','gc1','gc2','gc3','n_frac','start_atg','terminal_stop','codon_entropy','max_codon_frac','first60_gc','last60_gc','longest_homopolymer']+[f'codon_{c}' for c in CODONS]+[f'dinuc_{d}' for d in DINUCS]

def clean_protein(seq):
    lines=[x.strip() for x in str(seq).splitlines() if x.strip() and not x.strip().startswith('>')]
    return ''.join(lines).replace(' ','').upper().replace('*','')
def clean_dna(seq):
    lines=[x.strip() for x in str(seq).splitlines() if x.strip() and not x.strip().startswith('>')]
    return ''.join(lines).replace(' ','').upper().replace('U','T')
def protein_features(seq):
    s=clean_protein(seq)
    if not s: raise ValueError('Protein sequence is empty.')
    L=len(s); c=Counter(s); valid=''.join(a if a in AA_ORDER else 'G' for a in s)
    if ProteinAnalysis:
        pa=ProteinAnalysis(valid); mw=float(pa.molecular_weight()); gravy=float(pa.gravy())
    else:
        mw=sum(AA_MASS.get(a,75.0669) for a in valid)-18.01528*max(0,L-1); gravy=sum(KD.get(a,0) for a in valid)/L
    cur=mx=0
    for a in s:
        if a in HYDRO: cur+=1; mx=max(mx,cur)
        else: cur=0
    entropy=-sum((n/L)*math.log2(n/L) for n in c.values() if n); unknown=sum(n for a,n in c.items() if a not in AA_ORDER)/L
    x=[L,math.log1p(L),mw,gravy,(c['K']+c['R']+0.1*c['H']-c['D']-c['E'])/L,sum(c[a] for a in 'FWY')/L,sum(c[a] for a in 'AVILMFWY')/L,sum(c[a] for a in 'KRH')/L,sum(c[a] for a in 'DE')/L,sum(c[a] for a in 'KRHDE')/L,sum(c[a] for a in 'STNQC')/L,sum(c[a] for a in 'ACGSTV')/L,float(c['C']),c['C']/L,float(c['C']%2),c['G']/L,c['P']/L,(c['G']+c['P'])/L,unknown,entropy,max(c.values())/L,float(mx)]+[c[a]/L for a in AA_ORDER]
    return np.asarray(x,dtype=np.float32),dict(zip(PROTEIN_FEATURE_NAMES,x))
def _gc(s): return 0.0 if not s else (s.count('G')+s.count('C'))/len(s)
def dna_features(seq):
    s=clean_dna(seq)
    if not s: raise ValueError('DNA sequence is empty.')
    L=len(s); frame=[s[i:i+3] for i in range(0,L-2,3)]; valid=[x for x in frame if len(x)==3 and set(x)<=set('ACGT')]; cc=Counter(valid); nc=max(1,len(valid)); cent=-sum((n/nc)*math.log2(n/nc) for n in cc.values() if n); maxcf=max(cc.values())/nc if cc else 0.0
    gcpos=[]
    for off in range(3): gcpos.append(_gc(''.join(s[i] for i in range(off,L,3) if s[i] in 'ACGT')))
    longest=cur=0; prev=None
    for a in s:
        if a==prev: cur+=1
        else: prev=a; cur=1
        longest=max(longest,cur)
    x=[L,math.log1p(L),_gc(''.join(a for a in s if a in 'ACGT')),gcpos[0],gcpos[1],gcpos[2],sum(1 for a in s if a not in 'ACGT')/L,float(s.startswith('ATG')),float(s[-3:] in {'TAA','TAG','TGA'}),cent,maxcf,_gc(s[:60]),_gc(s[-60:]),float(longest)]
    x += [cc[codon]/nc for codon in CODONS]; dc=Counter(s[i:i+2] for i in range(L-1) if set(s[i:i+2])<=set('ACGT')); nd=max(1,sum(dc.values())); x += [dc[d]/nd for d in DINUCS]
    return np.asarray(x,dtype=np.float32),dict(zip(DNA_FEATURE_NAMES,x))
def _load(name): b=xgb.Booster(); b.load_model(str(MODEL_DIR/name)); return b
def _pred(model,x): return float(model.predict(xgb.DMatrix(np.asarray(x,dtype=np.float32).reshape(1,-1)))[0])
_MODELS=None
def models():
    global _MODELS
    if _MODELS is None:
        _MODELS={'general':[_load(f'general_solubility_fold{i}_xgb_fast.json') for i in range(5)],'cfps':[_load(f'cfps_esol_fold{i}_xgb_fast.json') for i in range(5)],'rp3p':_load('rp3_protein_final_xgb.json'),'rp3dna':_load('rp3_dna_protein_plus_dna_final_xgb.json'),'nesg_exp':_load('nesg_exp_final_xgb.json'),'nesg_sol':_load('nesg_sol_final_xgb.json')}
    return _MODELS

def diagnostics(f,dna_f,p_general,p_cfps,p_prod):
    out=[]
    def add(level,title,detail,consequence,solution,validation): out.append({'level':level,'title':title,'detail':detail,'consequence':consequence,'solution':solution,'validation':validation})
    if p_cfps<0.35: add('high','Low predicted CFPS soluble fraction',f'Model predicts {p_cfps:.1%} soluble fraction.','High probability that a large fraction of synthesized protein will remain insoluble/aggregated.','Test lower-temperature CFPS, reduced template concentration, and chaperone supplementation; consider construct/domain redesign if sequence risks agree.','Measure total versus soluble protein after standardized centrifugation.')
    elif p_cfps<0.60: add('medium','Moderate CFPS solubility',f'Model predicts {p_cfps:.1%} soluble fraction.','Soluble yield may be limited even when total synthesis is adequate.','Screen 20-30 C, chaperone condition, and template concentration before redesigning sequence.','Quantify total and soluble yield across a small condition matrix.')
    if p_prod<0.45: add('high','Low recombinant production probability',f'RP3 production probability is {p_prod:.1%}.','The construct may fail at expression, recovery, or production QC.','Review construct boundaries, N/C-terminal additions, codon/DNA risks, and compare CFPS versus cellular production.','Run small-scale expression plus total/soluble fraction readout.')
    elif p_prod<0.65: add('medium','Intermediate recombinant production probability',f'RP3 production probability is {p_prod:.1%}.','Outcome is uncertain and sensitive to construct/context.','Use a pilot expression screen before committing to scale-up.','Run replicate micro-scale expression under two temperatures.')
    if p_general<0.45: add('high','Low general solubility probability',f'General solubility probability is {p_general:.1%}.','Sequence-intrinsic properties are consistent with poor soluble recovery.','Consider removing nonessential low-complexity/hydrophobic regions or introducing conservative surface-polar substitutions outside functional sites.','Test a shorter/domain construct in parallel.')
    if f['hydrophobic_frac']>0.48 or f['longest_hydrophobic_run']>=8: add('high','Hydrophobic aggregation risk',f"Hydrophobic fraction {f['hydrophobic_frac']:.1%}; longest hydrophobic run {int(f['longest_hydrophobic_run'])} aa.",'Hydrophobic exposure can drive co-translational aggregation or membrane association.','Consider domain trimming or conservative surface-polar substitutions outside functional sites; lower temperature and chaperones may help.','Compare soluble fraction +/- GroEL/ES or DnaK/J/GrpE and at reduced temperature.')
    if f['odd_cys']>0: add('medium','Odd cysteine count',f"{int(f['cys_count'])} cysteines detected.",'Unpaired cysteine can increase disulfide mispairing or heterogeneous products in oxidative conditions.','Confirm intended disulfide pattern; use redox-controlled CFPS when disulfides are required.','Compare reducing versus oxidative/disulfide-supporting CFPS.')
    if f['cys_frac']>0.04: add('medium','High cysteine burden',f"Cysteine fraction {f['cys_frac']:.1%}.",'Disulfide formation may become a folding bottleneck.','Use oxidative/disulfide-supporting CFPS and consider PDI/DsbC-type support where appropriate.','Measure soluble product and non-reducing/reducing gel behavior.')
    if f['length']>800: add('medium','Long construct',f"Length {int(f['length'])} aa.",'Long proteins are less densely represented in the public training set and may impose translational/folding burden.','Consider domain-level constructs or validate with a matched long-protein control.','Compare full-length and domain constructs.')
    if abs(f['charge_proxy'])>0.12: add('medium','Extreme charge proxy',f"Charge proxy {f['charge_proxy']:+.3f}.",'Strong net charge may alter solubility, nonspecific interactions, or buffer sensitivity.','Review highly charged termini/patches and optimize salt/buffer conditions before sequence redesign.','Run a small salt/pH tolerance screen.')
    if f['entropy']<3.5 or f['max_aa_frac']>0.20: add('medium','Low-complexity bias',f"Entropy {f['entropy']:.2f}; max single-residue fraction {f['max_aa_frac']:.1%}.",'Low-complexity segments can contribute to poor folding or proteolysis.','Remove nonessential low-complexity tails or test shorter constructs.','Compare intact and trimmed variants.')
    if dna_f:
        if dna_f['gc']<0.35 or dna_f['gc']>0.70: add('medium','Extreme coding GC',f"Coding GC {dna_f['gc']:.1%}.",'Extremes may affect synthesis, transcription, or translation efficiency depending on system.','Use synonymous redesign toward a moderate host-appropriate GC range while preserving amino-acid sequence.','Compare original and synonymously optimized DNA with the same protein sequence.')
        if dna_f['longest_homopolymer']>=7: add('medium','Long DNA homopolymer',f"Longest homopolymer {int(dna_f['longest_homopolymer'])} nt.",'Long homopolymers can complicate synthesis or sequence stability.','Break the homopolymer with synonymous codon substitutions.','Sequence-verify the synthesized template and compare expression after synonymous redesign.')
        if dna_f['n_frac']>0: add('high','Ambiguous DNA bases',f"Ambiguous-base fraction {dna_f['n_frac']:.2%}.",'Model-based codon analysis is unreliable and synthesis template is not fully specified.','Resolve all ambiguous DNA bases before synthesis.','Re-run prediction using the finalized coding sequence.')
        if not dna_f['start_atg']: add('medium','Coding sequence does not start with ATG','No ATG start codon detected.','The supplied DNA may not represent the intended complete CDS.','Confirm vector architecture and whether the start codon is supplied by the vector.','Sequence-verify the final expression cassette.')
    if not out: add('low','No dominant sequence-level warning','No major first-pass red flag was detected.','This does not guarantee high yield; system-specific effects remain.','Proceed to experimental validation and retain model predictions as a prospective record.','Measure total yield, soluble yield, and replicate variability.')
    return out

def compatibility_release(p_general,p_cfps,p_prod,p_prod_dna,exp_score,dna_present):
    # Transparent compatibility layer using current v3.3 endpoints. This is NOT a retrained v2.5 model.
    protein_score=100*np.clip(0.40*p_cfps+0.30*p_general+0.20*p_prod+0.10*(exp_score/5),0,1)
    dna_score=None if not dna_present else 100*np.clip(p_prod_dna if p_prod_dna is not None else p_prod,0,1)
    combined=None if dna_score is None else 0.70*protein_score+0.30*dna_score
    protein_pass=protein_score>=72
    dna_pass=None if dna_score is None else dna_score>=78
    combined_pass=None if combined is None else combined>=85
    if combined is None: status='PROVISIONAL'
    elif combined_pass and protein_pass and dna_pass: status='GO'
    else: status='HOLD'
    return {'protein_score':float(protein_score),'dna_score':None if dna_score is None else float(dna_score),'combined_score':None if combined is None else float(combined),'protein_pass':bool(protein_pass),'dna_pass':None if dna_pass is None else bool(dna_pass),'combined_pass':None if combined_pass is None else bool(combined_pass),'status':status,'red_warning':bool(combined is not None and combined<85)}

def predict(protein_seq,dna_seq=None):
    x,f=protein_features(protein_seq); m=models(); p_general=float(np.mean([_pred(mm,x) for mm in m['general']])); p_cfps=float(np.clip(np.mean([_pred(mm,x) for mm in m['cfps']]),0,1)); p_prod=_pred(m['rp3p'],x); exp=float(np.clip(_pred(m['nesg_exp'],x),0,5)); sol=float(np.clip(_pred(m['nesg_sol'],x),0,5))
    df=None;p_prod_dna=None; dna_present=bool(dna_seq and clean_dna(dna_seq))
    if dna_present:
        dx,df=dna_features(dna_seq);p_prod_dna=_pred(m['rp3dna'],np.concatenate([x,dx]))
    prod=p_prod_dna if p_prod_dna is not None else p_prod; overall=float(np.clip(0.35*prod+0.35*p_cfps+0.20*p_general+0.10*(exp/5),0,1))
    scaler=np.load(MODEL_DIR/'public_engineered_scaler.npz',allow_pickle=True);z=np.abs((x-scaler['mean'])/scaler['std']);ood=float(np.mean(np.minimum(z,10)));confidence='higher' if ood<1.5 else ('moderate' if ood<2.5 else 'lower / OOD')
    diags=diagnostics(f,df,p_general,p_cfps,prod); validation=[]
    for d in diags:
        if d['validation'] not in validation: validation.append(d['validation'])
    release=compatibility_release(p_general,p_cfps,p_prod,p_prod_dna,exp,dna_present)
    recommendations=[]
    for d in diags:
        if d['solution'] not in recommendations: recommendations.append(d['solution'])
    limitations=['Absolute CFPS yield (mg/mL) is not trained in this release.','CHO secretion/yield is not yet model-backed in this release.','Aggregation is operationally defined as 1 - eSOL/PURE soluble fraction.','The v2.5 compatibility scores are a transparent release layer built from v3.3 endpoints; they are not a separately retrained v2.5 ML model.','The overall score is a prioritization score, not a directly measured endpoint.']
    return {'overall_score':overall,'general_solubility_probability':p_general,'cfps_soluble_fraction':p_cfps,'cfps_aggregation_fraction':1-p_cfps,'rp3_production_probability_protein':p_prod,'rp3_production_probability_with_dna':p_prod_dna,'nesg_expression_score_0_5':exp,'nesg_solubility_score_0_5':sol,'confidence':confidence,'ood_index':ood,'protein_features':f,'dna_features':df,'diagnostics':diags,'risks':[{'level':d['level'],'title':d['title'],'detail':d['detail']} for d in diags],'recommendations':recommendations,'validation_plan':validation,'release':release,'limitations':limitations}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--protein-file');ap.add_argument('--dna-file',default='');a=ap.parse_args();p=Path(a.protein_file).read_text();d=Path(a.dna_file).read_text() if a.dna_file else None;print(json.dumps(predict(p,d),indent=2))
