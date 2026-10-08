"""Independent calculations from source tables, without importing project code."""
from pathlib import Path
import csv, hashlib, json, math, sys, zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
D = ROOT / 'data'
R = {}
rng = np.random.default_rng(261008)

def get42(assay):
    m = pd.read_csv(D / f'03_sensitive_{assay}_sample_info.tsv', sep='\t').set_index('run')
    x = pd.read_csv(D / f'03_sensitive_{assay}_genus_pct.tsv', sep='\t', index_col=0)
    return m.loc[x.index], x

def auc_summary(x, y):
    x, y = np.asarray(x), np.asarray(y)
    p, n = x[y == 1], x[y == 0]
    pair = (p[:, None] > n[None, :]).astype(float) + .5*(p[:, None] == n[None, :])
    boot = pair[rng.integers(len(p), size=(10000,len(p),1)), rng.integers(len(n), size=(10000,1,len(n)))].mean(axis=(1,2))
    u = stats.mannwhitneyu(p, n, method='asymptotic', use_continuity=False)
    return {'auc': float(roc_auc_score(y,x)), 'ci_10000': np.quantile(boot,[.025,.975]).tolist(), 'median_sensitive': float(np.median(p)), 'median_non_sensitive': float(np.median(n)), 'p_asymptotic_tie_adjusted_no_continuity': float(u.pvalue), 'zero_count': int((x == 0).sum()), 'tie_count': int(len(x)-len(np.unique(x)))}

R['hashes'] = []
for line in (D/'SHA256SUMS').read_text().splitlines():
    expected, filename = line.split(maxsplit=1)
    p=D/filename.lstrip('*')
    actual=hashlib.sha256(p.read_bytes()).hexdigest()
    R['hashes'].append({'file': filename, 'ok': actual == expected})

mi, xi = get42('ITS')
mb, xb = get42('16S')
y = (mi['group'] == 'Sensitive').astype(int).to_numpy()
age = mi['age'].to_numpy(float)
mal = xi['g__Malassezia'].to_numpy()
R['n42'] = {'n': len(y), 'sensitive': int(y.sum()), 'unique_subjects': mi.subject.nunique()}
R['auc'] = {'Malassezia': auc_summary(mal,y)}
for genus in ['Cutibacterium','Staphylococcus','Corynebacterium']:
    R['auc'][genus] = auc_summary(xb[genus], (mb.group == 'Sensitive').astype(int))
R['shapiro'] = {g:{str(k):float(stats.shapiro((xi['g__Malassezia'] if g=='Malassezia' else xb[g]).to_numpy()[(y if g=='Malassezia' else (mb.group=='Sensitive').to_numpy())==k]).pvalue) for k in [0,1]} for g in R['auc']}
fungi = [c for c in xi.columns if c not in ['Unassigned','unidentified'] and (xi[c]>0).any()]
pvals = [stats.mannwhitneyu(xi.loc[y==1,g],xi.loc[y==0,g], method='asymptotic',use_continuity=False).pvalue for g in fungi]
R['all_fungal_bh'] = {'n':len(fungi), 'malassezia_q':float(multipletests(pvals,method='fdr_bh')[1][fungi.index('g__Malassezia')])}
ref = xi.drop(columns=['g__Malassezia','Unassigned','unidentified'],errors='ignore').sum(axis=1).to_numpy()
lr = np.log10((mal+.001)/(ref+.001))
R['auc']['Malassezia_logratio'] = auc_summary(lr,y)
R['auc']['sebum'] = auc_summary(mi.sebumeter_sebum_rigth.to_numpy(),y)
R['corynebacterium_depth_rho'] = float(stats.spearmanr(xb.Corynebacterium,mb.read_depth_final).statistic)
tax = pd.read_csv(D/'03_sensitive_ITS_asv_taxonomy.tsv',sep='\t',index_col=0).fillna('')
counts = pd.read_csv(D/'03_sensitive_ITS_asv_count.tsv',sep='\t',index_col=0).loc[:,xi.index]
mcounts = counts.loc[tax.index[tax.Genus=='g__Malassezia']]
restricta = counts.loc[tax.index[(tax.Genus=='g__Malassezia')&(tax.Species=='s__restricta')]].sum()
species_pct = 100*restricta/counts.sum()
within = restricta/mcounts.sum()
R['restricta'] = {'total_sample':auc_summary(species_pct,y), 'within_malassezia':auc_summary(within,y), 'rho_genus':float(stats.spearmanr(species_pct,mal).statistic), 'share_of_malassezia_reads':float(restricta.sum()/mcounts.to_numpy().sum()), 'species_labeled_malassezia_fraction':float(counts.loc[tax.index[(tax.Genus=='g__Malassezia')&(tax.Species!='')]].to_numpy().sum()/mcounts.to_numpy().sum())}
P=xi.to_numpy()/xi.sum(axis=1).to_numpy()[:,None]
H=stats.entropy(P,axis=1)
rest=xi.drop(columns='g__Malassezia').to_numpy()
R['shannon']={'rho_mal':float(stats.spearmanr(mal,H).statistic),'auc':float(roc_auc_score(y,H)),'rho_without_mal':float(stats.spearmanr(mal,stats.entropy(rest,axis=1)).statistic)}

def fit_predict(X, yy, Xtest=None):
    mod=sm.Logit(yy,sm.add_constant(X,has_constant='add')).fit(disp=0,maxiter=100)
    return mod, mod.predict(sm.add_constant(X if Xtest is None else Xtest,has_constant='add'))

base, bscores = fit_predict(age[:,None],y)
full, fscores = fit_predict(np.column_stack([age,mal]),y)
oof=[np.empty(len(y)),np.empty(len(y))]
for i in range(len(y)):
    tr=np.arange(len(y))!=i
    for j,X in enumerate([age[:,None],np.column_stack([age,mal])]):
        _,pred=fit_predict(X[tr],y[tr],X[~tr]); oof[j][i]=pred[0]
R['logistic']={'auc_age':float(roc_auc_score(y,bscores)),'auc_age_mal':float(roc_auc_score(y,fscores)), 'lr':float(2*(full.llf-base.llf)), 'lr_p':float(stats.chi2.sf(2*(full.llf-base.llf),1)), 'loocv_mle':[float(roc_auc_score(y,p)) for p in oof]}
lpos=[[],[]]
for a in np.where(y==1)[0]:
    for b in np.where(y==0)[0]:
        tr=np.ones(len(y),bool); tr[[a,b]]=False
        for j,X in enumerate([age[:,None],np.column_stack([age,mal])]):
            _,p=fit_predict(X[tr],y[tr],X[[a,b]])
            lpos[j].append(float(p[0]>p[1])+.5*float(p[0]==p[1]))
R['logistic']['leave_pair_out']=[float(np.mean(p)) for p in lpos]
_, lp=fit_predict(np.column_stack([age,lr]),y)
R['logistic']['logratio_mle_delta']=float(roc_auc_score(y,lp)-roc_auc_score(y,bscores))

def hm(A,n1,n2):
    return (A*(1-A)+(n1-1)*(A/(2-A)-A*A)+(n2-1)*(2*A*A/(1+A)-A*A))/(n1*n2)
def npower(A):
    return next(n for n in range(2,5001) if stats.norm.ppf(.975)*math.sqrt(hm(.5,n,n))+stats.norm.ppf(.8)*math.sqrt(hm(A,n,n)) <= A-.5)
R['sample_size']={str(a):{'analyzable_per_group':npower(a),'recruit_if_success_93_120':math.ceil(npower(a)/(93/120))} for a in [.65,.7,.75,.8]}
R['precision']={'current_halfwidth':float(stats.norm.ppf(.975)*math.sqrt(hm(.5,23,19))), '65_halfwidth':float(stats.norm.ppf(.975)*math.sqrt(hm(.5,65,65))), '65_after_dropout_required':math.ceil(65/(93/120)), '57_with_two_independent_successes':math.ceil(57/(93/120)**2)}

meta=pd.read_csv(D/'AAML_merged_data.csv')
selected=meta.loc[meta.Age.between(20,79)]
R['aaml_meta']={'rows':len(meta),'selected_n':len(selected),'sex':selected.Sex.value_counts().to_dict(),'unique_ids':selected['mb.selected_sample-id'].nunique()}
# Separate standard-library XLSX reader (does not use code/xlsx.py).
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def xrows(path):
    with zipfile.ZipFile(path) as z:
        ss=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            ss=[''.join(el.itertext()) for el in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS)]
        with z.open('xl/worksheets/sheet1.xml') as f:
            for _,el in ET.iterparse(f,events=['end']):
                if el.tag.endswith('}row'):
                    row={}
                    for c in el.findall('s:c',NS):
                        col=''.join(a for a in c.attrib['r'] if a.isalpha()); idx=0
                        for a in col: idx=idx*26+ord(a)-64
                        v=c.find('s:v',NS); value=v.text if v is not None else ''
                        if c.attrib.get('t')=='s': value=ss[int(value)]
                        elif c.attrib.get('t')=='inlineStr': value=''.join(c.find('s:is',NS).itertext())
                        row[idx-1]=value
                    yield row
                    el.clear()
it=xrows(D/'AAML_마이크로바이옴 분석 데이터_final.xlsx'); hdr=next(it)
ids=selected['mb.selected_sample-id'].tolist(); wanted={v:k for k,v in hdr.items() if v in set(ids)}
genera={}
for row in it:
    g=row.get(6,'')
    vals=np.array([float(row.get(wanted[s],'') or 0) for s in ids])
    if g not in genera: genera[g]=vals
    else: genera[g]+=vals
stage1=[]; stage2=[]; details={}
for g,v in genera.items():
    if g not in xb.columns: continue
    x=xb[g].to_numpy(); med=float(np.median(v)*100); med42=float(np.median(x))
    if med<=0 or med42<=0: continue
    d={'det595':float((v>0).mean()*100),'det42':float((x>0).mean()*100),'median595_pct':med,'fold':max(med/med42,med42/med)}
    details[g]=d
    if d['det595']>=90 and d['det42']>=90 and med>=.5:
        stage1.append(g)
        if d['fold']<5: stage2.append(g)
R['selection']={'n_genera595':len(genera),'stage1':stage1,'stage2':stage2,'details':{g:details[g] for g in set(stage1)|{'Neisseria','Haemophilus'}}}
R['detection595']={g:int((genera[g]>0).sum()) for g in ['Cutibacterium','Staphylococcus','Corynebacterium']}
R['cutibacterium_age_rho']=float(stats.spearmanr(genera['Cutibacterium'],selected.Age).statistic)
R['raw_table_row_totals']={a:sum(len(pd.read_csv(p,sep='\t')) for p in D.glob(f'*_{a}_genus_pct.tsv'))+(595 if a=='16S' else 0) for a in ['16S','ITS']}
cohorts=[('16S','03_sensitive',None),('16S','02_healthy','Forehead'),('16S','02_healthy','Cheek'),('16S','07_skincare','Cheek'),('16S','15_season','Forehead'),('16S','16_dandruff','Scalp'),('ITS','03_sensitive',None),('ITS','01_aging',None),('ITS','15_season','Forehead'),('ITS','16_dandruff','Scalp')]
R['analysis_cohorts']={}
for assay,tag,site in cohorts:
    x=pd.read_csv(D/f'{tag}_{assay}_genus_pct.tsv',sep='\t',index_col=0)
    m=pd.read_csv(D/f'{tag}_{assay}_sample_info.tsv',sep='\t').set_index('run').loc[x.index]
    if site:
        x=x.loc[(m.body_site==site)&(m.group!='음성대조')]
    R['analysis_cohorts'][f'{assay}:{tag}:{site}']={'n':len(x),'malassezia_detected':int((x['g__Malassezia']>0).sum()) if assay=='ITS' else None}
R['all_cohorts']={a:sum(v['n'] for k,v in R['analysis_cohorts'].items() if k.startswith(a+':'))+(595 if a=='16S' else 0) for a in ['16S','ITS']}
aging=pd.read_csv(D/'01_aging_ITS_genus_pct.tsv',sep='\t',index_col=0)
am=pd.read_csv(D/'01_aging_ITS_sample_info.tsv',sep='\t').set_index('run').loc[aging.index]
R['aging_cheek']={g:{'n':int(((am.body_site=='Cheek')&(am.group==g)).sum()), 'malassezia_median':float(aging.loc[(am.body_site=='Cheek')&(am.group==g),'g__Malassezia'].median())} for g in ['O','Y']}
R['environment']={'python':sys.version, 'seed':261008}
(OUT/'독립재계산_결과.json').write_text(json.dumps(R,ensure_ascii=False,indent=2,default=lambda x:int(x) if isinstance(x,np.integer) else str(x)))
print(json.dumps({k:v for k,v in R.items() if k not in ['hashes','selection']},ensure_ascii=False,indent=2))
