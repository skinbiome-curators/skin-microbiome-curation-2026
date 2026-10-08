"""Check the review's ICC and cross-validation claims independently."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
ROOT=Path(__file__).resolve().parents[2]
D=ROOT/'data'
m=pd.read_csv(D/'03_sensitive_ITS_sample_info.tsv',sep='\t').sort_values('subject').set_index('run')
x=pd.read_csv(D/'03_sensitive_ITS_genus_pct.tsv',sep='\t',index_col=0).loc[m.index]
y=(m.group=='Sensitive').to_numpy().astype(int)
Xs=[m[['age']].to_numpy(),np.column_stack([m.age,x['g__Malassezia']])]
splits=list(RepeatedStratifiedKFold(n_splits=5,n_repeats=50,random_state=20260918).split(Xs[0],y))
R={'cv':[],'icc':{}}
for X in Xs:
    pooled=[]; foldmean=[]
    for r in range(50):
        p=np.zeros(len(y)); a=[]
        for tr,te in splits[5*r:5*r+5]:
            mod=make_pipeline(StandardScaler(),LogisticRegression(penalty=None,solver='lbfgs',max_iter=10000)).fit(X[tr],y[tr])
            p[te]=mod.predict_proba(X[te])[:,1]; a.append(roc_auc_score(y[te],p[te]))
        pooled.append(roc_auc_score(y,p)); foldmean.append(np.mean(a))
    R['cv'].append({'pooled_mean':float(np.mean(pooled)), 'foldwise_mean':float(np.mean(foldmean))})
for tag,site,gs in [('15_season_16S','Forehead',['Cutibacterium','Staphylococcus','Corynebacterium']),('15_season_ITS','Forehead',['g__Malassezia']),('07_skincare_16S','Cheek',['Cutibacterium','Staphylococcus','Corynebacterium'])]:
    x=pd.read_csv(D/f'{tag}_genus_pct.tsv',sep='\t',index_col=0)
    m=pd.read_csv(D/f'{tag}_sample_info.tsv',sep='\t').set_index('run').loc[x.index]
    m=m[(m.body_site==site)&(m.group!='음성대조')]
    for g in gs:
        a=np.array([x.loc[ix,g].to_numpy() for ix in m.groupby('subject').groups.values()])
        n,k=a.shape
        msb=k*a.mean(axis=1).var(ddof=1)
        msw=((a-a.mean(axis=1)[:,None])**2).sum()/(n*(k-1))
        R['icc'][f'{tag}:{g}']={'n_subjects':n,'n_repeat':k,'icc':float((msb-msw)/(msb+(k-1)*msw))}
(Path(__file__).parent/'보충재계산_결과.json').write_text(json.dumps(R,ensure_ascii=False,indent=2))
print(json.dumps(R,ensure_ascii=False,indent=2))
