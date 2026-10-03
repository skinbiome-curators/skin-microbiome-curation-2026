import csv, json, math, os, random, statistics as st
random.seed(20260918)

# ── OS 이식성 (2026-09-18 추가) ──────────────────────────────────────────────
# 맥에서 만든 결과를 Windows 에서 그대로 재현하기 위한 두 가지.
#   1) stdout 을 UTF-8 로 고정 — 한글 Windows 기본 cp949 는 em-dash(—) 를 못 쓴다
#   2) save_json — encoding·newline 을 못 박아 OS 와 무관하게 같은 바이트를 낸다
import sys as _sys, json as _json
try:
    _sys.stdout.reconfigure(encoding='utf-8')
    _sys.stderr.reconfigure(encoding='utf-8')
except Exception:
    pass

def save_json(obj, path):
    with open(path, 'w', encoding='utf-8', newline=chr(10)) as _fh:
        _json.dump(obj, _fh, ensure_ascii=False, indent=2)
# ────────────────────────────────────────────────────────────────────────────

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(BASE, '분석표')
RAW  = os.path.join(BASE, '원자료')
CAND=['Cutibacterium','Staphylococcus','Corynebacterium']

def f(x):
    try: return float(x)
    except: return None

def rank(v):
    idx=sorted(range(len(v)),key=lambda i:v[i]); r=[0.0]*len(v); i=0
    while i<len(idx):
        j=i
        while j+1<len(idx) and v[idx[j+1]]==v[idx[i]]: j+=1
        rr=(i+j)/2+1
        for k in range(i,j+1): r[idx[k]]=rr
        i=j+1
    return r

def spearman(x,y):
    rx,ry=rank(x),rank(y); n=len(x)
    mx,my=st.mean(rx),st.mean(ry)
    num=sum((a-mx)*(b-my) for a,b in zip(rx,ry))
    den=math.sqrt(sum((a-mx)**2 for a in rx)*sum((b-my)**2 for b in ry))
    rho=num/den if den else 0.0
    t=rho*math.sqrt((n-2)/max(1e-12,1-rho**2))
    p=2*(1-0.5*(1+math.erf(abs(t)/math.sqrt(2))))
    return rho,p

def auc_p(a,b):            # a=양성군, b=음성군
    allv=sorted([(v,1) for v in a]+[(v,0) for v in b]); ranked=[]; i=0
    while i<len(allv):
        j=i
        while j+1<len(allv) and allv[j+1][0]==allv[i][0]: j+=1
        rr=(i+j)/2+1
        for k in range(i,j+1): ranked.append((rr,allv[k][1]))
        i=j+1
    n1,n2=len(a),len(b); R=sum(r for r,g in ranked if g==1); U=R-n1*(n1+1)/2
    mu=n1*n2/2; sd=math.sqrt(n1*n2*(n1+n2+1)/12); z=(U-mu)/sd
    return U/(n1*n2), 2*(1-0.5*(1+math.erf(abs(z)/math.sqrt(2))))

def auc_ci(a,b,B=4000):
    s=[]
    for _ in range(B):
        aa=[random.choice(a) for _ in a]; bb=[random.choice(b) for _ in b]
        s.append(auc_p(aa,bb)[0])
    s.sort(); return s[int(.025*B)], s[int(.975*B)]

# ---------- 595 ----------
gr=list(csv.DictReader(open(f'{OUT}/595_genus_relabund.tsv'),delimiter='\t'))
md={r['metadata_id']:r for r in csv.DictReader(open(f'{OUT}/595_metadata.tsv'),delimiter='\t')}
qc={r['metadata_id']:r for r in csv.DictReader(open(f'{OUT}/595_qc_flags.tsv'),delimiter='\t')}
mis=list(csv.DictReader(open(f'{OUT}/595_taxonomy_mismatch.tsv'),delimiter='\t'))
ids=[r['metadata_id'] for r in gr]
ab={g:[f(r[g]) or 0.0 for r in gr] for g in CAND}
age=[f(md[i]['Age']) for i in ids]; oil=[f(md[i]['Oil_avg']) for i in ids]
rds=[f(qc[i]['non_chimeric_reads']) for i in ids]
lowd=[qc[i]['flag_low_depth_p10']=='Y' for i in ids]
lowrich=[int(qc[i]['n_genus_detected'])<10 for i in ids]

res={'n':len(ids),'n_low_depth':sum(lowd),'n_low_richness':sum(lowrich),'genera':{}}
for g in CAND:
    v=ab[g]
    det=lambda sel:100*sum(1 for x,s in zip(v,sel) if s and x>0)/max(1,sum(sel))
    r_age,p_age=spearman(age,v); r_oil,p_oil=spearman(oil,v); r_dep,p_dep=spearman(rds,v)
    bands={}
    for lo,hi in [(20,39),(40,54),(55,79)]:
        s=[lo<=a<=hi for a in age]
        bands[f'{lo}-{hi}']=round(100*st.median([x for x,k in zip(v,s) if k]),3)
    keep=[not k for k in lowrich]
    mm=[m for m in mis if m['Genus']==g]
    mm_share=sum(float(m['sum_relabund_595']) for m in mm)/max(1e-12,sum(v))
    res['genera'][g]={
      'detect_pct':round(100*sum(1 for x in v if x>0)/len(v),2),
      'median_pct':round(100*st.median(v),3),
      'detect_pct_lowdepth':round(det(lowd),1),
      'detect_pct_rest':round(det([not k for k in lowd]),1),
      'rho_depth':round(r_dep,3),'p_depth':p_dep,
      'rho_age':round(r_age,3),'p_age':p_age,
      'rho_oil':round(r_oil,3),'p_oil':p_oil,
      'median_pct_by_ageband':bands,
      'median_pct_excl_lowrich':round(100*st.median([x for x,k in zip(v,keep) if k]),3),
      'mismatch_rows':len(mm),'mismatch_share_of_genus':round(mm_share,5),
      'mismatch_examples':[m['Species'] for m in mm[:3]],
    }

# 03의 나이차(민감군이 11세 젊음)로 기대되는 변화량: 595에서 30세 vs 41세 부근 중앙값 비교
for g in CAND:
    v=ab[g]
    a=[x for x,ag in zip(v,age) if 25<=ag<=34]; b=[x for x,ag in zip(v,age) if 36<=ag<=45]
    res['genera'][g]['expected_age_shift_pct']={'25-34세':round(100*st.median(a),3),'36-45세':round(100*st.median(b),3),
                                                'n':[len(a),len(b)]}

# ---------- 03 (④ 판별력) ----------
si={r['run']:r for r in csv.DictReader(open(f'{RAW}/03_sensitive_16S_sample_info.tsv'),delimiter='\t')}
rd=csv.reader(open(f'{RAW}/03_sensitive_16S_genus_pct.tsv'),delimiter='\t')
hdr=next(rd); col={g:hdr.index(g) for g in CAND if g in hdr}
rows03=[r for r in rd]
grp=[si[r[0]]['group'] for r in rows03]; ag03=[f(si[r[0]]['age']) for r in rows03]
res03={'n':len(rows03),'n_sensitive':grp.count('Sensitive'),'n_nonsensitive':grp.count('Non-sensitive'),
       'age_median':{'Sensitive':st.median([a for a,g_ in zip(ag03,grp) if g_=='Sensitive']),
                     'Non-sensitive':st.median([a for a,g_ in zip(ag03,grp) if g_=='Non-sensitive'])},'genera':{}}
aS=[a for a,g_ in zip(ag03,grp) if g_=='Sensitive']; aN=[a for a,g_ in zip(ag03,grp) if g_=='Non-sensitive']
au,pv=auc_p(aS,aN); lo,hi=auc_ci(aS,aN)
res03['auc_age']={'auc':round(au,3),'ci':[round(lo,3),round(hi,3)],'p':round(pv,4)}
for g in CAND:
    j=col[g]; v=[f(r[j]) or 0.0 for r in rows03]
    S=[x for x,g_ in zip(v,grp) if g_=='Sensitive']; N=[x for x,g_ in zip(v,grp) if g_=='Non-sensitive']
    au,pv=auc_p(S,N); lo,hi=auc_ci(S,N); r_a,p_a=spearman(ag03,v)
    res03['genera'][g]={'median_S':round(st.median(S),3),'median_N':round(st.median(N),3),
                        'auc':round(au,3),'ci':[round(lo,3),round(hi,3)],'p':round(pv,4),
                        'rho_age_in03':round(r_a,3),'detect_pct':round(100*sum(1 for x in v if x>0)/len(v),1)}
out={'gate123_595':res,'gate4_03':res03}
save_json(out, f'{OUT}/gate_results.json')
print(json.dumps(out,ensure_ascii=False,indent=2))
