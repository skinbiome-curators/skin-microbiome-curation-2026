# 재검증 계산과 근거 2026년 10월 8일

판단과 발표 수정안은 [클로드 검증 재검증 보고](../클로드_검증_재검증_261008.md)에 있다. 이 폴더는 계산을 다시 확인할 수 있도록 코드와 실행 기록을 보존한다.

| 파일 | 용도 |
|---|---|
| `independent.py` | 팀 코드를 import하지 않고 원자료에서 AUC·부트스트랩·로지스틱·후보 선정·시료 수 재계산 |
| `독립재계산_결과.json` | 독립 계산 값과 원자료 25개 해시 확인 |
| `supplementary.py` | 교차검증 집계 방식과 ICC 별도 계산 |
| `보충재계산_결과.json` | 교차검증·ICC 결과 |
| `rerun.py` | git archive 복사본에서 스크립트 1–23 실행. 16번은 저장된 BLAST 요약만 실행 |
| `재실행_결과.json` | 커밋·환경·실행 시간·종료 코드·결과 파일 SHA-256 |
| `재실행_차이.json` | 원본 결과와 다른 JSON 필드 |
| `실행로그/01.txt`–`23.txt` | 각 스크립트의 표준 출력과 오류 출력 |
| `저자저장소_대조.json` | 저자 GitHub와 로컬 AAML 3파일의 URL·해시 대조 |
| `p값_설정대조.json` | Mann–Whitney 동점·연속성 보정에 따른 수치 차이 |
| `Ahle2022.xml`, `Wang2023.xml`, `AAML2026.xml` | Europe PMC에서 받은 논문 원문 XML. 각각 PMC9452508, PMC10520257, PMC13595949 |

프로젝트 루트에서 다음과 같이 실행한다. 원본 분석 결과는 덮어쓰지 않으며, 이 폴더의 해당 재검증 결과와 로그는 갱신된다.

```bash
conda run -n microbiome python 검증/재검증_261008/independent.py
conda run -n microbiome python 검증/재검증_261008/supplementary.py
conda run -n microbiome python 검증/재검증_261008/rerun.py
```

사용 환경은 Python 3.10.20, numpy 2.2.6, scipy 1.15.2, pandas 2.3.3, statsmodels 0.14.6, scikit-learn 1.7.2다. 독립 부트스트랩 시드는 261008이고, 교차검증 비교는 프로젝트와 같은 분할을 확인하기 위해 20260918을 썼다.
