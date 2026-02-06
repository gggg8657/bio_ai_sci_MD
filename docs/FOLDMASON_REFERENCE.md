# FoldMason 완전 가이드

> **Multiple Protein Structure Alignment at Scale**  
> 저장소: https://github.com/steineggerlab/foldmason  
> 논문: Gilchrist CLM, Mirdita M, Steinegger M. *Science* (2026) doi:10.1126/science.ads6733

---

## 목차

1. [개요](#개요)
2. [설치 방법](#설치-방법)
3. [핵심 모듈](#핵심-모듈)
4. [easy-msa 워크플로우](#easy-msa-워크플로우)
5. [출력 형식](#출력-형식)
6. [주요 파라미터](#주요-파라미터)
7. [고급 기능](#고급-기능)
8. [사용 예제](#사용-예제)
9. [웹서버](#웹서버)

---

## 개요

FoldMason은 대규모 단백질 구조 집합으로부터 정확한 **다중 구조 정렬(MSA)**을 생성하는 도구입니다.

### 핵심 특징

| 특징 | 설명 |
|------|------|
| **구조 기반 정렬** | 서열이 아닌 3D 구조를 기반으로 정렬 |
| **3Di 알파벳** | Foldseek의 구조 알파벳 활용 |
| **대규모 처리** | 수천 개 구조를 효율적으로 정렬 |
| **LDDT 스코어링** | 정렬 품질을 LDDT로 평가 |
| **반복 정제** | 정렬 품질 향상을 위한 반복 최적화 |
| **대화형 시각화** | HTML 보고서로 정렬 결과 탐색 |

---

## 설치 방법

### Conda (권장)

```bash
conda install -c conda-forge -c bioconda foldmason
```

### 사전 컴파일 바이너리

```bash
# Linux AVX2
wget https://mmseqs.com/foldmason/foldmason-linux-avx2.tar.gz
tar xvzf foldmason-linux-avx2.tar.gz
export PATH=$(pwd)/foldmason/bin/:$PATH

# Linux SSE2
wget https://mmseqs.com/foldmason/foldmason-linux-sse2.tar.gz
tar xvzf foldmason-linux-sse2.tar.gz

# Linux ARM64
wget https://mmseqs.com/foldmason/foldmason-linux-arm64.tar.gz

# macOS (Universal)
wget https://mmseqs.com/foldmason/foldmason-osx-universal.tar.gz
```

### Docker

```bash
# 저장소의 Dockerfile로 빌드
git clone https://github.com/steineggerlab/foldmason.git
cd foldmason
docker build -t foldmason .
```

---

## 핵심 모듈

| 모듈 | 설명 | 용도 |
|------|------|------|
| **easy-msa** | 올인원 MSA 워크플로우 | 구조 파일 → 정렬 + 보고서 |
| **structuremsa** | 구조 DB 기반 MSA | 사전 생성된 DB에서 정렬 |
| **structuremsacluster** | 클러스터링 + MSA | 대규모 데이터셋 처리 |
| **msa2lddt** | LDDT 점수 계산 | 정렬 품질 평가 |
| **msa2lddtreport** | LDDT + HTML 보고서 | 시각화 보고서 생성 |
| **msa2lddtjson** | LDDT + JSON | 웹서버 연동용 |
| **refinemsa** | MSA 정제 | 반복적 품질 향상 |
| **createdb** | 구조 DB 생성 | 입력 전처리 |
| **convertalis** | 형식 변환 | BLAST-tab, SAM 등 출력 |

---

## easy-msa 워크플로우

가장 간단한 사용법:

```bash
foldmason easy-msa <PDB/mmCIF 파일들> <출력파일> <임시폴더> [옵션]
```

### 기본 예제

```bash
# 여러 PDB 파일 정렬
foldmason easy-msa *.pdb result.fasta tmpFolder

# HTML 보고서 포함
foldmason easy-msa *.pdb result tmpFolder --report-mode 1

# gzip 압축 파일 지원
foldmason easy-msa *.pdb.gz result.fasta tmpFolder
```

### 내부 동작

`easy-msa`는 다음 명령을 순차 실행:

```bash
# 1. 구조 DB 생성
foldmason createdb <입력파일들> myDb

# 2. 구조 MSA 수행
foldmason structuremsa myDb result

# 3. (--report-mode 1) LDDT 계산 + HTML
foldmason msa2lddtreport myDb result_aa.fa result.html --guide-tree result.nw
```

---

## 출력 형식

### 기본 출력 (--report-mode 0)

| 파일 | 설명 |
|------|------|
| `result_aa.fa` | 아미노산 서열 FASTA 정렬 |
| `result_3di.fa` | 3Di 구조 알파벳 FASTA 정렬 |
| `result.nw` | Newick 형식 가이드 트리 |

### 보고서 출력 (--report-mode 1)

| 파일 | 설명 |
|------|------|
| `result.html` | 대화형 HTML 시각화 |
| + 기본 출력 파일들 | |

### JSON 출력 (--report-mode 2)

| 파일 | 설명 |
|------|------|
| `result.json` | 웹서버 로드용 JSON |
| + 기본 출력 파일들 | |

---

## 주요 파라미터

### 정렬 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `--gap-open` | 25 | 갭 오픈 페널티 |
| `--gap-extend` | 2 | 갭 확장 페널티 |
| `--comp-bias-corr` | 1 | 아미노산 조성 편향 보정 (0-1) |

### 프로파일 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `--wg` | 1 | 전역 서열 가중치 사용 |
| `--match-ratio` | 0.9 | 잔기 비율 임계값 |
| `--filter-msa` | 1 | MSA 필터링 (0: 끔, 1: 켬) |
| `--diff` | 5 | 다양성 필터링 최소 서열 수 |
| `--qsc` | -20.0 | 쿼리 스코어 임계값 |
| `--mask-profile` | 1 | tantan 마스킹 |

### 스코어링 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `--pair-threshold` | 0.0 | LDDT 계산용 갭 비율 임계값 |
| `--bitfactor-aa` | 1.1 | AA 행렬 비트 팩터 |
| `--bitfactor-3di` | 2.1 | 3Di 행렬 비트 팩터 |

### 정제 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `--refine-iters` | 0 | 정제 반복 횟수 |
| `--refine-seed` | -1 | 랜덤 시드 (-1: 자동) |

### 이웃 스코어링 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `--nb-sigma` | 3.841 | 이웃 점수 감쇠 상수 |
| `--nb-multiplier` | 13.0 | 이웃 점수 승수 |
| `--nb-ang-cut` | 45.0 | 최대 거리 컷오프 (Å) |
| `--nb-low-cut` | 0.02 | 최소 이웃 점수 임계값 |
| `--fast` | 0 | 빠른 모드 (이웃 스코어링 비활성화) |

### 입력/출력 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `--input-format` | 0 | 입력 형식 (0: 자동, 1: PDB, 2: mmCIF, 3: mmJSON, 4: ChemComp, 5: Foldcomp) |
| `--file-include` | `.*` | 포함할 파일 정규식 |
| `--file-exclude` | `^$` | 제외할 파일 정규식 |
| `--report-mode` | 0 | 보고서 모드 |

### 성능 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `--threads` | 전체 | 사용 CPU 코어 수 |
| `--gpu` | 0 | GPU 사용 (CUDA) |
| `--precluster` | 0 | 사전 클러스터링 |
| `-v` | 3 | 상세 출력 레벨 (0-3) |

---

## 고급 기능

### 1. 대규모 데이터셋 처리

수천 개 구조를 정렬할 때 `--precluster` 옵션 사용:

```bash
foldmason easy-msa large_dataset/*.pdb result tmpFolder --precluster
```

### 2. MSA 정제 (Refinement)

초기 정렬 후 반복적으로 품질 향상:

```bash
# easy-msa에서 직접
foldmason easy-msa *.pdb result tmpFolder --refine-iters 100

# 별도 refinemsa 모듈
foldmason refinemsa myDb initial.fasta refined.fasta --refine-iters 1000
```

### 3. 외부 MSA의 LDDT 계산

다른 도구로 생성한 MSA의 품질 평가:

```bash
# DB 생성 (구조 파일 필요)
foldmason createdb structures/ myDb

# LDDT 계산
foldmason msa2lddt myDb external_msa.fa

# HTML 보고서 생성
foldmason msa2lddtreport myDb external_msa.fa report.html
```

### 4. 사용자 정의 가이드 트리

```bash
foldmason easy-msa *.pdb result tmpFolder --guide-tree custom_tree.nw
```

### 5. 구조 DB 사전 생성

동일 입력을 여러 번 정렬할 때:

```bash
# DB 한 번 생성
foldmason createdb structures/ structureDB

# 여러 번 사용
foldmason structuremsa structureDB result1
foldmason structuremsa structureDB result2 --gap-open 20
```

---

## 사용 예제

### 예제 1: 기본 MSA

```bash
# AlphaFold 모델들 정렬
foldmason easy-msa fold_test1_model_*.pdb alignment.fasta tmp
```

### 예제 2: 전체 보고서 생성

```bash
foldmason easy-msa \
  fold_test1_model_0.pdb \
  fold_test1_model_1.pdb \
  fold_test1_model_2.pdb \
  result tmp --report-mode 1

# 출력: result_aa.fa, result_3di.fa, result.nw, result.html
```

### 예제 3: 고품질 정렬 (정제 포함)

```bash
foldmason easy-msa *.pdb refined tmp \
  --refine-iters 500 \
  --gap-open 15 \
  --gap-extend 1 \
  --report-mode 1
```

### 예제 4: 대규모 구조 집합

```bash
foldmason easy-msa large_family/*.pdb family_msa tmp \
  --precluster \
  --threads 32 \
  --fast
```

### 예제 5: 형식 변환

```bash
# BLAST-tab 형식으로 변환
foldmason convertalis myDb result alnResult.m8 --format-mode 0
```

---

## 웹서버

온라인 웹서버: https://search.foldseek.com/foldmason

### 특징
- 구조 파일 업로드만으로 빠른 정렬
- 대화형 HTML 시각화
- 다운로드 가능한 결과

---

## LDDT (Local Distance Difference Test)

FoldMason은 정렬 품질을 **LDDT** 점수로 평가합니다.

### LDDT란?
- 구조 유사성을 측정하는 점수 (0-1)
- 거리 차이 기반 로컬 평가
- 전역 중첩 불필요

### 계산 방식
1. MSA의 각 컬럼에서 쌍별(pairwise) LDDT 계산
2. 전체 MSA 길이로 평균

### 사용

```bash
# LDDT만 출력
foldmason msa2lddt myDb result.fa

# LDDT + 시각화
foldmason msa2lddtreport myDb result.fa report.html
```

---

## 3Di 알파벳

FoldMason은 Foldseek의 **3Di (3D-interaction) 알파벳**을 사용합니다.

- 20자 구조 알파벳
- 잔기의 3D 이웃 관계 인코딩
- 서열 정렬 알고리즘으로 구조 비교 가능

### 출력 예

```
>structure1
VVDDLVVPLDDD...  # 3Di 서열
>structure2  
VADDLVVPLDDD...
```

---

## 참고 자료

- **GitHub**: https://github.com/steineggerlab/foldmason
- **논문**: Gilchrist et al. Science (2026)
- **웹서버**: https://search.foldseek.com/foldmason
- **Foldseek**: https://github.com/steineggerlab/foldseek
