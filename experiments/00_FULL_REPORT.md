# SSTR2 Drug Discovery Pipeline -- 전체 실험 보고서

> **프로젝트**: PRST_N_FM (PyRosetta + FoldMason + BioNeMo)
> **목표**: SSTR2 (Somatostatin Receptor Type 2)에 대해 기존 리간드(Somatostatin)보다 더 강하거나 오래 결합하는 분자 탐색
> **환경**: Ubuntu 22.04 (WSL2), Conda `bio-tools`, NVIDIA NIM API (GPU 불필요)
> **날짜**: 2026-02-09 ~ 2026-02-10

---

## 1. 프로젝트 개요

SSTR2는 신경내분비종양(NET) 치료의 핵심 표적 GPCR이다. 현재 승인된 약물(Octreotide, Lanreotide)은 모두 천연 리간드 Somatostatin-14의 펩타이드 유사체이다. 본 프로젝트는 3가지 경로를 병렬로 탐색하여 더 나은 바인더 후보를 발굴한다:

- **Arm 1**: 소분자 (Small Molecule) -- MolMIM 생성 + DiffDock 도킹
- **Arm 2**: 펩타이드 변이체 -- Somatostatin 돌연변이 분석
- **Arm 3**: De Novo 펩타이드 -- RFdiffusion + ProteinMPNN + ESMFold

```
                         AlphaFold3 복합체
                        (SSTR2 + Somatostatin)
                               |
                     바인딩 포켓 분석 (35잔기)
                        /      |       \
                Arm 1          Arm 2         Arm 3
              소분자         펩타이드       De Novo
             MolMIM         변이체        RFdiffusion
                |              |              |
            DiffDock       FlexPepDock    ProteinMPNN
              도킹            도킹            |
                |              |          ESMFold
                 \             |          /
                   통합 랭킹 & 비교
```

---

## 2. 환경 구축 (실험 01~04)

### 2.1 Conda 환경 (`bio-tools`)

| 도구 | 버전 | 용도 |
|------|------|------|
| PyRosetta | 2026.06 | 분자 모델링 (scoring, relax, docking) |
| FoldMason | 4.dd3c235 | 구조 기반 다중 정렬 + lDDT |
| PyMOL | 3.1.0 (OSS) | 분자 시각화 |
| Biopython | 1.86 | PDB/CIF 파싱 |
| RDKit | 2025.03.6 | 화학 정보학, SMILES-SDF 변환 |
| Meeko | 0.7.1 | AutoDock-GPU 전처리 |

```bash
conda env create -f environment-bio-tools.yml
conda activate bio-tools
python scripts/verify_bio_tools_env.py
```

### 2.2 NVIDIA NIM API

GPU 없이 NVIDIA 호스팅 API로 5개 AI 모델 사용:

| 모델 | 기능 | API 엔드포인트 |
|------|------|---------------|
| MolMIM | 소분자 생성/최적화 | `health.api.nvidia.com/.../nvidia/molmim` |
| DiffDock | 분자 도킹 (blind) | `health.api.nvidia.com/.../mit/diffdock` |
| RFdiffusion | 단백질 바인더 백본 설계 | `health.api.nvidia.com/.../ipd/rfdiffusion` |
| ProteinMPNN | 역접힘 (백본→서열) | `health.api.nvidia.com/.../ipd/proteinmpnn` |
| ESMFold | 서열→구조 예측 | `health.api.nvidia.com/.../nvidia/esmfold` |

인증: `nvapi-` 로 시작하는 API 키 (`molmim.key` 파일), 모든 모델 공통.

---

## 3. 데이터 준비 (실험 01~03)

### 3.1 AlphaFold3 구조 예측

SSTR2 + Somatostatin-14 복합체를 AlphaFold3 Server에서 예측:

| 항목 | 값 |
|------|-----|
| Chain A | Somatostatin-14: `AGCKNFFWKTFTSC` (14잔기 고리 펩타이드) |
| Chain B | SSTR2: 369잔기 GPCR |
| 모델 수 | 5개 (model_0 ~ model_4) |
| 최고 모델 | model_0 (ranking_score=0.83, ipTM=0.71, pTM=0.74) |
| 파일 위치 | `data/fold_test1/` |

### 3.2 CIF-PDB 변환

AlphaFold3 출력 mmCIF 13개를 PDB로 변환:

```bash
python scripts/cif_to_pdb.py data/fold_test1
# 결과: 13/13 변환 완료
```

### 3.3 FoldMason 구조 정렬

3개 모델(0,1,2)의 구조적 일관성 평가:

```bash
foldmason easy-msa model_0.pdb model_1.pdb model_2.pdb result --report-mode 1
```

- **Average MSA lDDT**: 0.664 (중간 수준 -- 일부 영역 유동적)
- HTML 리포트, AA/3Di MSA, Newick 트리 생성

### 3.4 PyMOL 시각화

PDB 로드, cartoon/surface 렌더링, B-factor 색상, 모델 중첩(align) 확인.

---

## 4. SSTR2 바인딩 포켓 분석 (Step 0)

**스크립트**: `bionemo/04_sstr2_pocket_analysis.py`

Somatostatin(Chain A) 기준 **5A** 이내 SSTR2(Chain B) 잔기를 Biopython NeighborSearch로 추출.

### 결과

- **바인딩 포켓 잔기**: 35개
- **핵심 접촉 잔기**: B122(ASP), B127(PHE), B184(ARG), B197(TRP), B205(TYR), B272(PHE), B294(PHE)

### Somatostatin 잔기별 접촉 요약

| 잔기 | 접촉 SSTR2 잔기 수 | 최소 거리 | 비고 |
|------|-------------------|----------|------|
| A1 Ala | 1 | 3.29A | 표면 접촉 |
| A6 Phe | 4 | 3.27A | 방향족 상호작용 |
| A7 Phe | 8 | 3.11A | 깊은 포켓 삽입 |
| **A8 Trp** | **13** | **2.91A** | **가장 깊은 삽입, 핵심 약효단** |
| **A9 Lys** | **9** | **2.58A** | **가장 가까운 접촉, 염기성 상호작용** |
| A10 Thr | 5 | 3.42A | 수소결합 |
| A12 Thr | 3 | 2.87A | 수소결합 |

> **핵심 발견**: Trp8과 Lys9이 SSTR2 결합의 핵심 잔기. 이 두 잔기의 상호작용을 보존하거나 강화하는 것이 새 바인더 설계의 핵심.

---

## 5. Arm 1: 소분자 스크리닝

**스크립트**: `bionemo/05_sstr2_smallmol_screen.py`

### 5.1 시드 분자

문헌 기반 SSTR2 관련 소분자 scaffold 4종:

| 시드 | SMILES | 근거 |
|------|--------|------|
| Paltusotine_core | `CC1=CC(=CC(=C1)OCC2=CC=NC=C2)C(=O)NC3CCCCC3` | SSTR2 소분자 작용제 |
| SSTR2_agonist_1 | `c1ccc2c(c1)c(=O)n(c(=O)[nH]2)CC(=O)O` | Quinazolinedione scaffold |
| Indole_scaffold | `c1ccc2c(c1)[nH]cc2CC(=O)NC` | 인돌 기반 scaffold |
| Benzimidazole_hit | `c1ccc2[nH]c(nc2c1)CCNC(=O)c1ccccc1` | 벤즈이미다졸 scaffold |

### 5.2 MolMIM 후보 생성

각 시드 → CMA-ES QED 최적화 → 시드당 10개 = **총 40개** 후보

최고 QED 후보:

| 순위 | SMILES | QED | 시드 |
|------|--------|-----|------|
| 1 | `Cn1c(=O)c(C(=O)O)c(CC(=O)C(C)(C)C)c2ccccc21` | 0.944 | SSTR2_agonist_1 |
| 2 | `O=C(COc1cc2ccccc2cc(=O)n1)NC1CCCCC1` | 0.941 | Paltusotine |
| 3 | `CNC(=O)Cc1c(OC(F)(F)F)ccc2ccccc12` | 0.940 | Indole |

### 5.3 DiffDock 도킹

QED 상위 15개를 SSTR2에 blind docking:

- **성공률**: 15/15 (100%)
- **DiffDock confidence**: -3.0 ~ -5.5 (higher = better)
- 각 분자당 5개 포즈 생성

### 5.4 해석

DiffDock은 blind docking으로 바인딩 포켓을 자동 탐색하며, confidence score가 높을수록 결합 가능성이 높다. 생성된 40개 후보 중 drug-likeness(QED > 0.9)를 가진 후보가 다수 확인되어, 후속 실험가치가 있다.

---

## 6. Arm 2: 펩타이드 변이체 분석

**스크립트**: `bionemo/06_sstr2_flexpep_dock.py`

### 6.1 변이체 설계

Somatostatin-14 (`AGCKNFFWKTFTSC`) 기반 13개 변이체:

| 종류 | 변이체 | 서열 | 목적 |
|------|--------|------|------|
| 야생형 | wildtype | `AGCKNFFWKTFTSC` | 기준 |
| Ala scan | F6A | `AGCKNAFWKTFTSC` | Phe6 기여도 평가 |
| Ala scan | F7A | `AGCKNFAWKTFTSC` | Phe7 기여도 평가 |
| Ala scan | W8A | `AGCKNFFAKTFTSC` | **Trp8 핵심 잔기** 확인 |
| Ala scan | K9A | `AGCKNFFWATFTSC` | **Lys9 핵심 잔기** 확인 |
| Ala scan | T10A | `AGCKNFFWKAFTSC` | Thr10 기여도 |
| Ala scan | F11A | `AGCKNFFWKTATSC` | Phe11 기여도 |
| 유사체 | octreotide | `FCFWKTCT` | 기존 승인 약물 코어 |
| 강화 | enhanced_2 | `AGCRNFFWKTFTSC` | K4R (양전하 강화) |
| 강화 | enhanced_3 | `AGCKNYFWKTFTSC` | F6Y (수소결합 추가) |
| 강화 | enhanced_4 | `AGCKNFFWRTFTSC` | K9R (양전하 강화) |
| 강화 | enhanced_5 | `AGCKNFFWKTYTSC` | F11Y (수소결합 추가) |

### 6.2 현재 상태

- 변이체 서열 분석 완료 (13개)
- PyRosetta FlexPepDock: Rosetta DB 초기화 이슈로 미수행
- **대안**: AlphaFold3 Server에서 각 변이체 + SSTR2 복합체 예측 → ipTM 비교

---

## 7. Arm 3: De Novo 펩타이드 바인더 설계

**스크립트**: `bionemo/07_sstr2_denovo_binder.py`

이 경로는 기존 리간드와 **완전히 무관한 새로운 펩타이드**를 AI로 설계한다.

### 7.1 파이프라인 (Baker Lab 프로토콜)

```
SSTR2 구조 + 핫스팟 잔기
        ↓
  RFdiffusion (백본 설계)
    - contigs: B1-369/0 10-30
    - hotspot: B50,B92,...B302
    - 50 diffusion steps
        ↓
  ProteinMPNN (서열 설계)
    - sampling_temp: 0.2
    - 4 서열/백본
        ↓
  ESMFold (폴딩 검증)
    - pLDDT ≥ 50 통과
        ↓
  검증된 바인더 후보
```

### 7.2 RFdiffusion 결과

| 백본 | 바인더 길이 | 소요 시간 | 상태 |
|------|-----------|----------|------|
| backbone_00 | 22잔기 | 38.4s | 성공 |
| backbone_01 | 14잔기 | 36.8s | 성공 |
| backbone_02 | 11잔기 | 36.1s | 성공 |
| backbone_03 | 17잔기 | 37.5s | 성공 |
| backbone_04 | - | - | 서버 오류 (500) |

### 7.3 ProteinMPNN 결과

4개 백본 x 4개 서열 = **16개** de novo 펩타이드 서열

### 7.4 ESMFold 검증 결과

**16/16 전부 통과** (구조 생성됨):

| 순위 | 서열 | 길이 | pLDDT | 백본 |
|------|------|------|-------|------|
| **1** | `AALARTIAARFRKELEA` | 17 | **81.4** | bb03 |
| **2** | `AALARTIRADFRAQQQA` | 17 | **81.2** | bb03 |
| **3** | `SGLTGGLLALRRYAELARRYLE` | 22 | **80.4** | bb00 |
| **4** | `AAALGLLLFEAAEQ` | 14 | **79.9** | bb01 |
| 5 | `AGLTGGLAAYREYCRLARRLLE` | 22 | 76.9 | bb00 |
| 6 | `AALWQTILTRFRRQQEE` | 17 | 74.7 | bb03 |
| 7 | `MAALGLLLFEYAEQ` | 14 | 73.6 | bb01 |
| 8 | `TPLTGGEAQLVRYASLARRYLE` | 22 | 73.3 | bb00 |

> **pLDDT > 70**: 높은 폴딩 신뢰도. 8/16 서열이 이 기준을 초과함.

### 7.5 서열 패턴 분석

bb03 유래 서열(1, 2위)에서 공통 모티프 발견:
- `AAL` 시작: 소수성 앵커
- `R...R` 패턴: 양전하 잔기 반복 → SSTR2 포켓의 음전하(D122) 보완
- `F/Y` 방향족: Trp8 위치의 방향족 상호작용 모방

---

## 8. 결과 비교 요약

| Arm | 방법 | 후보 수 | 주요 지표 | 상태 |
|-----|------|---------|----------|------|
| 1 | 소분자 (MolMIM+DiffDock) | 40 (15 도킹) | QED=0.94, confidence=-3.0 | 완료 |
| 2 | 펩타이드 변이체 | 13 | Ala scan + 강화 | 분석 완료, 도킹 미수행 |
| 3 | De Novo (RFdiff+MPNN+ESMFold) | 16 | pLDDT=81.4 (최고) | 완료 |

---

## 9. 다음 단계

### 즉시 가능
1. **AlphaFold3 Server 제출**: Top 4 de novo 펩타이드 + SSTR2 서열 → ipTM > 0.71 (Somatostatin 기준) 비교
2. **DiffDock으로 Arm 2 검증**: 변이체 서열을 ESMFold로 구조 예측 → DiffDock 도킹

### 중기
3. **PyRosetta 환경 정비**: FlexPepDock으로 Arm 2 에너지 비교
4. **Hit 확장**: Top de novo 서열을 시드로 2차 RFdiffusion 라운드
5. **Boltz-1/Chai-1**: 로컬 복합체 예측 (AlphaFold3 대안)

### 장기
6. **MD 시뮬레이션**: GROMACS로 바인딩 안정성 및 residence time 예측
7. **실험 검증**: 합성 가능성 평가 (synthetic accessibility), in vitro 바인딩 어세이

---

## 10. 사용 도구 전체 목록

### 로컬 도구 (Conda `bio-tools`)
| 도구 | 용도 |
|------|------|
| Biopython | PDB/CIF 파싱, 바인딩 포켓 분석 |
| RDKit | SMILES-SDF 변환, 화학 정보학 |
| PyMOL | 구조 시각화 |
| FoldMason | 구조 기반 MSA + lDDT |
| PyRosetta | 분자 모델링 (설치 완료, DB 이슈 미해결) |
| Meeko | AutoDock 전처리 |

### NVIDIA NIM API (GPU 불필요)
| 모델 | 파라미터 | 용도 |
|------|---------|------|
| MolMIM | 65.2M | 소분자 생성/최적화 |
| DiffDock | - | Blind molecular docking |
| RFdiffusion | - | De novo 단백질 바인더 설계 |
| ProteinMPNN | - | 역접힘 (backbone→sequence) |
| ESMFold | 650M | 서열→구조 예측 |

---

## 부록: 파일 구조

```
PRST_N_FM/
├── README.md
├── environment-bio-tools.yml
├── molmim.key / ngc.key          # API 키 (.gitignore)
│
├── data/fold_test1/              # AlphaFold3 입력 데이터
│   ├── fold_test1_model_{0-4}.pdb
│   ├── fold_test1_job_request.json
│   └── msas/, templates/
│
├── bionemo/                      # BioNeMo API 클라이언트 & 파이프라인
│   ├── api_base.py               # 공통 API 베이스 클래스
│   ├── molmim_client.py          # MolMIM 클라이언트
│   ├── diffdock_client.py        # DiffDock 클라이언트
│   ├── rfdiffusion_client.py     # RFdiffusion 클라이언트
│   ├── proteinmpnn_client.py     # ProteinMPNN 클라이언트
│   ├── esmfold_client.py         # ESMFold 클라이언트
│   ├── 01-03_*.py                # MolMIM 시나리오 스크립트
│   ├── 04_sstr2_pocket_analysis.py
│   ├── 05_sstr2_smallmol_screen.py
│   ├── 06_sstr2_flexpep_dock.py
│   └── 07_sstr2_denovo_binder.py
│
├── results/sstr2_docking/        # 실험 결과
│   ├── binding_pocket.json
│   ├── sstr2_receptor.pdb
│   ├── arm1_smallmol/
│   ├── arm2_flexpep/
│   └── arm3_denovo/
│
├── experiments/                  # 실험 기록
│   ├── 00_FULL_REPORT.md        # 이 보고서
│   ├── 01_cif_to_pdb.md
│   ├── 02_foldmason_msa.md
│   ├── 03_pymol_visualization.md
│   ├── 04_pyrosetta_setup.md
│   └── 05_sstr2_virtual_screening.md
│
├── scripts/                      # 실행 스크립트
└── docs/                         # 참조 문서
```
