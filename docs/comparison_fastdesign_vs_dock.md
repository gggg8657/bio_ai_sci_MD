# Peptide Binder Design: FastDesign vs Independent Mutation + Dock

## 1. 실험 개요

SSTR2(Somatostatin Receptor 2) + SST14(Somatostatin-14) 복합체를 대상으로,
펩타이드 바인더 설계에서 **두 가지 접근법**의 품질과 리소스 효율을 정량 비교한다.

### 대상 시스템

| 구성 요소 | 세부 정보 |
|-----------|----------|
| Receptor  | SSTR2, 369 residues (Chain A after standardization) |
| Peptide   | SST14, 14 residues: `AGCKNFFWKTFTSC` (Chain B after standardization) |
| Disulfide | Cys3 — Cys13 (고정, 변이 불가) |
| Design positions | 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14 (12개, Cys 제외) |

> **참고**: 원본 AlphaFold3 출력에서는 Chain A = Peptide(14 res), Chain B = Receptor(369 res)이나,
> `SSTR2_SST14_demo.ipynb`의 표준화 단계에서 **A = Receptor, B = Peptide**로 통일된다.
> 비교 노트북은 표준화된 `standardized_relaxed.pdb`를 입력으로 사용한다.

---

## 2. 비교 설계 원칙

### 공정한 비교를 위한 핵심 설계

| 원칙 | 설명 |
|------|------|
| **동일 출발점** | 두 접근법 모두 **분리된 receptor PDB + peptide 서열**에서 시작 |
| **동일 후보 수** | 각 20개 |
| **동일 스코어링** | InterfaceAnalyzerMover (dG, dSASA) + stability/PK proxy |
| **동일 모니터링** | ResourceMonitor (CPU, Memory, GPU, Wall-clock time) |

기존 `SSTR2_SST14_demo.ipynb`에서는 Approach A가 AlphaFold3로 이미 결합된 복합체를
"무료"로 사용하여 도킹 비용이 0이었다.
이번 비교에서는 **두 접근법 모두 조립 + 도킹 단계를 포함**하여 비교의 공정성을 확보한다.

---

## 3. 두 접근법 상세

### Approach A: Dock-then-Design (결합 후 변이)

```
receptor PDB + peptide seq
       │
       ▼
  make_complex_pose()     ← 조립 (receptor + peptide)
       │
       ▼
  FlexPepDockingProtocol  ← 초기 도킹 (binding pose 확보)
       │
       ▼
  FastDesign              ← 결합 상태에서 펩타이드 서열 변이
  (TaskFactory 제어:         - Receptor: 완전 고정
   receptor 고정,            - Peptide Cys: 고정 (이황화결합 보존)
   Cys 고정,                 - Design positions: 변이 허용
   design pos만 변이)        - 나머지: repack만
       │
       ▼
  InterfaceAnalyzerMover  ← 스코어링 (dG, dSASA)
```

**특징**:
- FastDesign이 receptor-peptide 상호작용을 **직접 고려**하며 서열 최적화
- 각 후보 생성 시 에너지 최소화 + Monte Carlo 사이클 반복 → **느림**
- 동일 시작점 근처의 local minimum으로 수렴하는 경향 → 서열 다양성 제한

**예상 후보당 소요 시간**: ~457s (FastDesign) + ~6s (FlexPepDock) = **~463s**

### Approach B: Mutate-then-Dock (변이 후 결합)

```
template complex (standardized_relaxed.pdb)
       │
       ▼
  generate_random_mutant()  ← 문자열 수준에서 랜덤 변이 생성
       │                       (design positions에서 Cys 제외 18종 AA 중 선택)
       ▼
  make_complex_pose()       ← template complex 로드 + MutateResidue로 서열 교체
       │                       (backbone은 binding pocket에 유지, sidechain만 변경)
       ▼
  FastRelax (peptide only)  ← receptor 고정, 펩타이드만 relax
       │                       (backbone + sidechain 자유)
       ▼
  FlexPepDockingProtocol    ← 도킹 리파인 (receptor-peptide 인터페이스 재최적화)
       │
       ▼
  InterfaceAnalyzerMover    ← 스코어링 (dG, dSASA)
```

**특징**:
- 변이 서열이 receptor context 없이 **순수 랜덤** 생성 (독립적)
- FastRelax가 peptide만 이동, receptor는 완전 고정 (독립적)
- FastDesign과 달리 서열 변이 시 에너지 피드백 없음
- 서열 다양성이 높음 (랜덤 변이이므로 중복 가능성 낮음)
- FastRelax + FlexPepDock은 FastDesign보다 **훨씬 빠름**

**Template-based positioning이 필요한 이유**:
- FlexPepDock은 **로컬 리파인** 프로토콜이라 펩타이드가 binding pocket 근처에 있어야 함
- `pose_from_sequence()`로 만든 확장 구조는 receptor와 너무 멀어 작동 불가
- 대신 원본 복합체의 backbone 위치를 활용하고 sidechain만 교체

**예상 후보당 소요 시간**: ~12s (FastRelax) + ~6s (FlexPepDock) = **~18s**

---

## 4. 스코어링 지표

### 결합 품질 (Quality)

| 지표 | 도구 | 해석 |
|------|------|------|
| **dG (REU)** | InterfaceAnalyzerMover | 인터페이스 결합 에너지. **낮을수록** 좋음 |
| **dSASA (Å²)** | InterfaceAnalyzerMover | 매몰 표면적. **높을수록** 좋음 |
| **rank_score** | 계산식 | `(-dG) - 0.5 * cleavage_risk - 1.0 * pk_penalty` |

### 안정성/PK Proxy

| 지표 | 계산 | 해석 |
|------|------|------|
| cleavage_risk | 2.0 * (K+R) + 1.0 * (F+Y+W) | 프로테아제 절단 위험 |
| hydrophobic_fraction | hydrophobic AA / length | 소수성 비율 |
| net_charge_proxy | (K+R+H) - (D+E) | 순 전하 |
| pk_penalty | 5.0 * max(0, hyd_frac - 0.5) + 0.5 * |net_charge| | PK 패널티 |

### 리소스 사용량 (Performance)

| 지표 | 도구 | 설명 |
|------|------|------|
| Wall-clock time | `time.perf_counter()` | 실제 경과 시간 |
| CPU % | `psutil.Process.cpu_percent()` | 프로세스 CPU 사용률 |
| RSS Memory (MB) | `psutil.Process.memory_info().rss` | 물리 메모리 사용량 |
| GPU Utilization % | `pynvml` | GPU 코어 사용률 (해당 시) |
| GPU Memory (MB) | `pynvml` | GPU VRAM 사용량 (해당 시) |

---

## 5. 시각화 출력

노트북 실행 시 `comparison_results/` 디렉토리에 다음 파일이 생성된다:

| 파일명 | 내용 |
|--------|------|
| `comparison_all_candidates.csv` | 40개 후보 전체 데이터 (A 20 + B 20) |
| `resource_timeseries.csv` | 0.5초 간격 리소스 샘플링 시계열 |
| `quality_comparison.png` | dG / dSASA / Time box plot (3패널) |
| `resource_timeseries.png` | CPU / Memory / GPU 시계열 (4패널) |
| `efficiency_frontier.png` | dG vs Wall Time scatter plot |
| `per_candidate_resources.png` | 후보별 Time / CPU / Memory bar chart |
| `diversity_analysis.png` | Hamming distance 분포 + Unique sequences |

---

## 6. 코드 출처 및 재활용

| 함수/클래스 | 출처 | 비고 |
|------------|------|------|
| `build_task_factory()` | `SSTR2_SST14_demo.ipynb` | Approach A용 |
| `make_complex_pose()` | 신규 구현 (template-based, `MutateResidue` 활용) | 양쪽 공유 |
| `analyze_interface()` | `SSTR2_SST14_demo.ipynb` | 양쪽 공유 |
| `stability_pk_proxy_scores()` | `SSTR2_SST14_demo.ipynb` | 양쪽 공유 |
| `ResourceMonitor` | 신규 구현 | `psutil` + `pynvml` 기반 context manager |
| `generate_random_mutant()` | 신규 구현 | Approach B용 랜덤 변이 생성 |
| `mutate_and_dock_single()` | 신규 구현 | Approach B 전체 파이프라인 |

---

## 7. 실행 환경

```bash
# 환경 활성화
conda activate bio-tools

# 추가 패키지 설치 (최초 1회)
pip install psutil pynvml

# 또는 environment 파일로 일괄 업데이트
conda env update -f environment-bio-tools.yml

# 실행
cd notebooks/
jupyter notebook comparison_fastdesign_vs_dock.ipynb
```

### 전제 조건
- `standardized_relaxed.pdb`가 `notebooks/` 디렉토리에 존재해야 함
  (→ `SSTR2_SST14_demo.ipynb`의 Section 1-6 실행 결과물)
- PyRosetta 라이선스가 유효해야 함

---

## 8. 검증 체크리스트

### Chain 일관성 검증

| 단계 | Chain 1 | Chain 2 | 확인 |
|------|---------|---------|------|
| `standardized_relaxed.pdb` 로드 | Receptor (369 res) | Peptide (14 res) | OK |
| `receptor_only.pdb` 추출 | Receptor (369 res) | — | OK (Chain A 필터링) |
| `make_complex_pose()` 결과 | Receptor (chain 1) | Peptide (chain 2) | OK (append_pose_by_jump) |
| `InterfaceAnalyzerMover(1)` | Jump 1 기준 분석 | — | OK (chain 1 vs 2) |
| `get_peptide_seq(pose, 2)` | — | chain 2 서열 추출 | OK |

### Approach A 검증

- [x] FlexPepDock 초기 도킹 후 docked_initial.pdb 저장 → 모든 후보가 동일 시작점
- [x] FastDesign에 TaskFactory 적용: receptor 고정, Cys 고정, design positions만 변이 허용
- [x] 각 후보별 고유 Rosetta random seed (SEED_BASE + k * 100)
- [x] ResourceMonitor로 전체 + 후보별 리소스 측정
- [x] 초기 도킹 시간이 total time에 포함됨

### Approach B 검증

- [x] `generate_random_mutant()`: Cys 위치(3, 13)는 design positions에 없으므로 변이 불가
- [x] `AA_NO_CYS`: 19종 (Cys 제외) — Cys로의 변이 자체를 차단
- [x] `mutate_and_dock_single()`: template 로드 → MutateResidue → FastRelax (peptide only) → FlexPepDock
- [x] Template-based positioning: 원본 backbone 유지 + sidechain만 교체 → FlexPepDock 작동 보장
- [x] FastRelax MoveMap: receptor 완전 고정, peptide만 backbone + sidechain 자유
- [x] 중복 서열 방지: `seen_seqs_b` set으로 최대 10회 재시도
- [x] ResourceMonitor로 전체 + 후보별 리소스 측정

### 공통 검증

- [x] Baseline 측정: 원본 서열로 assemble + FlexPepDock → 기준 dG/dSASA
- [x] 양쪽 모두 동일한 `InterfaceAnalyzerMover` 설정 (jump=1, pack_separated=True)
- [x] 양쪽 모두 동일한 `stability_pk_proxy_scores()` 적용
- [x] 양쪽 모두 동일한 `rank_score` 공식: `(-dG) - 0.5 * cleavage_risk - 1.0 * pk_penalty`

### FlexPepDock 제약 사항

- FlexPepDock은 **로컬 리파인** 프로토콜 — 펩타이드가 binding pocket 근처에 있어야 함
- `pose_from_sequence()`로 만든 확장 구조는 receptor와 너무 멀어 center-of-mass 에러 발생
- 해결: `make_complex_pose()`가 원본 복합체를 template로 사용, `MutateResidue`로 sidechain만 교체
- 이 방식은 Approach B의 '독립적 변이' 특성을 유지함 (서열 생성은 랜덤, relax는 peptide-only)

---

## 9. 예상 결과 및 가설

| 차원 | Approach A (Dock→Design) | Approach B (Mutate→Dock) |
|------|--------------------------|--------------------------|
| **속도** | 느림 (~150+ min) | 빠름 (~6-10 min) |
| **dG 품질** | 우수 (receptor context 고려) | 불확실 (blind mutation) |
| **서열 다양성** | 낮음 (local minimum 수렴) | 높음 (랜덤 변이) |
| **메모리** | 높음 (FastDesign 내부 상태) | 낮음 |
| **확장성** | 후보 수 증가 시 선형 증가 | 후보 수 증가 시 선형 증가 (기울기 낮음) |

### 핵심 가설
> Approach B는 Approach A 대비 **~25x 빠르지만**, dG 품질이 떨어질 수 있다.
> 만약 Approach B의 Top-5 후보 중 Approach A의 평균 dG 이하인 후보가 존재한다면,
> **B로 대량 탐색 → A로 정밀 최적화**하는 2단계 전략이 유효하다.
