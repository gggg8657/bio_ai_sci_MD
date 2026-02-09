# BioNeMo MolMIM API Client

NVIDIA BioNeMo 플랫폼의 **MolMIM** (Molecular Mutual Information Machine) 모델을 
NVIDIA 호스팅 API로 사용하는 클라이언트입니다. **로컬 GPU 불필요** -- HTTP 요청만 보냅니다.

## BioNeMo란?

NVIDIA BioNeMo는 생체분자 AI 모델 플랫폼입니다. MolMIM은 그 중 하나입니다:

| 모델 | 용도 | API 지원 |
|------|------|----------|
| **MolMIM** | 소분자 생성/최적화 | **구현 완료** |
| ESM-2 / ESMFold | 단백질 언어 모델 / 구조 예측 | 추후 추가 |
| DiffDock | 분자 도킹 | 추후 추가 |
| AlphaFold2 | 단백질 구조 예측 | 추후 추가 |
| MegaMolBART | 분자 생성 | 추후 추가 |

## 빠른 시작

### 1. API 키 발급

1. https://build.nvidia.com/nvidia/molmim-generate 접속
2. NVIDIA 계정으로 로그인
3. **"Get API Key"** 클릭
4. `nvapi-`로 시작하는 키 복사

### 2. 키 설정 (택 1)

```bash
# 방법 A: 프로젝트 루트에 키 파일
echo "nvapi-YOUR_KEY" > ../molmim.key

# 방법 B: .env 파일
cp .env.example .env
# .env에서 NGC_CLI_API_KEY=nvapi-YOUR_KEY 설정

# 방법 C: 환경변수
export NGC_CLI_API_KEY="nvapi-YOUR_KEY"
```

### 3. 의존성 설치

```bash
conda activate bio-tools
pip install requests python-dotenv
```

### 4. 테스트

```bash
cd bionemo/
python molmim_client.py
```

## 사용법

### Python 코드에서

```python
from molmim_client import get_client

client = get_client()

# CMA-ES로 QED 최적화 분자 생성
molecules = client.generate(
    smi="CCO",              # 시드 분자 (에탄올)
    num_molecules=5,         # 생성 수
    algorithm="CMA-ES",      # 최적화 알고리즘
    property_name="QED",     # 최적화 대상
    min_similarity=0.3,      # 최소 유사도
    particles=10,            # CMA-ES 파티클 수
    iterations=3,            # 반복 횟수
)

for mol in molecules:
    print(f"{mol['sample']:50s} QED={mol['score']:.4f}")
```

### 랜덤 샘플링

```python
# 시드 분자 주변에서 랜덤 분자 생성
samples = client.sampling(
    smi="CC(=O)Oc1ccccc1C(=O)O",  # 아스피린
    num_samples=10,
    scaled_radius=1.0,
)
```

## 엔드포인트

### 호스팅 API (build.nvidia.com)

| 엔드포인트 | 설명 | 지원 |
|-----------|------|------|
| `/generate` | CMA-ES 최적화 또는 랜덤 분자 생성 | O |
| `/embedding` | SMILES → 512차원 벡터 | Self-hosted만 |
| `/hidden` | SMILES → 숨은 상태 | Self-hosted만 |
| `/decode` | 숨은 상태 → SMILES 복원 | Self-hosted만 |
| `/sampling` | 잠재공간 샘플링 | Self-hosted만 |

> **참고**: 호스팅 API에서는 `/generate`만 사용 가능합니다.  
> 나머지 엔드포인트는 Docker NIM을 Self-host해야 사용 가능합니다 (GPU 필요: Ampere 이상).

### `/generate` 파라미터

**CMA-ES 최적화:**
```python
client.generate(
    smi="CCO",
    algorithm="CMA-ES",       # CMA-ES 알고리즘
    num_molecules=10,          # 1~100
    property_name="QED",       # "QED" 또는 "plogP"
    minimize=False,            # True: 최소화, False: 최대화
    min_similarity=0.3,        # 0~0.7
    particles=30,              # 2~1000
    iterations=5,              # 1~1000
)
```

**랜덤 샘플링:**
```python
client.generate(
    smi="CCO",
    algorithm="none",          # 랜덤 샘플링
    num_molecules=10,          # 1~100
    particles=10,              # 2~1000
    scaled_radius=1.0,         # 0~2
)
```

## 시나리오 스크립트

| 스크립트 | 설명 |
|---------|------|
| `01_embedding_similarity.py` | 4개 시드 분자별 CMA-ES 생성 + 랜덤 vs 최적화 비교 |
| `02_molecule_generation.py` | 시드 기반 랜덤/CMA-ES/plogP 3가지 모드 생성 |
| `03_property_optimization.py` | 다단계 최적화 (라운드별 최고 분자를 다음 시드로) |

```bash
# 전체 실행
bash ../scripts/run_scenarios.sh

# 개별 실행
python 01_embedding_similarity.py
python 02_molecule_generation.py --smi "c1ccccc1" --num-molecules 10
python 03_property_optimization.py --rounds 3 --num-molecules 10
```

## 모델 정보

- **모델명**: MolMIM-24.03
- **파라미터**: 65.2M
- **아키텍처**: Perceiver Encoder + Transformer Decoder
- **학습 데이터**: ZINC-15 (1.54B 분자)
- **학습 방식**: Mutual Information Machine (MIM) - 비지도 학습
- **최대 입력 길이**: 512 토큰
- **최대 출력 길이**: 128 토큰

## 파일 구조

```
bionemo/
├── README.md              ← 이 문서
├── .env.example           ← API 키 템플릿
├── requirements.txt       ← 의존성
├── __init__.py
├── molmim_client.py       ← MolMIM API 클라이언트 (핵심)
├── 01_embedding_similarity.py
├── 02_molecule_generation.py
└── 03_property_optimization.py
```

## 참고 자료

- [MolMIM NIM 문서](https://docs.nvidia.com/nim/bionemo/molmim/latest/index.html)
- [API Reference](https://docs.api.nvidia.com/nim/reference/nvidia-molmim)
- [build.nvidia.com (API 키 발급)](https://build.nvidia.com/nvidia/molmim-generate)
- [논문: Improving Small Molecule Generation using MIM](https://arxiv.org/abs/2208.09016)
