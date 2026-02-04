# Linear 이슈 초안 – bio-tools 환경 및 연동

아래 이슈를 Linear에 생성한 뒤, 작업 시 **브랜치명**과 **커밋 메시지**에 이슈 키를 붙여 사용하세요.

---

## 진행 현황

| # | 이슈 | 상태 | 비고 |
|---|------|------|------|
| 1 | PyRosetta 설치 | ✅ 완료 | `verify_bio_tools_env.py` 통과 |
| 2 | AutoDock-GPU 빌드 | ⏳ 스크립트+클론 완료, CUDA 설치 후 `./scripts/build_autodock_gpu.sh` 실행 | |
| 3 | AlphaFold3 env 정리 | 대기 | |

---

## 이슈 1: [Infra] bio-tools conda 환경에 PyRosetta 설치

**Type:** Engineering / Code  
**Project:** (기존 프로젝트 또는 Infra)

**Description:**
- **Background:** bio-tools env에 Biopython·Meeko는 설치 완료. PyRosetta(~1.4GB)는 미설치.
- **Goal:** `conda activate bio-tools` 후 PyRosetta import 및 `pyrosetta.init()` 동작 확인.
- **Definition of Done:**
  - `conda install -c https://conda.rosettacommons.org -c conda-forge pyrosetta` 완료
  - `python scripts/verify_bio_tools_env.py` 에서 PyRosetta OK 출력
- **Estimate:** 1h (다운로드 포함)

**Run:**
```bash
conda activate bio-tools
conda install -y -c https://conda.rosettacommons.org -c conda-forge pyrosetta
python scripts/verify_bio_tools_env.py
```

**Branch:** `feature/bio-tools-pyrosetta`  
**Commit:** `feat: add PyRosetta to bio-tools env [#ISSUE_KEY]`

---

## 이슈 2: [Infra] AutoDock-GPU 바이너리 빌드 및 PATH 설정 (WSL)

**Type:** Engineering / Code  
**Project:** (동일)

**Description:**
- **Background:** AutoDock-GPU는 C++ 빌드 바이너리. Meeko는 이미 bio-tools에 있음.
- **Goal:** WSL에서 AutoDock-GPU 소스 빌드 후, bin/autodock_gpu_*wi 를 PATH 또는 스크립트에서 호출 가능하게 함.
- **Definition of Done:**
  - CUDA/OpenCL 개발 패키지 설치, 저장소 클론, `make DEVICE=GPU` 성공
  - `bin/autodock_gpu_64wi` (또는 128wi) 실행 가능
  - 문서(docs/ENV_COMPATIBILITY.md 또는 ENVIRONMENT.md)에 경로/사용법 반영
- **Estimate:** 2h

**Branch:** `feature/autodock-gpu-build-wsl`  
**Commit:** `feat: build AutoDock-GPU binary and document path [#ISSUE_KEY]`

---

## 이슈 3: [Infra] AlphaFold3 전용 환경(Docker/env) 정리

**Type:** Engineering / Code  
**Project:** (동일)

**Description:**
- **Background:** AlphaFold2/3는 무거운 의존성으로 bio-tools와 같은 env에 두기 부적합. 전용 분리 권장.
- **Goal:** AlphaFold3(또는 AF2) 실행을 위한 전용 Docker 또는 conda env 선택 및 문서화.
- **Definition of Done:**
  - 선택: 공식 Docker 사용 vs 전용 conda env (environment-alphafold.yml)
  - ENVIRONMENT.md 또는 docs/ENV_COMPATIBILITY.md에 AlphaFold 섹션 추가
- **Estimate:** 2h (조사+문서)

**Branch:** `chore/alphafold-env-docs`  
**Commit:** `docs: add AlphaFold3 dedicated env guidance [#ISSUE_KEY]`

---

## Git 브랜치/커밋 규칙 (이번 계획 연동)

| 작업 | 브랜치 | 커밋 예시 |
|------|--------|------------|
| bio-tools env + 문서 일괄 | `feature/bio-tools-conda-env` | `feat: add bio-tools conda env and compatibility docs` |
| PyRosetta 설치 | (같은 브랜치 또는 `feature/bio-tools-pyrosetta`) | `feat: install PyRosetta in bio-tools [#XXX]` |
| AutoDock-GPU 빌드 | `feature/autodock-gpu-build-wsl` | `feat: build AutoDock-GPU and document path [#XXX]` |
| AlphaFold env 문서 | `chore/alphafold-env-docs` | `docs: AlphaFold3 env guidance [#XXX]` |

Linear 이슈 생성 후 `#ISSUE_KEY` 부분을 실제 키(예: INFRA-12)로 치환하면 됩니다.
