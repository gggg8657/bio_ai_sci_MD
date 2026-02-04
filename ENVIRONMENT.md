# ENVIRONMENT.md – Repo Conda 환경

## 기본 환경: `bio` (기존)

- **ENV_NAME:** bio
- **PYTHON:** 3.10–3.12
- **CUDA:** cpu
- **CREATE:** `conda env create -f environment.yml`
- **ACTIVATE:** `conda activate bio`
- **VERIFY:** `python -V` / `python -c "import sys; print(sys.version)"`
- **RUN:** (프로젝트별)
- **NOTES:** Windows·WSL 동기화용 단일 profile. 상세는 `CONDA_SETUP.md`.

---

## 신규 환경: `bio-tools` (PyRosetta + Biopython + AutoDock-GPU 워크플로)

- **ENV_NAME:** bio-tools
- **PYTHON:** 3.10–3.12
- **CUDA:** cpu (AutoDock-GPU 바이너리는 시스템 CUDA 사용)
- **CREATE:**
  ```bash
  conda env create -f environment-bio-tools.yml
  conda activate bio-tools
  # PyRosetta (선택, ~1.4GB, 학술 라이선스). conda 권장:
  conda install -y -c https://conda.rosettacommons.org -c conda-forge pyrosetta
  # 또는 pip quarterly wheel:
  # pip install pyrosetta --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release
  ```
  참고: Anaconda 채널 ToS 미수락 시 `conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main` 및 `.../pkgs/r` 실행 후 env 생성.
- **ACTIVATE:** `conda activate bio-tools`
- **VERIFY (최소 2개):**
  - `python -V`
  - `python -c "import pyrosetta; pyrosetta.init(); print('PyRosetta OK')"`
  - `python -c "import Bio; print('Biopython', Bio.__version__)"`
  - `python -c "import meeko; print('Meeko OK')"`
- **RUN:**
  - PyRosetta/Biopython 스크립트: `python your_script.py`
  - AutoDock-GPU: Meeko로 전처리 후 시스템에 빌드한 `autodock_gpu_*wi` 실행 (PATH 또는 경로 지정)
- **TEST/SMOKE:**
  - `python scripts/verify_bio_tools_env.py`
  - 또는 `python -c "import Bio, meeko; print('Bio+Meeko OK')"` (PyRosetta 없이)
- **NOTES:**
  - PyRosetta 학술/비영리 라이선스.
  - AutoDock-GPU 바이너리·AutoGrid는 별도 빌드 후 PATH 또는 스크립트 경로 지정. 호환성·빌드 요약은 `docs/ENV_COMPATIBILITY.md`.
  - AlphaFold3는 별도 Docker 또는 전용 env 권장.
