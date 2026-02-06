# PPT용 요약 자료 – Bio 레포 환경 및 구조 시각화

아래 내용을 슬라이드별로 복사해 사용할 수 있습니다. 표·불릿은 PPT에서 그대로 활용 가능합니다.

---

## 슬라이드 1: 제목 / 목차

**제목:** Bio 레포 개발 환경 및 단백질 구조 시각화 정리

**목차**
1. 통합 Conda 환경 (bio-tools)
2. 구조 파일 변환 (CIF → PDB) 및 시각화
3. 오픈소스 PDB 뷰어 · PyMOL
4. Linear/Git 연동
5. FoldMason 설치 가능 여부

---

## 슬라이드 2: 통합 Conda 환경 (bio-tools)

**목표**  
PyRosetta, Biopython, AutoDock-GPU(Meeko), AlphaFold3·MCP를 염두에 둔 단일/분리 환경 구성.

**구성 요약**

| 항목 | 내용 |
|------|------|
| 환경 이름 | `bio-tools` |
| Python | 3.10–3.12 |
| 주요 패키지 | PyRosetta, Biopython, Meeko(rdkit, scipy, gemmi), PyMOL |
| 생성 | `conda env create -f environment-bio-tools.yml` 후 PyRosetta conda 설치 |
| 검증 | `python scripts/verify_bio_tools_env.py` |

**환경 분리**
- **bio-tools 1개 env:** PyRosetta + Biopython + Meeko(AutoDock-GPU 전처리) ✅
- **AlphaFold3:** 전용 Docker 또는 별도 env 권장 (의존성 충돌 방지)

---

## 슬라이드 3: 구조 파일 변환 (CIF → PDB)

**대상**  
AlphaFold 등에서 받은 **fold_test1** 결과 (mmCIF).

**방법**
- **도구:** Biopython (`MMCIFParser` → `PDBIO`)
- **스크립트:** `scripts/cif_to_pdb.py`
- **실행 예:**  
  `python scripts/cif_to_pdb.py "fold_test1 (1)"`

**결과**
- 변환된 파일: **13개 PDB**
  - 루트: `fold_test1_model_0.pdb` ~ `fold_test1_model_4.pdb` (5개)
  - `templates/`: 템플릿 8개 → 동일 이름 `.pdb`

---

## 슬라이드 4: 오픈소스 PDB 시각화 도구

**조사한 오픈소스 뷰어 (요약)**

| 도구 | 유형 | 라이선스 | 비고 |
|------|------|----------|------|
| **Mol*** | 웹 | Apache 2.0 | PDBe/RCSB/AlphaFold DB 공식 뷰어 |
| **PyMOL (Open-Source)** | 데스크톱 | BSD-like | 업계 표준, 출판·스크립팅 |
| **NGL Viewer** | 웹/JS | MIT | Jupyter(NGLView), 단일 ngl.js |
| **3Dmol.js** | 웹/JS | BSD | URL/선언적 API, 임베드 용이 |
| **Jmol / JSmol** | 데스크톱/웹 | LGPL | Java + HTML5, 교육·다국어 |

상세: `docs/PDB_VISUALIZATION_TOOLS.md`

---

## 슬라이드 5: PyMOL 설치 및 사용

**설치 (bio-tools env)**  
`conda install -c conda-forge pymol-open-source`

**실행**
```bash
conda activate bio-tools
pymol "fold_test1 (1)/fold_test1_model_0.pdb"
```
- 여러 모델: `pymol model_0.pdb model_1.pdb ...`
- 스크립트: `./scripts/run_pymol_pdb.sh` (인자 없으면 model_0.pdb)

**환경**
- WSL: Windows 11 WSLg 또는 X 서버 필요 (GUI)
- 헤드리스: `pymol -c "file.pdb" -d "png out.png; quit"` 로 이미지만 저장 가능

---

## 슬라이드 6: Linear / Git 연동

**Git**
- 브랜치: `feature/bio-tools-conda-env`
- 주요 커밋: env 정의, PyRosetta 설치, AutoDock-GPU 스크립트·문서, AlphaFold 가이드, CIF→PDB·PyMOL

**Linear**
- 이슈 초안: `docs/LINEAR_ISSUES_BIO_TOOLS.md`
- 항목: (1) PyRosetta 설치 (2) AutoDock-GPU 빌드 (3) AlphaFold3 env 정리
- 브랜치/커밋 규칙: `feat: ... [#ISSUE_KEY]`

**AutoDock-GPU**
- 빌드 스크립트: `scripts/build_autodock_gpu.sh`
- 사전 조건: CUDA 툴킷 또는 `nvidia-cuda-dev` 설치 후 실행

---

## 슬라이드 7: FoldMason 설치 가능 여부

**저장소:** https://github.com/steineggerlab/foldmason  
**용도:** 대규모 단백질 구조 다중 정렬 (Multiple Protein Structure Alignment at Scale, *Science* 2026)

**설치 가능: ✅**

| 방법 | 명령/URL |
|------|----------|
| **Conda (권장)** | `conda install -c conda-forge -c bioconda foldmason` |
| **Linux 바이너리** | https://mmseqs.com/foldmason (AVX2/SSE2/ARM64 tar.gz) |
| **macOS** | foldmason-osx-universal.tar.gz |
| **Docker** | 저장소 내 Dockerfile로 빌드 |

**bio-tools와 함께:** Conda로 설치 시 `conda activate bio-tools` 후 `foldmason` 명령 사용 가능.  
**사용 예:** `foldmason easy-msa result.fasta tmpFolder --report-mode 1` (PDB/mmCIF → FASTA 정렬 + HTML 리포트)

---

## 슬라이드 8: 정리 / 다음 단계

**지금까지**
- bio-tools conda 환경 구축 (PyRosetta, Biopython, Meeko, PyMOL)
- fold_test1 mmCIF → PDB 변환 (Biopython)
- 오픈소스 PDB 뷰어 조사 및 PyMOL로 시각화
- AutoDock-GPU 빌드 스크립트·문서, AlphaFold 전용 env 가이드
- Linear 이슈 초안 및 Git 브랜치/커밋 정리

**다음 단계 제안**
- CUDA 설치 후 AutoDock-GPU 실제 빌드 및 Meeko 파이프라인 테스트
- AlphaFold3 전용 Docker/env 구축 (필요 시)
- FoldMason 공식 출처 확인 후 설치·통합 검토

---

## Mermaid 다이어그램 컴파일

- **소스:** `pipeline_orchestration.mermaid` (flowchart)
- **컴파일:** Node.js + `@mermaid-js/mermaid-cli` → SVG 생성
  ```bash
  npx -y @mermaid-js/mermaid-cli@latest -i pipeline_orchestration.mermaid -o docs/pipeline_orchestration.svg
  ```
  또는 `./scripts/compile_mermaid.sh`
- **출력:** `docs/pipeline_orchestration.svg` (PPT/웹에 삽입 가능)
- **문법:** 노드 라벨에 `\n` 또는 괄호만 쓰면 파서 오류 나는 경우 있음 → 큰따옴표로 감싸기: `A["Label text"]`

---

## 참고 문서 (레포 내)

| 문서 | 내용 |
|------|------|
| `ENVIRONMENT.md` | Conda env 생성/활성화/검증 |
| `docs/ENV_COMPATIBILITY.md` | PyRosetta·AutoDock·AlphaFold 호환성 |
| `docs/PDB_VISUALIZATION_TOOLS.md` | PDB 뷰어 상세 |
| `docs/LINEAR_ISSUES_BIO_TOOLS.md` | Linear 이슈 초안·진행 현황 |
