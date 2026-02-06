# FoldMason 설치 가능 여부

**저장소:** https://github.com/steineggerlab/foldmason  
**용도:** 대규모 단백질 구조 다중 정렬 (Multiple Protein Structure Alignment at Scale)  
**논문:** Gilchrist CLM, Mirdita M, Steinegger M. *Science* (2026), doi: 10.1126/science.ads6733

---

## 설치 가능: ✅

### 1. Conda (Linux / macOS, 권장)

```bash
conda install -c conda-forge -c bioconda foldmason
```

- **bio-tools와 같은 env:** `conda activate bio-tools` 후 위 명령 실행 가능.
- 의존성 충돌 시 전용 env 생성 후 설치.

### 2. 사전 빌드 바이너리 (Linux / macOS)

- **Linux AVX2:**  
  `wget https://mmseqs.com/foldmason/foldmason-linux-avx2.tar.gz; tar xvzf foldmason-linux-avx2.tar.gz; export PATH=$(pwd)/foldmason/bin/:$PATH`
- **Linux SSE2:**  
  `wget https://mmseqs.com/foldmason/foldmason-linux-sse2.tar.gz` 후 동일하게 압축 해제 및 PATH 추가.
- **Linux ARM64:**  
  `foldmason-linux-arm64.tar.gz`
- **macOS:**  
  `foldmason-osx-universal.tar.gz`

다운로드 목록: https://mmseqs.com/foldmason

### 3. Docker

- 저장소에 `Dockerfile` 포함. `docker build` 후 컨테이너에서 실행 가능.

---

## 사용 예 (Quick Start)

- **다중 정렬 (PDB/mmCIF):**
  ```bash
  foldmason easy-msa result.fasta tmpFolder --report-mode 1
  ```
- **예제 실행:**
  ```bash
  foldmason easy-msa ./lib/foldseek/example/d* example.fasta tmpFolder --report-mode 1
  ```
- **출력:** FASTA 정렬 (`_aa.fa`, `_3di.fa`), Newick 트리, `--report-mode 1` 시 HTML 시각화.

---

## bio-tools와 함께 쓰기

- Conda로 설치하면 `bio-tools` env에서 `foldmason` 명령 바로 사용 가능.
- 변환해 둔 PDB(`fold_test1 (1)/fold_test1_model_*.pdb`)로 예시 실행 가능:
  ```bash
  conda activate bio-tools
  foldmason easy-msa "fold_test1 (1)/fold_test1_model_"*.pdb result_foldmason tmpFolder --report-mode 1
  ```
  (실제 인자는 쉘에서 PDB 목록이 맞게 확장되는지 확인 후 사용.)

---

## 참고

- **웹 서버:** https://search.foldseek.com/foldmason
- **모듈:** `easy-msa`, `structuremsa`, `msa2lddt`, `refinemsa`, `msa2lddtreport` 등. 자세한 옵션은 README 및 `foldmason --help` 참고.
