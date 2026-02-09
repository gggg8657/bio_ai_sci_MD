"""
ESMFold API Client
===================
NVIDIA 호스팅 API를 통한 단백질 구조 예측 (MSA 불필요).

사용법:
    from esmfold_client import ESMFoldClient
    client = ESMFoldClient()
    pdb = client.predict("AGCKNFFWKTFTSC")
"""

from pathlib import Path
from typing import Optional
from api_base import NVIDIABaseClient


class ESMFoldClient(NVIDIABaseClient):
    """ESMFold NIM API 클라이언트"""

    BASE_URL = "https://health.api.nvidia.com/v1/biology/nvidia/esmfold"

    def predict(self, sequence: str) -> dict:
        """
        아미노산 서열 → 3D 구조 예측

        Args:
            sequence: 아미노산 서열 (예: "AGCKNFFWKTFTSC")

        Returns:
            {"pdbs": ["PDB content..."], "plddt_scores": [...]}
        """
        payload = {"sequence": sequence}
        return self._post("", payload, timeout=120)

    def predict_and_save(
        self,
        sequence: str,
        output_path: str | Path,
    ) -> dict:
        """구조 예측 후 PDB 파일 저장"""
        result = self.predict(sequence)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # PDB 내용 저장
        pdb_content = None
        if isinstance(result, dict):
            if "pdbs" in result:
                pdb_content = result["pdbs"][0] if result["pdbs"] else None
            elif "pdb" in result:
                pdb_content = result["pdb"]
            elif "output" in result:
                pdb_content = result["output"]

        if pdb_content:
            output_path.write_text(pdb_content)
            print(f"  PDB 저장: {output_path}")

        return result

    def batch_predict(
        self,
        sequences: list[str],
    ) -> list[dict]:
        """여러 서열 배치 예측"""
        results = []
        for i, seq in enumerate(sequences):
            print(f"  [{i+1}/{len(sequences)}] ESMFold 예측: {seq[:30]}...")
            try:
                result = self.predict(seq)
                result["input_sequence"] = seq
                results.append(result)
            except Exception as e:
                print(f"    오류: {e}")
                results.append({"input_sequence": seq, "error": str(e)})
        return results


def get_client(**kwargs) -> ESMFoldClient:
    return ESMFoldClient(**kwargs)


if __name__ == "__main__":
    client = get_client()
    print("ESMFold Client initialized")
    print(f"  Base URL: {client.base_url}")
    print(f"  API Key:  {'***' + client.api_key[-8:] if client.api_key else 'NOT SET'}")
