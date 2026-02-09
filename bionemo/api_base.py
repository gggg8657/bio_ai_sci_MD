"""
공통 API 베이스 클래스
=====================
모든 NVIDIA NIM API 클라이언트가 상속하는 베이스 클래스.
molmim.key / .env / 환경변수에서 키를 자동 로드한다.
"""

import os
import json
import requests
from pathlib import Path
from typing import Optional


class NVIDIABaseClient:
    """NVIDIA NIM API 공통 베이스 클라이언트"""

    BASE_URL = ""  # 서브클래스에서 오버라이드

    def __init__(
        self,
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
    ):
        self.api_key = api_key or self._load_api_key()
        self.base_url = (base_url or self.BASE_URL).rstrip("/")
        self.headers = {
            "Content-Type": "application/json",
            "Accept": "application/json",
        }
        if self.api_key:
            self.headers["Authorization"] = f"Bearer {self.api_key}"

    @staticmethod
    def _load_api_key() -> str:
        """환경변수 → .env → molmim.key/ngc.key 순서로 API 키 탐색"""
        # 1. 환경변수
        key = os.getenv("NGC_CLI_API_KEY") or os.getenv("NVIDIA_API_KEY")
        if key:
            return key

        # 2. .env 파일
        env_path = Path(__file__).parent / ".env"
        if env_path.exists():
            for line in env_path.read_text().splitlines():
                line = line.strip()
                if "API_KEY=" in line and not line.startswith("#"):
                    val = line.split("=", 1)[1].strip().strip('"').strip("'")
                    if val and "your-" not in val:
                        return val

        # 3. 키 파일 (프로젝트 루트)
        for search in [Path(__file__).parent.parent, Path.cwd()]:
            for name in ["molmim.key", "ngc.key"]:
                key_file = search / name
                if key_file.exists():
                    val = key_file.read_text().strip()
                    if val:
                        return val

        raise ValueError(
            "NVIDIA API 키를 찾을 수 없습니다.\n"
            "방법 1: 환경변수 NGC_CLI_API_KEY 설정\n"
            "방법 2: bionemo/.env 파일에 NGC_CLI_API_KEY=xxx\n"
            "방법 3: 프로젝트 루트에 molmim.key 파일"
        )

    def _post(self, endpoint: str, payload: dict, timeout: int = 120) -> dict:
        """API POST 요청"""
        url = f"{self.base_url}/{endpoint}" if endpoint else self.base_url
        resp = requests.post(url, headers=self.headers, json=payload, timeout=timeout)
        if resp.status_code != 200:
            raise RuntimeError(
                f"API 오류 [{resp.status_code}]: {resp.text[:500]}"
            )
        return resp.json()

    def _post_raw(self, endpoint: str, payload: dict, timeout: int = 120) -> requests.Response:
        """API POST 요청 (raw Response 반환)"""
        url = f"{self.base_url}/{endpoint}" if endpoint else self.base_url
        resp = requests.post(url, headers=self.headers, json=payload, timeout=timeout)
        if resp.status_code != 200:
            raise RuntimeError(
                f"API 오류 [{resp.status_code}]: {resp.text[:500]}"
            )
        return resp
