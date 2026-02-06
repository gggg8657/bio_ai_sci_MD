# Linear MCP 연결 확인 결과

## 1. 설정 상태 ✅

- **위치:** `~/.cursor/mcp.json`
- **서버명:** `linear`
- **엔드포인트:** `https://mcp.linear.app/mcp` (mcp-remote)
- **인증:** `LINEAR_API_KEY` 환경 변수로 설정됨
- **프로젝트 승인:** `mcp-approvals.json`에 `linear-c1daa4a20aaea8cf` 포함 → Linear MCP 사용 승인됨

## 2. 터미널에서 연결 테스트 ✅

아래 명령으로 **수동 실행 시** 정상 연결 확인됨:

```text
npx -y mcp-remote "https://mcp.linear.app/mcp"
# → "Connected to remote server using StreamableHTTPClientTransport"
# → "Proxy established successfully"
```

즉, Linear MCP 서버는 도달 가능하고, mcp-remote 프록시도 정상 동작합니다.

## 3. 에이전트에서 리소스 목록 조회 ⚠️

- **현재:** `list_mcp_resources(server: "linear")` 호출 시 **타임아웃** 발생 (`MCP error -32001: Request timed out`).
- **가능한 원인:**
  - Cursor → MCP → Linear 구간에서 응답 지연
  - Linear 쪽 리소스 목록 응답이 느림
  - mcp-remote/스트리밍 구간에서의 타임아웃 설정

## 4. 권장 사항

1. **Cursor UI에서 Linear 도구 사용**
   - 채팅/에이전트에서 Linear 관련 도구(이슈 목록, 이슈 생성 등)가 보이면, 그쪽으로 이슈 생성/조회를 시도해 보세요.
2. **재시도**
   - 네트워크나 Linear 상태에 따라 일시적일 수 있으므로, 같은 작업을 나중에 다시 시도해 보세요.
3. **보안**
   - `mcp.json`에 API 키가 평문으로 있으므로, 가능하면 Cursor의 env 설정으로 `LINEAR_API_KEY`를 넣고 `mcp.json`에는 `"env": {}`만 두거나, 환경 변수 참조 방식을 사용하는 편이 좋습니다.

## 5. 요약

| 항목           | 상태 |
|----------------|------|
| Linear MCP 설정 | ✅ 등록·승인됨 |
| 터미널 연결 테스트 | ✅ 성공 |
| 에이전트 `list_mcp_resources(linear)` | ⚠️ 타임아웃 (재시도 권장) |

Linear MCP는 **연결되어 있고**, 터미널에서는 정상 동작합니다. 에이전트에서만 리소스 목록 조회가 타임아웃되므로, Cursor에서 Linear 도구가 보이면 그 도구로 먼저 사용해 보시고, 필요하면 같은 호출을 다시 시도해 보시면 됩니다.
