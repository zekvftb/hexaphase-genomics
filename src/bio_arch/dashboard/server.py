"""Local Virtual Visualization Server.

Runs 100% locally on your laptop using Python's built-in HTTP server.
Provides:
- Interactive dashboard UI for inspecting runs, findings, and reports.
- Real-time Biological Disassembler IDE endpoint.
- DNA Linguistics and Zipf rank-frequency visualizer.
- Zero paid APIs, zero external server infrastructure.
"""

from __future__ import annotations

import argparse
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import sys
from urllib.parse import urlparse
import webbrowser

from bio_arch.logger import setup_logger
from bio_arch.modules.disassembler import disassemble_sequence
from bio_arch.modules.linguistics import analyze_linguistic_architecture

logger = setup_logger("bio_arch.dashboard")

STATIC_DIR = Path(__file__).resolve().parent / "static"
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent


class DashboardHandler(SimpleHTTPRequestHandler):
    """HTTP request handler supporting static asset delivery and local REST APIs."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, directory=str(STATIC_DIR), **kwargs)

    def log_message(self, format: str, *args: Any) -> None:
        # Suppress noisy standard HTTP logs on console
        pass

    def send_json(self, data: Any, status: int = 200) -> None:
        """Send JSON response with CORS allowed for local tools."""
        body = json.dumps(data).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)

        # 1. API: List available run directories
        if parsed.path == "/api/runs":
            runs_dir = REPO_ROOT / "outputs" / "runs"
            runs_list = []
            if runs_dir.is_dir():
                for d in sorted(runs_dir.iterdir(), reverse=True):
                    if d.is_dir() and d.name.startswith("run_"):
                        summary_file = d / "final_summary.json"
                        meta: dict[str, Any] = {"run_id": d.name}
                        if summary_file.is_file():
                            try:
                                content = json.loads(summary_file.read_text(encoding="utf-8"))
                                meta["status"] = content.get("status", "unknown")
                                meta["timestamp"] = content.get("timestamp", "")
                                meta["total_findings"] = content.get("total_findings", 0)
                                meta["total_interpretations"] = content.get("total_interpretations", 0)
                            except Exception:
                                pass
                        runs_list.append(meta)
            self.send_json({"runs": runs_list})
            return

        # 2. API: Fetch specific run summary and Markdown report
        if parsed.path.startswith("/api/run/"):
            run_id = parsed.path[len("/api/run/") :]
            run_dir = REPO_ROOT / "outputs" / "runs" / run_id
            if not run_dir.is_dir():
                self.send_json({"error": "Run not found"}, status=404)
                return

            summary_data = {}
            report_text = ""
            summary_file = run_dir / "final_summary.json"
            report_file = run_dir / "final_report.md"

            if summary_file.is_file():
                try:
                    summary_data = json.loads(summary_file.read_text(encoding="utf-8"))
                except Exception:
                    pass

            if report_file.is_file():
                try:
                    report_text = report_file.read_text(encoding="utf-8")
                except Exception:
                    pass

            self.send_json({"summary": summary_data, "report_markdown": report_text})
            return

        # 3. Default: Serve static assets
        super().do_GET()

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        content_length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(content_length).decode("utf-8")
        payload = json.loads(body) if body else {}

        # 1. API: Biological Disassembler
        if parsed.path == "/api/disassemble":
            sequence = payload.get("sequence", "").strip()
            name = payload.get("name", "interactive_routine")
            if not sequence:
                self.send_json({"error": "No sequence provided"}, status=400)
                return

            result = disassemble_sequence(sequence, routine_name=name)
            self.send_json(result.to_dict())
            return

        # 2. API: DNA Linguistics & Zipf analysis
        if parsed.path == "/api/linguistics":
            sequence = payload.get("sequence", "").strip()
            if not sequence:
                self.send_json({"error": "No sequence provided"}, status=400)
                return

            out, findings, interps = analyze_linguistic_architecture(sequence, record_id="user_input")
            self.send_json({
                "profile": out,
                "findings": [f.to_dict() for f in findings],
                "interpretations": [i.to_dict() for i in interps],
            })
            return

        self.send_json({"error": "Endpoint not found"}, status=404)


def start_server(port: int = 8080, open_browser: bool = False) -> None:
    """Start local visualization server on localhost."""
    server_address = ("127.0.0.1", port)
    httpd = ThreadingHTTPServer(server_address, DashboardHandler)
    url = f"http://localhost:{port}"

    print(f"\n=======================================================")
    print(f"  BIOLOGY AS INFORMATION ARCHITECTURE - LOCAL DASHBOARD")
    print(f"=======================================================")
    print(f"  Server running locally at: {url}")
    print(f"  100% Free & Private (Runs on local PC CPU)")
    print(f"  Press Ctrl+C to stop the server anytime.")
    print(f"=======================================================\n")

    if open_browser:
        webbrowser.open(url)

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down dashboard server...")
        httpd.shutdown()


def main() -> None:
    parser = argparse.ArgumentParser(description="Start bio_arch local visualization dashboard")
    parser.add_argument("--port", type=int, default=8080, help="Port to bind (default 8080)")
    parser.add_argument("--open", action="store_true", help="Automatically open browser")

    args = parser.parse_args()
    start_server(port=args.port, open_browser=args.open)


if __name__ == "__main__":
    main()
