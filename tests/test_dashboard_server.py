"""Unit and integration tests for Local Visualization Server."""

from http.server import ThreadingHTTPServer
import json
import threading
import time
import urllib.request
import pytest

from bio_arch.dashboard.server import DashboardHandler


@pytest.fixture(scope="module")
def local_server():
    """Start local test server on a free port in a background thread."""
    server = ThreadingHTTPServer(("127.0.0.1", 0), DashboardHandler)
    port = server.server_address[1]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    time.sleep(0.1)
    yield f"http://127.0.0.1:{port}"
    server.shutdown()


def test_server_index_html(local_server: str):
    """Verify that GET / returns the dashboard HTML."""
    req = urllib.request.urlopen(f"{local_server}/")
    assert req.status == 200
    content = req.read().decode("utf-8")
    assert "Bio-Information Architecture" in content


def test_server_api_runs(local_server: str):
    """Verify that GET /api/runs returns valid JSON list of runs."""
    req = urllib.request.urlopen(f"{local_server}/api/runs")
    assert req.status == 200
    data = json.loads(req.read().decode("utf-8"))
    assert "runs" in data
    assert isinstance(data["runs"], list)


def test_server_api_disassemble(local_server: str):
    """Verify that POST /api/disassemble correctly decompiles DNA sequence."""
    payload = json.dumps({
        "sequence": "TTGACACTGATATAATCGTAGGAGGCCATGGCCGCCTAATAA",
        "name": "api_test_operon",
    }).encode("utf-8")

    req = urllib.request.Request(
        f"{local_server}/api/disassemble",
        data=payload,
        headers={"Content-Type": "application/json"},
    )
    res = urllib.request.urlopen(req)
    assert res.status == 200
    data = json.loads(res.read().decode("utf-8"))
    assert data["routine_id"] == "api_test_operon"
    assert "def api_test_operon" in data["decompiled_pseudocode"]
    assert len(data["tokens"]) > 0
