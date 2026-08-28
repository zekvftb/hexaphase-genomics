// Biology as Information Architecture - Interactive Dashboard Client

let zipfChart = null;

document.addEventListener("DOMContentLoaded", () => {
  setupTabs();
  loadRunsList();
  setupDisassembler();
  setupLogicCircuit();
});

// ---------------------------------------------------------------------------
// Tab Navigation
// ---------------------------------------------------------------------------
function setupTabs() {
  const tabButtons = document.querySelectorAll(".tab-btn");
  const tabPanes = document.querySelectorAll(".tab-pane");

  tabButtons.forEach(btn => {
    btn.addEventListener("click", () => {
      tabButtons.forEach(b => b.classList.remove("active"));
      tabPanes.forEach(p => p.classList.remove("active"));

      btn.classList.add("active");
      const targetId = btn.getAttribute("data-tab");
      document.getElementById(targetId).classList.add("active");
    });
  });
}

// ---------------------------------------------------------------------------
// Pipeline Runs Inspector
// ---------------------------------------------------------------------------
async function loadRunsList() {
  const select = document.getElementById("run-select");
  const refreshBtn = document.getElementById("refresh-runs-btn");

  refreshBtn.addEventListener("click", loadRunsList);

  try {
    const res = await fetch("/api/runs");
    const data = await res.json();
    select.innerHTML = "";

    if (!data.runs || data.runs.length === 0) {
      select.innerHTML = "<option value=''>No runs found</option>";
      document.getElementById("report-content").innerText = "No runs executed yet.";
      return;
    }

    data.runs.forEach(r => {
      const opt = document.createElement("option");
      opt.value = r.run_id;
      opt.textContent = `${r.run_id} (${r.status || "done"})`;
      select.appendChild(opt);
    });

    select.addEventListener("change", () => loadRunDetails(select.value));
    // Load first run by default
    loadRunDetails(data.runs[0].run_id);
  } catch (err) {
    console.error("Error loading runs:", err);
  }
}

async function loadRunDetails(runId) {
  if (!runId) return;
  try {
    const res = await fetch(`/api/run/${runId}`);
    const data = await res.json();

    const summary = data.summary || {};
    document.getElementById("stat-status").innerText = (summary.status || "COMPLETED").toUpperCase();
    document.getElementById("stat-findings").innerText = summary.total_findings || 0;
    document.getElementById("stat-interpretations").innerText = summary.total_interpretations || 0;
    
    const timings = summary.timings_seconds || {};
    const totalSec = Object.values(timings).reduce((a, b) => a + b, 0);
    document.getElementById("stat-time").innerText = `${totalSec.toFixed(2)}s`;

    document.getElementById("report-content").innerText = data.report_markdown || "No report markdown found.";
  } catch (err) {
    console.error("Error loading run details:", err);
  }
}

// ---------------------------------------------------------------------------
// Biological Disassembler
// ---------------------------------------------------------------------------
const SAMPLE_LAC = 
  "TTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATT" +
  "TCACACAGGAAACAGCTATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTC" +
  "GTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTC" +
  "GCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGC";

const SAMPLE_PHIX = 
  "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTTGATAAAGCAGGA" +
  "ATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAAATTCGACCTATCCTTGCGCAGCTCG" +
  "AGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTGTCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAAT";

function setupDisassembler() {
  const input = document.getElementById("dna-input");
  const btn = document.getElementById("disassemble-btn");

  document.getElementById("btn-sample-lac").addEventListener("click", () => {
    input.value = SAMPLE_LAC;
  });

  document.getElementById("btn-sample-phix").addEventListener("click", () => {
    input.value = SAMPLE_PHIX;
  });

  // Default to Lac operon sample
  input.value = SAMPLE_LAC;

  btn.addEventListener("click", async () => {
    const seq = input.value.trim().replace(/\s+/g, "");
    if (!seq) return;

    btn.disabled = true;
    btn.innerText = "Analyzing...";

    try {
      // 1. Call Disassembler API
      const disRes = await fetch("/api/disassemble", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence: seq, name: "analyzed_routine" }),
      });
      const disData = await disRes.json();

      // Render tokens
      renderTokenStream(seq, disData.tokens || []);
      document.getElementById("pseudocode-output").textContent = disData.decompiled_pseudocode || "";
      document.getElementById("assembly-output").textContent = (disData.assembly_listing || []).join("\n");

      // 2. Call Linguistics API for Tab 3
      const lingRes = await fetch("/api/linguistics", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence: seq }),
      });
      const lingData = await lingRes.json();
      renderLinguistics(lingData.profile || {});

    } catch (err) {
      console.error("Disassembly error:", err);
    } finally {
      btn.disabled = false;
      btn.innerText = "⚡ Disassemble & Decompile";
    }
  });

  // Initial trigger
  btn.click();
}

function renderTokenStream(rawSeq, tokens) {
  const container = document.getElementById("token-stream");
  container.innerHTML = "";

  if (tokens.length === 0) {
    container.innerHTML = "<span>No canonical biological syntax tokens detected in this window.</span>";
    return;
  }

  tokens.forEach(tok => {
    const span = document.createElement("span");
    span.style.marginRight = "6px";
    span.style.padding = "2px 6px";
    span.style.borderRadius = "3px";
    span.style.display = "inline-block";

    if (tok.token_type.includes("PROMOTER")) {
      span.className = "bg-promoter";
    } else if (tok.token_type === "RBS_SHINE_DALGARNO") {
      span.className = "bg-rbs";
    } else if (tok.token_type === "START_CODON") {
      span.className = "bg-start";
    } else if (tok.token_type === "STOP_CODON") {
      span.className = "bg-stop";
    } else {
      span.className = "bg-cds";
    }

    span.textContent = `[${tok.start}:${tok.label} "${tok.sequence}"]`;
    span.title = JSON.stringify(tok.metadata);
    container.appendChild(span);
  });
}

// ---------------------------------------------------------------------------
// DNA Linguistics & Zipf's Law Chart
// ---------------------------------------------------------------------------
function renderLinguistics(profile) {
  document.getElementById("val-alpha").innerText = profile.codon_zipf_alpha ?? "--";
  document.getElementById("val-r2").innerText = profile.codon_zipf_r2 ?? "--";
  document.getElementById("val-null-alpha").innerText = profile.null_shuffled_alpha ?? "--";

  // Reusable subroutines list
  const list = document.getElementById("subroutine-list");
  list.innerHTML = "";
  const subroutines = profile.top_subroutines || [];
  if (subroutines.length === 0) {
    list.innerHTML = "<li>No recurring subroutines found.</li>";
  } else {
    subroutines.slice(0, 8).forEach(([phrase, count]) => {
      const li = document.createElement("li");
      li.innerHTML = `<strong>${phrase}</strong> (called ${count}x across sequence)`;
      list.appendChild(li);
    });
  }

  // Render Chart.js Zipf plot
  const ranked = profile.codon_rank_distribution || [];
  const labels = ranked.map(r => r.rank);
  const dataPoints = ranked.map(r => r.frequency);

  const ctx = document.getElementById("zipfChart").getContext("2d");
  if (zipfChart) {
    zipfChart.destroy();
  }

  zipfChart = new Chart(ctx, {
    type: "line",
    data: {
      labels: labels,
      datasets: [{
        label: "Codon Rank-Frequency (Observed)",
        data: dataPoints,
        borderColor: "#38bdf8",
        backgroundColor: "rgba(56, 189, 248, 0.1)",
        tension: 0.2,
        fill: true,
      }],
    },
    options: {
      responsive: true,
      scales: {
        x: {
          title: { display: true, text: "Rank (r)", color: "#94a3b8" },
          ticks: { color: "#94a3b8" },
          grid: { color: "#334155" },
        },
        y: {
          title: { display: true, text: "Frequency f(r)", color: "#94a3b8" },
          ticks: { color: "#94a3b8" },
          grid: { color: "#334155" },
        },
      },
      plugins: {
        legend: { labels: { color: "#f8fafc" } },
      },
    },
  });
}

// ---------------------------------------------------------------------------
// Regulatory Logic Gate Simulation
// ---------------------------------------------------------------------------
function setupLogicCircuit() {
  const btnLac = document.getElementById("toggle-lac");
  const btnGlu = document.getElementById("toggle-glu");
  const stateText = document.getElementById("operon-state-text");
  const outputBox = document.getElementById("operon-output-box");

  let lactose = false;
  let glucose = true;

  function updateGate() {
    btnLac.classList.toggle("active", lactose);
    btnLac.textContent = lactose ? "TRUE" : "FALSE";

    btnGlu.classList.toggle("active", glucose);
    btnGlu.textContent = glucose ? "TRUE" : "FALSE";

    // Biological logic: Operon is active ONLY if Lactose is TRUE and Glucose is FALSE
    const isActive = lactose && !glucose;

    if (isActive) {
      stateText.textContent = "TRANSCRIBING (ON)";
      stateText.className = "status-on";
      outputBox.style.borderColor = "var(--accent-green)";
    } else {
      stateText.textContent = "REPRESSED (OFF)";
      stateText.className = "status-off";
      outputBox.style.borderColor = "var(--accent-rose)";
    }
  }

  btnLac.addEventListener("click", () => {
    lactose = !lactose;
    updateGate();
  });

  btnGlu.addEventListener("click", () => {
    glucose = !glucose;
    updateGate();
  });

  updateGate();
}
