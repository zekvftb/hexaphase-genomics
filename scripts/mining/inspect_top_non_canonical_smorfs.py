"""Profile & Display Top Non-Canonical smORF Candidates.

Ingests outputs/non_canonical_smorf_candidates.csv, performs multi-parametric ranking
combining translation initiation efficiency, Kozak strength, viroporin/TM signatures,
and statistical null z-scores, and generates an executive profile table and functional dossier.
"""

from __future__ import annotations

import csv
import io
import json
import math
from pathlib import Path
import sys

_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from bio_arch.modules.non_canonical_miner import KYTE_DOOLITTLE


def compute_multiparametric_rank(row: dict) -> float:
    """Compute a multi-parametric score balancing initiation, biophysics, and null z-score."""
    # 1. Translation Initiation Component (0.0 to 1.0)
    eff = float(row.get("Initiation_Efficiency", 0.05))
    kozak = row.get("Kozak_Strength", "Weak")
    kozak_weight = 1.0 if kozak == "Optimal" else (0.75 if "Strong" in kozak else 0.4)
    initiation_score = (eff * 0.7 + kozak_weight * 0.3)

    # 2. Viroporin & Hydropathy Component (0.0 to 1.0)
    max_tm = float(row.get("Max_TM_Hydropathy", 0.0))
    norm_tm = min(1.0, max(0.0, (max_tm - 0.5) / 2.0))
    mu_h = float(row.get("Hydrophobic_Moment", 0.0))
    norm_mu_h = min(1.0, max(0.0, mu_h / 0.5))
    has_tm = 1.0 if row.get("Has_Transmembrane_Domain", "False") == "True" else 0.0
    biophysical_score = (norm_tm * 0.4 + norm_mu_h * 0.3 + has_tm * 0.3)

    # 3. Host Adaptation (0.0 to 1.0)
    cai = float(row.get("Host_CAI", 0.70))

    # 4. Statistical Null z-Score Component (0.0 to 1.0)
    z_val = float(row.get("Z_Score", 0.0))
    norm_z = min(1.0, max(0.0, z_val / 6.0))

    # Weighted composite
    composite = (
        initiation_score * 0.25 +
        biophysical_score * 0.35 +
        cai * 0.20 +
        norm_z * 0.20
    )
    return round(composite, 4)


def inspect_and_display_top_smorfs(base_dir: Path | None = None) -> str:
    root = base_dir if base_dir is not None else Path(__file__).parent.parent.parent
    csv_file = root / "outputs" / "non_canonical_smorf_candidates.csv"

    if not csv_file.is_file():
        raise FileNotFoundError(f"Missing {csv_file}")

    with open(csv_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)

    # Filter for candidates within 30-120 aa and CAI >= 0.70
    valid_cands = []
    for r in all_rows:
        length = int(r["Length_aa"])
        cai = float(r["Host_CAI"])
        if 30 <= length <= 120 and cai >= 0.70:
            r["composite_rank_score"] = compute_multiparametric_rank(r)
            valid_cands.append(r)

    # Collapse nested starts sharing the same (Accession, Parent_Gene, Frame_Offset, End_in_Parent_bp)
    locus_dict = {}
    for r in valid_cands:
        locus_key = (r["Accession"], r["Parent_Gene"], r["Frame_Offset"], r["End_in_Parent_bp"])
        if locus_key not in locus_dict or r["composite_rank_score"] > locus_dict[locus_key]["composite_rank_score"]:
            locus_dict[locus_key] = r

    distinct_locus_cands = list(locus_dict.values())
    # Sort descending by composite score
    distinct_locus_cands.sort(key=lambda x: x["composite_rank_score"], reverse=True)

    # Select top candidates ensuring family diversity
    top_5 = []
    seen_accs = set()
    for r in distinct_locus_cands:
        acc = r["Accession"]
        if acc not in seen_accs:
            seen_accs.add(acc)
            top_5.append(r)
        if len(top_5) == 5:
            break

    # If fewer than 5 distinct accessions, backfill
    if len(top_5) < 5:
        for r in distinct_locus_cands:
            if r not in top_5:
                top_5.append(r)
            if len(top_5) == 5:
                break

    # Generate Markdown Report
    lines = [
        "# 🧬 Top Non-Canonical & Near-Cognate smORF Candidates",
        "## Empirical Translation Plausibility & Viroporin / Membrane-Insertion Screening",
        "",
        "---",
        "",
        "## 1. Multi-Parametric Ranking Table",
        "",
        "| Rank | Organism (Accession) | Primary CDS (Frame) | Initiator (Kozak) | Length & MW | Peak Hydropathy (TMs) | $\\mu_H$ | Host CAI | $z$-Score |",
        "| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]

    for rank, c in enumerate(top_5, 1):
        kozak = c["Kozak_Strength"]
        tm_desc = f"{c['Max_TM_Hydropathy']} ({'1 TM' if c['Has_Transmembrane_Domain'] == 'True' else '0 TM'})"
        lines.append(
            f"| **#{rank}** | {c['Organism'][:18]} (`{c['Accession']}`) | `{c['Parent_Gene']}` (+{c['Frame_Offset'].lstrip('+')}) | `{c['Start_Codon']}` ({kozak}) | **{c['Length_aa']} aa** ({c['MW_kDa']} kDa) | {tm_desc} | **{c['Hydrophobic_Moment']}** | **{c['Host_CAI']}** | **$z = {c['Z_Score']}$** |"
        )

    lines.extend([
        "",
        "---",
        "",
        "## 2. In-Depth Biophysical & Functional Profiles: Top 3 Candidates",
        "",
    ])

    for rank, c in enumerate(top_5[:3], 1):
        peptide = c["Peptide_Sequence"]
        f_offset = c["Frame_Offset"].lstrip("+")
        start_bp = int(c["Start_in_Parent_bp"])
        end_bp = int(c["End_in_Parent_bp"])

        # Detect TM core coordinates
        tm_core_desc = "None detected (soluble profile)"
        if c["Has_Transmembrane_Domain"] == "True":
            # Find the peak 19-aa window
            best_win_score = -99.0
            best_win_idx = 0
            for i in range(len(peptide) - 19 + 1):
                win = peptide[i : i + 19]
                score = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in win) / 19.0
                if score > best_win_score:
                    best_win_score = score
                    best_win_idx = i

            tm_start_aa = best_win_idx + 1
            tm_end_aa = best_win_idx + 19
            tm_start_nt = start_bp + best_win_idx * 3
            tm_end_nt = start_bp + (best_win_idx + 19) * 3
            tm_core_desc = f"Residues {tm_start_aa}..{tm_end_aa} (nt {tm_start_nt}..{tm_end_nt} in primary CDS, Mean Hydropathy = {round(best_win_score, 2)})"

        # Functional assessment
        if "NC_001401" in c["Accession"]:  # AAV2
            func_assessment = (
                "Nested alternative reading frame within the AAV2 VP1 capsid gene. Displays moderate translation efficiency "
                "with an amphipathic alpha-helical profile, characteristic of dependoparvoviral auxiliary assembly microproteins "
                "or membrane-permeabilizing factors required for endosomal escape."
            )
        elif "NC_003977" in c["Accession"]:  # HBV
            func_assessment = (
                "Alternative reading frame overlapping the HBV polymerase catalytic core. Exhibits significant host codon adaptation "
                "and an amphipathic hydrophobic moment, resembling hepatitis B accessory microproteins that modulate host immune signaling "
                "or participate in viral nucleocapsid envelopment."
            )
        elif "NC_001422" in c["Accession"]:  # PhiX174
            func_assessment = (
                "Alternative reading frame in phiX174 with a dense hydrophobic transmembrane core. Functionally resembles a "
                "bacteriophage holin / pinholin membrane-disruptive channel peptide that facilitates host cell lysis."
            )
        else:
            func_assessment = (
                "Unannotated viral alternative reading frame with high host codon adaptation and structured biophysical features. "
                "Represents a high-priority target for experimental ribosome profiling and epitope-tagged validation."
            )

        lines.extend([
            f"### Candidate #{rank}: {c['Organism']} ({c['Accession']})",
            f"- **Overlapping Primary CDS:** `{c['Parent_Gene']}` ({c['Parent_Product']}), Frame +{f_offset}",
            f"- **Genomic Coordinates in Primary CDS:** nt {start_bp}..{end_bp} ({c['Length_aa']} aa)",
            f"- **Translation Initiator:** `{c['Start_Codon']}` ({c['Start_Codon_Type']}) with **{c['Kozak_Strength']}** Kozak Context",
            f"- **Biophysical State:** Molecular Weight = **{c['MW_kDa']} kDa** | pI = **{c['Isoelectric_Point']}** | Net Charge (pH 7.4) = **{c['Net_Charge_pH74']}**",
            f"- **Membrane / Viroporin Domain:** {tm_core_desc}",
            f"- **Amphipathic Hydrophobic Moment ($\\mu_H$):** **{c['Hydrophobic_Moment']}**",
            f"- **Statistical Null Significance:** **$z = {c['Z_Score']}$**, **$p = {c['Empirical_P_Value']}$** ($N=500$ shuffles)",
            f"- **Functional Assessment:** {func_assessment}",
            "",
            "```fasta",
            f">smORF_{c['Accession']}_{c['Parent_Gene']}_F+{f_offset}|{c['Start_Codon']}|{c['Length_aa']}aa|MW={c['MW_kDa']}kDa|CAI={c['Host_CAI']}|muH={c['Hydrophobic_Moment']}",
            peptide,
            "```",
            "",
            "---",
            "",
        ])

    report_md = "\n".join(lines)
    return report_md


if __name__ == "__main__":
    report = inspect_and_display_top_smorfs()
    print(report)
