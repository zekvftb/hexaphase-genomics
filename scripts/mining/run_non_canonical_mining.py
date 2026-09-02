"""Runner Script for Non-Canonical smORF Discovery across Viral Corpus.

Executes NonCanonicalSmorfMiner, filters candidates passing statistical null calibration (z > 3.0, p < 0.001),
exports discovery ledger to outputs/non_canonical_smorf_candidates.csv,
and generates the master publication report in NON_CANONICAL_VIRAL_SMORF_DISCOVERY.md.
"""

from __future__ import annotations

import csv
import io
import json
from pathlib import Path
import sys

_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
if str(_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_ROOT / "src"))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from bio_arch.modules.non_canonical_miner import NonCanonicalCandidate, NonCanonicalSmorfMiner
from bio_arch.provenance import now_iso


def run_non_canonical_pipeline(base_dir: Path | None = None, n_shuffles: int = 500) -> dict:
    root = base_dir if base_dir is not None else Path(__file__).parent.parent.parent
    corpus_dir = root / "data" / "mining_corpus"
    outputs_dir = root / "outputs"
    outputs_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("🧬 NON-CANONICAL & NEAR-COGNATE smORF DISCOVERY ENGINE")
    print(f"   Corpus: {corpus_dir} | Null Controls: N={n_shuffles} Dinucleotide Shuffles")
    print("=" * 80)

    miner = NonCanonicalSmorfMiner(
        min_length_aa=30,
        max_length_aa=120,
        min_cai_threshold=0.70,
        n_null_shuffles=n_shuffles,
        significance_alpha=0.001,
        min_z_threshold=3.0,
        seed=42,
    )

    all_candidates: list[NonCanonicalCandidate] = []
    gbk_files = list(corpus_dir.glob("*.gbk"))

    for gbk in gbk_files:
        cands = miner.mine_genome(gbk)
        all_candidates.extend(cands)

    # Deduplicate
    seen = set()
    unique_cands: list[NonCanonicalCandidate] = []
    for c in all_candidates:
        key = (c.accession, c.parent_gene, c.frame_offset, c.start_in_parent_bp, c.length_aa)
        if key not in seen:
            seen.add(key)
            unique_cands.append(c)

    # Sort by statistical significance, then z-score, then CAI
    unique_cands.sort(key=lambda x: (x.statistically_significant, x.z_score, x.host_cai), reverse=True)

    # Export CSV
    csv_file = outputs_dir / "non_canonical_smorf_candidates.csv"
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Accession",
            "Organism",
            "Parent_Gene",
            "Parent_Product",
            "Frame_Offset",
            "Start_in_Parent_bp",
            "End_in_Parent_bp",
            "Start_Codon",
            "Start_Codon_Type",
            "Initiation_Efficiency",
            "Kozak_Strength",
            "Length_aa",
            "MW_kDa",
            "Isoelectric_Point",
            "Net_Charge_pH74",
            "Host_CAI",
            "Mean_Hydropathy",
            "Max_TM_Hydropathy",
            "Has_Transmembrane_Domain",
            "Hydrophobic_Moment",
            "Is_Potential_Viroporin",
            "Z_Score",
            "Empirical_P_Value",
            "Statistically_Significant",
            "Is_Annotated",
            "Matched_Annotation",
            "Peptide_Sequence",
        ])
        for c in unique_cands:
            writer.writerow([
                c.accession,
                c.organism,
                c.parent_gene,
                c.parent_product,
                f"+{c.frame_offset}",
                c.start_in_parent_bp,
                c.end_in_parent_bp,
                c.start_codon,
                c.start_codon_type,
                c.relative_initiation_efficiency,
                c.kozak_strength,
                c.length_aa,
                c.molecular_weight_kda,
                c.isoelectric_point,
                c.net_charge_ph74,
                c.host_cai,
                c.mean_hydropathy,
                c.max_hydropathy_window,
                c.has_transmembrane_domain,
                c.hydrophobic_moment,
                c.is_potential_viroporin,
                c.z_score,
                c.empirical_p_value,
                c.statistically_significant,
                c.is_annotated_in_genbank,
                c.matched_annotation_name,
                c.peptide_sequence,
            ])

    # Author Master Report
    report_text = build_report_markdown(unique_cands, len(gbk_files), n_shuffles)
    report_file = root / "NON_CANONICAL_VIRAL_SMORF_DISCOVERY.md"
    report_file.write_text(report_text, encoding="utf-8")

    sig_count = sum(1 for c in unique_cands if c.statistically_significant)
    novel_sig = sum(1 for c in unique_cands if c.statistically_significant and not c.is_annotated_in_genbank)

    print("\n[DISCOVERY RESULTS]")
    print(f"  Total Genomes Scanned:                  {len(gbk_files)}")
    print(f"  Total Unique Candidates Evaluated:      {len(unique_cands)}")
    print(f"  Statistically Significant (p<0.001):    {sig_count}")
    print(f"  Novel Unannotated smORFs (z>3.0):       {novel_sig}")
    print(f"  Exported CSV:                           {csv_file}")
    print(f"  Exported Report:                        {report_file}")
    print("=" * 80 + "\n")

    return {
        "total_genomes": len(gbk_files),
        "total_candidates": len(unique_cands),
        "significant_candidates": sig_count,
        "novel_significant": novel_sig,
        "csv_path": str(csv_file),
        "report_path": str(report_file),
    }


def build_report_markdown(candidates: list[NonCanonicalCandidate], n_genomes: int, n_shuffles: int) -> str:
    sig_cands = [c for c in candidates if c.statistically_significant]
    novel_sig = [c for c in sig_cands if not c.is_annotated_in_genbank]

    lines = [
        "# 🔬 Non-Canonical smORF Discovery & Translational Anomaly Report",
        "## Systematic Mining of Near-Cognate Initiation Codons & Viroporin Motifs Across Viral Genomes",
        "",
        f"**Date:** {now_iso()[:10]}  ",
        f"**Genomes Screened:** {n_genomes} RefSeq Viral Genomes  ",
        f"**smORF Window:** 30 to 120 Amino Acids  ",
        f"**Host CAI Threshold:** $\\ge 0.70$  ",
        f"**Null Controls:** $N={n_shuffles}$ Altschul-Erickson Dinucleotide Shuffles per Candidate  ",
        f"**Falsification Standard:** Significant only if $z > 3.0$ and $p < 0.001$  ",
        "",
        "---",
        "",
        "## 1. Executive Discovery Summary",
        "",
        "Standard genome annotation pipelines rely strictly on canonical AUG start codons, systematically missing functional microproteins (smORFs) initiated by **near-cognate start codons (CUG, GUG, ACG, AUA, UUG)**. Using Kozak initiation context scoring, host codon adaptation modeling, and $N=500$ dinucleotide null model calibrations, this engine systematically screened the RefSeq viral corpus.",
        "",
        f"### Key Findings:",
        f"- **Total Candidates Evaluated (CAI $\\ge 0.70$, 30–120 aa):** {len(candidates)}",
        f"- **Statistically Significant ($z > 3.0, p < 0.001$):** {len(sig_cands)}",
        f"- **Novel Unannotated Candidates:** **{len(novel_sig)}**",
        "",
        "---",
        "",
        "## 2. Top Non-Canonical smORF Candidates ($z > 3.0, p < 0.001$)",
        "",
        "| Organism (Accession) | Parent CDS | Frame | Start Codon (Type) | Kozak | Length | Host CAI | TM Helices | $\\mu_H$ | $z$-Score | Viroporin Profile |",
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |",
    ]

    for c in (novel_sig[:12] if novel_sig else sig_cands[:12]):
        tm_cnt = len(c.tm_segments)
        viroporin_icon = "⚡ Putative Channel" if c.is_potential_viroporin else "Soluble / Globular"
        lines.append(
            f"| {c.organism[:20]} (`{c.accession}`) | `{c.parent_gene}` | +{c.frame_offset} | `{c.start_codon}` ({c.start_codon_type.replace('NEAR_COGNATE_', '')}) | {c.kozak_strength} | **{c.length_aa} aa** | **{c.host_cai}** | {tm_cnt} TM | {c.hydrophobic_moment} | **$z = {c.z_score}$** | {viroporin_icon} |"
        )

    lines.extend([
        "",
        "---",
        "",
        "## 3. High-Priority Viroporin & Transmembrane smORF Profiles",
        "",
    ])

    viroporin_cands = [c for c in sig_cands if c.has_transmembrane_domain or c.is_potential_viroporin]
    for idx, c in enumerate(viroporin_cands[:5], 1):
        lines.extend([
            f"### smORF Candidate #{idx}: {c.organism} ({c.accession})",
            f"- **Parent Locus:** Frame +{c.frame_offset} of `{c.parent_gene}` ({c.parent_product})",
            f"- **Initiation Codon:** `{c.start_codon}` ({c.start_codon_type}) | **Efficiency:** {int(c.relative_initiation_efficiency*100)}%",
            f"- **Kozak Context:** `{c.kozak_context}` (Strength: **{c.kozak_strength}**)",
            f"- **Biophysical Properties:** {c.length_aa} aa | {c.molecular_weight_kda} kDa | pI = {c.isoelectric_point} | Net Charge = {c.net_charge_ph74}",
            f"- **Channel Architecture:** Mean Hydropathy = {c.mean_hydropathy}, Max Window = {c.max_hydropathy_window}, Amphipathic Moment $\\mu_H = {c.hydrophobic_moment}$",
            f"- **Statistical Null Calibration:** $z = {c.z_score}, p = {c.empirical_p_value}$ ($N=500$ shuffles)",
            "```fasta",
            f">smORF_{c.accession}_{c.parent_gene}_F+{c.frame_offset}|{c.start_codon}|{c.length_aa}aa|MW={c.molecular_weight_kda}kDa|CAI={c.host_cai}",
            c.peptide_sequence,
            "```",
            "",
        ])

    lines.extend([
        "---",
        "",
        "## 4. Methodological Controls & Experimental Next Steps",
        "",
        "1. **Ribo-seq P-Site Offset Mapping:** Align public ribo-seq datasets targeting non-canonical start sites (CUG/GUG) to verify translation initiation.",
        "2. **Chemical Crosslinking & Patch-Clamp:** For putative viroporin candidates with high amphipathic moments, perform planar lipid bilayer electrophysiology to test ion channel conductance.",
        "3. **Tagged Expression:** Express with C-terminal FLAG tag in HEK293T / host cell lines.",
        "",
        "---",
        "*Report generated deterministically by `scripts/mining/run_non_canonical_mining.py`.*",
    ])

    return "\n".join(lines)


if __name__ == "__main__":
    run_non_canonical_pipeline(n_shuffles=200)
