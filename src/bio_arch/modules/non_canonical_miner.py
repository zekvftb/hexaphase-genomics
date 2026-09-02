"""Non-Canonical & Near-Cognate smORF Discovery Engine for Viral Genomes.

Systematically scans alternative reading frames (+1, +2) across viral CDS features
for unannotated microproteins (30-120 aa) initiated by near-cognate start codons (CUG, GUG, ACG, AUA, UUG),
evaluates Kozak initiation contexts, computes amphipathic viroporin/TM channel motifs,
and falsifies candidates against N=500 Altschul-Erickson dinucleotide-preserving null models.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import io
import json
import math
from pathlib import Path
import random
import sys
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq

from bio_arch.modules.information import shuffle_dinucleotide
from bio_arch.modules.recompiler import compute_cai, translate_sequence
from bio_arch.provenance import SeedManager, now_iso

# Near-cognate start codons and empirical relative initiation efficiency in eukaryotic/mammalian systems
START_CODON_EFFICIENCY = {
    "ATG": 1.00,  # Canonical AUG
    "CTG": 0.25,  # High-efficiency CUG (~15-30%)
    "GTG": 0.12,  # Moderate-efficiency GUG (~10-15%)
    "ACG": 0.08,  # Moderate-efficiency ACG (~5-10%)
    "ATA": 0.05,  # Low-moderate AUA (~5%)
    "TTG": 0.05,  # Low-moderate UUG (~5%)
}

STOP_CODONS = {"TAA", "TAG", "TGA"}

# Kyte-Doolittle Hydropathy Scale
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

AA_WEIGHTS = {
    "A": 71.08, "R": 156.20, "N": 114.11, "D": 115.09, "C": 103.14,
    "E": 129.12, "Q": 128.13, "G": 57.05, "H": 137.14, "I": 113.17,
    "L": 113.17, "K": 128.18, "M": 131.21, "F": 147.18, "P": 97.12,
    "S": 87.08, "T": 101.11, "W": 186.21, "Y": 163.18, "V": 99.13,
}

PKA_VALUES = {
    "N_term": 9.69, "C_term": 2.34,
    "C": 8.33, "D": 3.86, "E": 4.25,
    "H": 6.00, "K": 10.50, "R": 12.48, "Y": 10.07,
}


@dataclass
class NonCanonicalCandidate:
    """An unannotated smORF initiated by canonical or near-cognate start codons."""

    accession: str
    organism: str
    parent_gene: str
    parent_product: str
    parent_cds_length_bp: int
    frame_offset: int
    start_in_parent_bp: int
    end_in_parent_bp: int
    start_codon: str
    start_codon_type: str
    relative_initiation_efficiency: float
    kozak_strength: str
    kozak_context: str
    length_aa: int
    molecular_weight_kda: float
    isoelectric_point: float
    net_charge_ph74: float
    peptide_sequence: str
    host_cai: float
    gc_content: float
    shannon_entropy: float
    mean_hydropathy: float
    max_hydropathy_window: float
    has_transmembrane_domain: bool
    tm_segments: list[dict] = field(default_factory=list)
    hydrophobic_moment: float = 0.0
    is_potential_viroporin: bool = False
    null_mean_score: float = 0.0
    null_std_score: float = 0.0
    z_score: float = 0.0
    empirical_p_value: float = 1.0
    statistically_significant: bool = False
    is_annotated_in_genbank: bool = False
    matched_annotation_name: str = ""


class NonCanonicalSmorfMiner:
    """Mines viral genomes for near-cognate initiated smORFs and translational anomalies."""

    def __init__(
        self,
        min_length_aa: int = 30,
        max_length_aa: int = 120,
        min_cai_threshold: float = 0.70,
        n_null_shuffles: int = 500,
        significance_alpha: float = 0.001,
        min_z_threshold: float = 3.0,
        seed: int = 42,
    ) -> None:
        self.min_length_aa = min_length_aa
        self.max_length_aa = max_length_aa
        self.min_cai_threshold = min_cai_threshold
        self.n_null_shuffles = n_null_shuffles
        self.significance_alpha = significance_alpha
        self.min_z_threshold = min_z_threshold
        self.seed_mgr = SeedManager(seed)

    def scan_cds(
        self,
        cds_dna: str,
        parent_gene: str = "CDS",
        parent_product: str = "viral protein",
        accession: str = "UNKNOWN",
        organism: str = "Viral isolate",
        annotated_translations: list[dict] | None = None,
    ) -> list[NonCanonicalCandidate]:
        """Scan Frame +1 and Frame +2 of primary CDS for near-cognate smORFs."""
        cds_dna = cds_dna.upper()
        parent_len_bp = len(cds_dna)
        annotated_translations = annotated_translations or []
        candidates: list[NonCanonicalCandidate] = []

        for frame_offset in (1, 2):
            subseq = cds_dna[frame_offset:]
            num_codons = len(subseq) // 3

            i = 0
            while i < num_codons:
                codon = subseq[i * 3 : i * 3 + 3]
                if codon in START_CODON_EFFICIENCY:
                    # Scan for in-frame stop codon
                    j = i + 1
                    found_stop = False
                    while j < num_codons:
                        stop_cand = subseq[j * 3 : j * 3 + 3]
                        if stop_cand in STOP_CODONS:
                            found_stop = True
                            break
                        j += 1

                    orf_len = j - i
                    if self.min_length_aa <= orf_len <= self.max_length_aa:
                        orf_dna = subseq[i * 3 : j * 3]
                        peptide = str(Seq(orf_dna).translate())

                        if "*" not in peptide:
                            cai = round(compute_cai(orf_dna), 4)

                            # Only consider candidates meeting the strict host CAI threshold
                            if cai >= self.min_cai_threshold:
                                # Kozak analysis
                                abs_start = frame_offset + i * 3
                                kozak_res = self._evaluate_kozak(cds_dna, abs_start)

                                # Biophysical metrics
                                mw = self._calc_mw(peptide)
                                pi = self._calc_pi(peptide)
                                charge = self._calc_charge(peptide)
                                gc = round(sum(1 for b in orf_dna if b in "GC") / len(orf_dna), 4)
                                entropy = round(self._calc_entropy(orf_dna), 4)
                                mean_hydro, max_hydro, tm_segs = self._analyze_hydropathy(peptide)
                                h_moment = self._calc_hydrophobic_moment(peptide)
                                is_viroporin = (len(tm_segs) >= 1) and (h_moment >= 0.25) and (charge > 0)

                                # Null model calibration against N=500 Eulerian dinucleotide shuffles
                                null_mean, null_std, z_val, p_val = self._calibrate_against_null(
                                    parent_cds_dna=cds_dna,
                                    frame_offset=frame_offset,
                                    actual_len=orf_len,
                                    actual_cai=cai,
                                    context_name=f"{accession}_{parent_gene}_nc_f{frame_offset}_i{i}",
                                )

                                is_sig = (p_val <= self.significance_alpha) and (z_val >= self.min_z_threshold)

                                # Check existing annotation
                                is_annotated, matched_name = self._check_existing(peptide, annotated_translations)

                                codon_type = "CANONICAL_AUG" if codon == "ATG" else f"NEAR_COGNATE_{codon.replace('T', 'U')}"
                                eff = START_CODON_EFFICIENCY.get(codon, 0.05)

                                cand = NonCanonicalCandidate(
                                    accession=accession,
                                    organism=organism,
                                    parent_gene=parent_gene,
                                    parent_product=parent_product,
                                    parent_cds_length_bp=parent_len_bp,
                                    frame_offset=frame_offset,
                                    start_in_parent_bp=abs_start,
                                    end_in_parent_bp=frame_offset + j * 3 + (3 if found_stop else 0),
                                    start_codon=codon,
                                    start_codon_type=codon_type,
                                    relative_initiation_efficiency=eff,
                                    kozak_strength=kozak_res["strength"],
                                    kozak_context=kozak_res["context"],
                                    length_aa=orf_len,
                                    molecular_weight_kda=mw,
                                    isoelectric_point=pi,
                                    net_charge_ph74=charge,
                                    peptide_sequence=peptide,
                                    host_cai=cai,
                                    gc_content=gc,
                                    shannon_entropy=entropy,
                                    mean_hydropathy=mean_hydro,
                                    max_hydropathy_window=max_hydro,
                                    has_transmembrane_domain=len(tm_segs) > 0,
                                    tm_segments=tm_segs,
                                    hydrophobic_moment=h_moment,
                                    is_potential_viroporin=is_viroporin,
                                    null_mean_score=null_mean,
                                    null_std_score=null_std,
                                    z_score=z_val,
                                    empirical_p_value=p_val,
                                    statistically_significant=is_sig,
                                    is_annotated_in_genbank=is_annotated,
                                    matched_annotation_name=matched_name,
                                )
                                candidates.append(cand)

                    i += 1
                else:
                    i += 1

        return candidates

    def _evaluate_kozak(self, cds_dna: str, start_bp: int) -> dict[str, str]:
        """Evaluate Kozak initiation context."""
        upstream = cds_dna[max(0, start_bp - 6) : start_bp]
        codon = cds_dna[start_bp : start_bp + 3]
        downstream = cds_dna[start_bp + 3 : min(len(cds_dna), start_bp + 7)]

        pos_minus_3 = cds_dna[start_bp - 3] if start_bp >= 3 else "N"
        pos_plus_4 = cds_dna[start_bp + 3] if (start_bp + 3) < len(cds_dna) else "N"

        score = 0
        if pos_minus_3 in ("A", "G"):
            score += 3
        if pos_plus_4 == "G":
            score += 2

        if score == 5:
            strength = "Optimal"
        elif score in (2, 3):
            strength = "Strong"
        else:
            strength = "Weak"

        return {
            "strength": strength,
            "context": f"...{upstream}[{codon}]{downstream}...",
        }

    def _calibrate_against_null(
        self,
        parent_cds_dna: str,
        frame_offset: int,
        actual_len: int,
        actual_cai: float,
        context_name: str,
    ) -> tuple[float, float, float, float]:
        """Null model calibration via N=500 Eulerian dinucleotide shuffles."""
        actual_score = actual_len * actual_cai
        null_scores = []
        n_shuffles = self.n_null_shuffles

        for k in range(n_shuffles):
            sub_seed = self.seed_mgr.derive_seed(f"{context_name}_nc_null_{k}")
            rng = random.Random(sub_seed)
            shuffled_parent = shuffle_dinucleotide(parent_cds_dna, rng=rng)

            shuf_sub = shuffled_parent[frame_offset:]
            num_codons = len(shuf_sub) // 3
            max_shuf_score = 0.0

            idx = 0
            while idx < num_codons:
                if shuf_sub[idx * 3 : idx * 3 + 3] in START_CODON_EFFICIENCY:
                    jdx = idx + 1
                    while jdx < num_codons:
                        if shuf_sub[jdx * 3 : jdx * 3 + 3] in STOP_CODONS:
                            break
                        jdx += 1
                    orf_len = jdx - idx
                    if orf_len >= 15:
                        shuf_orf_dna = shuf_sub[idx * 3 : jdx * 3]
                        shuf_cai = compute_cai(shuf_orf_dna)
                        shuf_score = orf_len * shuf_cai
                        if shuf_score > max_shuf_score:
                            max_shuf_score = shuf_score
                    idx = jdx + 1
                else:
                    idx += 1

            null_scores.append(max_shuf_score)

        null_mean = sum(null_scores) / len(null_scores)
        null_var = sum((x - null_mean) ** 2 for x in null_scores) / max(1, len(null_scores) - 1)
        null_std = max(1e-4, null_var ** 0.5)

        z_val = round((actual_score - null_mean) / null_std, 2)
        empirical_p = sum(1 for x in null_scores if x >= actual_score) / n_shuffles

        return round(null_mean, 2), round(null_std, 2), z_val, round(empirical_p, 4)

    def _calc_mw(self, peptide: str) -> float:
        total = sum(AA_WEIGHTS.get(aa, 110.0) for aa in peptide) + 18.015
        return round(total / 1000.0, 2)

    def _calc_pi(self, peptide: str) -> float:
        def net_charge(pH: float) -> float:
            c = 10.0 ** (PKA_VALUES["N_term"] - pH) / (1.0 + 10.0 ** (PKA_VALUES["N_term"] - pH))
            for aa in peptide:
                if aa in ("K", "R", "H"):
                    c += 10.0 ** (PKA_VALUES[aa] - pH) / (1.0 + 10.0 ** (PKA_VALUES[aa] - pH))
            c -= 10.0 ** (pH - PKA_VALUES["C_term"]) / (1.0 + 10.0 ** (pH - PKA_VALUES["C_term"]))
            for aa in peptide:
                if aa in ("D", "E", "C", "Y"):
                    c -= 10.0 ** (pH - PKA_VALUES[aa]) / (1.0 + 10.0 ** (pH - PKA_VALUES[aa]))
            return c

        low, high = 0.0, 14.0
        for _ in range(50):
            mid = (low + high) / 2.0
            if net_charge(mid) > 0:
                low = mid
            else:
                high = mid
        return round((low + high) / 2.0, 2)

    def _calc_charge(self, peptide: str, pH: float = 7.4) -> float:
        pos = sum(1 for aa in peptide if aa in ("K", "R"))
        pos += 0.1 * sum(1 for aa in peptide if aa == "H")
        neg = sum(1 for aa in peptide if aa in ("D", "E"))
        return round(pos - neg, 2)

    def _calc_entropy(self, dna: str) -> float:
        counts = {}
        for b in dna:
            counts[b] = counts.get(b, 0) + 1
        n = len(dna)
        if n == 0:
            return 0.0
        return -sum((c / n) * math.log2(c / n) for c in counts.values())

    def _analyze_hydropathy(self, peptide: str, window: int = 19, threshold: float = 1.6) -> tuple[float, float, list[dict]]:
        n = len(peptide)
        if n == 0:
            return 0.0, 0.0, []

        scores = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in peptide]
        mean_hydro = round(sum(scores) / n, 2)

        if n < window:
            return mean_hydro, mean_hydro, []

        win_scores = []
        tm_segs = []
        for i in range(n - window + 1):
            w_score = sum(scores[i : i + window]) / window
            win_scores.append((i, w_score))

        max_window = round(max(s for _, s in win_scores), 2)

        # Detect TM spanning helices
        i = 0
        while i < len(win_scores):
            start_i, score = win_scores[i]
            if score >= threshold:
                j = i
                top_s = score
                while j < len(win_scores) and win_scores[j][1] >= threshold:
                    top_s = max(top_s, win_scores[j][1])
                    j += 1
                tm_segs.append({
                    "start_aa": start_i + 1,
                    "end_aa": j + window - 1,
                    "length_aa": (j + window - 1) - start_i,
                    "mean_hydropathy": round(top_s, 2),
                })
                i = j + 1
            else:
                i += 1

        return mean_hydro, max_window, tm_segs

    def _calc_hydrophobic_moment(self, peptide: str, delta_deg: float = 100.0) -> float:
        """Calculate mean amphipathic hydrophobic moment for alpha-helix (delta=100 degrees)."""
        n = len(peptide)
        if n == 0:
            return 0.0

        delta_rad = math.radians(delta_deg)
        sum_cos = sum(KYTE_DOOLITTLE.get(aa, 0.0) * math.cos(i * delta_rad) for i, aa in enumerate(peptide))
        sum_sin = sum(KYTE_DOOLITTLE.get(aa, 0.0) * math.sin(i * delta_rad) for i, aa in enumerate(peptide))

        moment = math.sqrt(sum_cos ** 2 + sum_sin ** 2) / n
        return round(moment, 3)

    def _check_existing(self, peptide: str, annotated_translations: list[dict]) -> tuple[bool, str]:
        for item in annotated_translations:
            ann_trans = item.get("translation", "")
            gene = item.get("gene_name", "Known Gene")
            if not ann_trans:
                continue
            if peptide in ann_trans or ann_trans in peptide:
                return True, gene
            if len(peptide) >= 25 and peptide[:25] in ann_trans:
                return True, gene
        return False, ""

    def mine_genome(self, gbk_path: Path) -> list[NonCanonicalCandidate]:
        record = SeqIO.read(str(gbk_path), "genbank")
        acc = record.id
        organism = record.annotations.get("organism", record.description)

        all_annotated = []
        for f in record.features:
            if f.type == "CDS":
                gene = f.qualifiers.get("gene", [""])[0] or f.qualifiers.get("product", [""])[0]
                trans = f.qualifiers.get("translation", [""])[0]
                all_annotated.append({"gene_name": gene, "translation": trans})

        candidates: list[NonCanonicalCandidate] = []
        for feat in record.features:
            if feat.type == "CDS":
                cds_dna = str(feat.extract(record.seq))
                if len(cds_dna) >= 90 and len(cds_dna) % 3 == 0:
                    gene = feat.qualifiers.get("gene", [""])[0] or "CDS"
                    prod = feat.qualifiers.get("product", [""])[0] or "Protein"

                    cands = self.scan_cds(
                        cds_dna=cds_dna,
                        parent_gene=gene,
                        parent_product=prod,
                        accession=acc,
                        organism=organism,
                        annotated_translations=all_annotated,
                    )
                    candidates.extend(cands)

        return candidates
