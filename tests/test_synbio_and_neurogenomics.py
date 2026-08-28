"""Unit tests for synthetic biology enzymes and neurodevelopmental synaptic genomics."""

import pytest
from bio_arch.modules.biological_circuits import scan_all_biological_circuits
from bio_arch.modules.information import compute_composition, shannon_entropy
from bio_arch.modules.logic_gates import scan_all_logic_gates
from bio_arch.modules.quantum_biology import calculate_wkb_proton_tunneling
from bio_arch.modules.recompiler import recompile_dual_protein_dna, RecompilationResult


def test_petase_dual_phase_recompilation():
    """Verify exact dual-protein recompilation of a PETase construct."""
    f0 = "MCLV"
    f1 = "CLWV"
    result = recompile_dual_protein_dna(f0, f1)
    assert isinstance(result, RecompilationResult)
    assert result.f0_identity_pct == 100.0
    assert result.compression_ratio >= 1.70
    assert len(result.synthesized_dna) == 13


def test_neurogenomics_shank3_synaptic_circuits():
    """Verify circuit architecture and CpG memory latches in SHANK3 scaffold."""
    # Synthetic SHANK3 GC-rich postsynaptic scaffold fragment
    shank3_fragment = "ATGGAGGGCCCGGGCCTGGGGCTCGGCCCGGGCGGCGGGGGCGGCTCCGGGGGCGGCTCC" * 10
    circuits = scan_all_biological_circuits(shank3_fragment, genome_id="SHANK3_SYNAPSE")

    assert circuits.genome_length_bp == len(shank3_fragment)
    assert len(circuits.cpg_memory_islands) >= 1
    # Check that high-GC synaptic density functions as promoter latch
    assert circuits.cpg_memory_islands[0].gc_content_pct >= 60.0


def test_drd4_dopamine_repeat_entropy():
    """Verify information entropy calculation across DRD4 dopamine receptor sequence."""
    drd4_seq = "ATGGGGAACCGCAGCGCCGCGGCCGCCGGGGGCGCCGACGCGGCTGCCGCTGCCGCCGCC" * 5
    comp = compute_composition(drd4_seq)
    entropy = shannon_entropy(drd4_seq)

    assert comp.gc_content >= 0.70
    assert 0.0 < entropy <= 2.0
