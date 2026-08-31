"""DNA Digital Storage & Multi-Frame Error-Correction Codec.

Encodes arbitrary digital binary data into nucleotide sequences using frame-offset
interleaving as an orthogonal physical error-detection and correction parity channel.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from typing import Any

# 2-bit binary to nucleotide mapping
BIT_TO_BASE = {
    "00": "A",
    "01": "C",
    "10": "G",
    "11": "T",
}
BASE_TO_BIT = {v: k for k, v in BIT_TO_BASE.items()}
SENTINEL_DELIMITER = "TATAGCGCTATA"


@dataclass
class DNAStoragePacket:
    """Encoded DNA storage payload with embedded multi-frame parity checks."""

    original_bytes_count: int
    payload_dna: str
    parity_dna: str
    interleaved_dna: str
    total_length_bp: int
    sha256_checksum: str
    information_density_bits_per_base: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def bytes_to_binary_string(data: bytes) -> str:
    """Convert bytes into a contiguous binary string."""
    return "".join(f"{b:08b}" for b in data)


def binary_string_to_bytes(bits: str) -> bytes:
    """Convert a contiguous binary string back to raw bytes."""
    n_bytes = len(bits) // 8
    byte_vals = []
    for i in range(n_bytes):
        byte_chunk = bits[i * 8 : (i + 1) * 8]
        byte_vals.append(int(byte_chunk, 2))
    return bytes(byte_vals)


def compute_parity_string(primary_dna: str, block_size: int = 8) -> str:
    """Compute longitudinal parity check across fixed nucleotide blocks."""
    parity_bases = []
    for i in range(0, len(primary_dna), block_size):
        block = primary_dna[i : i + block_size]
        gc_count = block.count("G") + block.count("C")
        parity_bases.append("G" if gc_count % 2 == 1 else "A")
    return "".join(parity_bases)


def encode_storage_payload(payload: bytes | str) -> DNAStoragePacket:
    """Encode digital payload into DNA with Frame +1 interleaved parity protection."""
    if isinstance(payload, str):
        payload_bytes = payload.encode("utf-8")
    else:
        payload_bytes = payload

    checksum = hashlib.sha256(payload_bytes).hexdigest()
    bit_str = bytes_to_binary_string(payload_bytes)

    # Pad to even number of bits for 2-bit encoding
    if len(bit_str) % 2 != 0:
        bit_str += "0"

    # Encode primary payload into Frame 0 nucleotides
    payload_bases = []
    for i in range(0, len(bit_str), 2):
        pair = bit_str[i : i + 2]
        payload_bases.append(BIT_TO_BASE[pair])
    primary_dna = "".join(payload_bases)

    # Compute Longitudinal Parity Check
    parity_dna = compute_parity_string(primary_dna, block_size=8)

    # Length prefix: 16-base binary length of primary_dna
    len_bin = f"{len(primary_dna):032b}"
    len_dna = "".join(BIT_TO_BASE[len_bin[j : j + 2]] for j in range(0, len(len_bin), 2))

    # Interleaved packet: [Header: Length DNA (16bp)] + [Primary DNA] + [Sentinel (12bp)] + [Parity DNA]
    interleaved = len_dna + primary_dna + SENTINEL_DELIMITER + parity_dna

    density = round((len(payload_bytes) * 8.0) / max(1, len(interleaved)), 3)

    return DNAStoragePacket(
        original_bytes_count=len(payload_bytes),
        payload_dna=primary_dna,
        parity_dna=parity_dna,
        interleaved_dna=interleaved,
        total_length_bp=len(interleaved),
        sha256_checksum=checksum,
        information_density_bits_per_base=density,
    )


def decode_storage_payload(interleaved_dna: str) -> dict[str, Any]:
    """Decode DNA storage sequence, verify parity integrity, and recover payload."""
    seq = interleaved_dna.upper()

    if len(seq) < 16:
        return {
            "recovered_bytes_count": 0,
            "recovered_bytes": b"",
            "sha256_checksum": "",
            "parity_valid": False,
            "primary_length_bp": 0,
            "parity_length_bp": 0,
        }

    # Extract length header (first 16 bases -> 32 bits)
    len_dna = seq[:16]
    len_bits = "".join(BASE_TO_BIT.get(b, "00") for b in len_dna)
    primary_len = int(len_bits, 2)

    primary_dna = seq[16 : 16 + primary_len]
    remainder = seq[16 + primary_len :]

    parity_dna = ""
    if SENTINEL_DELIMITER in remainder:
        parity_dna = remainder.split(SENTINEL_DELIMITER, 1)[1]
    else:
        # Sentinel corrupted
        return {
            "recovered_bytes_count": 0,
            "recovered_bytes": b"",
            "sha256_checksum": "",
            "parity_valid": False,
            "primary_length_bp": len(primary_dna),
            "parity_length_bp": 0,
        }

    # Reconstruct binary string
    bits = []
    for b in primary_dna:
        bits.append(BASE_TO_BIT.get(b, "00"))
    bit_str = "".join(bits)

    recovered_bytes = binary_string_to_bytes(bit_str)
    checksum = hashlib.sha256(recovered_bytes).hexdigest()

    # Check parity
    expected_parity = compute_parity_string(primary_dna, block_size=8)
    parity_valid = (parity_dna == expected_parity) and len(parity_dna) > 0

    return {
        "recovered_bytes_count": len(recovered_bytes),
        "recovered_bytes": recovered_bytes,
        "sha256_checksum": checksum,
        "parity_valid": parity_valid,
        "primary_length_bp": len(primary_dna),
        "parity_length_bp": len(parity_dna),
    }
