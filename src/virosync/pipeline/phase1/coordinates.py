"""
Coordinate helpers for 6-frame translated queries.

Shared between HHG seeding, marker validation, and legacy feature mapping.
"""

from typing import Optional


def parse_frame_id(seq_id: str) -> tuple[str, Optional[int], int]:
    """
    Parse a seqkit 6-frame translated header into (contig, frame, offset).

    Supports headers with `_frame=` or `_frame_` and optional `_offset=`.
    Returns (contig, frame, offset) where frame is 1..6 or negative for reverse.
    """
    offset = 0
    base = seq_id
    if "_offset=" in seq_id:
        base, offset_part = seq_id.rsplit("_offset=", 1)
        try:
            offset = int(offset_part)
        except ValueError:
            offset = 0
    if "_frame=" in base:
        contig, frame_part = base.rsplit("_frame=", 1)
        try:
            return contig, int(frame_part), offset
        except ValueError:
            return contig, None, offset
    if "_frame_" in base:
        contig, frame_part = base.rsplit("_frame_", 1)
        try:
            return contig, int(frame_part), offset
        except ValueError:
            return contig, None, offset
    return base, None, offset


def aa_to_nt_coords(
    aa_start: int,
    aa_end: int,
    contig_len: int,
    frame: int,
    offset: int = 0,
) -> tuple[int, int, str]:
    """
    Convert amino-acid coordinates to nucleotide coordinates for 6-frame translations.
    """
    if frame == 0:
        raise ValueError("Frame cannot be 0")
    aa_start = aa_start + offset
    aa_end = aa_end + offset
    if frame < 0:
        strand = "-"
        frame_offset = abs(frame) - 1
        nt_end = contig_len - (aa_start - 1) * 3 - frame_offset
        nt_start = contig_len - aa_end * 3 - frame_offset
    elif frame in (1, 2, 3):
        strand = "+"
        frame_offset = frame - 1
        nt_start = (aa_start - 1) * 3 + frame_offset
        nt_end = aa_end * 3 + frame_offset
    else:
        strand = "-"
        frame_offset = frame - 4
        nt_end = contig_len - (aa_start - 1) * 3 - frame_offset
        nt_start = contig_len - aa_end * 3 - frame_offset

    nt_start = max(0, nt_start)
    nt_end = min(contig_len, nt_end)
    if nt_end < nt_start:
        nt_start, nt_end = nt_end, nt_start
    return nt_start, nt_end, strand
