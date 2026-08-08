"""
pORF coordinate model.

Defines the PORF dataclass used to encode/decode pORF coordinates in FASTA
headers. The legacy six-frame translation helpers were retired; gene calling
is handled by prodigal-gv in the active workflow.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class PORF:
    """Represents a potential Open Reading Frame."""

    porf_id: str
    scaffold: str
    start: int  # Nucleotide start (0-based)
    end: int  # Nucleotide end (exclusive)
    strand: str  # '+' or '-'
    frame: int  # 1, 2, or 3
    sequence: str  # Amino acid sequence

    @classmethod
    def parse_header(cls, header: str) -> Optional["PORF"]:
        """
        Parse a pORF FASTA header to extract coordinates.

        Header format: pORF_1|scaffold:start-end_strand:+_frame:1
        Note: Scaffold names may contain underscores (e.g., NC_057019.1)

        Returns:
            PORF object with parsed coordinates, or None if parsing fails
        """
        try:
            parts = header.split("|")
            porf_id = parts[0]
            # Join remaining parts to handle scaffold names containing '|' (e.g., bin-7|contig_1001)
            coord_str = "|".join(parts[1:])

            # Parse from the right to handle scaffold names with underscores
            # Expected format: scaffold:start-end_strand:X_frame:N

            # Extract frame (rightmost)
            coord_str, frame_part = coord_str.rsplit("_frame:", 1)
            frame = int(frame_part)

            # Extract strand
            coord_str, strand_part = coord_str.rsplit("_strand:", 1)
            strand = strand_part

            # What remains is scaffold:start-end
            scaffold, range_str = coord_str.rsplit(":", 1)
            start, end = map(int, range_str.split("-"))

            return cls(
                porf_id=porf_id,
                scaffold=scaffold,
                start=start,
                end=end,
                strand=strand,
                frame=frame,
                sequence="",  # Not parsed from header
            )
        except (IndexError, ValueError):
            return None
