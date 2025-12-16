"""Shared type definitions for TUSCO selector."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

# Type alias for junction keys
JunctionKey = Tuple[str, int, int]  # (chromosome, start, end)


@dataclass
class Transcript:
    """Represents a transcript with its essential features.
    
    Attributes:
        transcript_id: Unique transcript identifier
        exons: List of exon coordinates as (start, end) tuples, sorted by start position
        strand: Strand orientation ('+' or '-')
        chrom: Chromosome name
    """
    transcript_id: str
    exons: List[Tuple[int, int]]
    strand: str
    chrom: str
    
    def __post_init__(self):
        """Ensure exons are sorted by start position."""
        self.exons = sorted(self.exons, key=lambda x: x[0])
    
    @property
    def num_exons(self) -> int:
        """Return the number of exons."""
        return len(self.exons)
    
    @property
    def is_single_exon(self) -> bool:
        """Check if transcript has only one exon."""
        return self.num_exons == 1
    
    @property
    def junctions(self) -> List[JunctionKey]:
        """Get all splice junctions as JunctionKey tuples."""
        junctions = []
        for i in range(len(self.exons) - 1):
            # Junction is from end of exon i to start of exon i+1
            junction_start = self.exons[i][1]
            junction_end = self.exons[i + 1][0]
            junctions.append((self.chrom, junction_start, junction_end))
        return junctions
    
    @property
    def start(self) -> int:
        """Get transcript start position (5' end considering strand)."""
        if self.strand == '+':
            return self.exons[0][0]
        else:
            return self.exons[-1][1]
    
    @property
    def end(self) -> int:
        """Get transcript end position (3' end considering strand)."""
        if self.strand == '+':
            return self.exons[-1][1]
        else:
            return self.exons[0][0]
    
    @property
    def tss(self) -> int:
        """Get transcription start site position."""
        return self.start
    
    @property
    def length(self) -> int:
        """Calculate total exonic length."""
        return sum(end - start for start, end in self.exons)


__all__ = [
    'Transcript',
    'JunctionKey',
]