"""
Phylogenetic Validation for Phase 3 Evidence Synthesis.

This module provides final validation of EVE predictions using:
1. GVClass: Phylogenetic classification against NCLDV/MIRUS reference genomes
2. Diamond BLASTp: Protein-level classification against multi-domain databases

These validations are run after Phase 2 boundary refinement to validate the
final predicted regions, not the rough seeds from Phase 1.

The key insight is that GVClass and Diamond provide independent phylogenetic
evidence that can confirm or reject candidate regions:
- High GVClass/Diamond viral scores increase confidence
- Non-viral classifications can override otherwise viral-looking candidates
- Mixed signals indicate chimeric or ambiguous regions
"""

import logging
import subprocess
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq

from virosync.utils.path_safety import require_strict_child, safe_filename_component

logger = logging.getLogger(__name__)


class EvidenceToolError(RuntimeError):
    """An evidence tool was run and failed.

    compute_verdict() scores an absent GVClass or DIAMOND result as the neutral
    0.5 that means "neither supports nor rejects". That is correct for a genuinely
    uncertain domain and wrong for a crash: it lets a broken tool read as an
    inconclusive one, and the caller cannot tell the difference. Failures raise
    so they cannot be scored.
    """


# Domain classifications
VIRAL_DOMAINS = {"NCLDV", "MIRUS", "NCLDV/MIRUS"}
NON_VIRAL_DOMAINS = {"EUK", "BAC", "ARC"}  # PHAGE can be considered viral-like
UNCERTAIN_DOMAINS = {"UNKNOWN", "PHAGE"}  # UNKNOWN = novel/uncertain, not non-viral


class PhylogeneticVerdict(Enum):
    """Verdict from phylogenetic validation."""
    CONFIRMED_VIRAL = "confirmed_viral"      # Strong viral signal
    LIKELY_VIRAL = "likely_viral"            # Moderate viral signal
    AMBIGUOUS = "ambiguous"                  # Mixed or weak signal
    LIKELY_NON_VIRAL = "likely_non_viral"    # Moderate non-viral signal
    CONFIRMED_NON_VIRAL = "confirmed_non_viral"  # Strong non-viral signal


@dataclass
class GVClassValidation:
    """
    GVClass validation result for a single EVE region.

    GVClass provides phylogenetic classification based on marker genes
    and k-nearest neighbors against NCLDV/MIRUS reference genomes.
    """
    eve_id: str
    domain: str  # Best/primary domain: NCLDV, MIRUS, EUK, BAC, etc.
    domain_percent: float  # Percentage of nearest neighbors in this domain
    mcp_count: int  # Major Capsid Protein markers
    mirus_count: int  # Mirus-specific markers
    gvog_count: int  # Giant Virus Orthologous Groups
    total_markers: int  # Total markers detected
    # All domains from nearest neighbors: {domain: (count, percent)}
    all_domains: dict[str, tuple[int, float]] = field(default_factory=dict)

    @property
    def is_viral(self) -> Optional[bool]:
        """
        Whether classified as viral.

        Returns:
            True: Domain is viral (NCLDV, MIRUS)
            False: Domain is definitively non-viral (BAC, ARC, EUK)
            None: Domain is uncertain/unknown (UNKNOWN, PHAGE) - don't reject!
        """
        if self.domain in VIRAL_DOMAINS:
            return True
        elif self.domain in NON_VIRAL_DOMAINS:
            return False
        else:
            # UNKNOWN, PHAGE, or anything else = uncertain
            return None

    @property
    def is_non_viral(self) -> bool:
        """Whether definitively classified as non-viral (BAC, ARC, EUK)."""
        return self.domain in NON_VIRAL_DOMAINS

    @property
    def is_uncertain(self) -> bool:
        """Whether classification is uncertain (UNKNOWN, PHAGE)."""
        return self.domain in UNCERTAIN_DOMAINS or self.domain not in VIRAL_DOMAINS.union(NON_VIRAL_DOMAINS)

    @property
    def has_viral_in_neighbors(self) -> bool:
        """
        Whether ANY nearest neighbor is viral (NCLDV or MIRUS).

        This is used to avoid rejecting regions where the best domain is cellular
        but viral domains are present among the top neighbors. Such cases may
        represent genuine EVEs that happen to have some cellular gene acquisitions.
        """
        if not self.all_domains:
            # Fall back to checking primary domain
            return self.domain in VIRAL_DOMAINS
        return any(dom in VIRAL_DOMAINS for dom in self.all_domains.keys())

    @property
    def confidence_score(self) -> float:
        """
        Confidence score (0-1) for viral classification.

        Based on domain percentage and marker count.

        Returns:
            For viral domains: score based on percentage and markers
            For uncertain domains: neutral 0.5 (don't reject or accept)
            For non-viral domains: low score but not 0 (allow other evidence)
        """
        # VIRAL: high confidence based on evidence
        if self.is_viral is True:
            # Base score from domain percentage
            base_score = self.domain_percent / 100.0
            # Marker bonus (up to 0.2)
            marker_bonus = min(0.2, self.total_markers * 0.02)
            # MCP bonus (key hallmark)
            mcp_bonus = 0.1 if self.mcp_count > 0 else 0.0
            return min(1.0, base_score + marker_bonus + mcp_bonus)

        # UNCERTAIN (UNKNOWN, PHAGE): neutral - let other evidence decide
        elif self.is_uncertain:
            # Still give some credit for having markers
            if self.total_markers > 0:
                return 0.5 + min(0.2, self.total_markers * 0.02)
            return 0.5  # Neutral - neither supports nor rejects

        # NON-VIRAL (BAC, ARC, EUK): low but not zero
        else:
            # Very low score but leave room for override by other evidence
            return 0.1


@dataclass
class DiamondValidation:
    """
    Diamond BLASTp validation result for a single EVE region.

    Diamond provides protein-level classification by searching ALL
    predicted proteins against reference proteomes.
    """
    eve_id: str
    total_proteins: int
    proteins_with_hits: int
    domain_counts: dict[str, int]  # Counts per domain (based on best hit per query)
    best_domain: str
    best_domain_percent: float
    viral_protein_count: int
    non_viral_protein_count: int
    # All top 5 hits per query: {query: [(domain, bitscore), ...]}
    all_hits_per_query: dict[str, list[tuple[str, float]]] = field(default_factory=dict)

    @property
    def is_viral(self) -> bool:
        """Whether majority of proteins classify as viral."""
        return self.best_domain in {"NCLDV", "PHAGE"} and self.best_domain_percent >= 50.0

    @property
    def confidence_score(self) -> float:
        """
        Confidence score (0-1) for viral classification.

        Based on proportion of viral vs non-viral hits.
        """
        if self.proteins_with_hits == 0:
            return 0.5  # No information = neutral

        viral_fraction = self.viral_protein_count / self.proteins_with_hits

        # Penalize if many non-viral hits
        if self.non_viral_protein_count > self.viral_protein_count:
            return max(0.0, 0.5 - (self.non_viral_protein_count / self.proteins_with_hits) * 0.5)

        return viral_fraction

    @property
    def contamination_score(self) -> float:
        """
        Score indicating potential contamination (0-1).

        High score means region contains mixed viral/host proteins,
        suggesting integration boundaries may be imprecise.
        """
        if self.proteins_with_hits == 0:
            return 0.0

        if sum(self.domain_counts.values()) == 0:
            return 0.0

        # Coarse banding on the dominant domain's share; not an entropy.
        if self.best_domain_percent >= 90:
            return 0.1  # Very clean
        elif self.best_domain_percent >= 70:
            return 0.3  # Mostly clean
        elif self.best_domain_percent >= 50:
            return 0.5  # Some mixing
        else:
            return 0.8  # High mixing

    @property
    def has_viral_in_neighbors(self) -> bool:
        """
        Whether ANY of the top 5 hits for any protein is viral (NCLDV or PHAGE).

        This is used to avoid rejecting regions where the best domain is cellular
        but viral domains are present among the top hits.
        """
        viral_domains = {"NCLDV", "PHAGE"}
        # Check all_hits_per_query if available
        if self.all_hits_per_query:
            for hits in self.all_hits_per_query.values():
                for domain, _ in hits:
                    if domain in viral_domains:
                        return True
            return False
        # Fall back to domain_counts
        return any(dom in viral_domains for dom in self.domain_counts.keys())


@dataclass
class PhylogeneticValidationResult:
    """
    Combined phylogenetic validation result.

    Integrates GVClass and Diamond results into a single verdict.
    """
    eve_id: str
    scaffold: str
    start: int
    end: int

    # Component results
    gvclass: Optional[GVClassValidation] = None
    diamond: Optional[DiamondValidation] = None

    # Combined verdict
    verdict: PhylogeneticVerdict = PhylogeneticVerdict.AMBIGUOUS
    combined_score: float = 0.5

    # Flags
    has_mcp: bool = False
    has_viral_markers: bool = False
    is_chimeric: bool = False  # Mixed viral/host signal

    def compute_verdict(self) -> None:
        """
        Compute combined verdict from GVClass and Diamond results.

        Logic:
        1. If both agree viral with high confidence -> CONFIRMED_VIRAL
        2. If one viral + other neutral/missing -> LIKELY_VIRAL
        3. If disagreement or weak signals -> AMBIGUOUS
        4. If one non-viral + other neutral -> LIKELY_NON_VIRAL
        5. If both agree non-viral -> CONFIRMED_NON_VIRAL
        """
        gv_score = self.gvclass.confidence_score if self.gvclass else 0.5
        dm_score = self.diamond.confidence_score if self.diamond else 0.5

        gv_viral = self.gvclass.is_viral if self.gvclass else None
        dm_viral = self.diamond.is_viral if self.diamond else None

        # Update flags
        if self.gvclass:
            self.has_mcp = self.gvclass.mcp_count > 0
            self.has_viral_markers = self.gvclass.total_markers > 0

        if self.diamond and self.diamond.contamination_score > 0.5:
            self.is_chimeric = True

        # Compute combined score (weighted average)
        if self.gvclass and self.diamond:
            # Both available: weighted combination
            self.combined_score = 0.6 * gv_score + 0.4 * dm_score
        elif self.gvclass:
            self.combined_score = gv_score
        elif self.diamond:
            self.combined_score = dm_score
        else:
            self.combined_score = 0.5

        # Determine verdict
        if gv_viral is True and dm_viral is True:
            if self.combined_score >= 0.8:
                self.verdict = PhylogeneticVerdict.CONFIRMED_VIRAL
            else:
                self.verdict = PhylogeneticVerdict.LIKELY_VIRAL
        elif gv_viral is True and dm_viral is None:
            self.verdict = PhylogeneticVerdict.LIKELY_VIRAL if gv_score >= 0.6 else PhylogeneticVerdict.AMBIGUOUS
        elif gv_viral is None and dm_viral is True:
            self.verdict = PhylogeneticVerdict.LIKELY_VIRAL if dm_score >= 0.6 else PhylogeneticVerdict.AMBIGUOUS
        elif gv_viral is False and dm_viral is False:
            if self.combined_score <= 0.2:
                self.verdict = PhylogeneticVerdict.CONFIRMED_NON_VIRAL
            else:
                self.verdict = PhylogeneticVerdict.LIKELY_NON_VIRAL
        elif gv_viral is False or dm_viral is False:
            # One says non-viral
            self.verdict = PhylogeneticVerdict.LIKELY_NON_VIRAL if self.combined_score < 0.4 else PhylogeneticVerdict.AMBIGUOUS
        else:
            self.verdict = PhylogeneticVerdict.AMBIGUOUS

    @property
    def supports_viral(self) -> bool:
        """Whether phylogenetic evidence supports viral classification."""
        return self.verdict in {
            PhylogeneticVerdict.CONFIRMED_VIRAL,
            PhylogeneticVerdict.LIKELY_VIRAL,
        }

    @property
    def has_viral_in_any_neighbor(self) -> bool:
        """
        Whether ANY neighbor (from GVClass or Diamond) is viral.

        This is used to avoid rejecting regions that have viral signal
        among their top 5 nearest neighbors, even if the best hit is cellular.
        Such regions may be genuine EVEs with some acquired cellular genes.
        """
        gv_has_viral = self.gvclass.has_viral_in_neighbors if self.gvclass else False
        dm_has_viral = self.diamond.has_viral_in_neighbors if self.diamond else False
        return gv_has_viral or dm_has_viral

    @property
    def rejects_viral(self) -> bool:
        """
        Whether phylogenetic evidence rejects viral classification.

        IMPORTANT: Does NOT reject if viral signal is found among any
        of the top 5 neighbors from GVClass or Diamond. This allows
        genuine EVEs with some cellular gene acquisitions to pass.
        """
        # If ANY neighbor is viral (NCLDV/MIRUS), don't reject
        if self.has_viral_in_any_neighbor:
            return False
        return self.verdict in {
            PhylogeneticVerdict.CONFIRMED_NON_VIRAL,
            PhylogeneticVerdict.LIKELY_NON_VIRAL,
        }


class PhylogeneticValidator:
    """
    Phylogenetic validation engine for Phase 3.

    Runs GVClass and Diamond on final EVE predictions to provide
    independent phylogenetic confirmation or rejection.
    """

    def __init__(
        self,
        genome_path: Path,
        work_dir: Path,
        gvclass_db: Optional[Path] = None,
        diamond_db: Optional[Path] = None,
        threads: int = 16,
        min_region_size: int = 20000,
    ):
        """
        Initialize phylogenetic validator.

        Args:
            genome_path: Path to input genome FASTA
            work_dir: Working directory for intermediate files
            gvclass_db: Path to GVClass database (optional)
            diamond_db: Path to Diamond database (optional)
            threads: Number of threads for external tools
            min_region_size: Minimum region size for GVClass (default 20kb)
        """
        self.genome_path = Path(genome_path)
        self.work_dir = Path(work_dir)
        self.gvclass_db = gvclass_db
        self.diamond_db = Path(diamond_db) if diamond_db else None
        self.threads = threads
        self.min_region_size = min_region_size

        # Load genome sequences
        self._genome_seqs = None

        # Create work directories
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.gvclass_dir = self.work_dir / "gvclass"
        self.diamond_dir = self.work_dir / "diamond"

    @property
    def genome_seqs(self) -> dict:
        """Lazy load genome sequences."""
        if self._genome_seqs is None:
            self._genome_seqs = {
                rec.id: rec for rec in SeqIO.parse(self.genome_path, "fasta")
            }
        return self._genome_seqs

    def validate_eve(
        self,
        eve_id: str,
        scaffold: str,
        start: int,
        end: int,
        run_gvclass: bool = True,
        run_diamond: bool = True,
    ) -> PhylogeneticValidationResult:
        """
        Validate a single EVE region.

        Args:
            eve_id: EVE identifier
            scaffold: Scaffold/chromosome name
            start: Start position (0-based)
            end: End position
            run_gvclass: Whether to run GVClass
            run_diamond: Whether to run Diamond

        Returns:
            PhylogeneticValidationResult with combined verdict
        """
        result = PhylogeneticValidationResult(
            eve_id=eve_id,
            scaffold=scaffold,
            start=start,
            end=end,
        )

        # Extract region sequence
        if scaffold not in self.genome_seqs:
            logger.warning(f"Scaffold {scaffold} not found in genome")
            return result

        seq_record = self.genome_seqs[scaffold]

        # Extend if too small for GVClass
        region_size = end - start
        actual_start, actual_end = start, end
        if region_size < self.min_region_size:
            extension = (self.min_region_size - region_size) // 2
            actual_start = max(0, start - extension)
            actual_end = min(len(seq_record), end + extension)

        region_seq = seq_record.seq[actual_start:actual_end]

        # Run GVClass
        if run_gvclass:
            try:
                result.gvclass = self._run_gvclass(
                    eve_id, str(region_seq), actual_start, actual_end
                )
            except EvidenceToolError:
                raise
            except Exception as e:
                raise EvidenceToolError(f"GVClass failed for {eve_id}: {e}") from e

        # Run Diamond
        if run_diamond and self.diamond_db:
            try:
                result.diamond = self._run_diamond(
                    eve_id, str(region_seq)
                )
            except EvidenceToolError:
                raise
            except Exception as e:
                raise EvidenceToolError(f"Diamond failed for {eve_id}: {e}") from e

        # Compute combined verdict
        result.compute_verdict()

        return result

    def _run_gvclass(
        self,
        eve_id: str,
        sequence: str,
        start: int,
        end: int,
    ) -> Optional[GVClassValidation]:
        """Run GVClass on a single region."""
        filename_component = safe_filename_component(eve_id)
        self.gvclass_dir.mkdir(parents=True, exist_ok=True)

        # Write sequence to temp file
        input_dir = self.gvclass_dir / filename_component
        require_strict_child(self.gvclass_dir, input_dir)
        input_dir.mkdir(exist_ok=True)
        input_fasta = input_dir / f"{filename_component}.fna"
        require_strict_child(input_dir, input_fasta)

        with open(input_fasta, "w") as f:
            f.write(f">{eve_id} start={start} end={end}\n{sequence}\n")

        output_dir = self.gvclass_dir / f"{filename_component}_output"
        require_strict_child(self.gvclass_dir, output_dir)

        # Build command
        cmd = [
            "gvclass",
            str(input_dir),
            "-o", str(output_dir),
            "-t", str(self.threads),
            "--mode-fast",
        ]

        if self.gvclass_db:
            cmd.extend(["-d", str(self.gvclass_db)])

        logger.debug(f"Running GVClass for {eve_id}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800,  # 30 min timeout
            )

            if result.returncode != 0:
                raise EvidenceToolError(
                    f"GVClass exited {result.returncode} for {eve_id}: "
                    f"{(result.stderr or '').strip()[-400:]}"
                )

        except subprocess.TimeoutExpired as exc:
            raise EvidenceToolError(f"GVClass timed out for {eve_id}") from exc

        # Parse results
        summary_files = list(output_dir.glob("*.summary.tab"))
        if not summary_files:
            summary_files = list(output_dir.glob("gvclass_summary.tsv"))

        if not summary_files:
            raise EvidenceToolError(
                f"GVClass produced no summary file for {eve_id}; it exited 0 but "
                "wrote nothing, so its verdict is unknown rather than uncertain"
            )

        return self._parse_gvclass_summary(eve_id, summary_files[0])

    def _parse_gvclass_summary(
        self,
        eve_id: str,
        summary_path: Path,
    ) -> Optional[GVClassValidation]:
        """Parse GVClass summary file, capturing ALL domains from nearest neighbors."""
        with open(summary_path) as f:
            f.readline()  # Skip header
            for line in f:
                if not line.strip():
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 17:
                    continue

                # Parse ALL domains from nearest neighbors breakdown
                # Format: "NCLDV:3(60.00%),BAC:2(40.00%)" or "BAC:2(100.00%)"
                domain_info = fields[8] if len(fields) > 8 else ""
                domain = "UNKNOWN"
                domain_percent = 0.0
                all_domains: dict[str, tuple[int, float]] = {}

                if domain_info:
                    parts = domain_info.split(",")
                    for i, part in enumerate(parts):
                        part = part.strip()
                        if ":" in part:
                            dom_name = part.split(":")[0]
                            # Parse count and percentage: "NCLDV:3(60.00%)"
                            count_pct = part.split(":")[1]
                            count = 1  # default
                            pct = 0.0
                            if "(" in count_pct:
                                count_str = count_pct.split("(")[0]
                                try:
                                    count = int(count_str)
                                except ValueError:
                                    count = 1
                                if "%" in count_pct:
                                    pct_str = count_pct.split("(")[1].split("%")[0]
                                    try:
                                        pct = float(pct_str)
                                    except ValueError:
                                        pct = 0.0
                            all_domains[dom_name] = (count, pct)
                            # First domain is the primary/best domain
                            if i == 0:
                                domain = dom_name
                                domain_percent = pct

                # Parse marker counts
                mcp_count = int(fields[16]) if len(fields) > 16 and fields[16] else 0
                mirus_count = int(fields[18]) if len(fields) > 18 and fields[18] else 0
                gvog_count = int(fields[14]) if len(fields) > 14 and fields[14] else 0

                return GVClassValidation(
                    eve_id=eve_id,
                    domain=domain,
                    domain_percent=domain_percent,
                    mcp_count=mcp_count,
                    mirus_count=mirus_count,
                    gvog_count=gvog_count,
                    total_markers=mcp_count + mirus_count + gvog_count,
                    all_domains=all_domains,
                )

        return None

    def _run_diamond(
        self,
        eve_id: str,
        sequence: str,
    ) -> Optional[DiamondValidation]:
        """Run Diamond BLASTp on proteins from region."""
        filename_component = safe_filename_component(eve_id)
        self.diamond_dir.mkdir(parents=True, exist_ok=True)

        # 6-frame translation
        proteins = []
        seq = Seq(sequence)

        for frame in [0, 1, 2]:
            # Forward
            trans = str(seq[frame:].translate(stop_symbol="*"))
            for j, prot in enumerate(trans.split("*")):
                if len(prot) >= 30:
                    proteins.append((f"{eve_id}_f{frame}_{j}", prot))
            # Reverse
            trans = str(seq.reverse_complement()[frame:].translate(stop_symbol="*"))
            for j, prot in enumerate(trans.split("*")):
                if len(prot) >= 30:
                    proteins.append((f"{eve_id}_r{frame}_{j}", prot))

        if not proteins:
            return DiamondValidation(
                eve_id=eve_id,
                total_proteins=0,
                proteins_with_hits=0,
                domain_counts={},
                best_domain="UNKNOWN",
                best_domain_percent=0.0,
                viral_protein_count=0,
                non_viral_protein_count=0,
            )

        # Write proteins
        proteins_fasta = self.diamond_dir / f"{filename_component}.faa"
        require_strict_child(self.diamond_dir, proteins_fasta)
        with open(proteins_fasta, "w") as f:
            for prot_id, prot_seq in proteins:
                f.write(f">{prot_id}\n{prot_seq}\n")

        # Run Diamond via the central hardened wrapper (per-call tempdir,
        # new session, env scrub, automatic retry on futex deadlock).
        from virosync.pipeline.search_backend import run_sequence_search

        output_file = self.diamond_dir / f"{filename_component}.diamond.tsv"
        require_strict_child(self.diamond_dir, output_file)
        logger.debug(f"Running Diamond for {eve_id}")
        try:
            run_sequence_search(
                query_fasta=proteins_fasta,
                db_path=self.diamond_db,
                output_tsv=output_file,
                threads=self.threads,
                backend="diamond",
                evalue=1e-5,
                max_target_seqs=5,
                output_columns=[
                    "qseqid", "sseqid", "pident", "length", "mismatch",
                    "gapopen", "qstart", "qend", "sstart", "send",
                    "evalue", "bitscore",
                ],
                extra_flags=["--sensitive"],
                timeout=600,
            )
        except (RuntimeError, subprocess.CalledProcessError) as exc:
            raise EvidenceToolError(f"Diamond failed for {eve_id}: {exc}") from exc

        # Parse results
        return self._parse_diamond_output(eve_id, output_file, len(proteins))

    def _parse_diamond_output(
        self,
        eve_id: str,
        output_file: Path,
        total_proteins: int,
    ) -> DiamondValidation:
        """Parse Diamond output and classify, tracking ALL top 5 hits per query."""
        domain_prefixes = {
            "NCLDV__": "NCLDV",
            "BAC__": "BAC",
            "ARC__": "ARC",
            "PHAGE__": "PHAGE",
            "CRESS__": "CRESS",
            "EUK__": "EUK",
        }

        # Track ALL hits per query (up to 5 from -k 5)
        all_hits_per_query: dict[str, list[tuple[str, float]]] = {}
        # Also track best hit per query for domain counting
        query_best_domain: dict[str, tuple[str, float]] = {}

        if output_file.exists():
            with open(output_file) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 12:
                        continue

                    query = parts[0]
                    subject = parts[1]
                    bitscore = float(parts[11])

                    # Infer domain from subject prefix
                    domain = "UNKNOWN"
                    for prefix, dom in domain_prefixes.items():
                        if subject.startswith(prefix):
                            domain = dom
                            break

                    # Store ALL hits per query (Diamond returns up to 5 with -k 5)
                    if query not in all_hits_per_query:
                        all_hits_per_query[query] = []
                    all_hits_per_query[query].append((domain, bitscore))

                    # Keep best hit per query for domain counting
                    if query not in query_best_domain or bitscore > query_best_domain[query][1]:
                        query_best_domain[query] = (domain, bitscore)

        # Count domains based on best hit per query
        domain_counts = Counter(dom for dom, _ in query_best_domain.values())
        proteins_with_hits = len(query_best_domain)

        # Determine best domain
        if domain_counts:
            best_domain, best_count = domain_counts.most_common(1)[0]
            best_domain_percent = (best_count / proteins_with_hits) * 100
        else:
            best_domain = "UNKNOWN"
            best_domain_percent = 0.0

        # Count viral vs non-viral
        viral_domains = {"NCLDV", "PHAGE"}
        non_viral_domains = {"BAC", "ARC", "EUK"}

        viral_count = sum(domain_counts.get(d, 0) for d in viral_domains)
        non_viral_count = sum(domain_counts.get(d, 0) for d in non_viral_domains)

        return DiamondValidation(
            eve_id=eve_id,
            total_proteins=total_proteins,
            proteins_with_hits=proteins_with_hits,
            domain_counts=dict(domain_counts),
            best_domain=best_domain,
            best_domain_percent=best_domain_percent,
            viral_protein_count=viral_count,
            non_viral_protein_count=non_viral_count,
            all_hits_per_query=all_hits_per_query,
        )
