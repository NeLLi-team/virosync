"""
Viral family marker definitions and classification.

This module defines marker gene profiles for different giant virus families,
enabling detection of diverse viral lineages beyond just NCLDVs.

Supported viral lineages:
- Nucleocytoviricota (NCLDVs): Mimivirus, Marseillevirus, Pandoravirus, etc.
- Mriyaviricetes: Yaravirus, Gamadviruses (small relatives of NCLDVs)
- Mirusviricota: Mirusviruses (discovered 2022, plankton-associated)
- Polintoviruses: DNA transposon-derived viruses

References:
- Schulz et al. 2020: Giant virus diversity and ecology
- Gaïa et al. 2023: Mirusviricota discovery
- Boratto et al. 2024: Mriyaviruses characterization
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ViralFamilyProfile:
    """Defines marker profile for a viral family/lineage."""

    name: str  # Full name
    short_name: str  # Code name

    # Markers highly specific to this lineage (rarely in host genomes)
    diagnostic_markers: set[str] = field(default_factory=set)

    # Markers shared with hosts or other lineages
    supporting_markers: set[str] = field(default_factory=set)

    # Minimum markers required for seed acceptance
    min_markers_default: int = 2  # For assembled genomes
    min_markers_fragmented: int = 1  # For fragmented assemblies/MAGs

    # Minimum bit score to accept single diagnostic marker
    single_marker_min_score: float = 100.0

    # Description
    description: str = ""

    @property
    def all_markers(self) -> set[str]:
        """All markers for this family."""
        return self.diagnostic_markers | self.supporting_markers



# =============================================================================
# VIRAL FAMILY DEFINITIONS
# =============================================================================

NCLDV_PROFILE = ViralFamilyProfile(
    name="Nucleocytoviricota (NCLDVs)",
    short_name="ncldv",
    diagnostic_markers={
        "mcp",      # Major Capsid Protein - most diagnostic
        "a32",      # ATPase A32
        "d5",       # D5 helicase-primase
        "vltf3",    # Virus late transcription factor 3
        "mrnac",    # mRNA capping enzyme
    },
    supporting_markers={
        "polb",     # DNA polymerase B - also in eukaryotes
        "rnapl",    # RNA polymerase large subunit - also in eukaryotes
        "rnaps",    # RNA polymerase small subunit - also in eukaryotes
        "rnr",      # Ribonucleotide reductase - also in eukaryotes
        "sfii",     # Superfamily II helicase - widespread
    },
    min_markers_default=2,
    min_markers_fragmented=1,
    single_marker_min_score=80.0,  # MCP at 80+ is quite diagnostic
    description="Giant viruses including Mimivirus, Marseillevirus, Pandoravirus, Pithovirus, Mollivirus"
)

MRIYAVIRUS_PROFILE = ViralFamilyProfile(
    name="Mriyaviricetes (Mriyaviruses)",
    short_name="mriyavirus",
    diagnostic_markers={
        "mcp",          # Major Capsid Protein (double jelly-roll)
        "vltf2",        # Virus late transcription factor 2 - very diagnostic
        "vltf3",        # Virus late transcription factor 3
        "atpase_pkg",   # DNA packaging ATPase
        "huh_endo",     # HUH endonuclease (rolling circle replication)
    },
    supporting_markers={
        "ruvc",         # RuvC Holliday junction resolvase
        "pddexk",       # PDDEXK endonuclease
        "sf3_hel",      # SF3 helicase
        "sf2_hel",      # SF2 helicase with primase domain
        "ssb",          # ssDNA binding protein
    },
    # Tightened 2026-04: min_markers=1 @ 60 bit was matching widespread
    # TE-borne HUH endonuclease / SF3 helicase homologs and calling them
    # Mriyavirus. Require ≥2 markers by default and in fragmented mode,
    # and raise the single-marker bit-score floor so any future config
    # override that re-enables min_markers=1 still rejects low-confidence
    # hits.
    min_markers_default=2,
    min_markers_fragmented=2,
    single_marker_min_score=100.0,
    description="Small relatives of NCLDVs including Yaravirus and Gamadviruses (35-45kb genomes)"
)

MIRUSVIRUS_PROFILE = ViralFamilyProfile(
    name="Mirusviricota (Mirusviruses)",
    short_name="mirusvirus",
    diagnostic_markers={
        "mcp_mirus",    # Mirusvirus MCP (distinct from NCLDV)
        "polb_mirus",   # Mirusvirus-specific PolB
        "hel_mirus",    # Mirusvirus helicase
    },
    supporting_markers={
        "polb",         # Generic DNA polymerase B
        "sfii",         # SF2 helicase
    },
    min_markers_default=1,  # Less characterized, accept single markers
    min_markers_fragmented=1,
    single_marker_min_score=80.0,
    description="Plankton-associated giant viruses discovered 2022, distinct from NCLDVs"
)

POLINTOVIRUS_PROFILE = ViralFamilyProfile(
    name="Polintoviruses",
    short_name="polintovirus",
    diagnostic_markers={
        "mcp_poli",     # Polintovirus MCP
        "ppolb",        # Protein-primed PolB
        "pro_c1",       # C1 cysteine protease
    },
    supporting_markers={
        "atpase",       # Packaging ATPase
        "int_tyr",      # Tyrosine integrase
    },
    min_markers_default=1,
    min_markers_fragmented=1,
    single_marker_min_score=70.0,
    description="DNA transposon-derived viruses, related to virophages"
)

VP_PLV_PROFILE = ViralFamilyProfile(
    name="Virophages and Polinton-like Viruses",
    short_name="vp_plv",
    diagnostic_markers={
        # VP MCP markers (>95% VP/PLV purity)
        "vp_mcp_1", "vp_mcp_2", "vp_mcp_3", "vp_mcp_4",
        "vp_mcp_5", "vp_mcp_6", "vp_mcp_7",
        # VP ATPase markers
        "vp_atpase_1", "vp_atpase_2", "vp_atpase_3", "vp_atpase_4",
        # VP Penton markers
        "vp_penton_1", "vp_penton_2", "vp_penton_3", "vp_penton_4",
        "vp_penton_6", "vp_penton_7",
        # VP Protease markers
        "vp_pro_1", "vp_pro_2",
        # PLV markers
        "plv_pc_054",
    },
    supporting_markers={
        "vp_penton_5",  # Lower purity penton marker
    },
    min_markers_default=1,  # Often fragmented or sparse
    min_markers_fragmented=1,
    single_marker_min_score=60.0,
    description="Virophages (VP) and Polinton-like Viruses (PLV) - small dsDNA viruses parasitizing giant viruses"
)

# Registry of all viral families
VIRAL_FAMILIES = {
    "ncldv": NCLDV_PROFILE,
    "mriyavirus": MRIYAVIRUS_PROFILE,
    "mirusvirus": MIRUSVIRUS_PROFILE,
    "polintovirus": POLINTOVIRUS_PROFILE,
    "vp_plv": VP_PLV_PROFILE,
}

# All known diagnostic markers across all families
ALL_DIAGNOSTIC_MARKERS = set()
for profile in VIRAL_FAMILIES.values():
    ALL_DIAGNOSTIC_MARKERS.update(m.lower() for m in profile.diagnostic_markers)

# All known supporting markers
ALL_SUPPORTING_MARKERS = set()
for profile in VIRAL_FAMILIES.values():
    ALL_SUPPORTING_MARKERS.update(m.lower() for m in profile.supporting_markers)




def get_family_for_markers(marker_set: set[str]) -> Optional[str]:
    """
    Determine the most likely viral family given a set of markers.

    Args:
        marker_set: Set of marker names found in a region

    Returns:
        Best matching family short_name, or None if ambiguous
    """
    markers_lower = {m.lower() for m in marker_set}

    best_family = None
    best_score = 0

    for family_name, profile in VIRAL_FAMILIES.items():
        diag_lower = {m.lower() for m in profile.diagnostic_markers}
        supp_lower = {m.lower() for m in profile.supporting_markers}

        # Score: 2 points for diagnostic, 1 for supporting
        diagnostic_hits = len(markers_lower & diag_lower)
        supporting_hits = len(markers_lower & supp_lower)
        score = diagnostic_hits * 2 + supporting_hits

        if score > best_score:
            best_score = score
            best_family = family_name

    return best_family if best_score > 0 else None


@dataclass
class AssemblyMode:
    """Configuration for different assembly types."""

    mode: str  # "default", "fragmented", "relaxed"

    # Diversity requirements
    min_marker_types: int = 2
    require_diagnostic: bool = True

    # Score thresholds
    single_marker_min_score: float = 100.0
    neighbor_score_threshold: float = 5.0

    # Acceptance rules
    accept_isolated_anchors: bool = True
    accept_single_diagnostic: bool = False

    @classmethod
    def default(cls) -> "AssemblyMode":
        """Standard mode for well-assembled genomes."""
        return cls(
            mode="default",
            min_marker_types=2,
            require_diagnostic=True,
            single_marker_min_score=100.0,
            neighbor_score_threshold=5.0,
            accept_isolated_anchors=True,
            accept_single_diagnostic=False,
        )

    @classmethod
    def fragmented(cls) -> "AssemblyMode":
        """
        Mode for fragmented assemblies (MAGs, highly fragmented genomes).

        Key differences:
        - Accept single diagnostic markers with high scores
        - Lower diversity requirements
        - Don't require markers to cluster
        """
        return cls(
            mode="fragmented",
            min_marker_types=1,
            require_diagnostic=False,  # Any marker accepted
            single_marker_min_score=70.0,  # Lower threshold
            neighbor_score_threshold=0.0,  # Don't require clustering
            accept_isolated_anchors=True,
            accept_single_diagnostic=True,  # Key difference
        )

    @classmethod
    def relaxed(cls) -> "AssemblyMode":
        """
        Relaxed mode for exploratory analysis.

        Use when sensitivity is more important than specificity.
        """
        return cls(
            mode="relaxed",
            min_marker_types=1,
            require_diagnostic=False,
            single_marker_min_score=50.0,
            neighbor_score_threshold=0.0,
            accept_isolated_anchors=True,
            accept_single_diagnostic=True,
        )

    @classmethod
    def strict(cls) -> "AssemblyMode":
        """
        Strict mode for high-confidence predictions only.

        Use when specificity is critical.
        """
        return cls(
            mode="strict",
            min_marker_types=3,
            require_diagnostic=True,
            single_marker_min_score=150.0,
            neighbor_score_threshold=10.0,
            accept_isolated_anchors=False,
            accept_single_diagnostic=False,
        )


# Assembly mode presets
ASSEMBLY_MODES = {
    "default": AssemblyMode.default,
    "fragmented": AssemblyMode.fragmented,
    "relaxed": AssemblyMode.relaxed,
    "strict": AssemblyMode.strict,
}


def get_assembly_mode(mode_name: str) -> AssemblyMode:
    """Get assembly mode configuration by name."""
    if mode_name not in ASSEMBLY_MODES:
        raise ValueError(
            f"Unknown assembly mode: {mode_name}. "
            f"Available modes: {list(ASSEMBLY_MODES.keys())}"
        )
    return ASSEMBLY_MODES[mode_name]()


# =============================================================================
# REGION CLASSIFICATION
# =============================================================================

# Marker prefix patterns for classification (based on >95% purity seed allowlist)
NCLDV_MARKER_PREFIXES = {"gvogm", "gamadvirusmcp"}
# Explicit NCLDV MCP models that do not follow the standard prefixes
NCLDV_MARKER_EXACT = {"og1352", "og484"}
VP_MARKER_PREFIXES = {"vp_"}
PLV_MARKER_PREFIXES = {"plv_"}
MIRUS_MARKER_PREFIXES = {"mirus_"}


def classify_region_by_markers(
    marker_names: set[str],
    seed_marker_allowlist: Optional[list[str]] = None,
) -> str:
    """
    Classify a region as NCLDV, VP, PLV, MIRUS, or MIXED based on seed markers.

    Uses the seed marker allowlist to determine region classification.
    Only markers from the allowlist (>95% purity) count for classification.

    Args:
        marker_names: Set of HMM marker target names found in the region
        seed_marker_allowlist: List of high-purity seed markers from config

    Returns:
        Classification string: "NCLDV", "VP", "PLV", "MIRUS", "MIXED", or "UNKNOWN"
    """
    markers_lower = {m.lower() for m in marker_names}

    # If we have an allowlist, only count markers that are in it
    if seed_marker_allowlist:
        allowlist_lower = {m.lower() for m in seed_marker_allowlist}
        markers_lower = markers_lower & allowlist_lower

    if not markers_lower:
        return "UNKNOWN"

    # Count markers by category
    ncldv_count = 0
    vp_count = 0
    plv_count = 0
    mirus_count = 0

    for marker in markers_lower:
        # Check NCLDV markers
        if marker in NCLDV_MARKER_EXACT or any(
            marker.startswith(prefix) for prefix in NCLDV_MARKER_PREFIXES
        ):
            ncldv_count += 1
        # Check VP markers
        elif any(marker.startswith(prefix) for prefix in VP_MARKER_PREFIXES):
            vp_count += 1
        # Check PLV markers
        elif any(marker.startswith(prefix) for prefix in PLV_MARKER_PREFIXES):
            plv_count += 1
        # Check MIRUS markers
        elif any(marker.startswith(prefix) for prefix in MIRUS_MARKER_PREFIXES):
            mirus_count += 1

    # Combine VP and PLV counts for threshold checks (they're related families)
    vp_plv_count = vp_count + plv_count

    # Determine classification
    total = ncldv_count + vp_plv_count + mirus_count
    if total == 0:
        return "UNKNOWN"

    # Check for mixed classification (allow dominant family if others are weak)
    categories_present = sum([ncldv_count > 0, vp_plv_count > 0, mirus_count > 0])
    if categories_present > 1:
        if vp_plv_count >= 2 and ncldv_count <= 1 and mirus_count == 0:
            # One Preplasmiviricota lineage. The vp_/plv_ marker split is kept in
            # the vp_plv_subclass column, not in the class token.
            return "PPV"
        if ncldv_count >= 2 and vp_plv_count <= 1 and mirus_count == 0:
            return "NCLDV"
        if mirus_count >= 2 and ncldv_count == 0 and vp_plv_count <= 1:
            return "MIRUS"
        return "MIXED"

    # Single category
    if ncldv_count > 0:
        return "NCLDV"
    elif vp_plv_count > 0:
        return "PPV"
    elif mirus_count > 0:
        return "MIRUS"

    return "UNKNOWN"


def get_region_classification_summary(
    marker_names: set[str],
    seed_marker_allowlist: Optional[list[str]] = None,
) -> dict:
    """
    Get detailed classification summary for a region.

    Args:
        marker_names: Set of HMM marker target names found in the region
        seed_marker_allowlist: List of high-purity seed markers from config

    Returns:
        Dictionary with classification details:
        - classification: Primary classification (NCLDV, VP, PLV, MIRUS, MIXED, UNKNOWN)
        - ncldv_markers: Count of NCLDV seed markers
        - vp_plv_markers: Combined count of VP and PLV seed markers
        - mirus_markers: Count of MIRUS seed markers
        - seed_markers: List of seed markers in region
    """
    markers_lower = {m.lower() for m in marker_names}

    # If we have an allowlist, only count markers that are in it
    seed_markers_found = []
    if seed_marker_allowlist:
        allowlist_lower = {m.lower() for m in seed_marker_allowlist}
        for marker in marker_names:
            if marker.lower() in allowlist_lower:
                seed_markers_found.append(marker)
        markers_lower = {m.lower() for m in seed_markers_found}

    ncldv_count = 0
    vp_count = 0
    plv_count = 0
    mirus_count = 0

    for marker in markers_lower:
        if marker in NCLDV_MARKER_EXACT or any(
            marker.startswith(prefix) for prefix in NCLDV_MARKER_PREFIXES
        ):
            ncldv_count += 1
        elif any(marker.startswith(prefix) for prefix in VP_MARKER_PREFIXES):
            vp_count += 1
        elif any(marker.startswith(prefix) for prefix in PLV_MARKER_PREFIXES):
            plv_count += 1
        elif any(marker.startswith(prefix) for prefix in MIRUS_MARKER_PREFIXES):
            mirus_count += 1

    classification = classify_region_by_markers(marker_names, seed_marker_allowlist)

    return {
        "classification": classification,
        "ncldv_markers": ncldv_count,
        "vp_plv_markers": vp_count + plv_count,  # Combined for backwards compatibility
        "mirus_markers": mirus_count,
        "seed_markers": seed_markers_found,
    }
