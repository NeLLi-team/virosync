"""Per-genome ANI clustering of accepted EVEs, with MCP class propagation.

Accepted regions are compared all-vs-all with skani. Regions that cluster are
the same integrated element seen more than once in one genome, so an MCP-bearing
member's taxonomy class names the class of MCP-free members that carry no
diagnostic capsid of their own.

Clustering itself only rewrites ``taxonomy_class``; ``clustering_bonus`` stays
0.0 because confidence scoring already ran. Publishing decisions live in
``unsupported_eve_ids``, which the caller applies separately.
"""

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO

from virosync.output_contract import (
    LINEAGE_EVE_CLASSES,
    normalize_effective_eve_class,
)
from virosync.utils.atomic_write import atomic_write
from virosync.utils.path_safety import (
    require_strict_child,
    safe_filename_components,
)

logger = logging.getLogger(__name__)

# One threshold pair governs clustering, propagation, and the notebook network.
# Module constants rather than config keys: nothing sets them, and the retired
# phase2.ani_threshold key is already rejected by the config loader.
MIN_CLUSTER_ANI = 95.0
MIN_CLUSTER_ALIGNED_FRACTION = 50.0

_EDGE_HEADER = "eve_a\teve_b\tani\taf_a\taf_b\n"
_UNSKETCHABLE_STDERR = "No genomes/sketches found"


def cluster_accepted_eves(
    accepted_results: list,
    *,
    genome_fasta: Path,
    output_dir: Path,
    threads: int = 4,
) -> tuple[Path, list[tuple[str, str, float, float, float]]]:
    """Cluster accepted EVEs by ANI and propagate MCP-backed taxonomy classes.

    Writes every pair skani reports to ``phase3_synthesis/eve_ani_edges.tsv``,
    including pairs below :data:`MIN_CLUSTER_ANI`, so the notebook can filter at
    95% or higher from one source of truth. Returns that path and the pairs, so
    a caller that drops EVEs afterwards can recluster the survivors.

    The table is the raw comparison over the gate-accepted set, written before
    any later removal, so it can name an EVE that :func:`unsupported_eve_ids`
    later drops. Readers should intersect it with the published predictions
    rather than assume every endpoint was published.

    Fewer than two accepted EVEs, a missing skani binary, and sequences skani
    cannot sketch all yield a header-only edge table with no class change. Any
    other skani failure raises, because a real error must not be reported as a
    genome without clusters.
    """
    edges_path = output_dir / "phase3_synthesis" / "eve_ani_edges.tsv"
    require_strict_child(output_dir, edges_path)

    if len(accepted_results) < 2:
        logger.info(
            "Phase 3 ANI clustering: %d accepted EVE(s), nothing to cluster",
            len(accepted_results),
        )
        _write_edges(edges_path, [])
        return edges_path, []

    skani_bin = shutil.which("skani")
    if not skani_bin:
        logger.warning(
            "Phase 3 ANI clustering: skani not found on PATH; keeping every "
            "accepted EVE as a singleton and publishing an empty edge table"
        )
        _write_edges(edges_path, [])
        return edges_path, []

    sequences = _accepted_region_sequences(accepted_results, genome_fasta)
    pairs = _run_skani_triangle(skani_bin, sequences, threads)
    if pairs is None:
        logger.info(
            "Phase 3 ANI clustering: skani could not sketch the accepted EVEs; "
            "keeping every accepted EVE as a singleton"
        )
        _write_edges(edges_path, [])
        return edges_path, []

    _write_edges(edges_path, pairs)
    _apply_clusters(accepted_results, pairs)
    return edges_path, pairs


def unsupported_eve_ids(accepted_results: list) -> set[str]:
    """Return accepted EVEs with no viral evidence and no clustered relative.

    An EVE reaches ``UNKNOWN`` only when it carries no validated marker and no
    gene of its own returned an identity-qualified viral hit in the all-gene
    search. Nothing in it is viral, so it is most likely host sequence that the
    length rule admitted.

    ANI rescues it: sharing a cluster with a marker-bearing EVE makes it another
    copy of a real element whose own copy has simply decayed past recognition.

    Call after :func:`cluster_accepted_eves`, which populates ``cluster_id``.
    Unlike everything else in this module this changes which EVEs are published,
    so it is deliberately a separate step the caller opts into.
    """
    marker_bearing_clusters = {
        result.cluster_id
        for result in accepted_results
        if result.cluster_id != -1 and int(getattr(result, "hallmark_count", 0) or 0) > 0
    }
    unsupported = {
        result.eve_id
        for result in accepted_results
        if normalize_effective_eve_class(result.taxonomy_class) == "UNKNOWN"
        and result.cluster_id not in marker_bearing_clusters
    }
    if unsupported:
        logger.info(
            "Phase 3: dropping %d accepted EVE(s) with no viral evidence and no "
            "marker-bearing ANI relative: %s",
            len(unsupported),
            ", ".join(sorted(unsupported)),
        )
    return unsupported


def recluster_survivors(
    accepted_results: list,
    pairs: list[tuple[str, str, float, float, float]],
) -> None:
    """Recompute cluster fields after unsupported EVEs were dropped.

    Without this a survivor keeps an ``ani_cluster_size`` that counts a member
    no longer published. Propagation is idempotent here: every surviving member
    already holds the class its donor gave it, so nothing changes class.
    """
    surviving = {result.eve_id for result in accepted_results}
    for result in accepted_results:
        result.cluster_id = -1
        result.cluster_size = 1
        result.max_cluster_ani = 0.0
    _apply_clusters(
        accepted_results,
        [pair for pair in pairs if pair[0] in surviving and pair[1] in surviving],
    )


def _accepted_region_sequences(
    accepted_results: list,
    genome_fasta: Path,
) -> dict[str, str]:
    """Extract each accepted region's nucleotides from the masked genome.

    Every accepted region must yield a sequence. A region whose scaffold is
    absent, whose slice is empty, or whose ID repeats is an inconsistency between
    the accepted set and the genome it came from, and clustering a subset of the
    accepted EVEs would read as a genome with fewer relatives than it has.
    """
    wanted = {result.scaffold for result in accepted_results}
    scaffolds = {
        record.id: str(record.seq)
        for record in SeqIO.parse(genome_fasta, "fasta")
        if record.id in wanted
    }
    sequences = {}
    for result in accepted_results:
        sequence = scaffolds.get(result.scaffold, "")[result.start:result.end]
        if not sequence:
            raise RuntimeError(
                "accepted EVE has no sequence in the genome used for ANI "
                f"clustering: {result.eve_id} ({result.scaffold}:"
                f"{result.start}-{result.end} in {genome_fasta})"
            )
        sequences[result.eve_id] = sequence
    if len(sequences) != len(accepted_results):
        raise RuntimeError(
            "accepted EVE IDs are not unique within one genome: "
            f"{len(accepted_results)} results, {len(sequences)} distinct IDs"
        )
    return sequences


def _run_skani_triangle(
    skani_bin: str,
    sequences: dict[str, str],
    threads: int,
) -> list[tuple[str, str, float, float, float]] | None:
    """Return every reported pair, or None when skani cannot sketch the input."""
    workdir = Path(tempfile.mkdtemp(prefix="virosync_eve_ani_"))
    try:
        components = safe_filename_components(sequences, label="EVE ID")
        list_path = workdir / "eve_fasta_list.txt"
        eve_by_file = {}
        with open(list_path, "w") as listing:
            for eve_id, sequence in sequences.items():
                fasta_path = require_strict_child(
                    workdir,
                    workdir / f"{components[eve_id]}.fna",
                )
                fasta_path.write_text(f">{eve_id}\n{sequence}\n")
                listing.write(f"{fasta_path}\n")
                eve_by_file[str(fasta_path)] = eve_id

        ani_path = workdir / "eve_ani_triangle.tsv"
        completed = subprocess.run(
            [
                skani_bin,
                "triangle",
                "-l", str(list_path),
                "-E",                       # sparse row-per-pair output
                "--medium",                 # accuracy mode for short sequences
                "-m", "200",                # low marker compression
                "-s", "80",                 # screen threshold
                "--min-af", str(int(MIN_CLUSTER_ALIGNED_FRACTION)),
                "-t", str(max(1, int(threads))),
                "-o", str(ani_path),
            ],
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            if _UNSKETCHABLE_STDERR in (completed.stderr or ""):
                return None
            raise RuntimeError(
                "skani triangle failed for accepted EVE clustering: "
                f"{completed.stderr}"
            )
        return _parse_skani_pairs(ani_path, eve_by_file)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _parse_skani_pairs(
    ani_path: Path,
    eve_by_file: dict[str, str],
) -> list[tuple[str, str, float, float, float]]:
    """Parse sparse skani output back into (eve_a, eve_b, ani, af_a, af_b).

    Strict on purpose: a row this cannot read is a reported pair that would
    otherwise vanish from the edge table and from clustering.
    """
    pairs = []
    with open(ani_path) as handle:
        next(handle, None)  # column header
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[0] not in eve_by_file or fields[1] not in eve_by_file:
                raise RuntimeError(
                    f"unreadable skani pair in {ani_path}: {line.rstrip()!r}"
                )
            pairs.append(
                (
                    eve_by_file[fields[0]],
                    eve_by_file[fields[1]],
                    float(fields[2]),
                    float(fields[3]),
                    float(fields[4]),
                )
            )
    return pairs


def _write_edges(
    edges_path: Path,
    pairs: list[tuple[str, str, float, float, float]],
) -> None:
    """Write the edge table, sorted so two runs of one genome agree byte for byte."""
    rows = [
        f"{eve_a}\t{eve_b}\t{ani:.4f}\t{af_a:.4f}\t{af_b:.4f}\n"
        for eve_a, eve_b, ani, af_a, af_b in sorted(pairs)
    ]
    atomic_write(edges_path, _EDGE_HEADER + "".join(rows))


def _apply_clusters(
    accepted_results: list,
    pairs: list[tuple[str, str, float, float, float]],
) -> None:
    """Populate the cluster fields and propagate within every multi-EVE cluster."""
    parent = {result.eve_id: result.eve_id for result in accepted_results}

    def find(eve_id: str) -> str:
        while parent[eve_id] != eve_id:
            parent[eve_id] = parent[parent[eve_id]]
            eve_id = parent[eve_id]
        return eve_id

    max_ani: dict[str, float] = {}
    # skani was run with --min-af MIN_CLUSTER_ALIGNED_FRACTION, so every pair it
    # reported already clears the aligned-fraction bar; only ANI is left to check.
    for eve_a, eve_b, ani, _af_a, _af_b in pairs:
        if ani < MIN_CLUSTER_ANI:
            continue
        root_a, root_b = find(eve_a), find(eve_b)
        if root_a != root_b:
            parent[root_a] = root_b
        for eve_id in (eve_a, eve_b):
            max_ani[eve_id] = max(max_ani.get(eve_id, 0.0), ani)

    components: dict[str, list] = {}
    for result in accepted_results:
        components.setdefault(find(result.eve_id), []).append(result)

    # Largest cluster first, ties broken by lowest member ID: cluster IDs are
    # then a function of the edge set alone, so a rerun reproduces them.
    clusters = sorted(
        (members for members in components.values() if len(members) > 1),
        key=lambda members: (
            -len(members),
            min(result.eve_id for result in members),
        ),
    )
    propagated = 0
    for cluster_id, members in enumerate(clusters):
        for result in members:
            result.cluster_id = cluster_id
            result.cluster_size = len(members)
            result.max_cluster_ani = max_ani.get(result.eve_id, 0.0)
        propagated += _propagate_mcp_class(members)

    logger.info(
        "Phase 3 ANI clustering: %d cluster(s) over %d accepted EVE(s) at "
        ">=%.0f%% ANI and >=%.0f%% aligned fraction; %d MCP-free EVE(s) "
        "relabelled from an MCP-bearing relative",
        len(clusters),
        len(accepted_results),
        MIN_CLUSTER_ANI,
        MIN_CLUSTER_ALIGNED_FRACTION,
        propagated,
    )


def _propagate_mcp_class(members: list) -> int:
    """Give members the lineage class their MCP-decided relatives agree on.

    A donor is an EVE whose own class an MCP marker's vote decided
    (``taxonomy_class_from_mcp``), not merely one with ``has_mcp`` set: that
    flag is also raised by a capsid annotation, a structural jelly-roll call,
    and phylogenetic evidence, none of which cast a taxonomy vote, so an EVE
    could otherwise donate a class no MCP ever voted for.

    Returns the number of members whose class actually changed. Disagreeing
    donors, and an agreed class that names no lineage, propagate nothing: each
    remaining member then keeps the class its own vote produced.

    VIRAL_UNKNOWN never propagates. Handing it on only destroys evidence: a
    member whose own weighted vote resolved a lineage would be pushed back to
    VIRAL_UNKNOWN by a relative whose capsid was itself ambiguous. A donated
    class is worth taking only when it names a lineage the member lacks.
    """
    with_mcp = sorted(
        (result for result in members if result.taxonomy_class_from_mcp),
        key=lambda result: result.eve_id,
    )
    without_mcp = [
        result for result in members if not result.taxonomy_class_from_mcp
    ]
    if not with_mcp or not without_mcp:
        return 0

    inherited_classes = {
        normalize_effective_eve_class(result.taxonomy_class)
        for result in with_mcp
    }
    if len(inherited_classes) != 1:
        return 0
    inherited = next(iter(inherited_classes))
    if inherited not in LINEAGE_EVE_CLASSES:
        return 0

    # Lowest-sorting MCP-bearing EVE, so the recorded source survives a rerun.
    source_eve_id = with_mcp[0].eve_id
    changed = 0
    for result in without_mcp:
        if normalize_effective_eve_class(result.taxonomy_class) == inherited:
            continue
        result.taxonomy_class_before_ani = result.taxonomy_class
        result.taxonomy_class_propagated_from = source_eve_id
        result.taxonomy_class = inherited
        changed += 1
    return changed
