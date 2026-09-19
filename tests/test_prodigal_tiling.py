import json
import subprocess
from pathlib import Path
from typing import TextIO

import pytest
from Bio import SeqIO

from virosync.orchestration._flows.single_genome.run_state import invalidate_from_phase
from virosync.pipeline.phase0 import prodigal


def _write_untiled_outputs(cmd: list[str]) -> tuple[Path, Path]:
    """Write matching Prodigal FAA/GFF output for each input record."""
    records = list(SeqIO.parse(cmd[cmd.index("-i") + 1], "fasta"))
    proteins = Path(cmd[cmd.index("-a") + 1])
    gff = Path(cmd[cmd.index("-o") + 1])
    proteins.write_text(
        "".join(
            f">{record.id}_1 # 1 # 6 # 1 # ID={index}_1;partial=00\nM*\n" for index, record in enumerate(records, 1)
        )
    )
    gff.write_text(
        "##gff-version 3\n"
        + "".join(
            f'# Sequence Data: seqnum={index};seqlen={len(record.seq)};seqhdr="{record.description}"\n'
            f"{record.id}\tProdigal\tCDS\t1\t6\t.\t+\t0\tID={index}_1\n"
            for index, record in enumerate(records, 1)
        )
    )
    return proteins, gff


@pytest.mark.parametrize("defect", ["partial_header", "missing_protein", "missing_sequence"])
def test_untiled_parallel_rejects_invalid_success_before_merge(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, defect: str
) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">first\nATGTAA\n>second\nATGTAA\n")
    output_dir = tmp_path / "phase0"
    output_dir.mkdir()
    merged = output_dir / "proteome.fasta"

    def invalid_output(cmd: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        proteins, gff = _write_untiled_outputs(cmd)
        if defect == "partial_header":
            proteins.write_text(proteins.read_text().split(">second")[0] + ">second_1 # 1 # ")
        elif defect == "missing_protein":
            proteins.write_text(proteins.read_text().split(">second")[0])
        else:
            gff.write_text(gff.read_text().split("# Sequence Data: seqnum=2")[0])
            proteins.write_text(proteins.read_text().split(">second")[0])
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(prodigal.subprocess, "run", invalid_output)
    with pytest.raises(RuntimeError):
        prodigal._run_prodigal_parallel(genome, output_dir, merged, output_dir / "genes.gff", threads=1)
    assert not merged.exists()
    assert not list(output_dir.glob("tmp*"))
    attempts = list((tmp_path / "prodigal_diagnostics" / "phase0").glob("*/attempt.json"))
    assert len(attempts) == 1
    assert json.loads(attempts[0].read_text())["returncode"] == 0
    assert (attempts[0].parent / "chunk_0.fasta").read_bytes() == genome.read_bytes()


@pytest.mark.parametrize(("retry_returncode", "valid"), [(0, True), (0, False), (1, True), (-6, True)])
def test_untiled_nonzero_requires_valid_record_retries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, retry_returncode: int, valid: bool
) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">first\nATGTAA\n>second\nATGTAA\n")
    output_dir = tmp_path / "phase0"
    output_dir.mkdir()
    calls = 0

    def retry_output(
        cmd: list[str], *, stderr: TextIO | int, check: bool, **_kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        proteins, _gff = _write_untiled_outputs(cmd)
        returncode = 23 if calls == 1 else retry_returncode
        if calls == 1:
            proteins.write_text(proteins.read_text().split(">second")[0])
        elif not valid:
            proteins.write_text("")
        if not isinstance(stderr, int):
            stderr.write("free(): invalid pointer\n")
        if check and returncode:
            raise subprocess.CalledProcessError(returncode, cmd)
        return subprocess.CompletedProcess(cmd, returncode)

    monkeypatch.setattr(prodigal.subprocess, "run", retry_output)
    merged = output_dir / "proteome.fasta"
    if retry_returncode == 0 and valid:
        _, genes = prodigal._run_prodigal_parallel(genome, output_dir, merged, output_dir / "genes.gff", threads=1)
        assert [gene.gene_id for gene in genes] == ["first_1", "second_1"]
        assert calls == 3
    else:
        with pytest.raises(RuntimeError):
            prodigal._run_prodigal_parallel(genome, output_dir, merged, output_dir / "genes.gff", threads=1)
        assert not merged.exists()
    attempts = list((tmp_path / "prodigal_diagnostics" / "phase0").glob("*/attempt.json"))
    assert any(json.loads(path.read_text())["returncode"] == 23 for path in attempts)
    assert any(
        (path.parent / "chunk_0.faa").read_text() == ">first_1 # 1 # 6 # 1 # ID=1_1;partial=00\nM*\n"
        for path in attempts
        if (path.parent / "chunk_0.faa").exists()
    )


def test_valid_untiled_parallel_preserves_bytes_and_order(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">first\nATGTAA\n>second\nATGTAA\n>third\nATGTAA\n>fourth\nATGTAA\n")
    output_dir = tmp_path / "phase0"
    output_dir.mkdir()

    def valid_output(cmd: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        _write_untiled_outputs(cmd)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(prodigal.subprocess, "run", valid_output)
    proteins, genes = prodigal._run_prodigal_parallel(
        genome, output_dir, output_dir / "proteome.fasta", output_dir / "genes.gff", threads=2
    )
    assert proteins.read_text() == (
        ">first_1 # 1 # 6 # 1 # ID=1_1;partial=00\nM*\n"
        ">third_1 # 1 # 6 # 1 # ID=2_1;partial=00\nM*\n"
        ">second_1 # 1 # 6 # 1 # ID=1_1;partial=00\nM*\n"
        ">fourth_1 # 1 # 6 # 1 # ID=2_1;partial=00\nM*\n"
    )
    assert [gene.gene_id for gene in genes] == ["first_1", "third_1", "second_1", "fourth_1"]


@pytest.mark.parametrize(("returncode", "valid"), [(0, True), (0, False), (-6, True)])
def test_serial_prodigal_requires_valid_zero_exit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, returncode: int, valid: bool
) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">first\nATGTAA\n>second\nATGTAA\n")
    output_dir = tmp_path / "phase0"
    output_dir.mkdir()

    def serial_output(
        cmd: list[str], *, stderr: TextIO | int, check: bool, **_kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        proteins, _gff = _write_untiled_outputs(cmd)
        if not valid:
            proteins.write_text(proteins.read_text().split(">second")[0])
        if not isinstance(stderr, int):
            stderr.write("free(): invalid pointer\n")
        if check and returncode:
            raise subprocess.CalledProcessError(returncode, cmd)
        return subprocess.CompletedProcess(cmd, returncode)

    monkeypatch.setattr(prodigal.subprocess, "run", serial_output)
    args = (genome, output_dir, output_dir / "proteome.fasta", output_dir / "genes.gff", "meta")
    if valid and returncode == 0:
        _, genes = prodigal._run_prodigal_single(*args)
        assert [gene.gene_id for gene in genes] == ["first_1", "second_1"]
        assert not (output_dir / "prodigal.stderr").exists()
    else:
        with pytest.raises(RuntimeError, match="diagnostics:"):
            prodigal._run_prodigal_single(*args)


@pytest.mark.parametrize("parallel", [False, True])
def test_prodigal_diagnostics_survive_phase0_invalidation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parallel: bool
) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">first\nATGTAA\n")
    output_dir = tmp_path / "phase0"

    def failed_output(
        cmd: list[str], *, stderr: TextIO | int, check: bool, **_kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        _write_untiled_outputs(cmd)
        if not isinstance(stderr, int):
            stderr.write("failed process\n")
        if check:
            raise subprocess.CalledProcessError(23, cmd)
        return subprocess.CompletedProcess(cmd, 23)

    monkeypatch.setattr(prodigal.subprocess, "run", failed_output)
    for _attempt in range(2):
        output_dir.mkdir()
        with pytest.raises(RuntimeError, match="diagnostics:"):
            if parallel:
                prodigal._run_prodigal_parallel(
                    genome, output_dir, output_dir / "proteome.fasta", output_dir / "genes.gff", threads=1
                )
            else:
                prodigal._run_prodigal_single(
                    genome, output_dir, output_dir / "proteome.fasta", output_dir / "genes.gff", "meta"
                )
        invalidate_from_phase(tmp_path, from_phase=0)
        assert not output_dir.exists()
    attempts = list((tmp_path / "prodigal_diagnostics" / "phase0").glob("*/attempt.json"))
    assert len(attempts) == (4 if parallel else 2)
    for path in attempts:
        metadata = json.loads(path.read_text())
        assert metadata["returncode"] == 23
        filenames = {
            "serial": ("genome.fasta", "proteome.fasta", "genes.gff", "prodigal.stderr"),
            "chunk": ("chunk_0.fasta", "chunk_0.faa", "chunk_0.gff", "chunk_0.stderr"),
            "record": ("record_0.fasta", "record_0.faa", "record_0.gff", "record_0.stderr"),
        }[metadata["stage"]]
        assert set(metadata["artifacts"]) == set(filenames)
        input_fasta, proteins, gff, stderr = [path.parent / name for name in filenames]
        assert input_fasta.read_bytes() == genome.read_bytes()
        assert proteins.read_text().startswith(">first_1")
        assert gff.is_file()
        assert stderr.read_text() == "failed process\n"


@pytest.mark.parametrize(
    ("start", "end", "expected"),
    [
        (8, 10, False),
        (9, 11, True),
        (9, 12, True),
        (18, 20, True),
        (19, 21, False),
        (19, 22, False),
    ],
)
def test_owned_midpoint_boundaries(start: int, end: int, expected: bool) -> None:
    assert prodigal._owns_midpoint(start, end, 10, 20) is expected


def test_long_scaffold_tiles_are_rebased_and_renumbered(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">long_scaffold\nACGTACGTACGTA\n>short_scaffold\nACGT\n")
    monkeypatch.setattr(prodigal, "_LONG_SCAFFOLD_BP", 6)
    monkeypatch.setattr(prodigal, "_TILE_CORE_BP", 4)
    monkeypatch.setattr(prodigal, "_TILE_OVERLAP_BP", 2)

    def fake_prodigal(
        chunk_fasta: str,
        chunk_out: str,
        *_args,
    ) -> str:
        with Path(chunk_out).open("w") as handle:
            for record in SeqIO.parse(chunk_fasta, "fasta"):
                handle.write(f">{record.id}_1 # 2 # 4 # 1 # ID=1_1;partial=00;genetic_code=11\nMKK\n")
        return chunk_out

    monkeypatch.setattr(prodigal, "_run_prodigal_on_chunk", fake_prodigal)
    proteins, genes = prodigal._run_prodigal_parallel(
        genome,
        tmp_path,
        tmp_path / "proteome.fasta",
        tmp_path / "genes.gff",
        threads=2,
    )

    assert [gene.gene_id for gene in genes] == [
        "long_scaffold_1",
        "long_scaffold_2",
        "long_scaffold_3",
        "long_scaffold_4",
        "short_scaffold_1",
    ]
    assert [(gene.scaffold, gene.start, gene.end) for gene in genes] == [
        ("long_scaffold", 1, 4),
        ("long_scaffold", 2, 5),
        ("long_scaffold", 5, 8),
        ("long_scaffold", 8, 11),
        ("short_scaffold", 1, 4),
    ]
    assert prodigal._TILE_ID_PREFIX not in proteins.read_text()
    loaded = prodigal.load_gene_predictions(proteins)
    assert sum(len(items) for items in loaded.values()) == len(genes)
    assert [
        (gene.scaffold, gene.start, gene.end, gene.strand) for scaffold in loaded.values() for gene in scaffold
    ] == [(gene.scaffold, gene.start, gene.end, gene.strand) for gene in genes]


def test_tiled_record_rejects_valid_output_after_nonzero_exit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    work_dir = tmp_path / "temporary"
    work_dir.mkdir()
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    chunk_fasta = work_dir / "chunk.fasta"
    chunk_out = work_dir / "chunk.faa"
    chunk_fasta.write_text(f">{record_id}\nACGT\n")
    calls = 0

    def abort_with_valid_output(
        cmd: list[str], *, stderr: TextIO | int, **_kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        nonlocal calls
        calls += 1
        if not isinstance(stderr, int):
            stderr.write("free(): invalid pointer\n")
        if calls == 2:
            output = Path(cmd[cmd.index("-a") + 1])
            gff = Path(cmd[cmd.index("-o") + 1])
            output.write_text(f">{record_id}_1 # 1 # 3 # 1 # ID=1_1;partial=00\nM\n")
            gff.write_text(
                "##gff-version 3\n"
                f'# Sequence Data: seqnum=1;seqlen=4;seqhdr="{record_id}"\n'
                f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
            )
        return subprocess.CompletedProcess(cmd, -6)

    monkeypatch.setattr(prodigal.subprocess, "run", abort_with_valid_output)

    with pytest.raises(RuntimeError, match="nonzero Prodigal-GV exit: -6; diagnostics:"):
        prodigal._run_prodigal_on_chunk(
            str(chunk_fasta),
            str(chunk_out),
            tmp_path / "diagnostics",
        )

    assert calls == 2
    attempts = {
        json.loads(path.read_text())["stage"]: path.parent for path in (tmp_path / "diagnostics").glob("*/attempt.json")
    }
    assert set(attempts) == {"chunk", "record"}
    assert json.loads((attempts["chunk"] / "attempt.json").read_text())["returncode"] == -6
    assert json.loads((attempts["record"] / "attempt.json").read_text())["returncode"] == -6
    assert (attempts["chunk"] / "chunk.fasta").read_bytes() == chunk_fasta.read_bytes()
    assert (attempts["chunk"] / "chunk.stderr").read_text() == "free(): invalid pointer\n"
    assert (attempts["record"] / "record_0.faa").read_text().startswith(f">{record_id}_1")
    assert (attempts["record"] / "record_0.gff").is_file()
    assert (attempts["record"] / "record_0.stderr").read_text() == "free(): invalid pointer\n"
    prodigal._validate_tiled_prodigal_output(
        attempts["record"] / "record_0.fasta",
        attempts["record"] / "record_0.faa",
        attempts["record"] / "record_0.gff",
    )


def test_strict_validation_rejects_no_delimiter_final_header(
    tmp_path: Path,
) -> None:
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    input_fasta = tmp_path / "input.fasta"
    proteins_faa = tmp_path / "proteins.faa"
    genes_gff = tmp_path / "genes.gff"
    input_fasta.write_text(f">{record_id}\n{'A' * 12}\n")
    proteins_faa.write_text(f">{record_id}_1 # 1 # 3 # 1 # ID=1_1\nM\n>__virosync_til\nM")
    genes_gff.write_text(
        "##gff-version 3\n"
        f'# Sequence Data: seqnum=1;seqlen=12;seqhdr="{record_id}"\n'
        f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
        f"{record_id}\tProdigal\tCDS\t7\t9\t.\t+\t0\tID=1_2\n"
    )

    with pytest.raises(RuntimeError, match="unparseable Prodigal protein header"):
        prodigal._validate_tiled_prodigal_output(
            input_fasta,
            proteins_faa,
            genes_gff,
        )


def test_tiled_chunk_retains_mismatched_faa_and_gff(
    tmp_path: Path,
    monkeypatch,
) -> None:
    work_dir = tmp_path / "temporary"
    work_dir.mkdir()
    chunk_fasta = work_dir / "chunk_2.fasta"
    chunk_out = work_dir / "chunk_2.faa"
    chunk_fasta.write_text(f">{prodigal._TILE_ID_PREFIX}0\nACGT\n")

    def mismatched_output(cmd, **kwargs):
        output = Path(cmd[cmd.index("-a") + 1])
        gff = Path(cmd[cmd.index("-o") + 1])
        output.write_text(f">{prodigal._TILE_ID_PREFIX}0_1 # 1 # 3 # 1 # ID=1_1\nM\n")
        gff.write_text(
            "##gff-version 3\n"
            f'# Sequence Data: seqnum=1;seqlen=4;seqhdr="{prodigal._TILE_ID_PREFIX}0"\n'
            f"{prodigal._TILE_ID_PREFIX}0\tProdigal\tCDS\t2\t4\t.\t+\t0\tID=1_1\n"
        )
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(prodigal.subprocess, "run", mismatched_output)

    with pytest.raises(RuntimeError, match="coordinates differ"):
        prodigal._run_prodigal_on_chunk(str(chunk_fasta), str(chunk_out), tmp_path / "diagnostics")

    failure_dir = next((tmp_path / "diagnostics").glob("chunk-*"))
    assert (failure_dir / "chunk_2.fasta").exists()
    assert (failure_dir / "chunk_2.faa").exists()
    assert (failure_dir / "chunk_2.gff").exists()


def test_tiled_chunk_rejects_truncated_protein(
    tmp_path: Path,
    monkeypatch,
) -> None:
    work_dir = tmp_path / "temporary"
    work_dir.mkdir()
    chunk_fasta = work_dir / "chunk_3.fasta"
    chunk_out = work_dir / "chunk_3.faa"
    chunk_fasta.write_text(f">{prodigal._TILE_ID_PREFIX}0\nACGTAC\n")

    def truncated_output(cmd, **kwargs):
        output = Path(cmd[cmd.index("-a") + 1])
        gff = Path(cmd[cmd.index("-o") + 1])
        output.write_text(f">{prodigal._TILE_ID_PREFIX}0_1 # 1 # 6 # 1 # ID=1_1\nM\n")
        gff.write_text(
            "##gff-version 3\n"
            f'# Sequence Data: seqnum=1;seqlen=6;seqhdr="{prodigal._TILE_ID_PREFIX}0"\n'
            f"{prodigal._TILE_ID_PREFIX}0\tProdigal\tCDS\t1\t6\t.\t+\t0\tID=1_1\n"
        )
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(prodigal.subprocess, "run", truncated_output)

    with pytest.raises(RuntimeError, match="protein length does not match"):
        prodigal._run_prodigal_on_chunk(str(chunk_fasta), str(chunk_out), tmp_path / "diagnostics")


def test_tiled_merge_rejects_unmapped_scaffold_ids(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">long_scaffold\nACGTACGT\n")
    monkeypatch.setattr(prodigal, "_LONG_SCAFFOLD_BP", 4)
    monkeypatch.setattr(prodigal, "_TILE_CORE_BP", 4)
    monkeypatch.setattr(prodigal, "_TILE_OVERLAP_BP", 1)

    def fake_prodigal(
        chunk_fasta: str,
        chunk_out: str,
        *_args,
    ) -> str:
        Path(chunk_out).write_text(">rogue_1 # 1 # 3 # 1 # ID=1_1;partial=00\nMKK\n")
        return chunk_out

    monkeypatch.setattr(prodigal, "_run_prodigal_on_chunk", fake_prodigal)

    with pytest.raises(RuntimeError, match="could not be mapped"):
        prodigal._run_prodigal_parallel(
            genome,
            tmp_path,
            tmp_path / "proteome.fasta",
            tmp_path / "genes.gff",
            threads=2,
        )


def test_tiled_genome_sets_diagnostic_root_for_every_chunk(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">long\nACGTACGT\n>short\nACGT\n")
    monkeypatch.setattr(prodigal, "_LONG_SCAFFOLD_BP", 4)
    monkeypatch.setattr(prodigal, "_TILE_CORE_BP", 4)
    monkeypatch.setattr(prodigal, "_TILE_OVERLAP_BP", 1)
    observed: list[Path] = []

    def fake_prodigal(
        _chunk_fasta: str,
        chunk_out: str,
        diagnostic_root: Path,
    ) -> str:
        observed.append(diagnostic_root)
        Path(chunk_out).write_text("")
        return chunk_out

    monkeypatch.setattr(prodigal, "_run_prodigal_on_chunk", fake_prodigal)

    prodigal._run_prodigal_parallel(
        genome,
        tmp_path,
        tmp_path / "proteome.fasta",
        tmp_path / "genes.gff",
        threads=3,
    )

    assert len(observed) == 3
    assert all(root == tmp_path.parent / "prodigal_diagnostics" / tmp_path.name for root in observed)


def test_long_scaffold_is_tiled_with_one_thread(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fasta"
    output_dir = tmp_path / "output"
    genome.write_text(">long\nACGTACG\n")
    monkeypatch.setattr(prodigal, "_LONG_SCAFFOLD_BP", 6)
    monkeypatch.setattr(prodigal.shutil, "which", lambda _name: "/bin/prodigal-gv")
    sentinel = output_dir / "proteome.fasta"

    def fake_parallel(*args):
        assert args[-1] == 1
        return sentinel, []

    monkeypatch.setattr(prodigal, "_run_prodigal_parallel", fake_parallel)

    assert prodigal.run_prodigal_genome(
        genome,
        output_dir,
        threads=1,
    ) == (sentinel, [])


def test_input_scaffold_rejects_reserved_tile_prefix(tmp_path: Path) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(f">{prodigal._TILE_ID_PREFIX}original\nACGT\n")

    with pytest.raises(RuntimeError, match="reserved tile prefix"):
        prodigal._run_prodigal_parallel(
            genome,
            tmp_path,
            tmp_path / "proteome.fasta",
            tmp_path / "genes.gff",
            threads=2,
        )
