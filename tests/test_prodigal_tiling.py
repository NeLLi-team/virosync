from pathlib import Path
import subprocess

from Bio import SeqIO
import pytest

from virosync.pipeline.phase0 import prodigal


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
    genome.write_text(
        ">long_scaffold\nACGTACGTACGTA\n"
        ">short_scaffold\nACGT\n"
    )
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
                handle.write(
                    f">{record.id}_1 # 2 # 4 # 1 # "
                    "ID=1_1;partial=00;genetic_code=11\nMKK\n"
                )
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
    assert [
        (gene.scaffold, gene.start, gene.end) for gene in genes
    ] == [
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
        (gene.scaffold, gene.start, gene.end, gene.strand)
        for scaffold in loaded.values()
        for gene in scaffold
    ] == [
        (gene.scaffold, gene.start, gene.end, gene.strand)
        for gene in genes
    ]


def test_tiled_chunk_retries_each_record_after_nonzero_exit(
    tmp_path: Path,
    monkeypatch,
) -> None:
    work_dir = tmp_path / "temporary"
    work_dir.mkdir()
    chunk_fasta = work_dir / "chunk.fasta"
    chunk_out = work_dir / "chunk.faa"
    chunk_fasta.write_text(f">{prodigal._TILE_ID_PREFIX}0\nACGT\n")
    calls = 0

    def fail_once(cmd, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            return subprocess.CompletedProcess(cmd, -6)
        output = Path(cmd[cmd.index("-a") + 1])
        gff = Path(cmd[cmd.index("-o") + 1])
        output.write_text(
            f">{prodigal._TILE_ID_PREFIX}0_1 # 1 # 3 # 1 # "
            "ID=1_1;partial=00\nM\n"
        )
        gff.write_text(
            "##gff-version 3\n"
            f"# Sequence Data: seqnum=1;seqlen=4;seqhdr=\"{prodigal._TILE_ID_PREFIX}0\"\n"
            f"{prodigal._TILE_ID_PREFIX}0\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
        )
        kwargs["stderr"].write("free(): invalid pointer\n")
        return subprocess.CompletedProcess(cmd, -6)

    monkeypatch.setattr(prodigal.subprocess, "run", fail_once)

    assert prodigal._run_prodigal_on_chunk(
        str(chunk_fasta),
        str(chunk_out),
        True,
        {f"{prodigal._TILE_ID_PREFIX}0": (0, 4)},
    ) == str(chunk_out)
    assert calls == 2
    assert chunk_out.read_text().startswith(
        f">{prodigal._TILE_ID_PREFIX}0_1"
    )
    assert (
        tmp_path
        / "accepted_cleanup_aborts"
        / f"{prodigal._TILE_ID_PREFIX}0.json"
    ).exists()


def test_tiled_chunk_accepts_only_incomplete_unowned_suffix(
    tmp_path: Path,
    monkeypatch,
) -> None:
    work_dir = tmp_path / "temporary"
    work_dir.mkdir()
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    chunk_fasta = work_dir / "chunk.fasta"
    chunk_out = work_dir / "chunk.faa"
    chunk_fasta.write_text(f">{record_id}\nACGTACGTACGTACGTAC\n")
    calls = 0

    def truncated_overlap(cmd, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            return subprocess.CompletedProcess(cmd, -6)
        output = Path(cmd[cmd.index("-a") + 1])
        gff = Path(cmd[cmd.index("-o") + 1])
        output.write_text(
            f">{record_id}_1 # 1 # 3 # 1 # ID=1_1\nM\n"
            f">{record_id}_2 # 7 # 12 # 1 # ID=1_2\nM"
        )
        gff.write_text(
            "##gff-version 3\n"
            f"# Sequence Data: seqnum=1;seqlen=18;seqhdr=\"{record_id}\"\n"
            f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
            f"{record_id}\tProdigal\tCDS\t7\t12\t.\t+\t0\tID=1_2\n"
            f"{record_id}\tProdigal\tCDS\t13\t18\t.\t+\t0\tID=1_3\n"
        )
        kwargs["stderr"].write("free(): invalid pointer\n")
        return subprocess.CompletedProcess(cmd, -6)

    monkeypatch.setattr(prodigal.subprocess, "run", truncated_overlap)

    assert prodigal._run_prodigal_on_chunk(
        str(chunk_fasta),
        str(chunk_out),
        True,
        {record_id: (0, 6)},
    ) == str(chunk_out)
    audit = (
        tmp_path / "accepted_cleanup_aborts" / f"{record_id}.json"
    ).read_text()
    assert '"start_0based": 6' in audit
    assert '"start_0based": 12' in audit
    assert chunk_out.read_text().endswith("M\n")
    assert f"{record_id}_2" not in chunk_out.read_text()


def test_tiled_chunk_rejects_owned_loss_without_gff_core_coverage(
    tmp_path: Path,
    monkeypatch,
) -> None:
    work_dir = tmp_path / "temporary"
    work_dir.mkdir()
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    chunk_fasta = work_dir / "chunk.fasta"
    chunk_out = work_dir / "chunk.faa"
    chunk_fasta.write_text(f">{record_id}\nACGTACGTAC\n")
    calls = 0

    def missing_owned_call(cmd, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            return subprocess.CompletedProcess(cmd, -6)
        output = Path(cmd[cmd.index("-a") + 1])
        gff = Path(cmd[cmd.index("-o") + 1])
        output.write_text(f">{record_id}_1 # 1 # 3 # 1 # ID=1_1\nM\n")
        gff.write_text(
            "##gff-version 3\n"
            f"# Sequence Data: seqnum=1;seqlen=10;seqhdr=\"{record_id}\"\n"
            f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
            f"{record_id}\tProdigal\tCDS\t7\t9\t.\t+\t0\tID=1_2\n"
        )
        kwargs["stderr"].write("free(): invalid pointer\n")
        return subprocess.CompletedProcess(cmd, -6)

    monkeypatch.setattr(prodigal.subprocess, "run", missing_owned_call)

    with pytest.raises(RuntimeError, match="GFF does not cover the owned core"):
        prodigal._run_prodigal_on_chunk(
            str(chunk_fasta),
            str(chunk_out),
            True,
            {record_id: (0, 10)},
        )


def _write_reconstruction_fixture(
    tmp_path: Path,
    first_protein: str = "M*",
) -> tuple[Path, Path, Path, str]:
    record_id = f"{prodigal._TILE_ID_PREFIX}reconstruct"
    input_fasta = tmp_path / "input.fasta"
    proteins_faa = tmp_path / "proteins.faa"
    genes_gff = tmp_path / "genes.gff"
    input_fasta.write_text(
        f">{record_id}\nGTGTAATTACATTTACACATGTAA\n"
    )
    proteins_faa.write_text(
        f">{record_id}_1 # 1 # 6 # 1 # "
        "ID=1_1;partial=00;start_type=GTG;genetic_code=11;gc_cont=0.5\n"
        f"{first_protein}\n"
        f">{record_id}_2 # 7 # 12 # -1 # "
        "ID=1_2;partial=00;start_type=ATG;genetic_code=11;gc_cont=0.5\n"
        "M*\n"
        ">__virosync_til\nM"
    )
    attributes = "partial=00;start_type=ATG;genetic_code=11;gc_cont=0.5"
    genes_gff.write_text(
        "##gff-version 3\n"
        f'# Sequence Data: seqnum=1;seqlen=24;seqhdr="{record_id}"\n'
        '# Model Data: version=Prodigal.v2.11.0-gv;transl_table=11;uses_sd=1\n'
        f"{record_id}\tProdigal\tCDS\t1\t6\t.\t+\t0\t"
        f"ID=1_1;partial=00;start_type=GTG;genetic_code=11;gc_cont=0.5;\n"
        f"{record_id}\tProdigal\tCDS\t7\t12\t.\t-\t0\t"
        f"ID=1_2;{attributes};\n"
        f"{record_id}\tProdigal\tCDS\t13\t18\t.\t-\t0\t"
        "ID=1_3;partial=00;start_type=GTG;genetic_code=11;gc_cont=0.5;\n"
        f"{record_id}\tProdigal\tCDS\t19\t24\t.\t+\t0\t"
        f"ID=1_4;{attributes};\n"
    )
    return input_fasta, proteins_faa, genes_gff, record_id


def test_cleanup_abort_reconstructs_owned_suffix_from_complete_gff(
    tmp_path: Path,
) -> None:
    input_fasta, proteins_faa, genes_gff, record_id = (
        _write_reconstruction_fixture(tmp_path)
    )
    validation = prodigal._validate_tiled_prodigal_output(
        input_fasta,
        proteins_faa,
        genes_gff,
        tile_cores={record_id: (0, 18)},
        allow_cleanup_recovery=True,
    )

    assert validation.reconstructed_coordinates == (
        (record_id, 12, 18, "-"),
    )
    assert validation.discarded_coordinates == (
        (record_id, 18, 24, "+"),
    )
    assert (
        prodigal._repair_cleanup_abort_proteins(
            input_fasta,
            proteins_faa,
            genes_gff,
            validation,
        )
        == 2
    )
    records = list(SeqIO.parse(proteins_faa, "fasta"))
    assert [record.id for record in records] == [
        f"{record_id}_1",
        f"{record_id}_2",
        f"{record_id}_3",
    ]
    assert [str(record.seq) for record in records] == ["M*", "M*", "M*"]
    assert all(
        prodigal.parse_prodigal_header(record.description, record.id)[0]
        == record_id
        for record in records
    )


def test_cleanup_abort_reconstruction_rejects_survivor_mismatch(
    tmp_path: Path,
) -> None:
    input_fasta, proteins_faa, genes_gff, record_id = (
        _write_reconstruction_fixture(tmp_path, first_protein="A*")
    )
    validation = prodigal._validate_tiled_prodigal_output(
        input_fasta,
        proteins_faa,
        genes_gff,
        tile_cores={record_id: (0, 18)},
        allow_cleanup_recovery=True,
    )

    with pytest.raises(RuntimeError, match="does not round-trip from GFF"):
        prodigal._repair_cleanup_abort_proteins(
            input_fasta,
            proteins_faa,
            genes_gff,
            validation,
        )


def test_cleanup_abort_reconstruction_rejects_unordered_gff(
    tmp_path: Path,
) -> None:
    input_fasta, proteins_faa, genes_gff, record_id = (
        _write_reconstruction_fixture(tmp_path)
    )
    lines = genes_gff.read_text().splitlines()
    lines[-2], lines[-1] = lines[-1], lines[-2]
    genes_gff.write_text("\n".join(lines) + "\n")

    with pytest.raises(RuntimeError, match="not strictly coordinate ordered"):
        prodigal._validate_tiled_prodigal_output(
            input_fasta,
            proteins_faa,
            genes_gff,
            tile_cores={record_id: (0, 18)},
            allow_cleanup_recovery=True,
        )


def test_cleanup_abort_reconstruction_requires_gff_final_newline(
    tmp_path: Path,
) -> None:
    input_fasta, proteins_faa, genes_gff, record_id = (
        _write_reconstruction_fixture(tmp_path)
    )
    genes_gff.write_text(genes_gff.read_text().rstrip("\n"))

    with pytest.raises(RuntimeError, match="GFF lacks a final newline"):
        prodigal._validate_tiled_prodigal_output(
            input_fasta,
            proteins_faa,
            genes_gff,
            tile_cores={record_id: (0, 18)},
            allow_cleanup_recovery=True,
        )


def test_cleanup_abort_reconstruction_requires_intact_survivor(
    tmp_path: Path,
) -> None:
    input_fasta, proteins_faa, genes_gff, record_id = (
        _write_reconstruction_fixture(tmp_path)
    )
    proteins_faa.write_text(">__virosync_til\nM")
    validation = prodigal._validate_tiled_prodigal_output(
        input_fasta,
        proteins_faa,
        genes_gff,
        tile_cores={record_id: (0, 18)},
        allow_cleanup_recovery=True,
    )

    with pytest.raises(RuntimeError, match="has no intact survivors"):
        prodigal._repair_cleanup_abort_proteins(
            input_fasta,
            proteins_faa,
            genes_gff,
            validation,
        )


def test_cleanup_abort_reconstruction_requires_complete_gff_metadata(
    tmp_path: Path,
) -> None:
    input_fasta, proteins_faa, genes_gff, record_id = (
        _write_reconstruction_fixture(tmp_path)
    )
    genes_gff.write_text(genes_gff.read_text().replace(";gc_cont=0.5", ""))
    validation = prodigal._validate_tiled_prodigal_output(
        input_fasta,
        proteins_faa,
        genes_gff,
        tile_cores={record_id: (0, 18)},
        allow_cleanup_recovery=True,
    )

    with pytest.raises(RuntimeError, match="lacks required attributes: gc_cont"):
        prodigal._repair_cleanup_abort_proteins(
            input_fasta,
            proteins_faa,
            genes_gff,
            validation,
        )


def test_tiled_chunk_reconstructs_owned_suffix_and_audits(
    tmp_path: Path,
    monkeypatch,
) -> None:
    fixture_dir = tmp_path / "fixture"
    fixture_dir.mkdir()
    source_input, source_faa, source_gff, record_id = (
        _write_reconstruction_fixture(fixture_dir)
    )
    work_dir = tmp_path / "temporary"
    work_dir.mkdir()
    chunk_fasta = work_dir / "chunk.fasta"
    chunk_out = work_dir / "chunk.faa"
    chunk_fasta.write_bytes(source_input.read_bytes())
    calls = 0

    def cleanup_abort(cmd, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            return subprocess.CompletedProcess(cmd, -6)
        output = Path(cmd[cmd.index("-a") + 1])
        gff = Path(cmd[cmd.index("-o") + 1])
        output.write_bytes(source_faa.read_bytes())
        gff.write_bytes(source_gff.read_bytes())
        kwargs["stderr"].write("free(): invalid pointer\n")
        return subprocess.CompletedProcess(cmd, -6)

    monkeypatch.setattr(prodigal.subprocess, "run", cleanup_abort)

    assert prodigal._run_prodigal_on_chunk(
        str(chunk_fasta),
        str(chunk_out),
        True,
        {record_id: (0, 18)},
    ) == str(chunk_out)
    assert calls == 2
    assert [record.id for record in SeqIO.parse(chunk_out, "fasta")] == [
        f"{record_id}_1",
        f"{record_id}_2",
        f"{record_id}_3",
    ]
    audit = (
        tmp_path / "accepted_cleanup_aborts" / f"{record_id}.json"
    ).read_text()
    assert '"survivor_check_count": 2' in audit
    assert '"start_0based": 12' in audit
    assert '"start_0based": 18' in audit


def test_tiled_chunk_rejects_cleanup_abort_for_untiled_record(
    tmp_path: Path,
    monkeypatch,
) -> None:
    work_dir = tmp_path / "temporary"
    work_dir.mkdir()
    tile_id = f"{prodigal._TILE_ID_PREFIX}0"
    short_id = "short_scaffold"
    chunk_fasta = work_dir / "chunk.fasta"
    chunk_out = work_dir / "chunk.faa"
    chunk_fasta.write_text(f">{tile_id}\nACGTAC\n>{short_id}\nACGTAC\n")
    calls = 0

    def cleanup_abort(cmd, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            return subprocess.CompletedProcess(cmd, -6)
        input_record = next(SeqIO.parse(cmd[cmd.index("-i") + 1], "fasta"))
        output = Path(cmd[cmd.index("-a") + 1])
        gff = Path(cmd[cmd.index("-o") + 1])
        output.write_text(
            f">{input_record.id}_1 # 1 # 3 # 1 # ID=1_1\nM\n"
        )
        gff.write_text(
            "##gff-version 3\n"
            f"# Sequence Data: seqnum=1;seqlen=6;seqhdr=\"{input_record.id}\"\n"
            f"{input_record.id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
        )
        if input_record.id == short_id:
            kwargs["stderr"].write("free(): invalid pointer\n")
            return subprocess.CompletedProcess(cmd, -6)
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(prodigal.subprocess, "run", cleanup_abort)

    with pytest.raises(RuntimeError, match="untiled record cannot be accepted safely"):
        prodigal._run_prodigal_on_chunk(
            str(chunk_fasta),
            str(chunk_out),
            True,
            {tile_id: (0, 6)},
        )


def test_cleanup_abort_rejects_noncontiguous_gff_loss(tmp_path: Path) -> None:
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    input_fasta = tmp_path / "input.fasta"
    proteins_faa = tmp_path / "proteins.faa"
    genes_gff = tmp_path / "genes.gff"
    input_fasta.write_text(f">{record_id}\n{'A' * 18}\n")
    proteins_faa.write_text(
        f">{record_id}_1 # 1 # 3 # 1 # ID=1_1\nM\n"
        f">{record_id}_3 # 13 # 15 # 1 # ID=1_3\nM\n"
    )
    genes_gff.write_text(
        "##gff-version 3\n"
        f"# Sequence Data: seqnum=1;seqlen=18;seqhdr=\"{record_id}\"\n"
        f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
        f"{record_id}\tProdigal\tCDS\t7\t9\t.\t+\t0\tID=1_2\n"
        f"{record_id}\tProdigal\tCDS\t13\t15\t.\t+\t0\tID=1_3\n"
    )

    with pytest.raises(RuntimeError, match="not a contiguous GFF suffix"):
        prodigal._validate_tiled_prodigal_output(
            input_fasta,
            proteins_faa,
            genes_gff,
            tile_cores={record_id: (0, 6)},
            allow_cleanup_recovery=True,
        )


@pytest.mark.parametrize(
    "truncated_header",
    [
        f">{prodigal._TILE_ID_PREFIX}0_2 # 7 # 1",
        ">__virosync_til",
    ],
)
def test_cleanup_abort_discards_malformed_final_header(
    tmp_path: Path,
    truncated_header: str,
) -> None:
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    input_fasta = tmp_path / "input.fasta"
    proteins_faa = tmp_path / "proteins.faa"
    genes_gff = tmp_path / "genes.gff"
    input_fasta.write_text(f">{record_id}\n{'A' * 12}\n")
    proteins_faa.write_text(
        f">{record_id}_1 # 1 # 3 # 1 # ID=1_1\nM\n"
        f"{truncated_header}\nM"
    )
    genes_gff.write_text(
        "##gff-version 3\n"
        f"# Sequence Data: seqnum=1;seqlen=12;seqhdr=\"{record_id}\"\n"
        f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
        f"{record_id}\tProdigal\tCDS\t7\t9\t.\t+\t0\tID=1_2\n"
    )

    validation = prodigal._validate_tiled_prodigal_output(
        input_fasta,
        proteins_faa,
        genes_gff,
        tile_cores={record_id: (0, 6)},
        allow_cleanup_recovery=True,
    )
    assert validation.discarded_coordinates == ((record_id, 6, 9, "+"),)

    prodigal._remove_discarded_proteins(
        proteins_faa,
        validation.discarded_coordinates,
    )
    assert proteins_faa.read_text() == (
        f">{record_id}_1 # 1 # 3 # 1 # ID=1_1\nM\n"
    )


def test_strict_validation_rejects_no_delimiter_final_header(
    tmp_path: Path,
) -> None:
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    input_fasta = tmp_path / "input.fasta"
    proteins_faa = tmp_path / "proteins.faa"
    genes_gff = tmp_path / "genes.gff"
    input_fasta.write_text(f">{record_id}\n{'A' * 12}\n")
    proteins_faa.write_text(
        f">{record_id}_1 # 1 # 3 # 1 # ID=1_1\nM\n>__virosync_til\nM"
    )
    genes_gff.write_text(
        "##gff-version 3\n"
        f"# Sequence Data: seqnum=1;seqlen=12;seqhdr=\"{record_id}\"\n"
        f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
        f"{record_id}\tProdigal\tCDS\t7\t9\t.\t+\t0\tID=1_2\n"
    )

    with pytest.raises(RuntimeError, match="unparseable Prodigal protein header"):
        prodigal._validate_tiled_prodigal_output(
            input_fasta,
            proteins_faa,
            genes_gff,
        )


def test_cleanup_abort_rejects_no_delimiter_nonfinal_header(
    tmp_path: Path,
) -> None:
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    input_fasta = tmp_path / "input.fasta"
    proteins_faa = tmp_path / "proteins.faa"
    genes_gff = tmp_path / "genes.gff"
    input_fasta.write_text(f">{record_id}\n{'A' * 12}\n")
    proteins_faa.write_text(
        ">__virosync_til\nM\n"
        f">{record_id}_2 # 7 # 9 # 1 # ID=1_2\nM\n"
    )
    genes_gff.write_text(
        "##gff-version 3\n"
        f"# Sequence Data: seqnum=1;seqlen=12;seqhdr=\"{record_id}\"\n"
        f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
        f"{record_id}\tProdigal\tCDS\t7\t9\t.\t+\t0\tID=1_2\n"
    )

    with pytest.raises(RuntimeError, match="unparseable Prodigal protein header"):
        prodigal._validate_tiled_prodigal_output(
            input_fasta,
            proteins_faa,
            genes_gff,
            tile_cores={record_id: (0, 6)},
            allow_cleanup_recovery=True,
        )


def test_cleanup_abort_rejects_malformed_final_without_missing_suffix(
    tmp_path: Path,
) -> None:
    record_id = f"{prodigal._TILE_ID_PREFIX}0"
    input_fasta = tmp_path / "input.fasta"
    proteins_faa = tmp_path / "proteins.faa"
    genes_gff = tmp_path / "genes.gff"
    input_fasta.write_text(f">{record_id}\n{'A' * 12}\n")
    proteins_faa.write_text(
        f">{record_id}_1 # 1 # 3 # 1 # ID=1_1\nM\n>__virosync_til\nM"
    )
    genes_gff.write_text(
        "##gff-version 3\n"
        f"# Sequence Data: seqnum=1;seqlen=12;seqhdr=\"{record_id}\"\n"
        f"{record_id}\tProdigal\tCDS\t1\t3\t.\t+\t0\tID=1_1\n"
    )

    with pytest.raises(
        RuntimeError,
        match="malformed final protein header has no matching GFF suffix",
    ):
        prodigal._validate_tiled_prodigal_output(
            input_fasta,
            proteins_faa,
            genes_gff,
            tile_cores={record_id: (0, 6)},
            allow_cleanup_recovery=True,
        )


@pytest.mark.parametrize(
    ("returncode", "stderr", "expected"),
    [
        (-6, "free(): invalid pointer\n", True),
        (1, "free(): invalid pointer\n", False),
        (-6, "unrelated error\n", False),
    ],
)
def test_known_cleanup_failure_is_narrow(
    tmp_path: Path,
    returncode: int,
    stderr: str,
    expected: bool,
) -> None:
    stderr_path = tmp_path / "prodigal.stderr"
    stderr_path.write_text(stderr)

    assert prodigal._known_cleanup_failure(returncode, stderr_path) is expected


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
        output.write_text(
            f">{prodigal._TILE_ID_PREFIX}0_1 # 1 # 3 # 1 # ID=1_1\nM\n"
        )
        gff.write_text(
            "##gff-version 3\n"
            f"# Sequence Data: seqnum=1;seqlen=4;seqhdr=\"{prodigal._TILE_ID_PREFIX}0\"\n"
            f"{prodigal._TILE_ID_PREFIX}0\tProdigal\tCDS\t2\t4\t.\t+\t0\tID=1_1\n"
        )
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(prodigal.subprocess, "run", mismatched_output)

    with pytest.raises(RuntimeError, match="coordinates differ"):
        prodigal._run_prodigal_on_chunk(str(chunk_fasta), str(chunk_out))

    failure_dir = tmp_path / "prodigal_failures" / "chunk_2"
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
        output.write_text(
            f">{prodigal._TILE_ID_PREFIX}0_1 # 1 # 6 # 1 # ID=1_1\nM\n"
        )
        gff.write_text(
            "##gff-version 3\n"
            f"# Sequence Data: seqnum=1;seqlen=6;seqhdr=\"{prodigal._TILE_ID_PREFIX}0\"\n"
            f"{prodigal._TILE_ID_PREFIX}0\tProdigal\tCDS\t1\t6\t.\t+\t0\tID=1_1\n"
        )
        return subprocess.CompletedProcess(cmd, 0)

    monkeypatch.setattr(prodigal.subprocess, "run", truncated_output)

    with pytest.raises(RuntimeError, match="protein length does not match"):
        prodigal._run_prodigal_on_chunk(str(chunk_fasta), str(chunk_out))


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
        Path(chunk_out).write_text(
            ">rogue_1 # 1 # 3 # 1 # ID=1_1;partial=00\nMKK\n"
        )
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


def test_tiled_genome_validates_every_chunk(
    tmp_path: Path,
    monkeypatch,
) -> None:
    genome = tmp_path / "genome.fasta"
    genome.write_text(">long\nACGTACGT\n>short\nACGT\n")
    monkeypatch.setattr(prodigal, "_LONG_SCAFFOLD_BP", 4)
    monkeypatch.setattr(prodigal, "_TILE_CORE_BP", 4)
    monkeypatch.setattr(prodigal, "_TILE_OVERLAP_BP", 1)
    observed: list[tuple[bool, dict[str, tuple[int, int]]]] = []

    def fake_prodigal(
        _chunk_fasta: str,
        chunk_out: str,
        has_tiles: bool,
        tile_cores: dict[str, tuple[int, int]],
    ) -> str:
        observed.append((has_tiles, tile_cores))
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
    assert all(has_tiles for has_tiles, _ in observed)
    assert any(not cores for _, cores in observed)


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
