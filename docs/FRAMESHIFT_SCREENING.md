# Run frameshift-sensitive VS marker screening

ViroSync screens the Phase 0 masked nucleotide assembly directly with the
shipped `VS######` marker HMMs. The Phase 1 screen reports marker homology
across frameshifts or in-frame stop codons that prevent
`prodigal-gv` from producing a protein call.

Raw BATH hits are candidates, not protein calls. ViroSync extracts each
model-conditioned aligned amino-acid domain, searches it against the Tier-1
viral marker database, and allows only confirmed loci to seed regions. A
confirmed domain that overlaps an accepted EVE appears in that EVE's FAA.
The domain does not enter the ordinary proteome, family clustering, or generic
gene-taxonomy searches.

## Install BATH

BATH is not available from ViroSync's configured conda-forge and bioconda
channels, so it cannot be a pixi dependency. Build it with:

```bash
pixi run setup-bath
```

This builds the tested revisions into `.bath/`. Later setup calls return
without a rebuild. Pass `--force` to rebuild:
`bash scripts/install_bath.sh --force`. `pixi.toml` puts `.bath/bin` on `PATH`
for every task. After setup, `pixi run example-frameshift` finds the tools
without a manual export.

Run it once before the frameshift example. It is not a dependency of
`example-frameshift`, which remains an explicit opt-in pinned by
`tests/test_ci_contracts.py::test_frameshift_example_task_is_explicit_opt_in`.

Check both entrypoints:

```bash
pixi run bathconvert -h
pixi run bathsearch -h
```

### Build by hand

`scripts/install_bath.sh` runs the commands in this section. Two commands
differ from the upstream instructions. A plain clone of easel does not contain
the pinned commit. Fetch the commit by SHA before checkout. The repository also
ships only `configure.ac`. Generate `configure` with `pixi exec`. This supplies
`autoconf` for one command without changing the project manifest.

```bash
git clone https://github.com/TravisWheelerLab/BATH.git
git -C BATH checkout 7842ebd58b96591b4b60863ee5c33e49eb79eccc
git clone https://github.com/EddyRivasLab/easel.git BATH/easel
git -C BATH/easel fetch origin 0f4e71832d6ba1e4c65039ba4b4663c546a041fa
git -C BATH/easel checkout 0f4e71832d6ba1e4c65039ba4b4663c546a041fa
cd BATH
pixi exec --spec autoconf --spec m4 --spec perl -- autoconf
./configure --prefix=/absolute/path/to/bath
make -j 16
make install
export PATH=/absolute/path/to/bath/bin:$PATH
```

The ViroSync CLI checks both commands before Phase 0 and stops if either is
missing.
Each BATH command has a one-hour timeout. A nonzero exit or timeout stops the
run.

## Configure the screen

The screen is off by default. Enable it for genomes where degraded viral
markers could be missed by protein prediction:

```yaml
phase1:
  frameshift_screening_enabled: true
```

The CLI flag overrides the run configuration:

```bash
pixi run virosync -i assembly.fna -o virosync_output \
  --config config/orchestration.yaml \
  --frameshift-screening
```

Use `--no-frameshift-screening` to disable a screen enabled in a custom
configuration.

## Run the shipped example

The repository includes three *Trichomonas vaginalis* G3 contigs that exercise
rescue seeding and boundary containment:

| NCBI accession | Length (bp) |
|---|---:|
| `DS113389.1` | 116,627 |
| `DS113495.1` | 92,297 |
| `DS113200.1` | 279,084 |

The 488,008-bp fixture is stored at
`example/frameshift/trichomonas-g3.fna`. The source records belong to
[BioProject PRJNA16084](https://www.ncbi.nlm.nih.gov/bioproject/16084) and the
[GCA_000002825.3 assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000002825.3/).
NCBI states that it places no restrictions on use of its molecular databases,
while noting that submitters may retain rights to submitted records. See the
[NCBI data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/).

After installing BATH and the v1.0.7 resource bundle, run:

```bash
pixi run example-frameshift
```

The command writes `results/example-frameshift/`. Four rescue domains pass
Tier-1 marker validation, and one accepted EVE FAA contains its confirmed
`VS000087` rescue domain. The expected counts depend on the resource schema:

| Resource | Pfam arbitration | Detailed | Canonical |
|---|---|---:|---:|
| Legacy schema-v1 v1.0.6 | No | 6 | 3 |
| Default schema-v2 v1.0.7 | Yes, 937 models | 5 | 2 |

The schema-v2 canonical EVE IDs are:

- `EVE_DS113200.1_129184-151689`
- `EVE_DS113495.1_58468-89417`

Its other detailed EVE IDs are `EVE_DS113495.1_2-471`,
`EVE_DS113389.1_48748-50484`, and `EVE_DS113200.1_26626-27034`.
`EVE_DS113495.1_18305-37386` is absent from both tables because `ResIII`
contradicts every `Pox_A32` marker assignment on its only marker protein.

The fixture checks the rescue and Pfam paths through the full pipeline. It is
not a sensitivity benchmark. BATH's raw event count can vary across repeated
multithreaded runs.

When enabled, ViroSync:

1. Streams the raw `models/combined.hmm` file and selects records whose names
   match `VS[0-9]{6}` exactly.
2. Writes the selected records to a fresh text HMM file. This prevents BATH
   from reading adjacent HMMER pressed files, whose binary model layout is not
   compatible with the tested BATH revision.
3. Converts the VS-only HMM file with `bathconvert`.
4. Runs `bathsearch --fs` against the masked Phase 0 assembly with reporting
   and inclusion E-value thresholds of `1e-5`.
5. Retains only hits that report at least one frameshift or stop codon in the
   normalized result table.
6. Extracts the model-conditioned aligned amino-acid domain from BATH's text
   report, removes gaps, and writes literal stop codons as `X` for protein
   search. The result is not a repaired full-length ORF or CDS.
7. Runs DIAMOND `blastp --sensitive` against the configured Tier-1 marker
   database with `--evalue 1e-5` and ten targets per query.
8. Confirms a candidate only when a validated viral reference aligns at 25%
   or greater identity, the BATH hit covers at least 50% of the VS model, and
   the qualifying DIAMOND HSP covers at least 50% of the candidate domain.
   Overlapping same-strand candidates collapse to the highest-scoring HMM hit.
   HMM-only `validated_novel` candidates do not qualify as rescued markers.
9. Adds confirmed loci to the normal marker seed set before region assembly.
   The normal region and EVE acceptance gates still apply.
10. Keeps rescue-derived and ordinary regions separate during Phase 2. A
    rescue-derived region can pass the marker floor with one marker only when
    the confirmed rescued marker remains inside its refined boundary.
11. Excludes a rescue-only candidate if refinement removes every confirmed
    rescued marker. A rescue-seeded candidate can remain on retained ordinary
    marker support. Phase 3 compares an active rescue branch with each ordinary
    branch it directly overlaps. Rescue replaces them only when it beats every
    overlap at a higher confidence tier, or at the same tier with a higher
    final score.

The resource bundle contains 808 VS profiles:

- `VS000001`-`VS000576`: 576 renamed OG HMMs
- `VS000577`-`VS000791`: 215 added OG-backed HMMs
- `VS000792`-`VS000808`: 17 CRESS HMMs

The separate 72-profile Bellas et al. set is not part of this screen.

## Runtime and reproducibility

On the *Trichomonas vaginalis* G3 assembly, screening took 10 minutes 50
seconds with eight threads. The Slurm step peaked at 2.9 GB of resident
memory, and the screen wrote 51 MB. At 16 threads, model conversion plus
assembly search took 8 minutes 41 seconds. DIAMOND validation of about 430
candidate domains against the 3,864,300-protein Tier-1 database took 20-30
seconds at 16 threads on the tested Dori nodes.

A matched eight-core Phase 1 run on G3 took 683.97 seconds (11:23.97) with
screening disabled and 1,384.50 seconds (23:04.50) with screening enabled.
The screen added 700.53 seconds (11:40.53), making Phase 1 2.02 times as long.
The enabled arm confirmed 97 rescued markers across 18 VS models. They
anchored 85 seed regions, including 28 that did not overlap a disabled-arm
seed. These 28 are Phase 1 EVE candidates, not accepted EVEs; they must still
pass the Phase 2 and Phase 3 gates.

The shipped three-contig fixture also has a matched schema-v1 v1.0.6
full-pipeline comparison at eight cores. Phase 1 took 33.4 seconds with
screening disabled and 137.4
seconds with screening enabled, a 104.0-second increase and a 4.11-fold Phase
1 runtime. Whole-command wall time changed from 4:21.07 to 5:06.46. The enabled
run added one canonical EVE, for three rather than two, and retained one
confirmed rescue domain in that EVE's FAA. The enabled arm ran second with warm
caches, so the whole-command and memory differences are not isolated flag
effects.

BATH's event counts changed across retained eight-thread runs despite
byte-identical models: 427, 428, and 433 event-bearing hits. Viral and coverage
filters retained 98-100 confirmed loci; each run pair shared 97-98 overlapping
loci with the same model. Three further runs with seed 42 and
a 2 Mb target block still returned 428, 438, and 441 candidates. Their 99,
99, and 100 confirmed loci had 97-99 same-model overlaps. A fixed target block
does not remove the limit. Keep the raw BATH and validation tables when
individual loci require audit.

## Read the output

The screen writes files under
`<genome_output>/phase1/frameshift_screening/`:

| File | Contents |
|---|---|
| `frameshift_hits.tsv` | Normalized table of raw event-bearing candidates |
| `frameshift_candidates.faa` | Model-conditioned candidate domains submitted to DIAMOND |
| `validation/diamond_top10.tsv` | Ranked Tier-1 protein matches for candidate domains |
| `validation/validated_marker_hits.tsv` | Validation status and reference evidence for every candidate |
| `confirmed_frameshift_markers.tsv` | Coverage-filtered, deduplicated markers allowed to seed |
| `confirmed_frameshift_proteins.faa` | Confirmed domains available for accepted per-EVE FAA output |
| `bathsearch.tblout` | Raw BATH per-domain table |
| `bathsearch.fstblout` | Raw BATH frameshift-event table |
| `bathsearch.txt` | Full BATH text report |

`frameshift_hits.tsv` always has a header, including when no event-bearing hit
passes the threshold. Its `ali_start` and `ali_end` fields use ViroSync's
0-based, half-open coordinate convention. `strand` records the orientation of
the original BATH alignment. `shifts` and `stops` are BATH event counts, and
`annotation_class` is always `frameshift_rescue_candidate`.

BATH reports a joined domain hit rather than separate translated HSPs, so
ViroSync does not add a second HSP-chaining step. It preserves the raw BATH
tables for review.

`frameshift_hits.tsv`, `confirmed_frameshift_markers.tsv`, and
`confirmed_frameshift_proteins.faa` are authenticated Phase 1 artifacts.
Deleting or changing one invalidates Phase 1 resume. Changing either BATH
executable changes the run fingerprint and invalidates resume from Phase 0.

The detailed Phase 3 TSV records canonical selection in
`canonical_selection_outcome`:

| Value | Meaning |
|---|---|
| `kept` | Accepted without ordinary/rescue overlap arbitration |
| `normal_gate_rejected` | Rejected by the standard EVE acceptance gate |
| `rescue_marker_excluded` | Rescue-only boundary lost its confirmed rescued marker |
| `overlap_selected` | Selected from directly overlapping ordinary and rescue branches |
| `overlap_suppressed_by:<candidate_id>` | Lost direct overlap arbitration to the named candidate |
| `unsupported_no_viral_evidence` | Removed after ANI clustering: neither it nor a marker-bearing relative retained viral evidence |

`total_proteins` counts ordinary predicted proteins. A canonical or detailed
EVE FAA can contain extra rescued domains.

## Interpret the result

Treat a confirmed domain as evidence for a degraded marker locus, not as a
coding-sequence annotation. The fixed `1e-5` BATH cutoff has not been
calibrated as a family-specific gathering threshold. Region assembly and EVE
acceptance provide further context, but confirmed-domain counts are not
full-length gene counts.

The screen uses only the shipped VS marker profiles. It does not run six-frame
translation and does not use the separate Bellas et al. hallmark HMM set.
