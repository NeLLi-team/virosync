# ViroSync Resource Bundle

This document records the core resource files expected by the public setup path
and the validation checks used before publishing a new bundle.

## Artifact layouts

The published v1.0.7 archive uses the schema-v2 runtime layout. It extracts into
`resources/virosync` and contains nine payloads:

- `DB_VERSION`
- `DATABASE_README.txt`
- `models/combined.hmm`
- `models/pfam_virosync_screening.hmm`
- `models/model_annotations_with_interpro.tsv`
- `models/og_marker_name_map.tsv`
- `marker/marker.dmnd`
- `genomes/combined_proteome.dmnd`
- `taxonomy/labels.tsv`

The archive also contains `RESOURCE_MANIFEST.json`, which authenticates the
nine payloads. Each entry records a canonical relative path, positive byte
size, SHA-256 digest, and file role. The manifest also records HMM, annotation,
marker, DIAMOND, and taxonomy counts. It does not list itself; code and both
shipped configs pin its SHA-256. Setup writes `DB_METADATA.json` after
validation as an install receipt, not as an archive member or trust root.

The builder can emit a bound source/repair companion with the five files
omitted from the runtime artifact:

- `models/combined.hmm.h3f`
- `models/combined.hmm.h3i`
- `models/combined.hmm.h3m`
- `models/combined.hmm.h3p`
- `marker/marker.faa`

Each artifact has its own `RESOURCE_MANIFEST.json`. The source/repair manifest
records the SHA-256 of the matching runtime manifest, so artifacts from
different builds cannot be paired. Both manifests carry the semantic counts
for the complete resource set. The runtime reads the raw HMM library through
pyhmmer and uses the prebuilt marker DIAMOND database. It does not read the
HMMER indices or marker FASTA. The v1.0.7 public setup publishes only the
runtime archive. The source/repair artifact supports inspection and can
rebuild the marker database. It cannot rebuild
`genomes/combined_proteome.dmnd`, so it is not a complete raw-source backup.

Core setup accepts schema-v1 legacy bundles and schema-v2 runtime bundles. A
schema-v1 bundle does not contain the Pfam screening HMM. ViroSync keeps its
HMM hits unchanged and logs a warning if it encounters ambiguous proteins.
Schema-v2 runtime manifests require the Pfam HMM and include `pfam_models` in
their semantic counts. Setup rejects a source/repair artifact at both fresh
install and repeat setup.

The [frameshift screening guide](FRAMESHIFT_SCREENING.md#run-the-shipped-example)
defines the resource-dependent output contract for the current schema-v2
bundle and the supported legacy schema-v1 bundle.

The split builder finishes and hashes both archives before replacing an output
path. The two replacements are separate filesystem operations, so they are not
atomic as a pair. Build each release under new versioned names and verify the
returned digests. Publish only the artifact selected for public setup. Do not
overwrite a published artifact in place.

## Release identity and installation

The v1.0.7 runtime archive has these identities:

- Archive: `resources_v1_0_7_runtime.tar.gz`
- Archive size: `5877324818` bytes
- Archive SHA-256: `57daed0b39bf2bc4c4f84ec3b612c6034a3d26ea38e7ec5fba4f4469da36e9a2`
- `RESOURCE_MANIFEST.json` SHA-256: `f3aeed77045f4728207c6997f5986ed155056e2b4b2a297574d57686982a18b3`

The exact 2,204-byte manifest is tracked at
`release-manifests/resources_v1_0_7/RESOURCE_MANIFEST.json`. The production
guard parses this file with the runtime manifest loader and compares its digest
with both shipped configurations and the database source record.

Setup authenticates the archive before extraction, preflights every tar member,
extracts into a same-filesystem sibling stage, validates the complete tree, and
only then activates it. The stable path is a relative symlink to a read-only
directory named from the resource version and manifest digest. A sibling lock
serializes installers. The first migration from a real directory records a
fsynced recovery journal and retains the prior directory. Payload data,
read-only modes, candidate directories, and the parent directory are fsynced
before the stable pointer changes. Setup always performs full semantic
validation before activation.

Normal pipeline startup checks the install receipt before hashing payload bytes.
It first resolves the configured source against the pinned version, archive
SHA-256, and manifest SHA-256, recovers any pending install journal, and
requires the stable path to select the exact relative versioned candidate. The
receipt must match a closed, read-only inventory and the device, inode, size,
modification time, and change time of every authenticated payload. A current
receipt returns without scanning or fsyncing the large payloads. A missing or
stale receipt triggers full manifest validation and refreshes the receipt only
after validation passes. The resolved versioned candidate, rather than the
mutable stable pointer, is passed to the run.

The directory containing the stable path, lock, journal, and versioned trees
must not be group- or world-writable. Do not run setup while a pipeline is
reading the active resources; the first migration from a real directory uses
two renames and briefly has no stable-path entry.

The installer does not prune retained version directories. Budget space for the
active release and at least one prior release. A hard process kill can also
leave a hidden `.virosync.stage-*` directory; it is never activated. An
unreadable recovery journal stops later setup attempts and remains in place for
operator inspection. A kill while preparing a recovery pointer may leave an
inert hidden `.virosync.recover-*` symlink; it does not replace the active
pointer.

Fast verification checks the pinned manifest, version, exact path set, regular
file types, and sizes without scanning the large payloads:

```bash
pixi run virosync orchestrate resources verify \
  --config config/orchestration.yaml
```

Full verification hashes every payload in the selected layout, checks the
available text-derived semantic counts and HMM-to-annotation identity, and runs
`diamond dbinfo` on each DIAMOND database present:

```bash
pixi run virosync orchestrate resources verify \
  --config config/orchestration.yaml --full
```

The stat receipt detects ordinary edits and replacement files but is not a
substitute for periodic full verification against silent storage corruption.
Scheduled and release smoke runs use `--full`.

## Database Contents

The public core bundle provides prebuilt search data for ViroSync's HMM-gated
GEVE workflow. It is not a raw-source database build. Version 1.0.7 retains the
v1.0.6 biological content and adds the authenticated Pfam screen in a
runtime-only package.

The bundle contains three linked biological resources:

- `models/combined.hmm`: the Phase 1 marker HMM library. In the v1.0.7
  bundle this contains 1,053 profiles, including existing gvclass/geNomad-style
  viral markers, Mirusviricota and Preplasmiviricota (PPV: virophage + PLV) markers,
  mitochondrial control markers, VS-renamed OG-backed markers, CRESS/ssDNA Rep/capsid
  markers, and capscan PLV major-capsid-protein markers (Bellas & Sommaruga 2026).
  `models/model_annotations_with_interpro.tsv` provides one annotation row per
  HMM profile, and `models/og_marker_name_map.tsv` maps VS marker IDs back to
  original source names where applicable. Current code requires
  `models/combined.hmm`; older `models/combined_ga.hmm` bundles are historical
  resources and are no longer accepted by the setup/preflight checks.
- `marker/marker.dmnd`: the small Tier 1 marker validation database used only
  for proteins that already hit marker HMMs. It contains 3,864,300 protein
  sequences. The bound source manifest records the same marker FASTA count,
  but the public runtime archive does not contain that FASTA.
- `genomes/combined_proteome.dmnd`: the large Tier 2 gene-taxonomy database
  used during Phase 2 boundary refinement for seed, interior, flanking, and
  sampled control genes. The current public lineage contains 46,765,751 protein
  sequences and 8,621,175,597 amino-acid letters. It is a broad host and
  virus/context proteome collection; the dominant prefixes in the current
  extracted FASTA are:

| Prefix | Role in taxonomy context | Protein sequences |
|---|---|---:|
| `EUK` | eukaryotic host/background proteomes | 42,369,796 |
| `BAC` | bacterial background proteomes | 2,051,841 |
| `ARC` | archaeal background proteomes | 167,784 |
| `NCLDV` | nucleocytoplasmic large DNA viruses | 1,331,717 |
| `MIRUS` | mirusvirus references | 177,872 |
| `PPV` | Preplasmiviricota references (virophage + PLV; lineage class rank distinguishes Virophaviricetes from Polintoviricetes/Aquintoviricetes) | 206,313 |
| `PHAGE` | bacteriophage references | 80,325 |
| `PLASTID` | plastid references | 160,715 |
| `MITO` | mitochondrial references | 219,388 |

Schema-v2 runtime bundles add `models/pfam_virosync_screening.hmm`. The
screening HMM prepared for v1.0.7 contains 937 Pfam 38.0 profiles, each with a
GA cutoff. The split builder rejects a profile without GA and requires
`model_name`, `pfam_signature`, and `source_scope` in the model annotation
table. The manifest authenticates the HMM and records its model count.

`taxonomy/labels.tsv` maps genome/source IDs to pipe-delimited lineage strings
used by marker validation, host-signature modeling, and gene-level taxonomy
trimming. Version 1.0.7 has 190,317 taxonomy-label rows.

Do not treat older private paths such as
`/path/to/legacy-db/virosync/genomes/combined_proteome.dmnd` as the
public bundle unless their `DB_METADATA.json` says so. That historical v1.0.0
database has 77,388,186 sequences and was used by earlier local/release
configurations; the current public lineage starts from the v1.0.2
46.8M-sequence combined proteome and refreshes the marker resources.

## May 2026 VS791 Refresh

The local May 2026 resource refresh starts from the public `v1.0.2` bundle and
adds the relaxed OG-to-VS marker migration:

- Resource bundle version: `v1.0.3`
- Renamed existing raw OG HMM models: 576 (`VS000001`-`VS000576`)
- Added OG-backed marker HMM models: 215 (`VS000577`-`VS000791`)
- Final HMM model count: 1,004
- Added pre-extracted marker proteins: 13,770
- Final `marker/marker.faa` record count: 3,817,304
- Final taxonomy label count: 167,836

Historical public `v1.0.3` release-check values:

- HMM models: `1004`
- `marker.faa` records: `3817304`
- marker annotation table lines: `1005`
- OG marker map lines: `792`
- taxonomy label lines: `167836`

The marker annotation table is `models/model_annotations_with_interpro.tsv`.
It must contain one row per HMM profile name in `models/combined.hmm`; in the
current v1.0.7 bundle that is 1,053 marker annotation rows plus the header.
`models/og_marker_name_map.tsv` maps each VS marker ID to its original OG model
where applicable.

## June 2026 CRESS v1.0.4 Refresh

The public June 2026 resource refresh starts from the public `v1.0.3` bundle and
adds CRESS/ssDNA Rep and capsid models from Oliver's HMM profiles:

- Resource bundle version: `v1.0.4`
- Candidate CRESS HMM models screened: 88
- Models passing the relaxed rule before redundancy consolidation: 21
- Hit-set Jaccard consolidation threshold: 0.75
- Added representative CRESS HMM models: 17 (`VS000792`-`VS000808`)
- Added pre-extracted marker proteins: 1,148
- Final HMM model count: 1,021
- Final `marker/marker.faa` record count: 3,818,452
- Final taxonomy label count: 167,836
- HMM library filename changed from `combined_ga.hmm` to `combined.hmm`

Detailed selected-model provenance is in `LOCAL_CRESS_SSDNA_INTEGRATION.md` and
`CRESS_V104_SELECTED_MODELS.tsv`.

## Public MetaVR CRESS Marker Enrichment

The public `v1.0.4` resource has additional Cressdnaviricota marker-taxonomy
enrichment layers on top of the v1.0.4 CRESS HMM refresh:

- Source checked first: `/path/to/metavr-mirror`
  (not present locally)
- Official FAA checked: `https://www.meta-virome.org/Data/Downloads/IMGVR5_UViG.faa.gz`
  (direct `curl` access was blocked by Cloudflare challenge)
- Local source used:
  `/path/to/gvclass-update-Apr24/IMGVR_all_Sequence_information-high_confidence.tsv.repsApr24.faa`
- Local lookup used:
  `/path/to/gvclass-update-Apr24/IMGVR_all_Sequence_information-high_confidence.tsv.repsApr24.lookup`
- CRESS HMMs searched: `VS000792`-`VS000808`
- Total local IMGVR representative proteins searched: 10,405
- HMM hit rows: 45
- Cressdnaviricota hit rows: 36
- Added deduplicated marker proteins: 16
- Added taxonomy-label rows: 10
- Intermediate `marker/marker.faa` record count: 3,818,468
- Intermediate `taxonomy/labels.tsv` line count: 167,846

Added records are rewritten from their local `PHAGE__IMGVR...` source IDs to
`CRESS__IMGVR...` IDs before appending.

The full MetaVR source is available locally at `/path/to/metavr` and was used
for the final enrichment:

- Full source FAA:
  `/path/to/metavr/data/source/viralproteins.faa`
- Full source metadata:
  `/path/to/metavr/data/source/metadata.tsv`
- Cressdnaviricota export: 31,961 UViGs, 136,449 proteins, 27,845,734 aa
- CRESS HMMs searched: `VS000792`-`VS000808`
- HMM hit rows: 83,649
- Unique Cressdnaviricota target proteins: 43,380
- Added deduplicated marker proteins: 34,237
- Skipped duplicate or missing records: 9,143
- Added taxonomy-label rows: 22,471
- Final v1.0.4 `marker/marker.faa` record count: 3,852,705
- Final v1.0.4 `marker/marker.dmnd` sequence count: 3,852,705
- Final installed `taxonomy/labels.tsv` line count: 190,317

Phase 1 marker validation now accepts a single `NCLDV__`, `MIRUS__`, `PPV__`
(Preplasmiviricota; legacy `PLV__`/`VP__` accepted transitionally), or `CRESS__`
Diamond top-10 hit when `pident >= 25.0`. `PHAGE__` does
not satisfy this marker-validation rule.

## Release Checks

Run these from the repository root against the installed resource directory:

```bash
rg -c '^NAME\s+' resources/virosync/models/combined.hmm
rg -n '^NAME\s+OG[0-9]+$' resources/virosync/models/combined.hmm
rg -c '^NAME\s+' resources/virosync/models/pfam_virosync_screening.hmm
pixi run diamond dbinfo --db resources/virosync/marker/marker.dmnd
pixi run diamond dbinfo --db resources/virosync/genomes/combined_proteome.dmnd
wc -l resources/virosync/models/model_annotations_with_interpro.tsv
wc -l resources/virosync/models/og_marker_name_map.tsv
wc -l resources/virosync/taxonomy/labels.tsv
```

Expected values for the public v1.0.7 runtime resource:

- HMM models: `1053`
- Raw `OG[0-9]+` HMM names: no matches
- Pfam models: `937`
- DIAMOND marker database sequences: `3864300`
- DIAMOND combined proteome sequences: `46765751`
- marker annotation table lines: `1054`
- OG marker map lines: `809`
- taxonomy label lines: `190317`

Build the bound schema-v2 artifacts from the repository root. The builder treats
`resources/virosync` as read-only input. It stages regenerated indices when
requested, synthesizes `DB_VERSION`, `DATABASE_README.txt`, and the manifest,
and writes a deterministic GNU-format tar stream with fixed gzip and member
metadata:

```bash
pixi run python scripts/build_resource_bundle.py \
  --version v1.0.7 \
  --split \
  --output resources_v1_0_7_runtime.tar.gz \
  --source-output resources_v1_0_7_source.tar.gz \
  --skip-hmmpress \
  --skip-marker-dmnd
```

The input tree must contain `models/pfam_virosync_screening.hmm` and the three
Pfam annotation columns. The runtime archive has ten regular-file members:
nine manifest-listed payloads plus its manifest. The source/repair archive has
six: five payloads plus its manifest. The builder prepares and authenticates
both archives before replacing either output. Publish only the runtime archive.
Manual `tar` packaging is not equivalent because it omits the authenticated
manifest and deterministic metadata contract. The two `--skip-*` flags reuse
the derived files in the input tree; omit them when a release build must
regenerate the HMMER indices and marker DIAMOND database from their sources.
