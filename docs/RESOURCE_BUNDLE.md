# ViroSync resource bundle

The current core resource version is `v1.0.7`, supplied as
`resources_v1_0_7_runtime.tar.gz`. It includes marker HMMs, Pfam profiles,
marker and broad reference-proteome DIAMOND databases, taxonomy labels, and
marker annotations used for EVE detection, boundary refinement, and
classification.

Setup downloads about 6 GB and installs about 13 GB. The archive and installed
data can use about 19 GB at the same time. Allow extra space for the Pixi
environment, results, filesystem overhead, optional resources, and retained
older versions.

Install the core resources:

```bash
pixi run setup-virosync-resources
```

To use another location, set the database root before setup and keep it set for
later runs:

```bash
export VIROSYNC_DB_ROOT=/data/virosync-db
pixi run setup-virosync-resources
```

Verify the complete installed bundle:

```bash
pixi run virosync orchestrate resources verify \
  --config config/orchestration.yaml \
  --full
```

`--full` hashes each resource file and checks the DIAMOND databases. Run it
after setup, after a storage error, or after copying the resource tree.

Do not run setup while ViroSync is using the databases.

## Gene-function profiles (bundle v1.2.0 and later)

A runtime bundle of v1.2.0 or later carries two more files next to
`combined.hmm`:

| File | Role in the manifest | Content |
| --- | --- | --- |
| `models/gene_function.hmm` | `gene_function_models` | HMMER3 text profiles from Pfam and VOGDB |
| `models/gene_function_labels.tsv` | `gene_function_labels` | One row per profile: `profile`, `source`, `label`, `module`, `max_evalue`, `description` |

The manifest counts both files (`gene_function_models` and
`gene_function_labels`). The counts must be equal, and the profile names in
the two files must match. `scripts/build_resource_bundle.py` parses both files
before it writes a manifest. It requires at least one profile, a valid HMMER3
library, a labels table that passes the runtime reader, and the same profile
names once each in both files. Installation checks hashes and counts only.
Bundles before v1.2.0 keep their file lists. A run on such a bundle has no
gene-function files, so it logs one INFO line and writes no
`gene_functions.tsv`. The [output reference](METHODS.md#gene-function-labels)
describes how a run uses the profiles.

The profiles label genes. They never change detection, boundaries, scores,
tiers or any existing output column. The label vocabulary follows Fig. 4 of
Bellas and Sommaruga (2026, bioRxiv, doi:10.64898/2026.06.19.733378). Each
label belongs to one module:

| Module | Labels |
| --- | --- |
| capsid | MCP, mCP, ATPase, Protease, Portal, Terminase, Triplex |
| replication | PolB, pPolB, PolA, SF1, SF3, AEP |
| integration | YR, RVE, Int |
| mriya | VLTF2, VLTF3, HUH |
| other | MTase, Tet, SLATT, Lipase |
| hallmark | RNAPL, TFIIB |

A profile enters the set in one of two ways:

- Pfam: the model NAME is one of the Pfam families of the vocabulary. A hit
  counts at an i-Evalue of 1e-5 or lower.
- VOGDB: the consensus description of the group, without its `REFSEQ` prefix,
  names a gene of the vocabulary, in any virus lineage. A hit counts at a
  full-sequence E-value of 1e-9 or lower. The label rules exclude look-alikes
  such as the PCNA "DNA polymerase sliding clamp", the lambda "DNA-packaging
  protein FI" and RecT recombinases.

Rebuild the two files from a Pfam-A and a VOGDB release:

```bash
pixi run python scripts/build_gene_function_profiles.py \
  --pfam Pfam-A.hmm --pfam-release 38.2 \
  --vog vogdb.hmm --vog-annotations vog.annotations.tsv \
  --vog-release <release> --output-dir <dir>
```

The script streams both HMM files, so it does not load the 4 GB VOGDB file
into memory. Copy the two outputs into `models/` of the resource tree before
`scripts/build_resource_bundle.py --split --version v1.2.0`.

### Credits

- VOGDB (Virus Orthologous Groups), CC BY 4.0: Trgovec-Greif et al. 2024,
  Viruses 16:1191, doi:10.3390/v16081191. The profiles in `gene_function.hmm`
  keep their VOGDB identifiers.
- Pfam, CC0: the Pfam profiles in `gene_function.hmm` keep their Pfam names.

The bundle README repeats both credits.
