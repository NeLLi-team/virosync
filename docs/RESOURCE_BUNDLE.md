# ViroSync resource bundle

The current core resource version is `v1.1.0`, supplied as
`resources_v1_1_0_runtime.tar.gz`. It includes marker HMMs, Pfam profiles,
marker and broad reference-proteome DIAMOND databases, taxonomy labels, and
marker annotations used for EVE detection, boundary refinement, and
classification.

The bundle contains 3,867,592 marker proteins, 1,086 marker profiles, and 937 Pfam
screening profiles. The marker profiles include 33
[structure-derived candidates](reference/capsid-models.md). This set comes
from a capsid search and includes candidates with adenain-like protease, tail-tube, or
unresolved protease-like evidence. Only `EVE_CAPSID_NOVEL_0029` has a putative
MCP role.

Marker references contain full proteins with identifiers linked to reference
genomes. The resource build tree stores these sequences in
`marker/marker.faa` and builds `marker/marker.dmnd` from that file. The runtime
archive contains the DIAMOND database; the matching source companion contains
the FASTA. Both encode the same marker reference set.

Update your ViroSync checkout before installing v1.1.0. The resource update
includes code changes for MCP roles and functional categories; the older
1.0.1 release lacks them.

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
