# ViroSync resource bundle

The current core resource version is `v1.1.0`, supplied as
`resources_v1_1_0_runtime.tar.gz`. It includes marker HMMs, Pfam profiles,
marker and broad reference-proteome DIAMOND databases, taxonomy labels, and
marker annotations used for EVE detection, boundary refinement, and
classification.

The bundle contains 1,086 marker profiles, including 33
[structure-derived candidates](reference/capsid-models.md), and 937 Pfam
screening profiles. The structure-derived set comes from a capsid search and
includes non-capsid candidates with adenain-like protease, tail-tube, or
unresolved protease-like evidence. Only `EVE_CAPSID_NOVEL_0029` has a putative
MCP role.

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
