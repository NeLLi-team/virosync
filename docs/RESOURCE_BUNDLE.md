# ViroSync resource bundle

Install the core databases from the repository root:

```bash
pixi run setup-virosync-resources
```

The released `v1.0.7` bundle, `resources_v1_0_7_runtime.tar.gz`, downloads
about 6 GB and installs about 13 GB. Setup can use about 19 GB while both are
on disk. Allow extra space for the Pixi environment, results, filesystem
overhead, optional resources, and retained older versions.

The bundle contains marker HMMs, Pfam profiles, marker and broad
reference-proteome DIAMOND databases, taxonomy labels, and marker annotations
for EVE detection, boundary refinement, and classification.

## Choose a database directory

Replace `/data/virosync-db` with your database directory and keep
`VIROSYNC_DB_ROOT` set for later runs:

```bash
export VIROSYNC_DB_ROOT=/data/virosync-db
pixi run setup-virosync-resources
```

Do not run setup while ViroSync is using the databases.

## Verify the installation

Verify the files after setup, a storage error, or copying the database directory:

```bash
pixi run virosync orchestrate resources verify \
  --config config/orchestration.yaml \
  --full
```

`--full` hashes each resource file and checks the DIAMOND databases.
