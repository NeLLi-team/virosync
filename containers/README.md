# ViroSync Apptainer image

`apptainer/virosync.def` builds an image with ViroSync, its locked Pixi
environment and the corrected Prodigal-GV. The resource bundle (v1.3.0, about
11.6 GB packed) is not in the image; users download it into the `virosync` subfolder of a bound `/db`
folder. The [README](../README.md#apptainer-image) shows how to run it.

## Build

From the repository root:

```bash
apptainer build virosync.sif containers/apptainer/virosync.def
```

On a cluster without root, add `--fakeroot`. The build needs network access
and about 2 GB of temporary space; the image is about 0.6 GB.

## Runscript commands

| Command | Action |
| --- | --- |
| `setup-resources` | Download and install the resource bundle into `$VIROSYNC_DB_ROOT` (default `/db/virosync`) |
| `verify-resources` | Hash every resource file and check the DIAMOND databases |
| `example OUTPUT_DIR` | Run the shipped example genome |
| any other arguments | Passed to the `virosync` command line; a run without `--config` gets the image's configuration |

## Publish

The `container` workflow builds the image, checks the version and the bundled
tools, and pushes it to `oras://ghcr.io/nelli-team/virosync:<version>` and
`:latest`. It publishes on a `v*` tag, or when started by hand with
`publish` set. It then pulls the image back without credentials and compares
its SHA-256 with the built file. A new GHCR package is private until an
organisation admin sets it to public once; the pull check fails until then.
