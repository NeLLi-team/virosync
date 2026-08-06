# Release ViroSync 1.0.0

ViroSync 1.0.0 has no GitHub Release, and the repository remains private.
Keep every software version field at `1.0.0` until this checklist passes.

## Release identity

- [x] `NeLLi-team/virosync` exists with `main` as its default branch.
- [x] The repository contains only the software, tests, setup files,
      documentation, release manifests, and shipped example.
- [ ] Merge the CLI/taxonomy v4 and repository-scope changes into `main`.
- [ ] Select the final release commit after all required checks pass.
- [ ] Resolve the existing `v1.0.0` tag. It points to `1b696b0` and does not
      contain the unreleased changes. Moving, deleting, or recreating the tag
      requires explicit approval. Do not publish a GitHub Release until the tag
      identifies the selected release commit.
- [ ] Confirm `1.0.0` in `pyproject.toml`, `pixi.toml`,
      `src/virosync/__init__.py`, the README badge, the production guard, and
      the latest release entry in `CHANGELOG.md`.
- [ ] Keep the current tunnel and hosting configuration unchanged during this
      merge.

## Validate the release commit

Run these commands from a clean checkout of the selected commit:

```bash
pixi install --locked
pixi run lint
pixi run python -m compileall -q src tests scripts
pixi run python -m pytest -q
pixi run check-production-ready
pixi run virosync --version
pixi run virosync orchestrate resources verify --full
pixi run virosync \
  -i example/ \
  -o results/release_example \
  --config config/orchestration.yaml \
  -w 1 \
  --threads-per-worker 16 \
  --clean-run
pixi run python scripts/ci/validate_example_run.py \
  results/release_example \
  --expect-predictions 6 \
  --expect-accepted 1
pixi run python scripts/ci/check_coordinate_outputs.py \
  results/release_example
pixi run example-frameshift
pixi run python scripts/ci/validate_example_run.py \
  results/example-frameshift \
  --expect-predictions 5 \
  --expect-accepted 2
pixi run python scripts/ci/check_coordinate_outputs.py \
  results/example-frameshift
```

The release checklist uses the public schema-v2 v1.0.7 runtime bundle. The
scheduled [`example-smoke` workflow](../.github/workflows/example-smoke.yml)
pins its exact EVE IDs. The
[frameshift screening guide](FRAMESHIFT_SCREENING.md#run-the-shipped-example)
defines the schema-v2 Pfam contract.

Then complete the repository checks:

- [ ] `tests` and `production-guards` pass on the release commit.
- [ ] `example-smoke` passes on the release commit with DB v1.0.7.
- [ ] `pixi run pre-commit run gitleaks --all-files` reports no secret.
- [ ] The tracked tree contains no private credentials, personal data,
      machine-specific paths, generated results, local resources, or agent
      working files.
- [ ] The configured resource URL and all pinned archive and manifest digests
      resolve to the intended v1.0.7 runtime bundle.
- [ ] A clean install reports `ViroSync 1.0.0`.

## Deferred release actions

These actions are outside the current merge and require explicit approval.

- [ ] Push the approved `v1.0.0` tag at the selected release commit.
- [ ] Create the GitHub Release from that tag.
- [ ] Verify the release assets and installed `--version`.
- [ ] Add the final release commit and tag to the manuscript provenance
      in `virosync-bench`.
