# Publishing

A release is a version bump on `main`:

```bash
# set the new version in Cargo.toml, then carry it into Cargo.lock
cargo build
git commit -am "1.0.3: what changed"
git push origin main
```

Nobody pushes a tag. The [Release workflow](.github/workflows/release.yml) sees
that `v1.0.3` has no GitHub release yet and

1. runs the tests and the byte-parity check against the Python tool,
2. builds the Linux x86_64 and macOS arm64 binaries,
3. publishes the crate to crates.io,
4. creates the `v1.0.3` tag and the GitHub release, with the binaries attached.

`gh run watch` follows the run. The version then shows on
[crates.io](https://crates.io/crates/gfa-to-pairwise-paf) and under
[releases](https://github.com/cmdcolin/gfa-to-pairwise-paf-rs/releases), whose
generated notes the release page lets you edit.

Only `Cargo.toml` and `Cargo.lock` carry the crate's version. The `v1.0.0` in
the README and in `push.yml` is the Python tool's version that CI checks parity
against.

## When a run fails

The tag and the release come last, so a failed run leaves neither. Fix the cause
on `main` and start the workflow again:

```bash
gh workflow run release.yml
```

The new run skips what already happened, so it never uploads a version that is
already on crates.io. If only the crates.io step fails, `cargo publish --locked`
from the release commit does that step by hand, and a new run does the rest.

crates.io never replaces a published version. A crate that went out broken
needs `cargo yank` and a new version.

## Rehearsing

From any branch other than `main`, the workflow does everything except publish:
it tests, builds, logs in to crates.io and runs `cargo publish --dry-run`.

```bash
git push origin HEAD:release-rehearsal
gh workflow run release.yml --ref release-rehearsal
gh run watch
git push origin --delete release-rehearsal
```

## crates.io login

The repository holds no crates.io token. The crate's Trusted Publishing setting
on crates.io names this repository and `release.yml` as its publisher, and each
run trades its GitHub identity for a short-lived upload token. Renaming the
workflow file or the repository breaks that match until the setting follows.
