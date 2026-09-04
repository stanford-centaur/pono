# Releasing Pono

Every version Pono reports derives from two variables at the top of the
top-level `CMakeLists.txt`:

- `project(pono VERSION <major>.<minor>.<patch>)`, the declared release version.
- `PONO_VERSION_SUFFIX`, a PEP 440 pre-release or development marker appended to
  it, empty for a final release.

Together they form `PONO_RELEASE_VERSION`, which feeds `pono --version`, the
Python wheel's version, and `pono.__version__`. No other file records the
version, so there is nothing else to edit.

`pono --version` also reports the commit the binary was built from, re-derived
on every build:

```
2.1.0.dev0 (v2.0.0-63-g1a2b3c4)   built from a git checkout
2.1.0.dev0 (1a2b3c4)              shallow checkout with no tags, as in CI
2.1.0.dev0                        release tarball, or no git available
```

The declared version always comes first, so anything parsing that output should
read the first whitespace-delimited token.

## Cutting a release

1. Clear the suffix in `CMakeLists.txt`:

   ```cmake
   set(PONO_VERSION_SUFFIX "")
   ```

   Check that `project(pono VERSION ...)` already names the version you intend
   to release. It is raised at the start of each development cycle, not here.

2. Commit, then confirm the version reads as expected. Editing `CMakeLists.txt`
   re-triggers configure on its own:

   ```
   cmake --build build --target pono-bin && ./build/pono --version
   ```

3. Tag the release commit. The tag must start with `v`:

   ```
   git tag -a v2.1.0 -m 'Pono 2.1.0'
   ```

   Prefer an annotated tag (`-a`). A lightweight tag is a plain ref that can be
   force-moved without leaving a tagger, date, or message behind, so nothing
   would record that a release tag had changed.

4. Push the commit and the tag.

5. Open the next development cycle in a follow-up commit: raise
   `project(pono VERSION ...)` and set `PONO_VERSION_SUFFIX` back to `".dev0"`,
   so mainline reports e.g. `2.2.0.dev0`.

## Why the tag needs the `v` prefix

The provenance suffix comes from

```
git describe --tags --always --dirty --match "v[0-9]*"
```

A tag without the `v` does not match, and `--always` then falls back to a bare
commit hash. Nothing fails: `pono --version` just quietly stops naming the
release a build descends from. Note the asymmetry is deliberate, in that the
tag carries the `v` and the declared version does not.

`--match` also stops a competition snapshot tag (`hwmcc*`) from winning on
distance, and `--tags` is what makes lightweight tags visible at all.

## Pre-releases

Set the suffix to a PEP 440 marker and tag to match:

| `PONO_VERSION_SUFFIX` | reported version | tag |
| --- | --- | --- |
| `".dev0"` | `2.1.0.dev0` | untagged mainline |
| `"b1"` | `2.1.0b1` | `v2.1.0b1` |
| `"rc1"` | `2.1.0rc1` | `v2.1.0rc1` |
| `""` | `2.1.0` | `v2.1.0` |

PEP 440 orders these `2.1.0.dev0` < `2.1.0b1` < `2.1.0rc1` < `2.1.0`.

Keep the suffix in PEP 440 spelling rather than something like `-beta.1`. The
same string becomes both the C++ and the Python version, and a version that
PEP 440 rejects breaks the wheel build. `project(VERSION)` accepts only numeric
components, which is why the marker lives in a separate variable.

## What a release does not do

- **Publish to PyPI.** `.github/workflows/wheels.yml` triggers on any tag push
  but is currently disabled, and the `contrib/wheels/` scripts behind it are
  unmaintained: they target Python 3.5 to 3.8 and carry their own hardcoded
  version. Publishing a release to PyPI needs that path replaced first.
- **Tag a container image.** `.github/workflows/docker.yml` publishes only
  `ghcr.io/stanford-centaur/pono:latest`, with no per-version tag.
- **Update a changelog.** There is none.
