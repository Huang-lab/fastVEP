# Releasing fastVEP

fastVEP follows [Semantic Versioning](https://semver.org). While the on-disk
formats (`.osa2`, cache, CSQ layout) are still evolving, we stay in the `0.x`
range: **`fix` → patch**, **`feat` → minor**, and a **breaking change → minor**
(reserving the first `1.0.0` for when the formats and CLI are stable enough to
promise backward compatibility).

The version is set once in `[workspace.package]` of the root `Cargo.toml`; every
crate inherits it via `version.workspace = true`.

## How the version is computed

We do **not** bump per commit. The next version is derived from the
[Conventional Commit](https://www.conventionalcommits.org) prefixes since the
last `vX.Y.Z` tag:

| Commit prefix                     | Bump   |
| --------------------------------- | ------ |
| `fix:` / `perf:`                  | patch  |
| `feat:`                           | minor  |
| `feat!:` / `BREAKING CHANGE:`     | minor (pre-1.0) / major (post-1.0) |
| `docs:` `chore:` `test:` `ci:` …  | none   |

`git-cliff --bumped-version` computes this for you.

> **Commit hygiene:** use real conventional prefixes with a type, e.g.
> `fix(acmg): …` — not `acmg: …`. Bare `scope:` messages can't be classified as
> feat vs fix and land in a generic "Changed" bucket in the draft.

## Cutting a release

```sh
# 1. Draft + bump (auto-computes the version, or pass one explicitly)
scripts/release.sh          # or: scripts/release.sh 0.4.0

# 2. Edit CHANGELOG.md — promote [Unreleased] to the new version heading and
#    polish the prose. The script prints a git-cliff draft to seed this; we keep
#    the changelog hand-curated (richer than raw commit subjects), so treat the
#    draft as a starting point, not the final text.

# 3. Commit and tag - but rehearse the build before the tag (next section)
git commit -am "release: v0.4.0"
git tag -a v0.4.0 -m "v0.4.0"
git push && git push --tags
```

## Before the tag: rehearse the build

The tag starts the release workflow, and the workflow publishes only if every build job succeeds.
A runner that never starts, or an archive that does not run, therefore stops the release at the moment it is cut.
Find that out first, on the commit you are about to tag:

```sh
gh workflow run release.yml --ref <branch>
gh run watch
```

A run started this way builds the archive for every platform, checks each one with `scripts/check-release-tarball.sh`, and publishes nothing.
Its archives are kept as workflow artifacts (`gh run download <run-id>`), and are named after the commit, not a tag.

The check unpacks the archive, verifies its `.sha256`, starts both binaries, and annotates the README's example.
On a tag it also requires `fastvep --version` to report the tag's version, so a tag pushed before the version was bumped fails there, and does not publish an archive that names the wrong release.

A Linux archive is run in Rocky Linux 9, whose glibc is 2.34.
That is what the v0.4.0 x86_64 binary needed, and the check is what keeps it so: a runner image whose toolchain links a newer glibc would otherwise build, pack and publish a binary that fails to start on RHEL 9.
If the check fails that way, build on an older runner image, or raise the floor on purpose and change the README's note on it with it.

## Release candidates

A tag `vX.Y.Z-rc.N` publishes a release candidate: the same archives, as a GitHub pre-release.
GitHub does not make a pre-release the "Latest" release, so `releases/latest` stays on the last final version.
Use one to let someone downstream test a build before the final one is cut.

```sh
scripts/release.sh 0.5.0-rc.1     # sets the version in Cargo.toml to 0.5.0-rc.1
git commit -am "release: v0.5.0-rc.1"
git tag -a v0.5.0-rc.1 -m "v0.5.0-rc.1"
git push && git push --tags
```

The version in `Cargo.toml` has to be the candidate's, because the check requires the binary to report the tag's version.
The release body is the one-line fallback: `CHANGELOG.md` gets its section when the final version is cut.
Skip the conda recipe below, which follows final releases only.
`cliff.toml` ignores candidate tags, so the draft changelog and the computed version of the final release still count from the last final one.
The final release sets the version back to `X.Y.Z`.

## After the tag: refresh the conda recipe

`conda/recipe/meta.yaml` pins the release tarball by checksum, so it must be
bumped once the tag is pushed and GitHub has generated the archive:

```sh
VERSION=0.4.0
curl -sL -o /tmp/fastvep.tar.gz \
  "https://github.com/Huang-lab/fastVEP/archive/refs/tags/v${VERSION}.tar.gz"
shasum -a 256 /tmp/fastvep.tar.gz
# → update `version` and `sha256` in conda/recipe/meta.yaml, reset `build.number` to 0
```

Once fastVEP is accepted into bioconda, the same two values also need updating
in `bioconda-recipes/recipes/fastvep/meta.yaml` (or left to bioconda's
autobump bot, which watches GitHub releases and opens the PR for you).

## Notes

`git-cliff` config lives in `cliff.toml`. It is a **drafting aid only** —
past release sections in `CHANGELOG.md` are written by hand and never
regenerated. `git describe` (e.g. `v0.3.0-12-gabcdef0`) gives the exact
commit offset past a tag for debugging.
