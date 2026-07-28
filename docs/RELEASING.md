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

# 3. Commit and tag
git commit -am "release: v0.4.0"
git tag -a v0.4.0 -m "v0.4.0"
git push && git push --tags
```

`git-cliff` config lives in `cliff.toml`. It is a **drafting aid only** —
past release sections in `CHANGELOG.md` are written by hand and never
regenerated. `git describe` (e.g. `v0.3.0-12-gabcdef0`) gives the exact
commit offset past a tag for debugging.
