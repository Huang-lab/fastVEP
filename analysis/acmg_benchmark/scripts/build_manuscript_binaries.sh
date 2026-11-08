#!/usr/bin/env bash
# Build the eight binaries the manuscript's runs are made by, each from a clean export of its source: the three of the family v38 (the
# latest classifier, from the tag of release 0.5.0), the two that Section 3.8 adds, and the three of the family v37.
#
#     build_manuscript_binaries.sh OUT_DIR
#     V38_REF=<commit> build_manuscript_binaries.sh OUT_DIR       # v38 as committed (the default is the tag v0.5.0-rc.1, release 0.5.0)
#     V38_TREE=. build_manuscript_binaries.sh OUT_DIR             # v38 as the checkout holds it, uncommitted files included
#     NAMES="v38 cut38 hook38" build_manuscript_binaries.sh OUT_DIR   # only the three binaries of the family v38 (the default family
#                                                                     # of 11_manuscript_runs.py), the others being left as they are
#     NAMES="v37 cut hook" build_manuscript_binaries.sh OUT_DIR   # only the three binaries of the family v37
#
#   v38   the latest classifier (V38_REF, default the tag v0.5.0-rc.1, release 0.5.0; or the files of the checkout named by V38_TREE)
#   cut38   v38 with replay_validation/combiner_cutoff_env.patch   the point-cutoff runs of the family v38 (FASTVEP_LP_CUT, FASTVEP_P_CUT)
#   hook38  v38 with replay_validation/bp_hooks_env.patch          the BP1 and BP7 runs of the family v38 (FASTVEP_BP1_LOO, ...)
#   base  the commit before the curated-inheritance change (BASE_REF, default 9d84be0, master before it), for Section 3.8
#   alt   v38 with replay_validation/inheritance_alternatives_env.patch   the alternatives of Section 3.8 (FASTVEP_ALT)
#   v37   the classifier at the close of the review rounds (V37_REF, default 5d8004a, the merge of PR #140), the before of what v38 changes
#   cut   v37 with replay_validation/combiner_cutoff_env.patch   the point-cutoff runs of the family v37
#   hook  v37 with replay_validation/bp_hooks_env.patch          the BP1 and BP7 runs of the family v37
#
# The two patches were written against v37 (5d8004a), and apply to the v38 tree of cut38 and hook38 (the tag v0.5.0-rc.1)
# as they are, with every context line matching (`patch --dry-run --verbose --fuzz=0` on an export of the tag): the
# cutoff patch lands on the line it was written for, and the hooks patch with offsets of 2 and 20 lines, because the
# v38 tree has documentation added above both of its hunks. The alternatives patch of `alt` was written against the working
# tree that became v38 and applies to the tag with one hunk at an offset of 30 lines and the rest where they were written,
# checked the same way. `apply_patch` refuses a hunk that needs fuzz, so a patch whose context has changed stops the build
# and is not quietly fitted to the tree by `patch`. (macOS's patch prints an offset only with --verbose; the build is
# quiet, so the offsets are not in its log.)
#
# Writes OUT_DIR/fastvep_<name> and OUT_DIR/SHA256SUMS, and prints the environment variables that
# `11_manuscript_runs.py` and `14_temporal_holdout.py` read their binaries from.
#
# Each binary is built from its own export of the source (`git archive`, so no uncommitted file and no .git: nothing
# here builds a tree that differs from the commit it names) in a target directory of its own. The target directory is
# not shared between trees, and that is the point of this script: cargo decides what is fresh by file times, and the
# trees of one commit and of the same commit with a patch have the same layout, so a patched tree built into the target
# directory of the unpatched one is judged fresh when its patched files are older than the artifacts, and the binary that
# comes out has no patch in it ("Finished `release` profile in 0.18s"). Nothing about such a binary looks wrong, and a
# run made by it is a run of the unpatched classifier under the patched one's name. The check at the end looks for what
# each patched build must contain and each plain build must not.
set -euo pipefail

OUT_DIR="${1:-}"
if [ -z "$OUT_DIR" ]; then
  echo "usage: $(basename "$0") OUT_DIR" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
PATCHES="$ROOT/manuscript/acmg/replay_validation"
V37_REF="${V37_REF:-5d8004a}"
BASE_REF="${BASE_REF:-9d84be0}"
V38_REF="${V38_REF:-v0.5.0-rc.1}"
V38_TREE="${V38_TREE:-}"
BUILD_JOBS="${BUILD_JOBS:-2}"
NAMES="${NAMES:-v38 cut38 hook38 base alt v37 cut hook}"
wants() { case " $NAMES " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

mkdir -p "$OUT_DIR"
OUT_DIR="$(cd "$OUT_DIR" && pwd -P)"
WORK="$OUT_DIR/build"
mkdir -p "$WORK"
for name in $NAMES; do rm -rf "$WORK/src_$name" "$WORK/target_$name"; done

export_ref() { # name ref
  mkdir -p "$WORK/src_$1"
  git -C "$ROOT" archive "$2" | tar -x -C "$WORK/src_$1"
}

# The files of the checkout as a commit of it would hold them: tracked files as they are now, new files that are not
# ignored, and not a tracked file that has been deleted.
export_tree() { # name directory
  mkdir -p "$WORK/src_$1"
  (
    cd "$2"
    git ls-files -z --cached --others --exclude-standard | while IFS= read -r -d '' f; do
      if [ -e "$f" ]; then printf '%s\0' "$f"; fi
    done | tar --null -T - -cf - | tar -x -C "$WORK/src_$1"
  )
}

apply_patch() { # name patch
  (cd "$WORK/src_$1" && patch -p1 --quiet --fuzz=0 < "$PATCHES/$2")
}

for name in $NAMES; do
  case "$name" in
    v37 | cut | hook) export_ref "$name" "$V37_REF" ;;
    base) export_ref base "$BASE_REF" ;;
    v38 | alt | cut38 | hook38)
      if [ -n "$V38_TREE" ]; then export_tree "$name" "$V38_TREE"; else export_ref "$name" "$V38_REF"; fi
      ;;
    *) echo "unknown name: $name (v38 cut38 hook38 base alt v37 cut hook)" >&2; exit 2 ;;
  esac
done
if wants cut; then apply_patch cut combiner_cutoff_env.patch; fi
if wants hook; then apply_patch hook bp_hooks_env.patch; fi
if wants alt; then apply_patch alt inheritance_alternatives_env.patch; fi
if wants cut38; then apply_patch cut38 combiner_cutoff_env.patch; fi
if wants hook38; then apply_patch hook38 bp_hooks_env.patch; fi

build() { # name
  local name=$1
  (cd "$WORK/src_$name" && CARGO_TARGET_DIR="$WORK/target_$name" cargo build --release -p fastvep-cli --locked) \
    > "$WORK/build_$name.log" 2>&1 || { echo "build of $name failed, see $WORK/build_$name.log" >&2; return 1; }
  cp -p "$WORK/target_$name/release/fastvep" "$OUT_DIR/fastvep_$name"
  echo "built $name"
}
export -f build
export WORK OUT_DIR
printf '%s\n' $NAMES | xargs -P "$BUILD_JOBS" -I{} bash -c 'build {}'

# What each build must carry. The environment hooks of a patched build are strings in its binary; a plain build has none.
carries() { grep -aq "$2" "$OUT_DIR/fastvep_$1"; }
fail=0
require() { # name marker
  if ! carries "$1" "$2"; then echo "fastvep_$1 does not contain $2: it was not built with its patch" >&2; fail=1; fi
}
forbid() { # name marker
  if carries "$1" "$2"; then echo "fastvep_$1 contains $2: it is not the plain build it is meant to be" >&2; fail=1; fi
}
forbid_other() { # name marker
  if carries "$1" "$2"; then echo "fastvep_$1 contains $2, which belongs to another patch: it was built from a tree that had both" >&2; fail=1; fi
}
for cutter in cut cut38; do
  if wants "$cutter"; then
    require "$cutter" FASTVEP_LP_CUT
    require "$cutter" FASTVEP_P_CUT
  fi
done
for hooker in hook hook38; do
  if wants "$hooker"; then
    require "$hooker" FASTVEP_BP1_LOO
    require "$hooker" FASTVEP_BP7_SPLICEAI_BENIGN
  fi
done
if wants alt; then require alt FASTVEP_ALT; fi
# Each patch adds its own variables and no other patch's: a build that carries a marker of the wrong patch was made from a
# tree that had both, and a run of it would be a run of something that no document describes.
for cutter in cut cut38; do
  if wants "$cutter"; then
    for marker in FASTVEP_BP1_LOO FASTVEP_BP7_SPLICEAI_BENIGN FASTVEP_ALT; do forbid_other "$cutter" "$marker"; done
  fi
done
for hooker in hook hook38; do
  if wants "$hooker"; then
    for marker in FASTVEP_LP_CUT FASTVEP_P_CUT FASTVEP_ALT; do forbid_other "$hooker" "$marker"; done
  fi
done
if wants alt; then
  for marker in FASTVEP_LP_CUT FASTVEP_P_CUT FASTVEP_BP1_LOO FASTVEP_BP7_SPLICEAI_BENIGN; do forbid_other alt "$marker"; done
fi
for plain in v37 base v38; do
  if wants "$plain"; then
    for marker in FASTVEP_LP_CUT FASTVEP_BP1_LOO FASTVEP_BP7_SPLICEAI_BENIGN FASTVEP_ALT; do
      forbid "$plain" "$marker"
    done
  fi
done
for patched in cut hook alt cut38 hook38; do
  if wants "$patched"; then
    for plain in v37 v38; do
      if [ -f "$OUT_DIR/fastvep_$plain" ] && cmp -s "$OUT_DIR/fastvep_$patched" "$OUT_DIR/fastvep_$plain"; then
        echo "fastvep_$patched is identical to fastvep_$plain" >&2
        fail=1
      fi
    done
  fi
done
[ "$fail" -eq 0 ] || exit 1

(cd "$OUT_DIR" && shasum -a 256 fastvep_* > SHA256SUMS)
cat "$OUT_DIR/SHA256SUMS"
# the variables of the binaries that were built, and of no other
exports() { # variable name
  if wants "$2"; then echo "export $1=$OUT_DIR/fastvep_$2"; fi
}
echo
exports FASTVEP_V38_BIN v38
exports FASTVEP_CUTOFF38_BIN cut38
exports FASTVEP_HOOKS38_BIN hook38
exports FASTVEP_V38_BASE_BIN base
exports FASTVEP_V38_ALT_BIN alt
exports FASTVEP_V37_BIN v37
exports FASTVEP_CUTOFF_BIN cut
exports FASTVEP_HOOKS_BIN hook
