#!/usr/bin/env bash
# Fail when personal genome data would be (or already is) in git.
#
#   scripts/check_no_genome_data.sh --staged          # what the next commit adds (pre-commit hook)
#   scripts/check_no_genome_data.sh --tracked         # every tracked file (CI)
#   scripts/check_no_genome_data.sh --range BASE..TIP # what each commit of a push adds (CI,
#                                                     # pre-push; any `git rev-list` arguments)
#
# Rules: no genomic or result paths (any case), no data file outside the data folders, no file
# over 5 MB, no path or text matching one of this clone's private IDs (a kit ID, a sample ID),
# and no VCF that names a sample other than a public one.
#
# The adversary is an honest mistake by the user or an agent: `dashboard --out x.html`,
# `query --json > x.json`, `cp -R reports reports.bak`, a run's output copied in as a fixture.
# A derived output carries no kit ID, so the data-file rule anchors on where it sits: a mistake
# names a file freely, but only a deliberate choice puts it under a data folder. There a VCF's
# sample columns are the anchor: the pipeline names every sample `sample`, the user's too,
# and an honest copy keeps that, while a fixture from a public sample's run renames it on
# purpose (docs/demo/ holds the README's demo, made from an invented genome, by the same
# deliberate choice). Bound: a Markdown or Python file, or a non-VCF file under a data folder,
# is checked for the private IDs and its size only; what such a file says is the author's to check.
# The private IDs are extended regular expressions, one per line, in the git directory's
# info/genome-private-ids, which is never committed: a public guard must not carry its user's
# kit ID. Without that file the ID rule is off (a fresh clone, CI), and the other rules hold.
# Paths are read NUL-delimited: git's default output quotes and escapes unusual names
# ("\320\264.vcf"), which would slip past the extension patterns. Every git call is checked:
# a failed listing must fail the guard, never read as an empty, clean list.

set -euo pipefail

mode=${1:-}
max_bytes=5242880 # 5 MB
id_file=$(git rev-parse --path-format=absolute --git-common-dir 2>/dev/null)/info/genome-private-ids
id_pattern=
if [ -f "$id_file" ]; then
  id_pattern=$(grep -v -e '^[[:space:]]*#' -e '^[[:space:]]*$' "$id_file" | paste -sd '|' -) || id_pattern=
fi
fail=0
vcf_header=$'^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'

die() { echo "genome-data guard: $*" >&2; exit 2; }
usage() { echo "usage: $0 --staged|--tracked|--range BASE..TIP [REV-LIST-ARG...]" >&2; exit 64; }

case "$mode" in
  --staged | --tracked) ;;
  --range) [ $# -ge 2 ] || usage ;;
  *) usage ;;
esac

tmp=$(mktemp) || die "cannot create a temporary file"
paths=$(mktemp) || die "cannot create a temporary file"
trap 'rm -f "$tmp" "$paths"' EXIT

# check_paths REV: the rules for each path in $paths, with its size from REV's blob (REV
# "": the index, "-": the working tree's file).
# Paths match whatever their case: macOS writes X.VCF.GZ as readily as x.vcf.gz, and git
# tracks either. The extensions are .gitignore's, every index (.tbi, .csi) included, and the
# result folders match at any depth, as .gitignore's do.
check_paths() {
  local f size
  shopt -s nocasematch
  while IFS= read -r -d '' f; do
    case "$f" in
      data/*|reports/*|results/*|work/*|*/data/*|*/reports/*|*/results/*|*/work/*)
        echo "genome-data guard: genomic/personal path: $f" >&2; fail=1; continue ;;
      *.cram|*.crai|*.bam|*.bai|*.sam|*.fq|*.fq.gz|*.fq.bz2|*.fastq|*.fastq.gz|*.fastq.bz2|*.vcf|*.vcf.gz|*.vcf.bgz|*.vcf.zst|*.tbi|*.csi|*.gvcf|*.gvcf.gz|*.bcf|*.bed|*.bed.gz|*.bim|*.fam|*.pgen|*.psam|*.pvar.zst|*.duckdb|*.duckdb.wal|*.sqlite|*.sqlite3|*.db|*.parquet)
        echo "genome-data guard: genomic/personal path: $f" >&2; fail=1; continue ;;
    esac
    # Data files (anything but source, docs and config) sit only where public data is kept.
    case "$f" in
      tests/fixtures/*|src/genokeeper/resources/*|docker/*|docs/demo/*|.github/*|.stayfixed/*) ;;
      .claude/settings.json|.githooks/*|Dockerfile|*/Dockerfile|LICENSE|NOTICE|CITATION.cff) ;;
      *.py|*.md|*.toml|*.yml|*.yaml|*.sh|*.lock|*.sql|*.js|*.css|.gitignore|*.gitkeep|.python-version) ;;
      *) echo "genome-data guard: a data file outside the data folders: $f" >&2; fail=1; continue ;;
    esac
    if [ -n "$id_pattern" ] && [[ "$f" =~ $id_pattern ]]; then
      echo "genome-data guard: path contains a private ID: $f" >&2; fail=1; continue
    fi
    if [ "$1" = - ]; then
      size=$(wc -c < "$f" 2>/dev/null | tr -d ' ') || size=
    else
      size=$(git cat-file -s "$1:$f" 2>/dev/null) || size=
    fi
    if ! [[ "$size" =~ ^[0-9]+$ ]]; then
      echo "genome-data guard: cannot read the size of: $f" >&2; fail=1
    elif [ "$size" -gt "$max_bytes" ]; then
      echo "genome-data guard: file over 5 MB: $f ($size bytes)" >&2; fail=1
    fi
  done < "$paths"
  shopt -u nocasematch
}

# grep_id ARG...: `git grep` for the private IDs in the blobs ARG names, binary ones included (a
# text diff shows a file with a NUL byte as "Binary files differ" and none of its lines).
# git grep: 0 = found, 1 = nothing found, anything else = error.
grep_id() {
  local status=0
  [ -n "$id_pattern" ] || return 0
  git grep -l -E -e "$id_pattern" "$@" > "$tmp" 2>/dev/null || status=$?
  case "$status" in
    0) echo "genome-data guard: text contains a private ID:" >&2
       cat "$tmp" >&2; fail=1 ;;
    1) ;;
    *) die "git grep failed (exit $status)" ;;
  esac
}

# grep_samples ARG...: every VCF column header in the blobs ARG names lists public samples
# only (HG002, NA12878). git grep -z ends each file name with a NUL, its line with a newline.
grep_samples() {
  local status=0 name line sample cols
  git grep -z -I -E -e "$vcf_header" "$@" > "$tmp" 2>/dev/null || status=$?
  case "$status" in
    0) ;;
    1) return ;;
    *) die "git grep failed (exit $status)" ;;
  esac
  while IFS= read -r -d '' name && IFS= read -r line; do
    IFS=$'\t' read -r -a cols <<< "$line"
    for sample in "${cols[@]:9}"; do
      if ! [[ "$sample" =~ ^(HG|NA)[0-9]+$ ]]; then
        echo "genome-data guard: a VCF names a sample that is not public (HG..., NA...): $name" >&2
        fail=1; break
      fi
    done
  done < "$tmp"
}

# check_rev REV: check_paths, and the private IDs in REV's blobs (REV "": the index) of the paths
# in $paths, each matched as written.
check_rev() {
  local f specs=()
  check_paths "$1"
  while IFS= read -r -d '' f; do specs+=(":(top,literal)$f"); done < "$paths"
  [ "${#specs[@]}" -eq 0 ] && return
  grep_id "${1:---cached}" -- "${specs[@]}"
  grep_samples "${1:---cached}" -- "${specs[@]}"
}

case "$mode" in
  --tracked)
    git ls-files -z > "$paths" || die "listing the tracked files failed"
    check_paths -
    grep_id
    grep_samples ;;
  --staged)
    git diff --cached --name-only --diff-filter=ACMR -z > "$paths" ||
      die "listing the staged files failed"
    check_rev "" ;;
  --range)
    # Each commit's own additions against its first parent (a root commit: its whole tree),
    # so a file a later commit deletes is still seen, and so is a merge's own change.
    commits=$(git rev-list "${@:2}") || die "listing the commits of ${*:2} failed"
    for c in $commits; do
      if git rev-parse -q --verify "$c^1" > /dev/null; then
        git diff-tree -r -z --name-only --diff-filter=ACMR "$c^1" "$c" > "$paths"
      else
        git diff-tree -r -z --name-only --no-commit-id --root --diff-filter=ACMR "$c" > "$paths"
      fi || die "listing the files of commit $c failed"
      check_rev "$c"
    done ;;
esac

[ "$fail" -eq 0 ] || { echo "Personal genome data stays out of git." >&2; exit 1; }
