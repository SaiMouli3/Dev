#!/bin/bash
# Build the submission bundle in research_paper/final/ from the current manuscript.
#   bash research_paper/experiments/scripts/finalize.sh
# Writes ACM_Research_Paper.pdf, .tex, _Source.zip and COMMIT.txt, and fails loudly if the
# manuscript still has unresolved references or LaTeX errors.
set -eu
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
M=$ROOT/research_paper/manuscript
F=$ROOT/research_paper/final
cd "$M"

# ACM journal format (acmart, acmsmall); appendices are part of main.tex
pdflatex -interaction=nonstopmode main.tex >/dev/null
bibtex main >/dev/null 2>&1 || true
pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode main.tex >/dev/null
# anonymous version for double-blind review: authors, affiliations, acknowledgements and the repository
# link are hidden by acmart's anonymous option; review adds line numbers
sed 's/\\documentclass\[acmsmall,screen\]{acmart}/\\documentclass[acmsmall,screen,review,anonymous]{acmart}/' main.tex > main_anon.tex
grep -q "review,anonymous" main_anon.tex || { echo "finalize: anonymous class option not set"; exit 1; }
pdflatex -interaction=nonstopmode main_anon.tex >/dev/null
bibtex main_anon >/dev/null 2>&1 || true
pdflatex -interaction=nonstopmode main_anon.tex >/dev/null
pdflatex -interaction=nonstopmode main_anon.tex >/dev/null

fail=0
for d in main main_anon; do
  if grep -qE '^! ' $d.log; then echo "LaTeX errors in $d:"; grep -E '^! ' $d.log | head; fail=1; fi
  if grep -q "undefined on input" $d.log; then
    echo "Unresolved references/citations in $d:"
    grep "undefined on input" $d.log | sed "s/.*\`\([^']*\)'.*/  \1/" | sort -u
    fail=1
  fi
  if pdftotext $d.pdf - | grep -qE "TBD|@@"; then echo "TBD or unfilled @@ markers remain in $d.pdf"; fail=1; fi
done
[ "$fail" = 0 ] || { echo "finalize: refusing to publish an inconsistent build"; exit 1; }

mkdir -p "$F"
cp main.pdf "$F/ACM_Research_Paper.pdf"
cp main.tex "$F/ACM_Research_Paper.tex"
cp main_anon.pdf "$F/ACM_Research_Paper_Anonymous.pdf"
if pdftotext main_anon.pdf - | grep -qiE "Mouli|Geetha|Savanur|SaiMouli3|msrit|Ramaiah"; then echo "finalize: anonymous PDF still identifies the authors"; exit 1; fi
# The paper's content commit. finalize.sh runs before the bundle is committed, so record the commit
# that the build came from and note that the bundle itself is committed on top of it.
git -C "$ROOT" rev-parse HEAD > "$F/COMMIT.txt"
{ echo "Repository: https://github.com/SaiMouli3/Dev (directory d3ra-dna)"
  echo "Built: $(date -u '+%Y-%m-%d %H:%M UTC')"
  echo "Environment: $(python3 -V 2>&1), DT4DDS $(pip show dt4dds 2>/dev/null | awk '/^Version/{print $2}')"
  echo "Note: this is the commit the PDF was built from; the built bundle is committed immediately after,"
  echo "so the repository tip is normally one commit ahead and differs only by these files."
} >> "$F/COMMIT.txt"

rm -f "$F/ACM_Research_Paper_Source.zip"
zip -qr "$F/ACM_Research_Paper_Source.zip" main.tex references.bib tables figures ACM-Reference-Format.bst
echo "finalize: $(pdfinfo main.pdf | awk '/Pages/{print $2}') pages, $(unzip -l "$F/ACM_Research_Paper_Source.zip" | tail -1 | awk '{print $2}') files in the source archive"
