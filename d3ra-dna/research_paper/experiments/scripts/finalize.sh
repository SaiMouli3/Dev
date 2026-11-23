#!/bin/bash
# Build the submission bundle in research_paper/final/ from the current manuscript.
#   bash research_paper/experiments/scripts/finalize.sh
# Writes IEEE_Research_Paper.pdf, .tex, _Source.zip and COMMIT.txt, and fails loudly if the
# manuscript still has unresolved references or LaTeX errors.
set -eu
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
M=$ROOT/research_paper/manuscript
F=$ROOT/research_paper/final
cd "$M"

# main.tex and supplement.tex reference each other through the xr package
pdflatex -interaction=nonstopmode main.tex >/dev/null
bibtex main >/dev/null 2>&1 || true
pdflatex -interaction=nonstopmode supplement.tex >/dev/null
bibtex supplement >/dev/null 2>&1 || true
for i in 1 2; do
  pdflatex -interaction=nonstopmode main.tex >/dev/null
  pdflatex -interaction=nonstopmode supplement.tex >/dev/null
done

fail=0
for d in main supplement; do
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
cp main.pdf "$F/IEEE_Research_Paper.pdf"
cp main.tex "$F/IEEE_Research_Paper.tex"
cp supplement.pdf "$F/IEEE_Research_Paper_Supplement.pdf"
# The paper's content commit. finalize.sh runs before the bundle is committed, so record the commit
# that the build came from and note that the bundle itself is committed on top of it.
git -C "$ROOT" rev-parse HEAD > "$F/COMMIT.txt"
{ echo "Repository: https://github.com/SaiMouli3/Dev (directory d3ra-dna)"
  echo "Built: $(date -u '+%Y-%m-%d %H:%M UTC')"
  echo "Environment: $(python3 -V 2>&1), DT4DDS $(pip show dt4dds 2>/dev/null | awk '/^Version/{print $2}')"
  echo "Note: this is the commit the PDF was built from; the built bundle is committed immediately after,"
  echo "so the repository tip is normally one commit ahead and differs only by these files."
} >> "$F/COMMIT.txt"

rm -f "$F/IEEE_Research_Paper_Source.zip"
zip -qr "$F/IEEE_Research_Paper_Source.zip" main.tex supplement.tex references.bib tables figures IEEEtran.cls IEEEtran.bst
echo "finalize: $(pdfinfo main.pdf | awk '/Pages/{print $2}') pages, $(unzip -l "$F/IEEE_Research_Paper_Source.zip" | tail -1 | awk '{print $2}') files in the source archive"
