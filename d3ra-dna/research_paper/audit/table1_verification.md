# Table I (literature comparison): cell-by-cell verification

Source of truth: `supplementary/literature_review.md`, Sec. 13, the comparison matrix built from the papers (full text, model sections or abstract, as noted per entry in Secs. 1–9). Table: `manuscript/tables/literature_comparison.tex`.

Mapping of symbols: n/s → "--"; "Implicit" → I.

| Row | Multi | Fail | Heterog. | Corr. err. | Adaptive | Stopping | Real reads | DT | Matches S1 §13 |
|---|---|---|---|---|---|---|---|---|---|
| Lopez 2019 | N | N | N | -- | N | P | Y | N | yes |
| Gimpel 2023 | N | N | P | N | N | N | Y | Y (DT4DDS itself) | yes |
| Bar-Lev 2025 | N | N | N | N | N | I | N | N | yes |
| Gruica 2024 | N | N | N | N | N | I | N | N | yes |
| Cao & Chen 2025 | N | N | N | N | N | I | P | N | yes |
| Kas Hanna 2026 | N | N | N | N | N | N | N | N | yes |
| Levy 2026 | Y | Y | N | N | N | I | N | N | yes |
| Kong 2026 | Y | N | P | N | N | I | N | N | yes (abstract only; now stated in the table notes) |
| Leong 2012 | Y | Y | N | N | N | N | N | N | yes |
| Weilguny 2023 | N | N | -- | N | Y | -- | Y | N | yes |
| Bubeck 2013 | Y | N | Y | N | Y | N | N | N | yes |
| This work (rev. 3) | Y | P | Y | P | Y | Y | P | Y | yes, with revision-3 qualifiers (notes below the table) |

Changes against the S1 §13 "This work" row, which describes revision 2:
- **Container failures Y → P:** failures are simulated or injected, not observed.
- **Correlated errors Y → P:** they are validated on nanopore reads but do not enter the controller's likelihood.
- **Real-data evaluation Y → P:** the replay uses oracle decodability and same-batch aliquots.
- **New column "End-to-end digital twin":** Y only for Gimpel 2023 (the DT itself) and this work.

The "1 1 1" corruption reported for an earlier PDF could not be reproduced. The table now uses plain letters with no superscripts in its cells. The rendered table in the compiled PDF was inspected after every build.
