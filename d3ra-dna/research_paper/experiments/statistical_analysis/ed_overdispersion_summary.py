"""Summarise the exploratory overdispersed-likelihood experiment (ED-OD) and its held-out check.

ICC 0.01 was chosen after inspecting the Twist + simulation run (paper_ed_overdispersion.csv); the
Genscript run (paper_ed_holdout_genscript.csv) uses that value unchanged. The Genscript PCR_30c
container (54% dropout) is reported separately: its 'true relative efficiency' (ratio of median
*finite* recovery times) is not meaningful when most molecules are absent."""
from pathlib import Path
import pandas as pd

RES = Path(__file__).resolve().parents[3] / "results"
rows = []
for f, tag in [("paper_ed_overdispersion.csv", "Twist+sim (selection)"), ("paper_ed_holdout_genscript.csv", "Genscript (held-out)")]:
    d = pd.read_csv(RES / f)
    for (setting, od), g in d.groupby([d.setting.str.split(" ").str[:2].str.join(" "), "overdispersion"]):
        ok = g[g.container != "PCR_30c"]
        healthy = ok[ok.true_relative_efficiency >= 0.95]
        moderate = ok[(ok.true_relative_efficiency >= 0.7) & (ok.true_relative_efficiency < 0.95)]
        degr = g[g.truly_degraded | (g.container == "PCR_30c")]
        rows.append(dict(source=tag, setting=setting, overdispersion=od,
                         ci90_cov_mean=ok.ci90_coverage.mean(), ci90_cov_min=ok.ci90_coverage.min(), ci90_cov_max=ok.ci90_coverage.max(),
                         flag_healthy_max=healthy.flag_rate.max(), flag_moderate_max=moderate.flag_rate.max() if len(moderate) else float("nan"),
                         flag_degraded_min=degr.flag_rate.min() if len(degr) else float("nan"),
                         abs_err_mean=(ok.posterior_mean_median - ok.true_relative_efficiency).abs().mean()))
out = pd.DataFrame(rows)
out.to_csv(RES / "paper_ed_overdispersion_summary.csv", index=False)
print(out.round(3).to_string())
