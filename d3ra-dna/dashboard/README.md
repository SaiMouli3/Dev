# D3RA benchmarking dashboard

Streamlit app for browsing the experiment results in `results/*.csv` and running
small interactive simulations.

## Run

From the repository root:

```bash
pip install streamlit matplotlib pandas   # plus the d3ra dependencies (numpy, scipy)
streamlit run dashboard/app.py
```

The app adds `../src` to `sys.path`, so the package does not need to be installed.
Pick a different results directory in the sidebar; press **Reload** after the
experiments write new CSVs.

## Tabs

| Tab | Source |
|---|---|
| Overview | `baselines.csv` - headline table per scenario, cheapest non-oracle design |
| Baselines | `baselines.csv` - spend / cap / P(exact) bars with CIs, comparison forest plot |
| Hypotheses | `h1_model_accuracy.csv`, `h5_robustness.csv` |
| Ablations | `ablations.csv` - change relative to `full (E)` |
| Redundancy & Replication | `redundancy.csv`, `replication_vs_parity.csv` |
| Placement optimisers | `placement_opt.csv` - gap vs ILP and runtime |
| Scalability | `scalability.csv` - log-log runtime vs fragments |
| Interactive simulator | live run of `AdaptiveController` vs a static uniform allocation |

Missing CSVs are reported in place; the rest of the app still works.

The simulator uses `calibrated=False` scenarios (analytic consensus curves) for
speed. The adaptive controller's cost cap is a multiple of the static policy's
0.99 budget (or its 0.97 quantile when 0.99 is unreachable, e.g. under container
failures); trials that cannot complete stop at the cap. A few hundred trials
take seconds to a minute, depending on the scenario.

Costs are in cost units = reads x cost per read.

## Static figures

The same figures can be written without Streamlit:

```bash
PYTHONPATH=src python -m d3ra.visualization.plots results figures
```
