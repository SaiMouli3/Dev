# D3RA results tables
Costs are in cost units (reads x cost per read); `cap@t` is the budget cap reaching reliability t, `mean_spend@t` the mean spend under that cap. Brackets: 95% CI.

## Baselines A-F per scenario

_source: `baselines.csv`_

| scenario         | design                 | variant     |   B_ref | cap@0.99                | mean_spend@0.99         | cap@0.999               | mean_spend@0.999        | p_exact@ref          |   mean_cost_success |   mean_rounds |   mean_containers_accessed |   redundancy_overhead |   runtime_s |
|:-----------------|:-----------------------|:------------|--------:|:------------------------|:------------------------|:------------------------|:------------------------|:---------------------|--------------------:|--------------:|---------------------------:|----------------------:|------------:|
| s1_homogeneous   | A_uniform_fixed        |             |  11,281 | 11,281 [10,434, 12,009] | 11,281 [10,434, 12,009] | 19,010 [17,168, 21,387] | 19,010 [17,168, 21,387] | 0.990 [0.989, 0.991] |               6,406 |               |                            |                   0.5 |    5.25e-06 |
| s1_homogeneous   | B_reliability_fixed    |             |  11,281 | 11,315 [10,501, 12,154] | 11,315 [10,501, 12,154] | 17,645 [16,306, 21,074] | 17,645 [16,306, 21,074] | 0.990 [0.988, 0.991] |               6,401 |               |                            |                   0.5 |    6.468    |
| s1_homogeneous   | C_uniform_adaptive     |             |  11,281 | 7,950 [7,750, 8,650]    | 5,894 [5,869, 5,921]    | 11,800 [10,100, 17,150] | 5,911 [5,881, 5,942]    | 0.998 [0.994, 0.999] |               5,915 |         118.3 |                          8 |                   0.5 |   63.1      |
| s1_homogeneous   | D_reliability_adaptive |             |  11,281 | 8,100 [7,850, 8,550]    | 5,909 [5,882, 5,935]    | 11,250 [9,400, 13,100]  | 5,921 [5,891, 5,951]    | 0.999 [0.996, 1.000] |               5,922 |         118.4 |                          8 |                   0.5 |   58.7      |
| s1_homogeneous   | E_joint_adaptive       | reliability |  11,281 | 8,100 [7,850, 8,550]    | 5,909 [5,882, 5,935]    | 11,250 [9,400, 13,100]  | 5,921 [5,891, 5,951]    | 0.999 [0.996, 1.000] |               5,922 |         118.4 |                          8 |                   0.5 |   58.1      |
| s1_homogeneous   | F_oracle               | reliability |  11,281 | 8,150 [7,850, 8,650]    | 5,906 [5,879, 5,932]    | 10,350 [9,550, 12,250]  | 5,916 [5,887, 5,944]    | 1.000 [0.997, 1.000] |               5,917 |         118.3 |                          8 |                   0.5 |   56.1      |
| s2_heterogeneous | A_uniform_fixed        |             |  20,755 | 20,755 [18,986, 23,023] | 20,755 [18,986, 23,023] | inf [inf, inf]          | inf [inf, inf]          | 0.990 [0.989, 0.991] |               8,115 |               |                            |                   0.5 |    4.77e-06 |
| s2_heterogeneous | B_reliability_fixed    |             |  20,755 | 9,893 [9,763, 10,103]   | 9,893 [9,763, 10,103]   | inf [inf, inf]          | inf [inf, inf]          | 0.999 [0.998, 0.999] |               5,789 |               |                            |                   0.5 |    6.945    |
| s2_heterogeneous | C_uniform_adaptive     |             |  20,755 | 12,960 [12,225, 13,505] | 7,194 [7,134, 7,254]    | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.995, 0.999] |               7,189 |         190   |                          8 |                   0.5 |  102.4      |
| s2_heterogeneous | D_reliability_adaptive |             |  20,755 | 7,620 [7,345, 7,823]    | 5,529 [5,504, 5,555]    | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.995, 0.999] |               5,530 |         130   |                          8 |                   0.5 |   65.3      |
| s2_heterogeneous | E_joint_adaptive       | failure1    |  20,755 | 5,450 [5,350, 5,700]    | 4,330 [4,313, 4,346]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.992, 0.998] |               4,327 |          86.2 |                          1 |                   0.5 |   47.6      |
| s2_heterogeneous | F_oracle               | failure1    |  20,755 | 5,450 [5,350, 5,700]    | 4,330 [4,313, 4,346]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.992, 0.998] |               4,327 |          86.2 |                          1 |                   0.5 |   45.8      |
| s6_tight_budget  | A_uniform_fixed        |             |  12,453 | 20,755 [18,986, 23,023] | 20,755 [18,986, 23,023] | inf [inf, inf]          | inf [inf, inf]          | 0.965 [0.963, 0.968] |               8,115 |               |                            |                   0.5 |    5.25e-06 |
| s6_tight_budget  | B_reliability_fixed    |             |  12,453 | 9,893 [9,763, 10,103]   | 9,893 [9,763, 10,103]   | inf [inf, inf]          | inf [inf, inf]          | 0.997 [0.996, 0.998] |               5,789 |               |                            |                   0.5 |    6.257    |
| s6_tight_budget  | C_uniform_adaptive     |             |  12,453 | 12,960 [12,225, 13,505] | 7,194 [7,134, 7,254]    | inf [inf, inf]          | inf [inf, inf]          | 0.986 [0.980, 0.990] |               7,189 |         190   |                          8 |                   0.5 |   95.9      |
| s6_tight_budget  | D_reliability_adaptive |             |  12,453 | 7,620 [7,345, 7,823]    | 5,529 [5,504, 5,555]    | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.995, 0.999] |               5,530 |         130   |                          8 |                   0.5 |   65.6      |
| s6_tight_budget  | E_joint_adaptive       | failure1    |  12,453 | 5,450 [5,350, 5,700]    | 4,330 [4,313, 4,346]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.992, 0.998] |               4,327 |          86.2 |                          1 |                   0.5 |   48.6      |
| s6_tight_budget  | F_oracle               | failure1    |  12,453 | 5,450 [5,350, 5,700]    | 4,330 [4,313, 4,346]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.992, 0.998] |               4,327 |          86.2 |                          1 |                   0.5 |   45.9      |

## Design comparisons (bootstrap mean cost difference, Cliff's delta, Holm-adjusted McNemar p)

_source: `baselines.csv`_

| scenario         | design                                            |   mean_cost_diff |   mean_cost_diff_lo |   mean_cost_diff_hi |   cliffs_delta |   p_holm |
|:-----------------|:--------------------------------------------------|-----------------:|--------------------:|--------------------:|---------------:|---------:|
| s1_homogeneous   | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |          -486.1  |              -550   |            -424.5   |       -0.306   | 0.000837 |
| s1_homogeneous   | cmp:D_reliability_adaptive-vs-B_reliability_fixed |          -464.4  |              -525.7 |            -405.1   |       -0.296   | 8.53e-06 |
| s1_homogeneous   | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |             6.9  |               -36.7 |              49.1   |        0.00757 | 1        |
| s1_homogeneous   | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |             0    |               -43.2 |              43     |        0       | 1        |
| s1_homogeneous   | cmp:F_oracle-vs-E_joint_adaptive                  |            -5.25 |               -15   |               4.651 |       -0.00438 | 1        |
| s2_heterogeneous | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |          -937.8  |             -1188   |            -675.7   |       -0.394   | 0.00033  |
| s2_heterogeneous | cmp:D_reliability_adaptive-vs-B_reliability_fixed |          -283.2  |              -509.7 |             -66.2   |       -0.155   | 1        |
| s2_heterogeneous | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |         -2698    |             -2960   |           -2416     |       -0.991   | 1        |
| s2_heterogeneous | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |         -1043    |             -1308   |            -772.8   |       -0.94    | 1        |
| s2_heterogeneous | cmp:F_oracle-vs-E_joint_adaptive                  |             0    |                 0   |               0     |        0       | 1        |
| s6_tight_budget  | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |          -937.8  |             -1188   |            -675.7   |       -0.394   | 2.71e-07 |
| s6_tight_budget  | cmp:D_reliability_adaptive-vs-B_reliability_fixed |          -283.2  |              -509.7 |             -66.2   |       -0.155   | 1        |
| s6_tight_budget  | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |         -2698    |             -2960   |           -2416     |       -0.991   | 0.00454  |
| s6_tight_budget  | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |         -1043    |             -1308   |            -772.8   |       -0.94    | 1        |
| s6_tight_budget  | cmp:F_oracle-vs-E_joint_adaptive                  |             0    |                 0   |               0     |        0       | 1        |

## H1 - predictive accuracy of channel models

_source: `h1_model_accuracy.csv`_

| model                     |   mae_p |   max_abs_err |   B99_pred |   B99_mc |   B99_rel_err |
|:--------------------------|--------:|--------------:|-----------:|---------:|--------------:|
| heterogeneous(full)       | 0.00674 |         0.015 |     73,304 |   95,451 |        -0.232 |
| homogeneous-avg           | 0.073   |         0.168 |     79,569 |   95,451 |        -0.166 |
| poisson-no-overdispersion | 0.36    |         0.892 |     11,900 |   95,451 |        -0.875 |
| no-failure-model          | 0.067   |         0.148 |     24,059 |   95,451 |        -0.748 |

## Redundancy allocation (uniform vs risk-adaptive parity)

_source: `redundancy.csv`_

| experiment    |   reads_per_fragment | method        |   parity_total |   overhead |   parity_min |   parity_max |   p_file_pred |   p_file_mc |
|:--------------|---------------------:|:--------------|---------------:|-----------:|-------------:|-------------:|--------------:|------------:|
| redundancy    |                    6 | uniform       |          2,080 |      1.625 |              |              |         0.999 |             |
| redundancy    |                    6 | risk-adaptive |          1,740 |      1.359 |           35 |           56 |         0.999 |             |
| redundancy-mc |                    6 | uniform       |                |            |              |              |               |       0.987 |
| redundancy-mc |                    6 | risk-adaptive |                |            |              |              |               |       0.984 |
| redundancy    |                   10 | uniform       |          1,760 |      1.375 |              |              |         0.999 |             |
| redundancy    |                   10 | risk-adaptive |          1,200 |      0.938 |           24 |           45 |         0.999 |             |
| redundancy-mc |                   10 | uniform       |                |            |              |              |               |       0.991 |
| redundancy-mc |                   10 | risk-adaptive |                |            |              |              |               |       0.983 |

## Replication vs RS parity under shared bias

_source: `replication_vs_parity.csv`_

|   shared_variance_fraction | scheme         |   stored_molecules | cap@0.99                |   mean_cost_success |   p_success_uncapped |
|---------------------------:|:---------------|-------------------:|:------------------------|--------------------:|---------------------:|
|                       0    | replication-2x |              1,280 | inf [inf, inf]          |              79,299 |                0.537 |
|                       0    | rs-parity-2x   |              1,280 | 10,742 [10,576, 10,909] |               6,943 |                1     |
|                       0.25 | replication-2x |              1,280 | inf [inf, inf]          |              60,857 |                0.541 |
|                       0.25 | rs-parity-2x   |              1,280 | 10,596 [10,431, 10,797] |               6,853 |                1     |
|                       0.5  | replication-2x |              1,280 | inf [inf, inf]          |              53,900 |                0.541 |
|                       0.5  | rs-parity-2x   |              1,280 | 10,200 [10,066, 10,426] |               6,777 |                1     |
|                       0.75 | replication-2x |              1,280 | inf [inf, inf]          |              51,717 |                0.539 |
|                       0.75 | rs-parity-2x   |              1,280 | 9,961 [9,809, 10,113]   |               6,709 |                1     |
|                       1    | replication-2x |              1,280 | inf [inf, inf]          |              53,829 |                0.539 |
|                       1    | rs-parity-2x   |              1,280 | 10,152 [9,974, 10,361]  |               6,709 |                1     |

## Placement optimisers

_source: `placement_opt.csv`_

| scenario               |   stripes | method              |   objective |   gap_vs_ilp |      time |
|:-----------------------|----------:|:--------------------|------------:|-------------:|----------:|
| s2_heterogeneous       |         4 | random              |       2.049 |        0.961 |  0.00021  |
| s2_heterogeneous       |         4 | round_robin         |       1.298 |        0.242 |  5.15e-05 |
| s2_heterogeneous       |         4 | greedy_tiered       |       1.06  |        0.014 |  0.000396 |
| s2_heterogeneous       |         4 | local_search        |       1.06  |        0.014 |  0.091    |
| s2_heterogeneous       |         4 | simulated_annealing |       1.06  |        0.014 |  0.747    |
| s2_heterogeneous       |         4 | genetic             |       1.06  |        0.014 |  0.116    |
| s2_heterogeneous       |         4 | ilp                 |       1.045 |        0     |  0.381    |
| s4_correlated_failures |         4 | random              |       2.99  |        1.002 |  0.000193 |
| s4_correlated_failures |         4 | round_robin         |       2.329 |        0.559 |  6.77e-05 |
| s4_correlated_failures |         4 | greedy_tiered       |       1.493 |        0     |  0.00084  |
| s4_correlated_failures |         4 | local_search        |       1.493 |        0     |  0.099    |
| s4_correlated_failures |         4 | simulated_annealing |       1.493 |        0     |  0.78     |
| s4_correlated_failures |         4 | genetic             |       1.493 |        0     |  0.169    |
| s4_correlated_failures |         4 | ilp                 |       1.493 |        0     |  3.261    |
| s10_mixed              |         6 | random              |       6.325 |        3.067 |  0.000165 |
| s10_mixed              |         6 | round_robin         |       2.036 |        0.309 |  5.03e-05 |
| s10_mixed              |         6 | greedy_tiered       |       1.679 |        0.08  |  0.000729 |
| s10_mixed              |         6 | local_search        |       1.679 |        0.08  |  0.101    |
| s10_mixed              |         6 | simulated_annealing |       1.679 |        0.08  |  0.887    |
| s10_mixed              |         6 | genetic             |       1.679 |        0.08  |  0.111    |
| s10_mixed              |         6 | ilp                 |       1.555 |        0     | 11.3      |

## Scalability

_source: `scalability.csv`_

|   containers |   stripes |   fragments |   placement_s |   predictor_bisection_s |   mc_static_2000_s |   adaptive_s |   peak_mem_mb |   B99_pred |   B99_mc |
|-------------:|----------:|------------:|--------------:|------------------------:|-------------------:|-------------:|--------------:|-----------:|---------:|
|           10 |        20 |         880 |         0.023 |                   2.102 |              0.462 |        2.502 |         171.7 |        inf |      inf |
|           50 |       100 |        4400 |         0.147 |                   0.161 |              4.732 |       45.6   |         858.7 |        inf |      inf |
|          100 |       200 |        8800 |         0.243 |                   2.487 |              8.152 |      110     |        1721   |      55008 |      inf |
|          100 |       600 |       26400 |         1.328 |                   4.942 |             28.7   |      473     |        5158   |     179702 |      inf |
