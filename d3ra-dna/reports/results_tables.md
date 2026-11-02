# D3RA results tables
Costs are in cost units (reads x cost per read); `cap@t` is the budget cap reaching reliability t, `mean_spend@t` the mean spend under that cap. Brackets: 95% CI.

## Baselines A-F per scenario

_source: `baselines.csv`_

| scenario                | design                 | variant     |   B_ref | cap@0.99                | mean_spend@0.99         | cap@0.999               | mean_spend@0.999        | p_exact@ref          |   mean_cost_success |   mean_rounds |   mean_containers_accessed |   redundancy_overhead |   runtime_s |
|:------------------------|:-----------------------|:------------|--------:|:------------------------|:------------------------|:------------------------|:------------------------|:---------------------|--------------------:|--------------:|---------------------------:|----------------------:|------------:|
| s1_homogeneous          | A_uniform_fixed        |             |  11,281 | 11,281 [10,434, 12,009] | 11,281 [10,434, 12,009] | 19,010 [17,168, 21,387] | 19,010 [17,168, 21,387] | 0.990 [0.989, 0.991] |               6,406 |               |                            |                 0.5   |    5.25e-06 |
| s1_homogeneous          | B_reliability_fixed    |             |  11,281 | 11,315 [10,501, 12,154] | 11,315 [10,501, 12,154] | 17,645 [16,306, 21,074] | 17,645 [16,306, 21,074] | 0.990 [0.988, 0.991] |               6,401 |               |                            |                 0.5   |    6.468    |
| s1_homogeneous          | C_uniform_adaptive     |             |  11,281 | 7,950 [7,750, 8,650]    | 5,894 [5,869, 5,921]    | 11,800 [10,100, 17,150] | 5,911 [5,881, 5,942]    | 0.998 [0.994, 0.999] |               5,915 |         118.3 |                          8 |                 0.5   |   63.1      |
| s1_homogeneous          | D_reliability_adaptive |             |  11,281 | 8,100 [7,850, 8,550]    | 5,909 [5,882, 5,935]    | 11,250 [9,400, 13,100]  | 5,921 [5,891, 5,951]    | 0.999 [0.996, 1.000] |               5,922 |         118.4 |                          8 |                 0.5   |   58.7      |
| s1_homogeneous          | E_joint_adaptive       | reliability |  11,281 | 8,100 [7,850, 8,550]    | 5,909 [5,882, 5,935]    | 11,250 [9,400, 13,100]  | 5,921 [5,891, 5,951]    | 0.999 [0.996, 1.000] |               5,922 |         118.4 |                          8 |                 0.5   |   58.1      |
| s1_homogeneous          | F_oracle               | reliability |  11,281 | 8,150 [7,850, 8,650]    | 5,906 [5,879, 5,932]    | 10,350 [9,550, 12,250]  | 5,916 [5,887, 5,944]    | 1.000 [0.997, 1.000] |               5,917 |         118.3 |                          8 |                 0.5   |   56.1      |
| s2_heterogeneous        | A_uniform_fixed        |             |  20,755 | 20,755 [18,986, 23,023] | 20,755 [18,986, 23,023] | inf [inf, inf]          | inf [inf, inf]          | 0.990 [0.989, 0.991] |               8,115 |               |                            |                 0.5   |    4.77e-06 |
| s2_heterogeneous        | B_reliability_fixed    |             |  20,755 | 9,893 [9,763, 10,103]   | 9,893 [9,763, 10,103]   | inf [inf, inf]          | inf [inf, inf]          | 0.999 [0.998, 0.999] |               5,789 |               |                            |                 0.5   |    6.945    |
| s2_heterogeneous        | C_uniform_adaptive     |             |  20,755 | 12,960 [12,225, 13,505] | 7,194 [7,134, 7,254]    | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.995, 0.999] |               7,189 |         190   |                          8 |                 0.5   |  102.4      |
| s2_heterogeneous        | D_reliability_adaptive |             |  20,755 | 7,620 [7,345, 7,823]    | 5,529 [5,504, 5,555]    | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.995, 0.999] |               5,530 |         130   |                          8 |                 0.5   |   65.3      |
| s2_heterogeneous        | E_joint_adaptive       | failure1    |  20,755 | 5,450 [5,350, 5,700]    | 4,330 [4,313, 4,346]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.992, 0.998] |               4,327 |          86.2 |                          1 |                 0.5   |   47.6      |
| s2_heterogeneous        | F_oracle               | failure1    |  20,755 | 5,450 [5,350, 5,700]    | 4,330 [4,313, 4,346]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.992, 0.998] |               4,327 |          86.2 |                          1 |                 0.5   |   45.8      |
| s3_nonuniform_abundance | A_uniform_fixed        |             |  73,974 | 73,974 [68,729, 79,776] | 73,974 [68,729, 79,776] | inf [inf, inf]          | inf [inf, inf]          | 0.990 [0.989, 0.991] |              19,201 |               |                            |                 0.5   |    5.01e-06 |
| s3_nonuniform_abundance | B_reliability_fixed    |             |  73,974 | 29,659 [28,759, 31,173] | 29,659 [28,759, 31,173] | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.998, 0.999] |              12,127 |               |                            |                 0.5   |    6.147    |
| s3_nonuniform_abundance | C_uniform_adaptive     |             |  73,974 | 39,005 [32,189, 49,500] | 15,559 [15,379, 15,748] | inf [inf, inf]          | inf [inf, inf]          | 0.997 [0.993, 0.999] |              15,624 |         439.9 |                          8 |                 0.5   |  530.4      |
| s3_nonuniform_abundance | D_reliability_adaptive |             |  73,974 | 23,235 [20,075, 28,505] | 10,857 [10,753, 10,962] | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.995, 0.999] |              10,886 |         262.9 |                          8 |                 0.5   |  153.3      |
| s3_nonuniform_abundance | E_joint_adaptive       | failure1    |  73,974 | 11,100 [10,724, 11,950] | 7,760 [7,715, 7,806]    | inf [inf, inf]          | inf [inf, inf]          | 0.997 [0.993, 0.998] |               7,754 |         154.5 |                          1 |                 0.5   |   83.4      |
| s3_nonuniform_abundance | F_oracle               | failure1    |  73,974 | 11,100 [10,724, 11,950] | 7,760 [7,715, 7,806]    | inf [inf, inf]          | inf [inf, inf]          | 0.997 [0.993, 0.998] |               7,754 |         154.5 |                          1 |                 0.5   |   79.2      |
| s4_correlated_failures  | A_uniform_fixed        |             |  27,637 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.972 [0.970, 0.974] |               8,239 |               |                            |                 0.625 |    5.96e-06 |
| s4_correlated_failures  | B_reliability_fixed    |             |  27,637 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.951 [0.948, 0.954] |               5,792 |               |                            |                 0.625 |    6.921    |
| s4_correlated_failures  | C_uniform_adaptive     |             |  27,637 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.978 [0.970, 0.983] |               6,897 |         177.8 |                          8 |                 0.625 |  109.8      |
| s4_correlated_failures  | D_reliability_adaptive |             |  27,637 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.955 [0.944, 0.963] |               5,375 |         121.6 |                          8 |                 0.625 |   72.9      |
| s4_correlated_failures  | E_joint_adaptive       | failure1    |  27,637 | 20,055 [inf, inf]       | 5,217 [5,120, 5,317]    | inf [inf, inf]          | inf [inf, inf]          | 0.991 [0.985, 0.994] |               5,077 |         129.5 |                          5 |                 0.625 |   87.5      |
| s4_correlated_failures  | F_oracle               | failure1    |  27,637 | 20,145 [inf, inf]       | 5,214 [5,117, 5,314]    | inf [inf, inf]          | inf [inf, inf]          | 0.991 [0.985, 0.994] |               5,071 |         129.2 |                          5 |                 0.625 |   77.4      |
| s5_high_errors          | A_uniform_fixed        |             |  55,573 | 55,573 [50,237, 62,375] | 55,573 [50,237, 62,375] | inf [inf, inf]          | inf [inf, inf]          | 0.990 [0.989, 0.991] |              15,856 |               |                            |                 0.5   |    8.11e-06 |
| s5_high_errors          | B_reliability_fixed    |             |  55,573 | 13,904 [13,223, 14,531] | 13,904 [13,223, 14,531] | inf [inf, inf]          | inf [inf, inf]          | 0.994 [0.993, 0.995] |               8,148 |               |                            |                 0.5   |    6.99     |
| s5_high_errors          | C_uniform_adaptive     |             |  55,573 | 33,440 [28,345, 37,360] | 15,019 [14,875, 15,174] | inf [inf, inf]          | inf [inf, inf]          | 0.997 [0.993, 0.999] |              15,002 |         373.8 |                          8 |                 0.5   |  204.8      |
| s5_high_errors          | D_reliability_adaptive |             |  55,573 | 12,980 [12,265, 13,720] | 8,894 [8,846, 8,942]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.991, 0.998] |               8,879 |         188.3 |                          7 |                 0.5   |   98.9      |
| s5_high_errors          | E_joint_adaptive       | failure1    |  55,573 | 6,950 [6,750, 7,700]    | 5,595 [5,577, 5,614]    | inf [inf, inf]          | inf [inf, inf]          | 0.994 [0.990, 0.997] |               5,589 |         111.1 |                          1 |                 0.5   |   64.9      |
| s5_high_errors          | F_oracle               | failure1    |  55,573 | 6,950 [6,750, 7,700]    | 5,595 [5,577, 5,614]    | inf [inf, inf]          | inf [inf, inf]          | 0.994 [0.990, 0.997] |               5,589 |         111.1 |                          1 |                 0.5   |   59.2      |
| s6_tight_budget         | A_uniform_fixed        |             |  12,453 | 20,755 [18,986, 23,023] | 20,755 [18,986, 23,023] | inf [inf, inf]          | inf [inf, inf]          | 0.965 [0.963, 0.968] |               8,115 |               |                            |                 0.5   |    5.25e-06 |
| s6_tight_budget         | B_reliability_fixed    |             |  12,453 | 9,893 [9,763, 10,103]   | 9,893 [9,763, 10,103]   | inf [inf, inf]          | inf [inf, inf]          | 0.997 [0.996, 0.998] |               5,789 |               |                            |                 0.5   |    6.257    |
| s6_tight_budget         | C_uniform_adaptive     |             |  12,453 | 12,960 [12,225, 13,505] | 7,194 [7,134, 7,254]    | inf [inf, inf]          | inf [inf, inf]          | 0.986 [0.980, 0.990] |               7,189 |         190   |                          8 |                 0.5   |   95.9      |
| s6_tight_budget         | D_reliability_adaptive |             |  12,453 | 7,620 [7,345, 7,823]    | 5,529 [5,504, 5,555]    | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.995, 0.999] |               5,530 |         130   |                          8 |                 0.5   |   65.6      |
| s6_tight_budget         | E_joint_adaptive       | failure1    |  12,453 | 5,450 [5,350, 5,700]    | 4,330 [4,313, 4,346]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.992, 0.998] |               4,327 |          86.2 |                          1 |                 0.5   |   48.6      |
| s6_tight_budget         | F_oracle               | failure1    |  12,453 | 5,450 [5,350, 5,700]    | 4,330 [4,313, 4,346]    | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.992, 0.998] |               4,327 |          86.2 |                          1 |                 0.5   |   45.9      |
| s7_channel_change       | A_uniform_fixed        |             |  36,433 | 36,433 [34,859, 37,959] | 36,433 [34,859, 37,959] | 92,098 [inf, inf]       | 92,098 [inf, inf]       | 0.990 [0.989, 0.991] |              14,784 |               |                            |                 0.5   |    5.72e-06 |
| s7_channel_change       | B_reliability_fixed    |             |  36,433 | 36,443 [35,557, 37,280] | 36,443 [35,557, 37,280] | inf [inf, inf]          | inf [inf, inf]          | 0.990 [0.988, 0.991] |              16,258 |               |                            |                 0.5   |    6.501    |
| s7_channel_change       | C_uniform_adaptive     |             |  36,433 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.898 [0.884, 0.911] |              20,288 |         752.2 |                          8 |                 0.5   |  499.6      |
| s7_channel_change       | D_reliability_adaptive |             |  36,433 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.756 [0.737, 0.774] |              25,865 |         755.2 |                          8 |                 0.5   |  440        |
| s7_channel_change       | E_joint_adaptive       | failure1    |  36,433 | 30,600 [30,300, 32,300] | 25,102 [25,023, 25,181] | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.991, 0.998] |              25,084 |         499.4 |                          1 |                 0.5   |  270.9      |
| s7_channel_change       | F_oracle               | failure1    |  36,433 | 30,600 [30,300, 32,300] | 25,102 [25,023, 25,181] | inf [inf, inf]          | inf [inf, inf]          | 0.996 [0.991, 0.998] |              25,084 |         499.4 |                          1 |                 0.5   |  265        |
| s9_uncertain_estimates  | A_uniform_fixed        |             |  20,755 | 20,755 [18,986, 23,023] | 20,755 [18,986, 23,023] | inf [inf, inf]          | inf [inf, inf]          | 0.990 [0.989, 0.991] |               8,115 |               |                            |                 0.5   |    7.63e-06 |
| s9_uncertain_estimates  | B_reliability_fixed    |             |  20,755 | 12,454 [12,119, 12,776] | 12,454 [12,119, 12,776] | inf [inf, inf]          | inf [inf, inf]          | 0.998 [0.997, 0.998] |               5,984 |               |                            |                 0.5   |    7.062    |
| s9_uncertain_estimates  | C_uniform_adaptive     |             |  20,755 | 28,120 [16,470, 44,825] | 9,023 [8,917, 9,138]    | inf [inf, inf]          | inf [inf, inf]          | 0.987 [0.980, 0.991] |               9,141 |         225.3 |                          8 |                 0.5   |  174.4      |
| s9_uncertain_estimates  | D_reliability_adaptive |             |  20,755 | 8,030 [7,555, 8,860]    | 5,363 [5,337, 5,390]    | inf [inf, inf]          | inf [inf, inf]          | 0.999 [0.996, 0.999] |               5,370 |         124.4 |                          6 |                 0.5   |   62.6      |
| s9_uncertain_estimates  | E_joint_adaptive       | failure1    |  20,755 | 7,455 [6,915, 8,725]    | 4,753 [4,730, 4,778]    | inf [inf, inf]          | inf [inf, inf]          | 0.997 [0.993, 0.998] |               4,751 |         100.7 |                          3 |                 0.5   |   54.6      |
| s9_uncertain_estimates  | F_oracle               | failure1    |  20,755 | 7,495 [6,915, 8,645]    | 4,771 [4,747, 4,797]    | inf [inf, inf]          | inf [inf, inf]          | 0.997 [0.993, 0.998] |               4,768 |         101.9 |                          3 |                 0.5   |   52.8      |
| s10_mixed               | A_uniform_fixed        |             | 155,477 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.964 [0.962, 0.967] |              25,219 |               |                            |                 0.625 |    6.2e-06  |
| s10_mixed               | B_reliability_fixed    |             | 155,477 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.937 [0.934, 0.941] |              10,356 |               |                            |                 0.625 |    8.068    |
| s10_mixed               | C_uniform_adaptive     |             | 155,477 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.972 [0.963, 0.978] |              25,234 |         626.3 |                          8 |                 0.625 |  663.8      |
| s10_mixed               | D_reliability_adaptive |             | 155,477 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.940 [0.929, 0.950] |              10,033 |         212.8 |                          6 |                 0.625 |  143.5      |
| s10_mixed               | E_joint_adaptive       | domain      | 155,477 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.969 [0.960, 0.975] |              25,634 |         633.2 |                          8 |                 0.625 |  762.3      |
| s10_mixed               | F_oracle               | domain      | 155,477 | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | inf [inf, inf]          | 0.969 [0.960, 0.975] |              26,231 |         643.2 |                          8 |                 0.625 |  721.3      |

## Design comparisons (bootstrap mean cost difference, Cliff's delta, Holm-adjusted McNemar p)

_source: `baselines.csv`_

| scenario                | design                                            |   mean_cost_diff |   mean_cost_diff_lo |   mean_cost_diff_hi |   cliffs_delta |    p_holm |
|:------------------------|:--------------------------------------------------|-----------------:|--------------------:|--------------------:|---------------:|----------:|
| s1_homogeneous          | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |         -486.1   |            -550     |            -424.5   |       -0.306   | 0.000837  |
| s1_homogeneous          | cmp:D_reliability_adaptive-vs-B_reliability_fixed |         -464.4   |            -525.7   |            -405.1   |       -0.296   | 8.53e-06  |
| s1_homogeneous          | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |            6.9   |             -36.7   |              49.1   |        0.00757 | 1         |
| s1_homogeneous          | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |            0     |             -43.2   |              43     |        0       | 1         |
| s1_homogeneous          | cmp:F_oracle-vs-E_joint_adaptive                  |           -5.25  |             -15     |               4.651 |       -0.00438 | 1         |
| s2_heterogeneous        | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |         -937.8   |           -1188     |            -675.7   |       -0.394   | 0.00033   |
| s2_heterogeneous        | cmp:D_reliability_adaptive-vs-B_reliability_fixed |         -283.2   |            -509.7   |             -66.2   |       -0.155   | 1         |
| s2_heterogeneous        | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |        -2698     |           -2960     |           -2416     |       -0.991   | 1         |
| s2_heterogeneous        | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |        -1043     |           -1308     |            -772.8   |       -0.94    | 1         |
| s2_heterogeneous        | cmp:F_oracle-vs-E_joint_adaptive                  |            0     |               0     |               0     |        0       | 1         |
| s3_nonuniform_abundance | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |        -6130     |          -13049     |           -2235     |       -0.277   | 0.00432   |
| s3_nonuniform_abundance | cmp:D_reliability_adaptive-vs-B_reliability_fixed |        -1178     |           -2004     |            -384.2   |       -0.301   | 1         |
| s3_nonuniform_abundance | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |        -7562     |           -8548     |           -6570     |       -0.992   | 1         |
| s3_nonuniform_abundance | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |        -2694     |           -3639     |           -1705     |       -0.895   | 1         |
| s3_nonuniform_abundance | cmp:F_oracle-vs-E_joint_adaptive                  |            0     |               0     |               0     |        0       | 1         |
| s4_correlated_failures  | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |        -1282     |           -2230     |            -290.1   |       -0.265   | 0.583     |
| s4_correlated_failures  | cmp:D_reliability_adaptive-vs-B_reliability_fixed |         -378.4   |           -1797     |             923.4   |       -0.11    | 1         |
| s4_correlated_failures  | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |        -3151     |           -3956     |           -2347     |       -0.805   | 0.00573   |
| s4_correlated_failures  | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |        -4082     |           -5175     |           -3094     |       -0.609   | 4.82e-12  |
| s4_correlated_failures  | cmp:F_oracle-vs-E_joint_adaptive                  |           -6.315 |             -17.8   |               4.873 |       -0.00188 | 1         |
| s5_high_errors          | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |         -571.2   |           -1354     |             226.8   |        0.037   | 0.00432   |
| s5_high_errors          | cmp:D_reliability_adaptive-vs-B_reliability_fixed |          719.1   |            -147.7   |            1605     |        0.456   | 1         |
| s5_high_errors          | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |        -8734     |           -9565     |           -7799     |       -0.988   | 0.947     |
| s5_high_errors          | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |        -2950     |           -3907     |           -2021     |       -0.988   | 1         |
| s5_high_errors          | cmp:F_oracle-vs-E_joint_adaptive                  |            0     |               0     |               0     |        0       | 1         |
| s6_tight_budget         | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |         -937.8   |           -1188     |            -675.7   |       -0.394   | 2.71e-07  |
| s6_tight_budget         | cmp:D_reliability_adaptive-vs-B_reliability_fixed |         -283.2   |            -509.7   |             -66.2   |       -0.155   | 1         |
| s6_tight_budget         | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |        -2698     |           -2960     |           -2416     |       -0.991   | 0.00454   |
| s6_tight_budget         | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |        -1043     |           -1308     |            -772.8   |       -0.94    | 1         |
| s6_tight_budget         | cmp:F_oracle-vs-E_joint_adaptive                  |            0     |               0     |               0     |        0       | 1         |
| s7_channel_change       | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |         9621     |            8457     |           10765     |        0.387   | 2.04e-103 |
| s7_channel_change       | cmp:D_reliability_adaptive-vs-B_reliability_fixed |        11727     |           10414     |           13168     |        0.0026  | 0         |
| s7_channel_change       | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |         1451     |             184.8   |            2650     |        0.68    | 9.07e-51  |
| s7_channel_change       | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |        -2635     |           -3974     |           -1322     |        0.436   | 5.48e-145 |
| s7_channel_change       | cmp:F_oracle-vs-E_joint_adaptive                  |            0     |               0     |               0     |        0       | 1         |
| s9_uncertain_estimates  | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |         1011     |             724     |            1328     |        0.47    | 0.648     |
| s9_uncertain_estimates  | cmp:D_reliability_adaptive-vs-B_reliability_fixed |         -586     |            -796.4   |            -387.7   |       -0.474   | 1         |
| s9_uncertain_estimates  | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |        -4264     |           -4536     |           -3927     |       -0.986   | 0.00391   |
| s9_uncertain_estimates  | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |         -461.7   |            -702.6   |            -198.2   |       -0.669   | 1         |
| s9_uncertain_estimates  | cmp:F_oracle-vs-E_joint_adaptive                  |           17     |               5.288 |              28.4   |        0.024   | 1         |
| s10_mixed               | cmp:C_uniform_adaptive-vs-A_uniform_fixed         |         -213.6   |           -6335     |            6165     |        0.461   | 0.39      |
| s10_mixed               | cmp:D_reliability_adaptive-vs-B_reliability_fixed |         -377.8   |           -9159     |            8835     |       -0.018   | 1         |
| s10_mixed               | cmp:E_joint_adaptive-vs-C_uniform_adaptive        |         2178     |           -4083     |            8548     |        0.044   | 1         |
| s10_mixed               | cmp:E_joint_adaptive-vs-D_reliability_adaptive    |        -2925     |          -10783     |            4664     |        0.867   | 9.84e-05  |
| s10_mixed               | cmp:F_oracle-vs-E_joint_adaptive                  |          579.2   |             464.6   |             702.2   |        0.029   | 1         |

## H1 - predictive accuracy of channel models

_source: `h1_model_accuracy.csv`_

| model                     |   mae_p |   max_abs_err |   B99_pred |   B99_mc |   B99_rel_err |
|:--------------------------|--------:|--------------:|-----------:|---------:|--------------:|
| heterogeneous(full)       | 0.00674 |         0.015 |     73,304 |   95,451 |        -0.232 |
| homogeneous-avg           | 0.073   |         0.168 |     79,569 |   95,451 |        -0.166 |
| poisson-no-overdispersion | 0.36    |         0.892 |     11,900 |   95,451 |        -0.875 |
| no-failure-model          | 0.067   |         0.148 |     24,059 |   95,451 |        -0.748 |

## H5 - robustness to mis-estimated channel parameters

_source: `h5_robustness.csv`_

| scenario               | controller                     | cap@0.99       | mean_spend@0.99   | cap@0.999      | mean_spend@0.999   |   mean_cost_success |   p_success_uncapped |
|:-----------------------|:-------------------------------|:---------------|:------------------|:---------------|:-------------------|--------------------:|---------------------:|
| s7_channel_change      | point-estimate (prior_sd=0.02) | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               39002 |                0.743 |
| s7_channel_change      | bayes (prior_sd=0.5)           | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               22483 |                0.974 |
| s7_channel_change      | thompson (prior_sd=0.5)        | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               16625 |                0.78  |
| s7_channel_change      | bayes (prior_sd=1.0)           | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               12428 |                0.977 |
| s7_channel_change      | oracle-params                  | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               11072 |                0.977 |
| s7_channel_change      | static-optimised-on-belief     | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                 inf |                0     |
| s9_uncertain_estimates | point-estimate (prior_sd=0.02) | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                5150 |                0.965 |
| s9_uncertain_estimates | bayes (prior_sd=0.5)           | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                5074 |                0.965 |
| s9_uncertain_estimates | thompson (prior_sd=0.5)        | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                5056 |                0.965 |
| s9_uncertain_estimates | bayes (prior_sd=1.0)           | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                5070 |                0.965 |
| s9_uncertain_estimates | oracle-params                  | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                5033 |                0.965 |
| s9_uncertain_estimates | static-optimised-on-belief     | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                 inf |                0     |
| s10_mixed              | point-estimate (prior_sd=0.02) | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               10073 |                0.92  |
| s10_mixed              | bayes (prior_sd=0.5)           | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               10257 |                0.92  |
| s10_mixed              | thompson (prior_sd=0.5)        | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                9885 |                0.92  |
| s10_mixed              | bayes (prior_sd=1.0)           | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               10790 |                0.92  |
| s10_mixed              | oracle-params                  | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |               10066 |                0.92  |
| s10_mixed              | static-optimised-on-belief     | inf [inf, inf] | inf [inf, inf]    | inf [inf, inf] | inf [inf, inf]     |                 inf |                0     |

## Ablations

_source: `ablations.csv`_

| scenario               | ablation                                                            | variant     | cap@0.99                | mean_spend@0.99      | cap@0.999         | mean_spend@0.999     |   mean_cost_success |
|:-----------------------|:--------------------------------------------------------------------|:------------|:------------------------|:---------------------|:------------------|:---------------------|--------------------:|
| s2_heterogeneous       | full (E)                                                            | failure1    | 6,000 [inf, inf]        | 4,347 [4,328, 4,364] | inf [inf, inf]    | inf [inf, inf]       |               4,331 |
| s2_heterogeneous       | A: -reliability-aware placement                                     |             | 13,560 [12,590, 15,455] | 7,232 [7,169, 7,294] | inf [inf, inf]    | inf [inf, inf]       |               7,222 |
| s2_heterogeneous       | B: -adaptive sequencing (static optimised weights)                  |             | 5,739 [5,613, 5,835]    | 5,739 [5,613, 5,835] | inf [inf, inf]    | inf [inf, inf]       |               4,301 |
| s2_heterogeneous       | E: fixed stopping (controller's mean allocation, committed)         |             | 5,739 [5,613, 5,835]    | 5,739 [5,613, 5,835] | inf [inf, inf]    | inf [inf, inf]       |               4,301 |
| s2_heterogeneous       | F: perfect channel information (oracle)                             |             | 6,000 [inf, inf]        | 4,347 [4,328, 4,364] | inf [inf, inf]    | inf [inf, inf]       |               4,331 |
| s2_heterogeneous       | G: independent instead of joint optimisation                        |             | 12,895 [11,990, 13,805] | 7,176 [7,118, 7,234] | inf [inf, inf]    | inf [inf, inf]       |               7,177 |
| s2_heterogeneous       | joint with static-budget proxy (v1) instead of adaptive-in-the-loop | failure0    | 8,910 [8,330, 9,605]    | 6,195 [6,165, 6,226] | 15,490 [inf, inf] | 6,214 [6,180, 6,249] |               6,205 |
| s2_heterogeneous       | adaptive stopping only (uniform round-robin allocation)             |             | 6,000 [inf, inf]        | 4,347 [4,328, 4,364] | inf [inf, inf]    | inf [inf, inf]       |               4,331 |
| s2_heterogeneous       | coverage-greedy instead of pivotal                                  |             | 6,000 [inf, inf]        | 4,347 [4,328, 4,364] | inf [inf, inf]    | inf [inf, inf]       |               4,331 |
| s4_correlated_failures | full (E)                                                            | failure1    | 20,635 [inf, inf]       | 5,202 [5,106, 5,309] | inf [inf, inf]    | inf [inf, inf]       |               5,046 |
| s4_correlated_failures | A: -reliability-aware placement                                     |             | inf [inf, inf]          | inf [inf, inf]       | inf [inf, inf]    | inf [inf, inf]       |               6,952 |
| s4_correlated_failures | B: -adaptive sequencing (static optimised weights)                  |             | inf [inf, inf]          | inf [inf, inf]       | inf [inf, inf]    | inf [inf, inf]       |               8,426 |
| s4_correlated_failures | D: -correlated failure modelling                                    | reliability | inf [inf, inf]          | inf [inf, inf]       | inf [inf, inf]    | inf [inf, inf]       |               5,360 |
| s4_correlated_failures | E: fixed stopping (controller's mean allocation, committed)         |             | inf [inf, inf]          | inf [inf, inf]       | inf [inf, inf]    | inf [inf, inf]       |               5,834 |
| s4_correlated_failures | F: perfect channel information (oracle)                             |             | 20,185 [inf, inf]       | 5,184 [5,091, 5,288] | inf [inf, inf]    | inf [inf, inf]       |               5,033 |
| s4_correlated_failures | G: independent instead of joint optimisation                        |             | inf [inf, inf]          | inf [inf, inf]       | inf [inf, inf]    | inf [inf, inf]       |               6,931 |
| s4_correlated_failures | joint with static-budget proxy (v1) instead of adaptive-in-the-loop |             |                         |                      |                   |                      |                     |
| s4_correlated_failures | adaptive stopping only (uniform round-robin allocation)             |             | 84,920 [inf, inf]       | 9,323 [8,876, 9,795] | inf [inf, inf]    | inf [inf, inf]       |               8,559 |
| s4_correlated_failures | coverage-greedy instead of pivotal                                  |             | 20,725 [inf, inf]       | 5,124 [5,027, 5,232] | inf [inf, inf]    | inf [inf, inf]       |               4,966 |

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
