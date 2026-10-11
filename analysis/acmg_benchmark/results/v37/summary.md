# Run comparison

Baseline `v36`.
Labels: ClinVar 2-star+ classifications (`clinvar_2star_labels.tsv`), a reference and not ground truth.

| run | classified | exact | same-direction | opposite | P same | LP same | VUS exact | LB same | B same | P/LP called B/LB | B/LB called P/LP | reviewed rows opposite |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| v36 | 673,278 | 425,328 | 522,542 | 59 | 51,796 | 6,750 | 285,985 | 66,046 | 111,965 | 30 | 29 | 59 |
| v37 | 673,278 | 425,869 | 522,886 | 49 | 51,855 | 6,755 | 287,070 | 65,461 | 111,745 | 24 | 25 | 49 |
| v37_bs2off | 673,278 | 425,289 | 522,388 | 49 | 51,793 | 6,747 | 286,168 | 65,825 | 111,855 | 24 | 25 | 49 |
| v37_bs2prev_3e-5 | 673,278 | 425,773 | 522,873 | 49 | 51,839 | 6,753 | 286,834 | 65,610 | 111,837 | 24 | 25 | 49 |
| v37_bs2prev_1e-4 | 673,278 | 425,877 | 522,938 | 49 | 51,852 | 6,755 | 287,028 | 65,500 | 111,803 | 24 | 25 | 49 |
| v37_bs2prev_1e-3 | 673,278 | 425,753 | 522,819 | 49 | 51,855 | 6,755 | 287,079 | 65,437 | 111,693 | 24 | 25 | 49 |
| v37_bs2prev_3e-3 | 673,278 | 425,653 | 522,775 | 49 | 51,855 | 6,755 | 287,086 | 65,427 | 111,652 | 24 | 25 | 49 |
| v37_bs2ac_20 | 673,278 | 425,721 | 522,827 | 49 | 51,833 | 6,753 | 286,741 | 65,656 | 111,844 | 24 | 25 | 49 |
| v37_bs2ac_100 | 673,278 | 425,872 | 522,918 | 49 | 51,853 | 6,755 | 287,036 | 65,488 | 111,786 | 24 | 25 | 49 |
| v37_bs2ac_1000 | 673,278 | 425,726 | 522,798 | 49 | 51,855 | 6,755 | 287,084 | 65,434 | 111,670 | 24 | 25 | 49 |
| v37_mpc | 673,278 | 425,645 | 522,484 | 47 | 51,689 | 6,655 | 287,063 | 65,355 | 111,722 | 23 | 24 | 47 |
| v37_rank | 673,278 | 423,481 | 519,120 | 452 | 52,821 | 6,895 | 286,479 | 63,171 | 109,754 | 12 | 440 | 37 |
| v37_adtag | 673,278 | 429,678 | 525,698 | 50 | 51,851 | 6,755 | 286,881 | 66,176 | 114,035 | 26 | 24 | 48 |
| v37_moi_global | 673,278 | 429,693 | 525,695 | 50 | 51,852 | 6,755 | 286,877 | 66,176 | 114,035 | 26 | 24 | 48 |

## `v36` to `v37`

2,052 calls changed: 1,202 toward ClinVar's label, 776 away, 74 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 1,508 | 1,003 | 432 | 73 | 0 | 70 | 963 | 475 | 0 | 0 |
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 16 | 7 | 9 | 0 | 0 | 0 | 7 | 9 | 0 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| BP4,BS2 -> - | 13 | 13 | 0 | 0 | 0 | 1 | 12 | 0 | 1 | 0 |
| BP7,BS2 -> - | 10 | 4 | 6 | 0 | 0 | 0 | 4 | 6 | 0 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_bs2off`

534 calls changed: 195 toward ClinVar's label, 338 away, 1 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 29 | 20 | 9 | 0 | 0 | 1 | 19 | 9 | 1 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_bs2prev_3e-5`

1,535 calls changed: 938 toward ClinVar's label, 559 away, 38 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 998 | 741 | 220 | 37 | 0 | 52 | 714 | 232 | 0 | 0 |
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 18 | 9 | 9 | 0 | 0 | 0 | 9 | 9 | 0 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| BP4,BS2 -> - | 11 | 11 | 0 | 0 | 0 | 1 | 10 | 0 | 1 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| BP7,BS2 -> - | 3 | 2 | 1 | 0 | 0 | 0 | 2 | 1 | 0 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_bs2prev_1e-4`

1,901 calls changed: 1,154 toward ClinVar's label, 682 away, 65 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 1,359 | 955 | 340 | 64 | 0 | 67 | 918 | 374 | 0 | 0 |
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 17 | 8 | 9 | 0 | 0 | 0 | 8 | 9 | 0 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| BP4,BS2 -> - | 12 | 12 | 0 | 0 | 0 | 1 | 11 | 0 | 1 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| BP7,BS2 -> - | 8 | 4 | 4 | 0 | 0 | 0 | 4 | 4 | 0 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_bs2prev_1e-3`

2,240 calls changed: 1,225 toward ClinVar's label, 933 away, 82 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 1,694 | 1,026 | 587 | 81 | 0 | 70 | 973 | 651 | 0 | 0 |
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 16 | 7 | 9 | 0 | 0 | 0 | 7 | 9 | 0 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| BP4,BS2 -> - | 13 | 13 | 0 | 0 | 0 | 1 | 12 | 0 | 1 | 0 |
| BP7,BS2 -> - | 11 | 4 | 7 | 0 | 0 | 0 | 4 | 7 | 0 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| BS2 -> BP4 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_bs2prev_3e-3`

2,456 calls changed: 1,271 toward ClinVar's label, 1,095 away, 90 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 1,909 | 1,072 | 748 | 89 | 0 | 70 | 980 | 859 | 0 | 0 |
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 16 | 7 | 9 | 0 | 0 | 0 | 7 | 9 | 0 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| BP4,BS2 -> - | 13 | 13 | 0 | 0 | 0 | 1 | 12 | 0 | 1 | 0 |
| BP7,BS2 -> - | 12 | 4 | 8 | 0 | 0 | 0 | 4 | 8 | 0 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| BS2 -> BP4 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_bs2ac_20`

1,368 calls changed: 831 toward ClinVar's label, 509 away, 28 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 831 | 634 | 170 | 27 | 0 | 46 | 605 | 180 | 0 | 0 |
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 20 | 11 | 9 | 0 | 0 | 0 | 11 | 9 | 0 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| BP4,BS2 -> - | 9 | 9 | 0 | 0 | 0 | 1 | 8 | 0 | 1 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| BP7,BS2 -> - | 3 | 2 | 1 | 0 | 0 | 0 | 2 | 1 | 0 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_bs2ac_100`

1,944 calls changed: 1,165 toward ClinVar's label, 711 away, 68 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 1,401 | 966 | 368 | 67 | 0 | 68 | 929 | 404 | 0 | 0 |
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 17 | 8 | 9 | 0 | 0 | 0 | 8 | 9 | 0 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| BP4,BS2 -> - | 12 | 12 | 0 | 0 | 0 | 1 | 11 | 0 | 1 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| BP7,BS2 -> - | 9 | 4 | 5 | 0 | 0 | 0 | 4 | 5 | 0 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_bs2ac_1000`

2,331 calls changed: 1,248 toward ClinVar's label, 998 away, 85 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 10 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 1,785 | 1,049 | 652 | 84 | 0 | 70 | 978 | 737 | 0 | 0 |
| BP7 -> - | 465 | 141 | 324 | 0 | 0 | 2 | 139 | 324 | 2 | 0 |
| BP4 -> - | 16 | 7 | 9 | 0 | 0 | 0 | 7 | 9 | 0 | 0 |
| PVS1 -> - | 14 | 10 | 4 | 0 | 0 | 4 | 9 | 1 | 1 | 0 |
| BP4,BS2 -> - | 13 | 13 | 0 | 0 | 0 | 1 | 12 | 0 | 1 | 0 |
| BP7,BS2 -> - | 11 | 4 | 7 | 0 | 0 | 0 | 4 | 7 | 0 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| PVS1 -> PP3 | 8 | 8 | 0 | 0 | 0 | 0 | 6 | 2 | 1 | 0 |
| PS1,PVS1 -> - | 2 | 1 | 0 | 1 | 0 | 1 | 1 | 0 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| BS2 -> BP4 | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_mpc`

2,890 calls changed: 1,418 toward ClinVar's label, 1,386 away, 86 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 12 before, 0 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BS2 -> - | 1,508 | 1,003 | 432 | 73 | 0 | 70 | 963 | 475 | 0 | 0 |
| BP7 -> - | 673 | 150 | 523 | 0 | 0 | 2 | 148 | 523 | 2 | 0 |
| PVS1 -> - | 137 | 14 | 123 | 0 | 0 | 123 | 13 | 1 | 1 | 0 |
| PVS1 -> BP4 | 86 | 11 | 75 | 0 | 0 | 75 | 9 | 2 | 2 | 0 |
| - -> BP7 | 70 | 68 | 2 | 0 | 0 | 0 | 2 | 68 | 0 | 0 |
| - -> PVS1 | 39 | 25 | 14 | 0 | 0 | 25 | 12 | 2 | 0 | 0 |
| PVS1 -> PP3 | 32 | 10 | 22 | 0 | 0 | 22 | 8 | 2 | 1 | 0 |
| PM1,PM5,PP2 -> - | 21 | 1 | 20 | 0 | 0 | 20 | 1 | 0 | 0 | 0 |
| BP4 -> - | 18 | 8 | 10 | 0 | 0 | 0 | 8 | 10 | 0 | 0 |
| PM1 -> - | 18 | 9 | 9 | 0 | 0 | 9 | 8 | 1 | 1 | 0 |
| - -> - | 17 | 5 | 12 | 0 | 0 | 2 | 0 | 15 | 0 | 0 |
| PM5 -> PP3 | 17 | 12 | 5 | 0 | 0 | 11 | 6 | 0 | 0 | 0 |
| BP4 -> PVS1 | 16 | 15 | 1 | 0 | 0 | 15 | 0 | 1 | 0 | 0 |
| PM1,PM5 -> PP3 | 15 | 1 | 14 | 0 | 0 | 14 | 1 | 0 | 0 | 0 |
| BP4,BS2 -> - | 13 | 13 | 0 | 0 | 0 | 1 | 12 | 0 | 1 | 0 |
| PS1,PVS1 -> - | 13 | 2 | 6 | 5 | 0 | 11 | 2 | 0 | 0 | 0 |
| PS1,PVS1 -> PP3 | 13 | 1 | 6 | 6 | 0 | 12 | 1 | 0 | 0 | 0 |
| - -> BP4 | 12 | 10 | 2 | 0 | 0 | 0 | 2 | 10 | 0 | 0 |
| PP2,PP3 -> - | 11 | 8 | 3 | 0 | 0 | 3 | 8 | 0 | 0 | 0 |
| PM1,PM5,PS1 -> PP3 | 11 | 0 | 11 | 0 | 0 | 11 | 0 | 0 | 0 | 0 |
| - -> PP2 | 11 | 9 | 2 | 0 | 0 | 0 | 0 | 11 | 0 | 0 |
| BP7,BS2 -> - | 10 | 4 | 6 | 0 | 0 | 0 | 4 | 6 | 0 | 0 |
| PP3 -> PM5 | 8 | 5 | 3 | 0 | 0 | 3 | 5 | 0 | 0 | 0 |
| PP3 -> PVS1 | 7 | 6 | 1 | 0 | 0 | 6 | 1 | 0 | 0 | 0 |
| - -> PM1 | 6 | 1 | 5 | 0 | 0 | 1 | 5 | 0 | 0 | 0 |

The 25 largest of 82 groups; all of them, and the totals above, are in `attribution_v36_to_v37_mpc.tsv`.

## `v36` to `v37_rank`

13,303 calls changed: 4,951 toward ClinVar's label, 8,244 away, 108 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 22 before, 415 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| BP7 -> - | 3,538 | 336 | 3,202 | 0 | 0 | 4 | 247 | 3,287 | 4 | 0 |
| BP4,BP7 -> - | 2,052 | 11 | 2,041 | 0 | 0 | 4 | 7 | 2,041 | 4 | 0 |
| - -> BP7 | 1,921 | 1,911 | 10 | 0 | 0 | 0 | 10 | 1,911 | 0 | 0 |
| BS2 -> - | 1,497 | 998 | 425 | 74 | 0 | 70 | 958 | 469 | 0 | 0 |
| - -> PVS1 | 863 | 519 | 344 | 0 | 0 | 519 | 183 | 161 | 0 | 144 |
| BP4 -> PVS1 | 625 | 299 | 326 | 0 | 0 | 297 | 189 | 139 | 0 | 78 |
| BP4,BP7 -> PP2 | 322 | 0 | 322 | 0 | 0 | 0 | 0 | 322 | 0 | 0 |
| PP3 -> PVS1 | 283 | 137 | 146 | 0 | 0 | 137 | 123 | 23 | 0 | 12 |
| - -> - | 259 | 227 | 31 | 1 | 0 | 3 | 7 | 249 | 0 | 0 |
| BP4 -> - | 245 | 10 | 233 | 2 | 0 | 0 | 10 | 235 | 0 | 0 |
| BP4,BP7 -> PVS1 | 219 | 12 | 202 | 5 | 0 | 6 | 11 | 202 | 6 | 114 |
| - -> BP4 | 91 | 77 | 3 | 11 | 0 | 0 | 3 | 88 | 0 | 0 |
| BP7 -> PP2 | 71 | 15 | 56 | 0 | 0 | 0 | 3 | 68 | 0 | 0 |
| - -> PS1,PVS1 | 69 | 66 | 1 | 2 | 0 | 68 | 1 | 0 | 0 | 0 |
| BP4,BP7 -> PM1 | 65 | 0 | 65 | 0 | 0 | 0 | 0 | 65 | 0 | 0 |
| BP4,BP7 -> PM4 | 59 | 4 | 55 | 0 | 0 | 0 | 4 | 55 | 0 | 0 |
| PVS1 -> - | 58 | 12 | 46 | 0 | 0 | 46 | 11 | 1 | 1 | 0 |
| - -> PM1,PM5,PP2 | 55 | 27 | 28 | 0 | 0 | 27 | 20 | 8 | 0 | 8 |
| - -> BP1 | 50 | 13 | 37 | 0 | 0 | 0 | 36 | 14 | 0 | 0 |
| BP7,BS2 -> - | 48 | 8 | 40 | 0 | 0 | 0 | 7 | 41 | 0 | 0 |
| BP4,BP7 -> PM1,PM5 | 39 | 0 | 39 | 0 | 0 | 0 | 0 | 39 | 0 | 0 |
| PP2 -> PVS1 | 36 | 0 | 36 | 0 | 0 | 0 | 36 | 0 | 0 | 0 |
| BP7 -> PM2 | 34 | 0 | 34 | 0 | 0 | 0 | 0 | 34 | 0 | 0 |
| - -> PP2,PP3 | 34 | 13 | 21 | 0 | 0 | 13 | 21 | 0 | 0 | 0 |
| - -> PP2 | 32 | 26 | 6 | 0 | 0 | 0 | 0 | 32 | 0 | 0 |

The 25 largest of 160 groups; all of them, and the totals above, are in `attribution_v36_to_v37_rank.tsv`.

## `v36` to `v37_adtag`

10,040 calls changed: 7,353 toward ClinVar's label, 2,195 away, 492 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 11 before, 2 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| - -> BS2 | 7,960 | 6,125 | 1,417 | 418 | 0 | 6 | 212 | 7,742 | 1 | 2 |
| BS2 -> - | 1,530 | 1,008 | 449 | 73 | 0 | 70 | 967 | 493 | 0 | 0 |
| BP7 -> - | 451 | 141 | 310 | 0 | 0 | 2 | 139 | 310 | 2 | 0 |
| BP7 -> BS2 | 19 | 19 | 0 | 0 | 0 | 0 | 0 | 19 | 0 | 0 |
| PVS1 -> - | 16 | 11 | 4 | 1 | 0 | 5 | 10 | 1 | 1 | 0 |
| BP4 -> - | 15 | 7 | 8 | 0 | 0 | 0 | 7 | 8 | 0 | 0 |
| BP4,BS2 -> - | 13 | 13 | 0 | 0 | 0 | 1 | 12 | 0 | 1 | 0 |
| BP7,BS2 -> - | 10 | 4 | 6 | 0 | 0 | 0 | 4 | 6 | 0 | 0 |
| PVS1 -> BP4 | 10 | 9 | 1 | 0 | 0 | 1 | 7 | 2 | 2 | 0 |
| PVS1 -> PP3 | 9 | 9 | 0 | 0 | 0 | 0 | 7 | 2 | 1 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| BP4 -> BS2 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BP3,PM4 -> BP4,BS2 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |

## `v36` to `v37_moi_global`

10,066 calls changed: 7,374 toward ClinVar's label, 2,200 away, 492 the same distance from it, 0 not scored (a NoCall on either side, or no ClinVar label).
Opposite-direction 11 before, 2 after, among the changed calls.

| criteria lost to gained | changed | toward | away | level | not scored | ClinVar P/LP | ClinVar VUS | ClinVar B/LB | opposite before | after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| - -> BS2 | 7,960 | 6,125 | 1,417 | 418 | 0 | 6 | 212 | 7,742 | 1 | 2 |
| BS2 -> - | 1,530 | 1,008 | 449 | 73 | 0 | 70 | 967 | 493 | 0 | 0 |
| BP7 -> - | 443 | 135 | 308 | 0 | 0 | 2 | 133 | 308 | 2 | 0 |
| - -> PM2 | 26 | 21 | 5 | 0 | 0 | 22 | 4 | 0 | 0 | 0 |
| BP7 -> BS2 | 19 | 19 | 0 | 0 | 0 | 0 | 0 | 19 | 0 | 0 |
| PVS1 -> - | 16 | 11 | 4 | 1 | 0 | 5 | 10 | 1 | 1 | 0 |
| BP4 -> - | 15 | 7 | 8 | 0 | 0 | 0 | 7 | 8 | 0 | 0 |
| BP4,BS2 -> - | 13 | 13 | 0 | 0 | 0 | 1 | 12 | 0 | 1 | 0 |
| BP7,BS2 -> - | 10 | 4 | 6 | 0 | 0 | 0 | 4 | 6 | 0 | 0 |
| PVS1 -> PP3 | 9 | 9 | 0 | 0 | 0 | 0 | 7 | 2 | 1 | 0 |
| PVS1 -> BP4 | 9 | 8 | 1 | 0 | 0 | 1 | 6 | 2 | 2 | 0 |
| BP7 -> PM2 | 8 | 6 | 2 | 0 | 0 | 0 | 6 | 2 | 0 | 0 |
| BS1 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| PM4,PVS1 -> - | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| BP4 -> BS2 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| - -> BP4 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| PVS1 -> BP4,PM2 | 1 | 1 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| BP3,PM4 -> BP4,BS2 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| BS1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
| BA1,BS2 -> - | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 |
