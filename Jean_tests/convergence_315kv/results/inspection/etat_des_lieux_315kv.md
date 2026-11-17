# État des lieux — réseau 315 kV « convergé » d'Etienne

Fichier inspecté : `Jean_tests\convergence_315kv\results\elec_315kv_shunt_pf.nc` (reproduction de la séquence d'Etienne, commit 2ef0c21d, PyPSA 1.3.0). Référence avant corrections : `network/networks_current/elec_solved.nc`. Généré par `inspect_315kv.py`.

Chaîne de traitement : `elec_solved.nc` → `correct_sending_end_circuits.py` → `apply_series_compensation.py` (70 %, 15 lignes de Nathan + 4 ajoutées par Etienne) → `add_shunt_impedance.py` (réacteurs ×1,0, condensateurs ×0,7) → AC PF (`run_pf.py`, `--pv-min-capacity 0 --pv-max-load-ratio 1`).

## 0. Synthèse

- Snapshots : 168 ; Newton converge sur **162** ; divergés : 6 (01-03 09h, 01-03 16h, 01-03 19h, 01-03 22h, 01-03 23h, 01-04 16h).
- Heures convergées avec tous les bus AC dans 0.9-1.1 pu : **0** ; dans 0.95-1.05 pu : **0**.
- Tension sur les heures convergées : min 0.910 pu, max 1.615 pu ; bus AC qui sortent au moins une fois de 0.9-1.1 : 96 / 200.
- Demande : moyenne 28871 MW, pointe 34959 MW (01-03 17h) ; moyenne des heures divergées 33020 MW, des heures convergées 28717 MW.

## 1. Inventaire du réseau

**Bus**

| carrier | v_nom | bus |
|---|---|---|
| AC | 315 | 140 |
| AC | 345 | 2 |
| AC | 735 | 56 |
| AC | 765 | 2 |
| DC | 450 | 5 |

**Lignes AC**

| v_nom | lignes | km | km_max | circuits |
|---|---|---|---|---|
| 315 | 167 | 8216 | 410 | 287 |
| 345 | 1 | 10 | 10 | 2 |
| 735 | 107 | 20164 | 623 | 117 |
| 765 | 1 | 79 | 79 | 1 |

**Transformateurs** : 22, s_nom = [100000.0] MVA (valeur plate de `run_lopf_main_island.py`), x = [5.0] (pu sur s_nom), r = x/30 (générique, run_pf).

**Liens (HVDC / B2B)**

| name | bus0 | bus1 | carrier | p_nom |
|---|---|---|---|---|
| 1386 | 65 | 3549 | DC | 737 |
| 1387 | 161 | 1015 | DC | 737 |
| 1390 | 3549 | 161 | DC | 737 |
| 2396 | 2206 | 1015 | DC | 737 |
| convert_25_65 | 65 | 64 | B2B | 2000 |

**Générateurs** (contrôle après `prepare_for_ac_pf`)

| carrier | control | nb | p_nom_MW | p_moy_conv_MW |
|---|---|---|---|---|
| AC | Slack | 2 | 7722 | 1940 |
| OCGT | PV | 1 | 411 | 16 |
| hydro | PV | 14 | 20528 | 0 |
| load_shedding | PQ | 185 | 34451 | 0 |
| load_shedding | Slack | 4 | 516 | 0 |
| onwind | PQ | 21 | 1436 | 231 |
| onwind | PV | 12 | 2059 | 458 |
| ror | PQ | 19 | 2492 | 1203 |
| ror | PV | 16 | 11078 | 3463 |
| solar | PQ | 2 | 4 | 0 |
| solar | PV | 1 | 6 | 0 |

**Stockage (réservoirs hydro)**

| carrier | nb | p_nom_MW | p_moy_conv_MW |
|---|---|---|---|
| hydro | 22 | 28583 | 23345 |

**Charges** : 647 sur 189 bus ; P moyen 28871 MW, Q = P·tan(acos 0,95) (générique, run_pf).

**Slack par sous-réseau**

| name | bus | carrier | p_nom |
|---|---|---|---|
| 1015 load_shedding | 1015 | load_shedding | 207 |
| 161 load_shedding | 161 | load_shedding | 23 |
| 3549 load_shedding | 3549 | load_shedding | 88 |
| 2206 load_shedding | 2206 | load_shedding | 198 |
| 339 slack-placeholder | 339 | AC | 7722 |
| 65 pf-slack-placeholder | 65 | AC | 0 |

## 2. Ce que la chaîne a modifié par rapport à `elec_solved.nc`

20 lignes modifiées (circuits, x ou b). Rapport x_après / x_avant : 0,3 = compensation série seule (70 %) ; 1/3 = passage de 1 à 3 circuits ; ≈ 0,1 = les deux.

| name | bus0 | bus1 | v_nom | km | circ_avant | circ_apres | x_ratio | s_nom_avant | s_nom_apres |
|---|---|---|---|---|---|---|---|---|---|
| 162 | 340 | 457 | 735 | 457 | 1 | 3 | 0.1 | 5.68e+03 | 1.7e+04 |
| 2408 | 114 | 1291 | 735 | 501 | 1 | 3 | 0.1 | 5.24e+03 | 1.57e+04 |
| 1938 | 66 | 457 | 735 | 433 | 1 | 1.5 | 0.2 | 5.96e+03 | 8.94e+03 |
| 156 | 66 | 457 | 735 | 433 | 1 | 1.5 | 0.2 | 5.96e+03 | 8.94e+03 |
| 1049 | 114 | 148 | 735 | 563 | 1 | 1.5 | 0.2 | 4.75e+03 | 7.12e+03 |
| 1448 | 114 | 148 | 735 | 568 | 1 | 1.5 | 0.2 | 4.72e+03 | 7.07e+03 |
| 880 | 312 | 308 | 735 | 433 | 1 | 1 | 0.3 | 5.95e+03 | 5.95e+03 |
| 163 | 457 | 312 | 735 | 345 | 1 | 1 | 0.3 | 7.34e+03 | 7.34e+03 |
| 1536 | 308 | 2944 | 735 | 256 | 1 | 1 | 0.3 | 9.74e+03 | 9.74e+03 |
| 1538 | 308 | 1081 | 735 | 266 | 1 | 1 | 0.3 | 9.39e+03 | 9.39e+03 |
| 1935 | 312 | 308 | 735 | 410 | 1 | 1 | 0.3 | 6.26e+03 | 6.26e+03 |
| 1653 | 345 | 310 | 735 | 372 | 1 | 1 | 0.3 | 6.83e+03 | 6.83e+03 |
| 2386 | 457 | 312 | 735 | 378 | 1 | 1 | 0.3 | 6.74e+03 | 6.74e+03 |
| 1936 | 308 | 136 | 735 | 348 | 1 | 1 | 0.3 | 7.28e+03 | 7.28e+03 |
| 2387 | 345 | 310 | 735 | 397 | 1 | 1 | 0.3 | 6.45e+03 | 6.45e+03 |
| 161 | 451 | 345 | 735 | 437 | 1 | 1 | 0.3 | 5.91e+03 | 5.91e+03 |
| 3535 | 2944 | 286 | 735 | 471 | 1 | 1 | 0.3 | 5.52e+03 | 5.52e+03 |
| 2410 | 312 | 308 | 735 | 410 | 1 | 1 | 0.3 | 6.26e+03 | 6.26e+03 |
| 1655 | 457 | 312 | 735 | 345 | 1 | 1 | 0.3 | 7.34e+03 | 7.34e+03 |
| 3756 | 651 | 25 | 735 | 623 | 1 | 3 | 0.333 | 4.38e+03 | 1.31e+04 |

**Shunts fixes ajoutés** : 6 condensateurs (+3423 MVAr à 1 pu), 4 réacteurs (-255 MVAr). Le dimensionnement a ignoré 120 bus qui sortent de la bande dans les deux sens au cours de la semaine (un shunt fixe ne peut corriger qu'un sens).

| name | bus | v_nom | Q_MVAr_1pu | Q_moy_MVAr | V_min | V_max |
|---|---|---|---|---|---|---|
| 1638 shunt-reactor | 1638 | 315 | -33.3 | -33.8 | 0.989 | 1.03 |
| 1762 shunt-reactor | 1762 | 315 | -10.9 | -11.5 | 1.02 | 1.03 |
| 1772 shunt-reactor | 1772 | 315 | -11.9 | -12.6 | 1.02 | 1.03 |
| 2913 shunt-reactor | 2913 | 315 | -199 | -202 | 0.991 | 1.03 |
| 1092 shunt-capacitor | 1092 | 315 | 5.03 | 5.11 | 1.01 | 1.01 |
| 2892 shunt-capacitor | 2892 | 315 | 11 | 12.2 | 0.963 | 1.11 |
| 3836 shunt-capacitor | 3836 | 315 | 154 | 152 | 0.95 | 1.02 |
| 3859 shunt-capacitor | 3859 | 315 | 324 | 356 | 0.963 | 1.11 |
| 45 shunt-capacitor | 45 | 315 | 852 | 937 | 0.954 | 1.12 |
| 800 shunt-capacitor | 800 | 315 | 2.08e+03 | 2.29e+03 | 0.916 | 1.15 |

## 3. Tensions sur les heures convergées

**Par niveau de tension**

| v_nom | bus | V_min | V_max | bus_hors_bande |
|---|---|---|---|---|
| 315.000 | 140 | 0.910 | 1.615 | 47 |
| 345.000 | 2 | 1.085 | 1.141 | 2 |
| 735.000 | 56 | 1.000 | 1.338 | 47 |
| 765.000 | 2 | 1.001 | 1.093 | 0 |

**Bus les plus hauts** (V max)

| name | v_nom | controle | V_min | V_max | h_sur_110 | h_sous_090 | degre | charge_moy_MW | gen_p_nom_MW | Q_fixe_genPQ_MVAr | shunt_MVAr_1pu |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1810 | 315 | PQ | 1.05 | 1.62 | 161 | 0 | 1 | 579 | 411 | 199 | 0 |
| 1806 | 315 | PQ | 1.07 | 1.58 | 161 | 0 | 2 | 78.7 | 0 | 0 | 0 |
| 1804 | 315 | PQ | 1.12 | 1.47 | 162 | 0 | 2 | 0.423 | 0 | 0 | 0 |
| 1803 | 315 | PQ | 1.14 | 1.43 | 162 | 0 | 2 | 0.278 | 0 | 0 | 0 |
| 1800 | 315 | PQ | 1.17 | 1.39 | 162 | 0 | 2 | 0.302 | 0 | 0 | 0 |
| 311 | 315 | PQ | 1.19 | 1.35 | 162 | 0 | 2 | 0.382 | 0 | 0 | 0 |
| 312 | 735 | PQ | 1.18 | 1.34 | 162 | 0 | 8 | 0 | 0 | 0 | 0 |
| 457 | 735 | PQ | 1.2 | 1.31 | 162 | 0 | 8 | 0 | 0 | 0 | 0 |
| 469 | 735 | PQ | 1.26 | 1.3 | 162 | 0 | 4 | 0 | 0 | 0 | 0 |
| 468 | 315 | PQ | 1.26 | 1.29 | 162 | 0 | 1 | 165 | 0 | 0 | 0 |
| 308 | 735 | PQ | 1.11 | 1.29 | 162 | 0 | 6 | 1.44e+03 | 0 | 0 | 0 |
| 3319 | 735 | PQ | 1.18 | 1.28 | 162 | 0 | 2 | 1.67 | 0 | 0 | 0 |
| 345 | 735 | PQ | 1.19 | 1.28 | 162 | 0 | 7 | 0.694 | 0 | 0 | 0 |
| 1291 | 735 | PQ | 1.2 | 1.28 | 162 | 0 | 2 | 129 | 0 | 0 | 0 |
| 310 | 735 | PQ | 1.18 | 1.28 | 162 | 0 | 7 | 12.7 | 0 | 0 | 0 |
| 3582 | 735 | PQ | 1.18 | 1.26 | 162 | 0 | 2 | 20.1 | 0 | 0 | 0 |
| 66 | 735 | PQ | 1.16 | 1.25 | 162 | 0 | 5 | 0 | 0 | 0 | 0 |
| 25 | 735 | PQ | 1.16 | 1.24 | 162 | 0 | 4 | 0 | 0 | 0 | 0 |
| 1081 | 735 | PQ | 1.06 | 1.24 | 150 | 0 | 2 | 1.09e+03 | 0 | 0 | 0 |
| 456 | 315 | PQ | 1.15 | 1.24 | 162 | 0 | 2 | 0.584 | 0 | 0 | 0 |

*(180 lignes de plus dans le CSV)*

**Bus les plus bas** (V min)

| name | v_nom | controle | V_min | V_max | h_sur_110 | h_sous_090 | degre | charge_moy_MW | gen_p_nom_MW | Q_fixe_genPQ_MVAr | shunt_MVAr_1pu |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 263 | 315 | PQ | 0.91 | 1.13 | 17 | 0 | 3 | 3.62e+03 | 0 | 0 | 0 |
| 800 | 315 | PQ | 0.916 | 1.15 | 35 | 0 | 1 | 2.5e+03 | 0 | 0 | 2.08e+03 |
| 271 | 315 | PQ | 0.935 | 1.12 | 14 | 0 | 4 | 970 | 0 | 0 | 0 |
| 3836 | 315 | PQ | 0.95 | 1.02 | 0 | 0 | 3 | 486 | 0 | 0 | 154 |
| 45 | 315 | PQ | 0.954 | 1.12 | 17 | 0 | 3 | 2.81e+03 | 0 | 0 | 852 |
| 1664 | 315 | PQ | 0.956 | 1.12 | 15 | 0 | 4 | 930 | 0 | 0 | 0 |
| 194 | 315 | PQ | 0.958 | 1.12 | 23 | 0 | 10 | 7.16 | 0 | 0 | 0 |
| 1715 | 315 | PQ | 0.96 | 1.12 | 25 | 0 | 5 | 8.02 | 0 | 0 | 0 |
| 1706 | 315 | PQ | 0.963 | 1.1 | 4 | 0 | 1 | 10.3 | 0 | 0 | 0 |
| 1669 | 315 | PQ | 0.963 | 1.1 | 4 | 0 | 3 | 3.61 | 0 | 0 | 0 |
| 2892 | 315 | PQ | 0.963 | 1.11 | 12 | 0 | 1 | 35.1 | 0 | 0 | 11 |
| 3859 | 315 | PQ | 0.963 | 1.11 | 12 | 0 | 3 | 3.63 | 0 | 0 | 324 |
| 676 | 315 | PQ | 0.964 | 1.12 | 18 | 0 | 5 | 558 | 0 | 0 | 0 |
| 979 | 315 | PQ | 0.969 | 1.11 | 9 | 0 | 2 | 0.571 | 0 | 0 | 0 |
| 1716 | 315 | PQ | 0.976 | 1.13 | 32 | 0 | 3 | 10 | 0 | 0 | 0 |

*(185 lignes de plus dans le CSV)*

Par heure : le nombre de bus > 1.1 pu varie de 35 à 96 ; corrélation avec la demande : V_max -0.72, V_min -0.79, nb de bus > 1,10 -0.97 (CSV `tensions_par_heure.csv`).

**Zoom : lignes autour des 6 bus les plus hauts**. Une ligne qui transporte moins que sa puissance naturelle (SIL) produit plus de réactif qu'elle n'en consomme, et la tension monte en bout de ligne (effet Ferranti).

| name | bus0 | bus1 | v_nom | length | num_parallel | Q_charge_MVAr | SIL_MW | P0_moy_MW | Q0_moy_MVAr |
|---|---|---|---|---|---|---|---|---|---|
| 1656 | 1800 | 311 | 315 | 29.2 | 1 | 12.5 | 343 | -446 | 250 |
| 1662 | 1803 | 1800 | 315 | 40.2 | 1 | 17.1 | 343 | -444 | 248 |
| 1663 | 1804 | 1803 | 315 | 30.2 | 1 | 12.9 | 343 | -442 | 244 |
| 1664 | 1806 | 1804 | 315 | 125 | 1 | 53.5 | 343 | -436 | 219 |
| 3341 | 1806 | 1810 | 315 | 154 | 2 | 132 | 686 | 358 | -244 |

Générateurs sur ces bus (q_set = Q imposé si PQ) :

| name | bus | carrier | control | p_nom | q_set |
|---|---|---|---|---|---|
| 1788 ror-Premiere-Chute | 1810 | ror | PQ | 124 | 60.1 |
| 1783 ror-Rapides-des-Iles | 1810 | ror | PQ | 176 | 85.2 |
| 2186 ror-Rapides-des-Quinze | 1810 | ror | PQ | 111 | 53.8 |

## 4. Puissance réactive

Limite générique utilisée pour juger : |Q| ≤ 0.484 × p_nom (facteur de puissance 0,9, la même hypothèse que run_pf.py). PyPSA n'impose **aucune** limite Q en PF : les PV et le slack fournissent ce qu'il faut.

- 46 générateurs PV/Slack ; **29** dépassent la limite au moins une heure (ou ont p_nom = 0).
- 42 générateurs PQ injectent un Q **fixe** de 0,484 × p_nom, même à P = 0 : 1904 MVAr au total, toute la semaine.
- Slack `1015 load_shedding` (bus 1015, p_nom 207 MW) : P -0 à 0 MW, Q 43 à 68 MVAr.
- Slack `161 load_shedding` (bus 161, p_nom 23 MW) : P -0 à 0 MW, Q 5 à 8 MVAr.
- Slack `3549 load_shedding` (bus 3549, p_nom 88 MW) : P -0 à 0 MW, Q 18 à 29 MVAr.
- Slack `2206 load_shedding` (bus 2206, p_nom 198 MW) : P -0 à 0 MW, Q 41 à 65 MVAr.
- Slack `339 slack-placeholder` (bus 339, p_nom 7722 MW) : P 1470 à 2568 MW, Q -8028 à -4344 MVAr.
- Slack `65 pf-slack-placeholder` (bus 65, p_nom 0 MW) : P nan à nan MW, Q 0 à 0 MVAr.

**PV / Slack les plus sollicités en réactif**

| name | bus | carrier | controle | p_nom | P_moy | Q_min | Q_max | Q_max_sur_pnom |
|---|---|---|---|---|---|---|---|---|
| 0 solar-Roussillon-Solar | 1 | solar | PV | 6.02 | 0 | -4.04e+03 | -29.9 | 671 |
| 3602 ror-Trenche | 286 | ror | PV | 302 | 246 | -1.35e+04 | -8.37e+03 | 44.5 |
| 2167 onwind-Gros-Morne | 85 | onwind | PV | 224 | 73.3 | -8.06e+03 | -4.98e+03 | 36.1 |
| 330 onwind-Des Moulins | 330 | onwind | PV | 166 | 166 | -4.99e+03 | -2.82e+03 | 30.1 |
| 1121 onwind-Témiscouata 2 | 2973 | onwind | PV | 54.6 | 7.77 | -807 | -475 | 14.8 |
| 470 onwind-Mont-Louis | 84 | onwind | PV | 106 | 21.8 | 131 | 703 | 6.62 |
| 1693 onwind-Belle-Rivière | 1693 | onwind | PV | 25.4 | 1.85 | -146 | 28.3 | 5.75 |
| 2811 onwind-Seigneurie de Beaupré 3 | 2811 | onwind | PV | 149 | 11.1 | -761 | -230 | 5.12 |
| 331 ror-La-Tuque | 1237 | ror | PV | 287 | 169 | -1.42e+03 | -690 | 4.95 |
| 1187 ror-La-Gabelle | 1188 | ror | PV | 124 | 8.48 | -607 | 51.2 | 4.89 |
| 1078 ror-Carillon | 3373 | ror | PV | 762 | 486 | -2.81e+03 | 260 | 3.69 |
| 451 storage-pv-placeholder | 451 | hydro | PV | 2.42e+03 | nan | -7.77e+03 | -5.45e+03 | 3.22 |
| 559 ror-Rocher-de-Grand-Mere | 166 | ror | PV | 230 | 35.6 | -732 | -425 | 3.18 |
| 1674 ror-Jean-Lesage | 1674 | ror | PV | 1.23e+03 | 170 | -3.17e+03 | -2.22e+03 | 2.58 |
| 107 ror-Outardes-2 | 109 | ror | PV | 523 | 12.9 | -1.15e+03 | -861 | 2.2 |
| 1665 ror-Rene-Levesque | 1665 | ror | PV | 1.33e+03 | 438 | -2.14e+03 | -1.07e+03 | 1.61 |
| 2812 onwind-Rivière-du-Moulin | 2812 | onwind | PV | 370 | 38.8 | -25.1 | 554 | 1.5 |
| 1005 storage-pv-placeholder | 1005 | hydro | PV | 785 | nan | -1.03e+03 | -515 | 1.31 |
| 459 storage-pv-placeholder | 459 | hydro | PV | 1.25e+03 | nan | -1.49e+03 | -958 | 1.2 |
| 339 slack-placeholder | 339 | AC | Slack | 7.72e+03 | 1.94e+03 | -8.03e+03 | -4.34e+03 | 1.04 |

*(9 lignes de plus dans le CSV)*

**Bilan réactif moyen (heures convergées)**, + = fourni au réseau. Les lignes produisent (capacité, b) et consomment (I²x) ; leur solde net est la somme q0 + q1, changée de signe.

| index | MVAr_moy |
|---|---|
| Q_generateurs | -55170 |
| Q_charges | -9439 |
| Q_shunts | 3493 |
| Q_net_lignes | 68437 |
| Q_net_transfos | -7460 |

**Réactif produit par la capacité des lignes** (b·V²). Une ligne chargée sous sa puissance naturelle (SIL) produit plus de réactif qu'elle n'en consomme. Le modèle n'a pas de réactances shunt de ligne sur le 735 kV ; ce surplus doit donc être absorbé par les PV et le slack.

| v_nom | lignes | Q_capacite_1pu_MVAr | Q_capacite_V_reel_MVAr | P_sur_SIL_median | part_lignes_sous_SIL |
|---|---|---|---|---|---|
| 315 | 167 | 5.69e+03 | 6.44e+03 | 0.478 | 0.754 |
| 345 | 1 | 9.93 | 12.4 | 0.00963 | 1 |
| 735 | 107 | 6.45e+04 | 8.9e+04 | 0.501 | 0.953 |
| 765 | 1 | 226 | 251 | 0.295 | 1 |

## 5. Lignes : charge et angle de transport (heures convergées, AC)

- Charge max (S/s_nom) : 75.9 % ; s_nom = St Clair × 3, donc gonflé (voir CLAUDE.md) : la charge thermique n'est pas le bon indicateur ici.
- Lignes avec un angle > 30° au moins une heure : 2 / 276.

| name | bus0 | bus1 | v_nom | km | circ | delta_max_deg | charge_max_pct | P_max_MW |
|---|---|---|---|---|---|---|---|---|
| 3360 | 710 | 84 | 315 | 410 | 1 | 78.1 | 75.9 | 696 |
| 3352 | 3836 | 1281 | 315 | 182 | 1 | 51.5 | 55.2 | 1.17e+03 |
| 3375 | 2858 | 305 | 315 | 186 | 1 | 22.3 | 28.1 | 617 |
| 3808 | 85 | 114 | 735 | 259 | 1 | 21.2 | 30.1 | 2.67e+03 |
| 1469 | 85 | 114 | 735 | 260 | 1 | 21.2 | 30.1 | 2.66e+03 |
| 170 | 297 | 453 | 315 | 141 | 2 | 20.8 | 23.7 | 1.36e+03 |
| 2406 | 306 | 85 | 735 | 349 | 1 | 20.4 | 26 | 1.86e+03 |
| 340 | 306 | 85 | 735 | 349 | 1 | 20.4 | 26 | 1.86e+03 |
| 2405 | 306 | 85 | 735 | 349 | 1 | 20.4 | 26 | 1.86e+03 |
| 2382 | 2202 | 84 | 315 | 111 | 1 | 20.3 | 24.3 | 872 |
| 1664 | 1806 | 1804 | 315 | 125 | 1 | 16.3 | 22.2 | 732 |
| 3809 | 85 | 112 | 735 | 278 | 1 | 15.8 | 24.5 | 1.86e+03 |
| 3756 | 651 | 25 | 735 | 623 | 3 | 15.2 | 32 | 2.78e+03 |
| 159 | 454 | 345 | 735 | 358 | 1 | 15.1 | 25.8 | 1.58e+03 |
| 2398 | 454 | 345 | 735 | 358 | 1 | 15.1 | 25.8 | 1.58e+03 |

*(261 lignes de plus dans le CSV)*

## 6. Les heures divergées, vues en DC power flow

Pour chaque heure : demande, rang (1 = pointe), plus grand angle DC et rapport P / Pmax le plus élevé (Pmax = V²/x). En DC, sin δ ≈ δ : un angle DC de plus de ~60-70° indique qu'une solution AC n'existe sans doute pas sans soutien de tension local.

**Heures divergées**

| index | charge_MW | rang_demande | diverge | delta_DC_max_deg | ligne_delta_max | P_sur_Pmax_max | delta_3360_deg | P_3360_MW | delta_3352_deg | P_3352_MW |
|---|---|---|---|---|---|---|---|---|---|---|
| 2022-01-03 09:00:00 | 3.26e+04 | 21 | True | 57.7 | 3360 | 1.01 | 57.7 | 677 | 39.3 | -1.04e+03 |
| 2022-01-03 16:00:00 | 3.37e+04 | 10 | True | 55.6 | 3360 | 0.97 | 55.6 | 653 | 42.2 | -1.11e+03 |
| 2022-01-03 19:00:00 | 3.41e+04 | 7 | True | 54.7 | 3360 | 0.954 | 54.7 | 642 | 43.1 | -1.14e+03 |
| 2022-01-03 22:00:00 | 3.23e+04 | 25 | True | 55.7 | 3360 | 0.972 | 55.7 | 654 | 40.8 | -1.08e+03 |
| 2022-01-03 23:00:00 | 3.19e+04 | 29 | True | 54.5 | 3360 | 0.951 | 54.5 | 640 | 40.2 | -1.06e+03 |
| 2022-01-04 16:00:00 | 3.35e+04 | 14 | True | 50.2 | 3360 | 0.876 | 50.2 | 590 | 39.2 | -1.04e+03 |

**Comparaison convergées / divergées (moyennes)**

| diverge | charge_MW | rang_demande | delta_DC_max_deg | P_sur_Pmax_max |
|---|---|---|---|---|
| False | 2.87e+04 | 87 | 46.5 | 0.811 |
| True | 3.3e+04 | 17.7 | 54.7 | 0.955 |

**Lignes les plus tendues aux heures divergées** (angle DC max), comparées aux heures convergées

| name | bus0 | bus1 | v_nom | km | circ | x_ohm | delta_DC_max_div | delta_DC_max_conv | P_sur_Pmax_div | P_sur_Pmax_conv |
|---|---|---|---|---|---|---|---|---|---|---|
| 3360 | 710 | 84 | 315 | 410 | 1 | 147 | 57.7 | 59.2 | 1.01 | 1.03 |
| 3352 | 3836 | 1281 | 315 | 182 | 1 | 65.6 | 43.1 | 44.2 | 0.753 | 0.771 |
| 3808 | 85 | 114 | 735 | 259 | 1 | 85.6 | 24.6 | 25.3 | 0.429 | 0.442 |
| 1469 | 85 | 114 | 735 | 260 | 1 | 85.9 | 24.6 | 25.3 | 0.429 | 0.442 |
| 3375 | 2858 | 305 | 315 | 186 | 1 | 67.1 | 23.9 | 24.4 | 0.418 | 0.426 |
| 340 | 306 | 85 | 735 | 349 | 1 | 115 | 22.5 | 23.1 | 0.393 | 0.403 |
| 2406 | 306 | 85 | 735 | 349 | 1 | 115 | 22.5 | 23.1 | 0.393 | 0.403 |
| 2405 | 306 | 85 | 735 | 349 | 1 | 115 | 22.5 | 23.1 | 0.393 | 0.403 |
| 3756 | 651 | 25 | 735 | 623 | 3 | 68.5 | 20.4 | 20.5 | 0.355 | 0.358 |
| 1664 | 1806 | 1804 | 315 | 125 | 1 | 45.1 | 20.1 | 18.9 | 0.351 | 0.329 |
| 2382 | 2202 | 84 | 315 | 111 | 1 | 39.9 | 19.7 | 20.1 | 0.344 | 0.35 |
| 170 | 297 | 453 | 315 | 141 | 2 | 25.4 | 18.8 | 20.7 | 0.329 | 0.361 |
| 2398 | 454 | 345 | 735 | 358 | 1 | 118 | 18.4 | 18.7 | 0.321 | 0.327 |
| 159 | 454 | 345 | 735 | 358 | 1 | 118 | 18.4 | 18.7 | 0.321 | 0.327 |
| 3809 | 85 | 112 | 735 | 278 | 1 | 91.9 | 18.2 | 18.4 | 0.317 | 0.321 |

*(261 lignes de plus dans le CSV)*

## 6b. Ponts (lignes uniques dont la perte isole une poche d'au moins 2 bus)

64 ponts à circuit unique au total (lignes + transfos), dont 30 isolent une poche d'au moins 2 bus. Pmax = V²/x à 1 pu ; s_nom est la limite que le LOPF a vue (St Clair × 3), le LOPF ne voit **pas** Pmax.

| ligne | bus0 | bus1 | v_nom | km | bus_poche | charge_moy_MW | prod_installee_MW | prod_moy_MW | Pmax_MW | P_DC_max_MW | P_sur_Pmax_max | s_nom_MVA | bus |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3360 | 84 | 710 | 315 | 410 | 2 | 29.1 | 1.66e+03 | 578 | 673 | 696 | 1.03 | 1.05e+03 | 653, 710 |
| 3352 | 1281 | 3836 | 315 | 182 | 9 | 963 | 780 | 49.7 | 1.51e+03 | 1.17e+03 | 0.771 | 2.28e+03 | 1548, 1693, 1762, 1772, 2207, 240, 2812, 3719, 3836 |
| 3375 | 305 | 2858 | 315 | 186 | 4 | 3.72 | 1.28e+03 | 541 | 1.48e+03 | 629 | 0.426 | 2.23e+03 | 2858, 3227, 582, 684 |
| 170 | 453 | 297 | 315 | 141 | 4 | 22.1 | 3.01e+03 | 1.14e+03 | 3.9e+03 | 1.41e+03 | 0.361 | 5.87e+03 | 452, 453, 686, 801 |
| 1664 | 1804 | 1806 | 315 | 125 | 2 | 658 | 411 | 213 | 2.2e+03 | 773 | 0.351 | 3.3e+03 | 1806, 1810 |
| 168 | 452 | 453 | 315 | 183 | 2 | 16 | 1.26e+03 | 398 | 3.01e+03 | 661 | 0.219 | 4.55e+03 | 452, 801 |
| 3350 | 240 | 1693 | 315 | 196 | 6 | 97.2 | 385 | 9.63 | 1.41e+03 | 248 | 0.176 | 2.13e+03 | 1548, 1762, 1772, 2207, 240, 3719 |
| 4107 | 684 | 3227 | 315 | 69.4 | 2 | 1.84 | 1.28e+03 | 541 | 3.97e+03 | 632 | 0.159 | 3.59e+03 | 582, 684 |
| 2397 | 459 | 456 | 315 | 95.1 | 2 | 1.61 | 2.65e+03 | 652 | 5.8e+03 | 811 | 0.14 | 7.19e+03 | 459, 460 |
| 1662 | 1800 | 1803 | 315 | 40.2 | 4 | 659 | 411 | 213 | 6.86e+03 | 774 | 0.113 | 3.59e+03 | 1803, 1804, 1806, 1810 |
| 3351 | 3836 | 1693 | 315 | 114 | 7 | 107 | 410 | 11.4 | 2.41e+03 | 237 | 0.0984 | 3.59e+03 | 1548, 1693, 1762, 1772, 2207, 240, 3719 |
| 1663 | 1803 | 1804 | 315 | 30.2 | 3 | 659 | 411 | 213 | 9.13e+03 | 773 | 0.0847 | 3.59e+03 | 1804, 1806, 1810 |
| 600 | 156 | 2973 | 315 | 96.1 | 16 | 444 | 1.48e+03 | 214 | 5.73e+03 | 476 | 0.083 | 7.19e+03 | 1092, 141, 1519, 1520, 1523, 1538, 1539, 1542, 1543, 233, 2813, 2973 … |
| 1656 | 1800 | 311 | 315 | 29.2 | 5 | 659 | 411 | 213 | 9.42e+03 | 774 | 0.0821 | 3.59e+03 | 1800, 1803, 1804, 1806, 1810 |
| 603 | 3067 | 141 | 315 | 96.6 | 13 | 336 | 1.14e+03 | 151 | 5.7e+03 | 373 | 0.0653 | 7.19e+03 | 1092, 141, 1519, 1520, 1523, 1538, 1539, 1542, 1543, 233, 2813, 3078 … |
| 4106 | 2858 | 3227 | 315 | 27.9 | 3 | 2.37 | 1.28e+03 | 541 | 9.9e+03 | 631 | 0.0638 | 3.59e+03 | 3227, 582, 684 |
| 601 | 2973 | 3066 | 315 | 48.2 | 15 | 408 | 1.37e+03 | 198 | 1.14e+04 | 444 | 0.0388 | 7.19e+03 | 1092, 141, 1519, 1520, 1523, 1538, 1539, 1542, 1543, 233, 2813, 3066 … |
| 602 | 3066 | 3067 | 315 | 46.3 | 14 | 352 | 1.14e+03 | 151 | 1.19e+04 | 391 | 0.0329 | 7.19e+03 | 1092, 141, 1519, 1520, 1523, 1538, 1539, 1542, 1543, 233, 2813, 3067 … |
| 609 | 141 | 1092 | 315 | 50.2 | 12 | 320 | 1.14e+03 | 151 | 1.1e+04 | 355 | 0.0323 | 7.19e+03 | 1092, 1519, 1520, 1523, 1538, 1539, 1542, 1543, 233, 2813, 3078, 3837 |
| 1359 | 2913 | 305 | 315 | 248 | 2 | 24.8 | 0 | 0 | 1.11e+03 | 30.1 | 0.0271 | 1.69e+03 | 1638, 2913 |

*(10 lignes de plus dans le CSV)*

## 7. Hypothèses génériques actives dans ce résultat (rappel)

- Q des charges : fp 0,95 partout. Q des générateurs PQ : 0,484 × p_nom fixe (pas lié au dispatch).
- Aucune limite Q sur PV/Slack (PyPSA PF). Consigne de tension des PV : [1.0] pu.
- Classement PV : `--pv-min-capacity 0 --pv-max-load-ratio 1` (réglage 735 kV, pas 315 kV).
- Compensation série 70 % (générique) ; 4 lignes 735 kV ajoutées par Etienne (3535, 161, 1653, 2387), choisies empiriquement sur la convergence.
- Shunts fixes dimensionnés à l'heure la pire, réacteurs ×1,0 / condensateurs ×0,7.
- s_nom = St Clair × 3 ; transformateurs plats à 100 000 MVA ; dispatch = LOPF du réseau **non corrigé** (elec_solved), jamais ré-optimisé après les corrections.
