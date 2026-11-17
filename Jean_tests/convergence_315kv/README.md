# Convergence AC du réseau 315 kV — dossier de travail de Jean

Créé le 16 nov. 2026. Ici vit tout ce qu'on code avec Claude Code pour faire converger
l'AC power flow du réseau 315 kV de Nathan (`network/networks_current/elec_solved.nc`),
sur la branche `jean`. On ne modifie jamais `network/networks_current/` : les fichiers
intermédiaires vont dans `results/`.

## Où en est Etienne (branche `origin/etienne`, commits du 11 au 14 nov.)

Notebook : `etienne_tests/315_convergence.ipynb`. Il applique au 315 kV les corrections
que Nathan avait faites pour le 735 kV, puis ajoute de la compensation série :

| Étape | Snapshots convergés |
|---|---|
| `elec_solved.nc` brut | 3/168 |
| + `correct_sending_end_circuits.py` | 10/168 |
| + `apply_series_compensation.py` (70 %) | 106/168 |
| + compensation des lignes 735 kV 3535, 161, 1653, 2387 | 157/168 |
| + `add_shunt_impedance.py` (SCALE_GRID réduit) | **161/168** |

Il a aussi modifié deux scripts de Nathan (`network/add_shunt_impedance.py` : accepte le
meilleur couple réacteurs/condensateurs même sans 168/168 ; `network/apply_series_compensation.py` :
4 lignes ajoutées) et écrit `etienne_tests/diagnose_315kv.py`.
Jana a ajouté `jana_tests/find_bad_hours.py` (liste les heures divergées d'un `*_pf.nc`).

### Points à vérifier avant de dire « ça converge »
1. **Convergé ≠ physiquement acceptable.** Sur les heures « convergées », Etienne voit des
   tensions jusqu'à ~1,5 pu (bus 1806, 1804, 1803, 1800, 311, une chaîne radiale 315 kV)
   et des minima ~0,75 pu. Un critère honnête : Newton converge **et** 0,90 ≤ V ≤ 1,10 pu
   (idéalement 0,95-1,05) **et** réactif des générateurs dans des limites plausibles.
2. **Les 7 heures restantes** : probablement les heures de pointe et les deux « ponts »
   identifiés dans le cours (ligne 3360 = 710-84, 410 km ; ligne 3352 = 3836-1281,
   182 km, qui alimente seule une antenne de 9 bus). Le flux demandé y dépasse la limite
   V²/x : aucun réglage de Newton ne les fera converger. À confirmer en premier.
3. **Les lignes compensées par Etienne sont des 735 kV** (161, 1653, 2387, 3535 : 370-470 km),
   ce qui est physiquement plausible (HQ compense en série son 735 kV). À documenter.
4. **Version de PyPSA** : Etienne tourne avec PyPSA 0.30.x (pixi, comme Nathan). Le venv
   `pypsa-env` de Jean est en PyPSA 1.3.0. Pour reproduire ses chiffres à l'identique,
   utiliser un venv en 0.30.3, ou vérifier que 1.3 donne les mêmes résultats.
   → **Vérifié le 16 nov.** : en 1.3.0, on retrouve 3 / 10 / 157 à l'identique et 162 au lieu de
   161 pour la dernière étape. On reste donc en 1.3.0.

## Plan de la séance (objectif : un 315 kV convergé sur `jean` pour le 17 nov.)
1. Synchroniser : `git fetch origin`, mettre de côté la modif en cours du notebook
   `Jean_tests/Comprendre_PyPSA_et_le_code_de_Nathan.ipynb`, récupérer le travail d'Etienne.
2. Environnement : venv PyPSA 0.30.3 (ou valider 1.3.0) + HiGHS.
3. Reproduire les 161/168 d'Etienne, avec un script à nous (pas de notebook à `!python`).
4. Diagnostiquer les 7 heures divergées (ponts 3360 / 3352 ?) et les surtensions.
5. Corriger, hypothèse par hypothèse (chaque correction = une ligne dans le journal ci-dessous).
6. Commit sur `jean` avec un message clair.

## Journal des hypothèses / corrections
| Date | Changement | Justification physique | Effet (convergés, V min-max) |
|---|---|---|---|
| 2026-11-16 | Reproduction de la séquence d'Etienne (2ef0c21d) avec `reproduce_etienne.py`, PyPSA 1.3.0, scripts de Nathan non modifiés (versions d'Etienne) | Aucune nouvelle hypothèse : on vérifie seulement que les chiffres d'Etienne se reproduisent | Brut 3/168 (0 dans la bande 0,90-1,10) ; + circuits 10 (0) ; + compensation 157 (0) ; + shunts **162** (0), V des heures convergées 0,910-1,615 pu (max : bus 1810, 1806, 1804). Etienne avait 161 : son propre dimensionnement des shunts annonçait déjà 162. Heures divergées : 01-03 09h, 16h, 19h, 22h, 23h et 01-04 16h |
| 2026-11-16 | Inspection du réseau final (`inspect_315kv.py`, lecture seule + DC PF) → `results/inspection/` | Aucun changement de modèle | 0/162 dans la bande ; lignes 735 kV : +65 GVAr de capacité non compensée, absorbés par des PV au-delà de toute capacité réelle ; les 6 heures divergées : ligne 3360 (pont, 410 km) à P/Pmax ≈ 1 |
