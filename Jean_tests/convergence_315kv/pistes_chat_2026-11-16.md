> **Mise à jour du 17 nov.** Ces pistes ont été menées jusqu'au bout : voir `reseau_315kv/` (README + notebook `Reseau_315kV_converge.ipynb`). Deux corrections par rapport à ce qui suit : (1) la règle des longueurs a été affinée (longueur géodésique du tracé OSM quand il est fiable, sinon vol d'oiseau × 1,08 : la méthode « stocké × vrai/Mercator » rendait certaines lignes plus courtes que le vol d'oiseau) ; (2) les 19 dépassements de limites Q sont réglés avec du réactif réel (SVC et compensateurs synchrones d'HQ, inductances et condensateurs commutables). Résultat : 168/168 h convergées avec limites Q, 165/168 dans 0,95-1,05 pu, V 0,922-1,035 pu.

# Diagnostic du 315 kV — tests faits depuis la conversation claude.ai (16 nov. 2026, soir)

Ces tests ont été faits en dehors de Claude Code, sur une copie du dépôt (PyPSA 1.3.0 + HiGHS),
à partir de `network/networks_current/elec_solved.nc` et des scripts de la branche `etienne`.
**Rien n'est encore intégré au pipeline** : le but est de refaire chaque étape avec Claude Code,
une hypothèse à la fois, en notant chacune dans le journal du README.
Seul fichier fourni : `fix_line_lengths_mercator.py` (testé, lecture seule sur `networks_current/`).

## 1. La cause principale : des longueurs de lignes gonflées de ~50 % (Web Mercator)

PyPSA-Earth mesure la longueur des lignes OSM dans `distance_crs: EPSG:3857` (config.default.yaml,
ligne 77), dans `scripts/build_osm_network.py` :

    df_l["length"] = df_l.to_crs(distance_crs).geometry.length

Web Mercator (la projection de Google Maps) n'est pas équidistante : elle étire les distances
d'un facteur 1/cos(latitude), soit ×1,41 à Montréal (45,5° N) et ×1,66 à la Baie-James (53° N).
Comme r, x et b sont proportionnels à la longueur (`apply_hq_line_characteristics.py`), **toutes
les lignes ont environ 50 % de trop en réactance série ET en capacité shunt.**

Preuves (sur `elec_full.nc`, 273 lignes 315/735 kV avec géométrie OSM) :
- longueur stockée / longueur de la même géométrie en EPSG:3857 = 1,00 à 1 % près pour la plupart
  des lignes ; longueur stockée / longueur géodésique = 1/cos(latitude) (ex. 3360 : 1,563 contre
  1/cos(50,3°) = 1,567 ; ligne 161 : 1,647 contre 1,650) ;
- les géométries sont complètes : extrémités à < 1 km des bus, sinuosité médiane 1,10-1,14 ;
- 735 kV : **19 486 km stockés, contre 12 895 km mesurés sur la géométrie**. Hydro-Québec
  déclare **12 581 km** de lignes à 765/735 kV (rapport annuel 2023), Labrador non compris.
  Après correction, le réseau 315 kV a 13 337 circuit-km de 735 kV.

Ça avait été pressenti le 3 nov. (Obsidian, « Hypothèses du modèle », A1.3 : « probablement la
plus grosse erreur du modèle »), sans en trouver la cause. Le test « longueurs × 0,7 » de la
nuit du 3 nov. (B8) allait dans le bon sens ; la correction ci-dessous est ligne par ligne
(facteur 0,58 à 0,70 selon la latitude), et non un facteur unique.

Conséquences :
- **Le pont 3360 n'est plus un problème** : 262 km au lieu de 410 km, donc V²/x passe de 673 MW
  à ~1 050 MW, au-dessus des 696 MW exportés par la Romaine.
- La capacité des lignes 735 kV passe de 64,5 à ~43 GVAr : beaucoup moins de surtension.
- La compensation série à 70 % et le passage à 3 circuits de `correct_sending_end_circuits.py`
  ont été calibrés sur des x trop grands : à réexaminer (HQ réel : 17-40 % sur 37 lignes).
- **Tous les réseaux du dépôt sont touchés** (735 kV de Nathan compris), y compris `s_nom` (St Clair
  dépend de la longueur) : il faudra relancer le LOPF après correction.

## 2. Résultats des tests (AC PF, 168 h, réglages PV d'Etienne : `--pv-min-capacity 0 --pv-max-load-ratio 1`)

| Variante | Convergés | 0,90-1,10 pu | V (h convergées) |
|---|---|---|---|
| Séquence d'Etienne (reproduction de Claude Code) | 162 | **0** | 0,91-1,62 |
| Longueurs corrigées seules | 56 | 0 | 0,78-1,19 |
| Longueurs + circuits | 100 | 0 | 0,75-1,23 |
| Longueurs + circuits + compensation 70 % (19 lignes) | 161 | 0 | 0,79-1,20 |
| **+ shunts** (`add_shunt_impedance.py`, inchangé) = séquence d'Etienne précédée des longueurs | 162 | **56** | 0,85-1,22 |
| + bus 1810 en PV (seule règle changée, voir § 3) | **168** | 92 | 0,94-1,12 |
| + réactances shunt fixes = 30 % de la capacité des lignes 735 kV | **168** | **168** | 0,90-1,07 |
| idem avec 50 % | 168 | 157 (93 dans 0,95-1,05) | 0,88-1,05 |

Pertes dans les lignes de la dernière variante : 3,6 % de la charge.

Ce qui n'a PAS marché (à ne pas refaire) : mettre toutes les centrales synchrones en PV et l'éolien
en PQ à q = 0 (112/168) ; Q = 0 sur les PQ (114/168) ; retirer les shunts d'Etienne (Montréal
descend à 0,73-0,76 pu) ; compensation série 40 % sans autre changement (moins bon que 70 %).

## 3. Les 6 heures qui divergeaient encore (01-03 9h, 16h, 19h, 22h, 23h, 01-04 16h)

Avec les longueurs corrigées, ce n'est plus la 3360. En réduisant la charge de cette heure jusqu'à
la convergence (λ = 0,94 à 19h), le bus le plus bas est **1810** (0,86 pu), au bout de la chaîne
radiale 311-1800-1803-1804-1806-1810 (Abitibi-Témiscamingue). Ce bus a 579 MW de charge et trois
centrales au fil de l'eau (Première-Chute, Rapides-des-Îles, Rapides-des-Quinze, 411 MW) laissées
en PQ par la règle de `run_pf.py` (charge locale / production = 1,4 > 1).
Correction testée : tout bus avec ≥ 100 MW d'hydro synchrone devient PV, quelle que soit la charge
locale (une centrale hydro tient la tension de son poste). Seul le bus 1810 change, et les 6 heures
convergent.

## 4. Ce qui reste faux (le critère « réactif plausible » n'est PAS encore tenu)

19 des 47 générateurs PV/slack dépassent |Q| ≤ 0,484·p_nom dans la meilleure variante. Les pires :
solaire Roussillon (6 MW) de -1,4 à +1,9 GVAr, Trenche (302 MW) jusqu'à -4,8 GVAr, éolien Gros-Morne
(224 MW) jusqu'à -2,7 GVAr. Ce sont de petites centrales placées en PV sur des bus 735 kV : elles
jouent le rôle des SVC et compensateurs synchrones qu'HQ a vraiment dans ces postes.
→ Prochaines étapes : limites Q avec bascule PV → PQ (déjà codée pour le 735 kV la nuit du 3 nov. :
`AC_convergence_nuit/ac_tools.py`, B1), puis remplacer ces « faux régulateurs » par des SVC
explicites (P = 0, Q limité) aux postes où HQ en a.

## 5. Plan proposé pour Claude Code (une ligne du journal par étape)

1. Relire `fix_line_lengths_mercator.py` avec Claude Code (comprendre EPSG:3857 et le facteur
   1/cos φ), puis l'exécuter :
   `python Jean_tests/convergence_315kv/fix_line_lengths_mercator.py --network network/networks_current/elec_solved.nc --output Jean_tests/convergence_315kv/results/elec_solved_len.nc`
2. Relancer `reproduce_etienne.py` avec `SOURCE` = ce fichier (attendu : 56 / 100 / 161 / 162
   convergés, 56 dans la bande à la dernière étape).
3. Ajouter la règle « hydro ≥ 100 MW → PV » (attendu : 168/168, 92 dans la bande).
4. Ajouter les réactances shunt de ligne sur le 735 kV **avant** le dimensionnement des shunts
   (dans les tests ci-dessus, elles ont été ajoutées après : à refaire dans le bon ordre).
   Chercher le niveau réel chez HQ pour justifier la valeur.
5. Relancer le LOPF sur le réseau corrigé (s_nom et x ont changé), puis l'AC PF.
6. Limites Q (étape 4 ci-dessus), puis revoir la compensation série vers les niveaux HQ.
7. Prévenir Etienne, Jana (et Nathan / la Prof. Wang) : la correction des longueurs touche tous les réseaux.
