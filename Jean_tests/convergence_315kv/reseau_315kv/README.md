# Réseau 315 kV convergé (16-17 nov. 2026)

**À ouvrir en premier : `Reseau_315kV_converge.ipynb`** (déjà exécuté, sorties et figures incluses).
Travail fait par Claude dans la conversation claude.ai, en autonomie, à partir du modèle de Nathan et du
travail d'Etienne. **À relire et refaire étape par étape avec Claude Code** ; à déclarer comme usage d'IA.

## Résultat
Semaine test 1-7 janv. 2022 (168 h) : **168/168 h convergées en AC avec les limites Q** de tous les
générateurs, SVC et compensateurs synchrones ; tensions **0,922-1,035 pu** (168/168 h dans 0,90-1,10,
165/168 dans 0,95-1,05 : les 3 heures les plus chargées, à Montréal-Ouest) ; LOPF sans délestage ; même
solution avec trois solveurs (moteur maison, PyPSA, pandapower). Etienne : 162/168 h, 0 h dans la bande,
réactif illimité.

## Les étapes (détail et physique dans le notebook)
| | Correction | Ce que ça règle |
|---|---|---|
| E1 | Longueurs recalculées sur le tracé OSM (ou vol d'oiseau × 1,08 si le tracé n'est pas fiable) : PyPSA-Earth les mesure en Web Mercator (×1,41 à ×1,66) | x et b trop grands de ~50 % partout |
| E2 | Churchill Falls : centrale rattachée au vrai poste, 3 lignes en double supprimées | 6 circuits au lieu de 3 |
| E2b | Équivalents de Kron du réseau < 315 kV (Saguenay, Montréal) | la réduction a supprimé les chemins < 315 kV (20 postes 735 kV sans transformation) |
| E3 | Circuits de Nathan (inchangé) | sans elle, 2 h hors limites Q (§12 du notebook) |
| E4 | Compensation série 40 % sur 8 corridors réels d'HQ | au lieu de 70 % générique |
| E5 | LOPF relancé, plafond de stabilité (V²/x) sin 40° | Romaine (ligne 3360) |
| E6 | Régulateurs à limites Q, SVC/CS réels d'HQ, inductances et condensateurs commutables, consigne 1,02 pu | réactif réaliste |
| E7 | Newton-Raphson heure par heure + bascule PV→PQ + manœuvre des shunts (réserve des SVC) | — |

## Fichiers
| Fichier | Rôle |
|---|---|
| `reseau315.py` | toutes les étapes (une fonction chacune), moteur Newton-Raphson, export MATPOWER, validations |
| `construire.py` | reconstruit tout et écrit `resultats/` (`--valider` : + PyPSA et pandapower, ≈ 10 min) |
| `sensibilites.py` | alternatives de conception : une hypothèse changée à la fois, + marge de charge à la pointe (16 variantes, ≈ 50 min avec `python sensibilites.py 0 2` et `1 2` en parallèle, puis `--fusion`) |
| `Reseau_315kV_converge.ipynb` | présentation complète (exécuté) |
| `resultats/reseau_315kv_final.nc` | réseau prêt pour l'AC (E1-E6) |
| `resultats/*.csv` | résultats heure par heure (tensions, Q, états des régulateurs et des shunts), validations, sensibilités |
| `resultats/quebec_315kv_pointe.m` | cas MATPOWER de la pointe (3 janv. 17 h), état final, pour AC-CFM |
| `resultats/figures/` | figures du notebook |
| `_faire_notebook.py`, `_reseau315_exploration.py` | génération du notebook ; version de travail du module (avec les essais abandonnés) |

## Relancer (venv `pypsa-env`, PyPSA 1.3.0)
```
cd quebec-pypsa-earth
pypsa-env\Scripts\python.exe -m pip install matplotlib pandapower
pypsa-env\Scripts\python.exe Jean_tests\convergence_315kv\reseau_315kv\construire.py --valider
```
Rien n'est écrit dans `network/networks_current/`. Testé avec PyPSA 1.3.0 et pandas 2.3.3 : installer
pandapower a fait passer pandas de 3.x à 2.3.3 dans l'environnement de test. Si on veut garder pandas 3 dans
`pypsa-env`, lancer `construire.py` sans `--valider` (pandapower n'est utilisé que pour la validation).

## Hypothèses les plus fragiles (détail : §13 du notebook)
1. Bancs de condensateurs aux postes de charge = Q de la charge à la pointe (10,8 GVAr, aucune donnée HQ).
2. Capacité réactive des centrales sur la puissance installée (optimiste).
3. Puissance des compensateurs synchrones de Duvernay, Abitibi, Manicouagan supposée (= Lévis).
4. **Pas de marge de charge à la pointe** (moins de +0,5 % ; le vrai réseau a tenu 43,1 GW, ×1,23).
   Toutes les sources de réactif du sud sont en butée. La marge monte à ×1,115 avec 50 % de condensateurs en
   plus et à ×1,06 avec 70 % de compensation série : elle dépend des données de réactif les moins sourcées.
   Avant AC-CFM à la pointe, il faut ces données (ou étudier une heure moins chargée).
5. Le LOPF a plusieurs optimums de même coût : le code fixe l'ordre des composants (résultats reproductibles),
   mais un autre optimum donnerait des tensions un peu différentes (~0,01-0,04 pu à la pointe).
