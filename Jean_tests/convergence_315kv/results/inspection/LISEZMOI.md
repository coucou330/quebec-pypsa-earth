# Inspection du réseau 315 kV d'Etienne — synthèse (16 nov. 2026)

Fichiers de ce dossier : `etat_des_lieux_315kv.md` (rapport chiffré complet, généré par
`Jean_tests/convergence_315kv/inspect_315kv.py`) et les CSV derrière chaque section.
Réseau inspecté : `results/elec_315kv_shunt_pf.nc`, la reproduction de la séquence d'Etienne
(commit 2ef0c21d) faite avec `reproduce_etienne.py` en PyPSA 1.3.0.

## Faits mesurés

| Étape | Newton converge | Converge et 0,90 ≤ V ≤ 1,10 pu |
|---|---|---|
| `elec_solved.nc` brut | 3/168 | 0 |
| + correction des circuits (5 corridors 735 kV) | 10/168 | 0 |
| + compensation série 70 % (19 lignes 735 kV) | 157/168 | 0 |
| + 10 shunts fixes | 162/168 (Etienne : 161) | **0** |

1. **Aucune heure n'est acceptable en tension.** Sur les 162 heures convergées, V va de 0,91 à
   1,62 pu ; 96 bus AC sur 200 sortent de 0,90-1,10 pu, et 47 des 56 bus 735 kV sont au-dessus
   de 1,10 pu toute la semaine (jusqu'à 1,34 pu au bus 312). Plus la demande baisse, plus les
   tensions montent : la corrélation entre la demande et le nombre de bus au-dessus de 1,10 pu est de −0,97.
2. **Le réseau a un énorme surplus de réactif.** La capacité des lignes 735 kV produit
   64,5 GVAr à 1 pu (89 GVAr aux tensions calculées) ; 95 % de ces lignes transitent moins que
   leur puissance naturelle (médiane P/SIL = 0,50). Bilan moyen : les lignes fournissent
   +68 GVAr net, et les générateurs absorbent −55 GVAr. Seuls 0,26 GVAr de réacteurs ont été posés.
3. **Le réactif demandé aux générateurs n'est pas physique.** 29 des 46 générateurs PV/slack
   dépassent |Q| ≤ 0,484·p_nom. Exemples : la centrale au fil de l'eau Trenche (302 MW) absorbe
   jusqu'à 13,5 GVAr, l'éolien Gros-Morne (224 MW) 8 GVAr, le solaire Roussillon (6 MW) 4 GVAr,
   et le slack 339 entre 4,3 et 8 GVAr. PyPSA n'impose aucune limite de Q en PF : c'est ce qui
   rend la « convergence » possible.
4. **Les 6 heures divergées** (le 3 janvier à 9h, 16h, 19h, 22h et 23h, le 4 janvier à 16h) sont
   des heures de forte demande (rang moyen 18/168, 33,0 GW contre 28,7 GW). Dans chacune, la ligne
   la plus tendue en DC PF est la **3360** (710 → 84, 315 kV, 410 km, un seul circuit, pont).
   Elle exporte la production de la poche 710/653 (Romaine, ~580 MW en moyenne, jusqu'à 696 MW)
   alors que sa limite statique vaut V²/x = 673 MW à 1 pu : P/Pmax atteint 0,88 à 1,01. Même
   sur les heures convergées, son angle AC monte à 78°.
   La **3352** (1281 → 3836, 182 km, pont) alimente seule une poche de 9 bus : 963 MW de charge
   en moyenne, 780 MW installés localement mais seulement 50 MW dispatchés par le LOPF ;
   P/Pmax ≈ 0,77.
5. **Le LOPF ne voit pas ces limites.** Pour la 3360, s_nom = 1052 MVA (St Clair × 3), au-dessus
   de Pmax = 673 MW. De plus, le dispatch vient du LOPF du réseau non corrigé et n'est jamais
   ré-optimisé après les corrections.
6. **Chaîne radiale 311-1800-1803-1804-1806-1810** (V jusqu'à 1,62 pu au bus 1810) : elle
   importe ~440 MW vers le bus 1810 (579 MW de charge, 3 centrales au fil de l'eau en PQ avec un
   Q fixe de 199 MVAr) et part déjà d'un bus 311 à 1,35 pu, collé au 735 kV surélevé (312 à 1,34 pu).

## Interprétation (à vérifier, pas encore testée)

- La surtension généralisée vient d'abord de la **capacité des lignes 735 kV non compensée en
  shunt**. Sur un vrai réseau 735 kV, des réactances shunt sont raccordées aux lignes en
  permanence (à confirmer et chiffrer avec des données HQ). Le modèle n'en a pas ; les PV et le
  slack absorbent donc un réactif impossible, et les bus PQ éloignés montent.
- Le Q fixe de 0,484·p_nom imposé à chaque générateur PQ (même à P = 0, hypothèse de `run_pf.py`)
  ajoute 1,9 GVAr dans le même sens.
- La divergence des 6 heures est un **problème de transfert** (3360 à sa limite V²/x), pas de
  réglage de Newton. Pistes : vérifier dans les données réelles le nombre de circuits ou le
  niveau de tension du corridor Romaine–Arnaud, ou contraindre le LOPF par la limite de stabilité.
- Le « 162/168 » repose sur trois hypothèses sans fondement physique : Q illimité aux PV, 4 lignes
  compensées choisies empiriquement, et des PV sans seuil de taille (`--pv-min-capacity 0`).

## Pistes, une hypothèse à la fois

1. Réactances shunt de ligne sur le 735 kV (fraction de b·V² par ligne, à documenter) : on
   attend une baisse générale des tensions et du réactif absorbé par les PV.
2. Limites Q sur les PV (contrôle a posteriori, ou PV → PQ quand la limite est atteinte).
3. Corridor 3360 : vérification des données, puis limite de stabilité dans le LOPF.
4. Revoir le Q fixe des générateurs PQ (le lier au dispatch, ou à 0 quand P = 0).
