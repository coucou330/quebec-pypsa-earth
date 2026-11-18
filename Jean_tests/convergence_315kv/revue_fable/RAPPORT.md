# Revue critique du réseau 315 kV convergé — nuit du 16 au 17 nov. 2026

Revue faite par Claude (Fable) sur le travail de la nuit (`../reseau_315kv/`). Rien n'a été modifié dans
`reseau315.py`, `resultats/` ni `network/networks_current/`. Aucun commit git. Tous les scripts de cette
revue sont dans ce dossier ; ils importent `reseau315` depuis le dossier voisin (`_commun.py`).

« Convergé » = Newton converge **et** tensions dans 0,90-1,10 pu **et** régulateurs dans leurs limites Q.

## Résumé en 10 lignes

1. **Solide** : règle des longueurs (aucune ligne aberrante), algèbre de Kron (exacte à 10⁻¹⁴), cohérence
   état PV/PQ ↔ tensions sur 168 h, code compatible pandas 3, SVC conformes à la Régie (sauf Qmin de Chénier).
2. **Faux dans le cas final** : la solution de la pointe est sur la **branche instable** de la courbe P-V
   (dV/dQ < 0 à 59 bus). Le bon chiffre est 167/168 heures stables, et la marge de la pointe est négative.
3. **Cause** : le LOPF (hydraulique toute à coût nul, optimum dégénéré) arrête **3,4 GW de fil de l'eau du sud**
   (Beauharnois, Carillon…) à la pointe, et leur réactif est quand même compté (1,7 GVAr de machines arrêtées).
4. **Une ligne de code** (fil de l'eau prioritaire, coût −0,2 $/MWh) : **168/168 Newton, limites Q, 0,95-1,05
   pu et stabilité V-Q**, V 0,965-1,037, **marge stable +7 %** (+6 % avec le réactif des seules machines en marche).
5. Les pertes sont comptées deux fois si la demande d'HQ les inclut déjà, et toutes produites par LG-2 :
   corriger les deux, en plus du point 4, porte la marge à ~+11 % (stabilité non vérifiée ; à confirmer :
   définition du jeu de données de demande).
6. **Douteux** : E3 (circuits de Nathan) ajoute 2 700 circuit-km au 735 kV (+17 % par rapport à HQ) et
   ~7 GVAr de capacité ; il n'est plus nécessaire une fois le dispatch corrigé.
7. Kron : deux défauts mineurs (r des transformateurs équivalents écrasé par la préparation de Nathan ;
   longueurs par l'ancienne méthode). **Étendu à tout le sud (bus du Québec seulement), il ajoute ~3 % de
   marge** : avec le point 4, marge stable **+10 % (38,5 GW)**, 168/168 sur les cinq critères. Attention : une
   boîte rectangulaire au sud contient ~500 bus de l'Ontario et des États-Unis, à exclure.
   Duvernay : 900 MVAr supposés, au point le plus faible du réseau.
8. **Sous 315 kV** : converge avec limites Q, mais seulement en laissant 7,2 GW sur 35 au 315 kV, avec 8 à
   11 bus hors 0,90-1,10 et sans marge. Le blocage principal : des charges impossibles (851 MW au bout
   d'une ligne 69 kV) issues de la répartition régionale de la demande. Garder Kron pour l'instant.
9. **Pour la Prof. Wang** : présenter le bug Web Mercator comme découverte principale ; le 315 kV avec le fil
   de l'eau prioritaire, Kron étendu au sud et le critère de stabilité ; la marge (+8-10 % contre +23 % en
   réalité) comme limite connue, avec ses causes probables (pertes comptées deux fois, réactif du sud non
   sourcé, répartition des charges).
10. Ne pas présenter « 168/168 convergées » sans dire que la pointe du cas final était instable : c'est
    précisément ce que la revue a trouvé, et c'est réglé par le point 4.

---

## Partie A — revue critique

### A1. Règle des longueurs : appliquée correctement — confiance élevée
- **Constat.** Dans le réseau final, les 273 lignes réelles ont toutes une longueur entre 1,000 et 1,599 fois
  le vol d'oiseau entre leurs bus (médiane 1,084). Aucune ligne plus courte que le vol d'oiseau, aucune au-dessus
  de 1,7. Aucune ligne « sautée » : 271 corrigées, 2 déjà bonnes (les 3 « sans géométrie » sont les doublons de
  Churchill Falls, supprimés ensuite).
- **Preuve.** `a1_longueurs.py` → `a1_longueurs.csv`.

### A1b. Trop de circuits à 735 kV après la correction des circuits (E3) — confiance moyenne
- **Constat.** Les longueurs sont bonnes, mais le **nombre de circuits** ne l'est probablement pas. Le réseau
  final a **15 418 circuit-km** à 735 kV. Sans E3, il en a 12 696 (dont 689 km vers Churchill Falls, au
  Labrador), soit 12 007 km au Québec pour 12 581 km déclarés par HQ (−5 %). E3 (script de Nathan,
  `correct_sending_end_circuits`) ajoute **2 684 circuit-km** sur 7 lignes (ex. 651-25 Micoua-Laurentides,
  409 km, ×3 ; 114-1291, 334 km, ×3 ; 340-457, 277 km, ×3) : le total dépasse alors celui d'HQ de **17 %**.
- **Conséquence.** Ces circuits en trop ajoutent ~7 GVAr de capacité de ligne (2,65 MVAr/km) : c'est à peu
  près l'écart entre les inductances que le modèle doit enclencher la nuit (jusqu'à 32,4 GVAr) et tout ce
  qu'HQ a d'installé (25,5 GVAr en 1996). Vérifié : sans E3, le maximum enclenché tombe à 24,6 GVAr (variante
  « A1b » de la partie B). Sans E3 seul, 2 heures sortent des limites Q ; avec le dispatch corrigé (B4 + B5),
  E3 n'est plus nécessaire (168/168 partout).
- **Preuve.** `a1b_km735.py`.
- **Réserve.** OSM peut avoir oublié des lignes parallèles, et c'est ce qu'E3 voulait corriger. Mais le total
  d'HQ dit qu'il en manquait ~600 km, pas 2 700. À vérifier ligne par ligne avec la carte d'HQ.

### A2. Réduction de Kron : l'algèbre est exacte, deux défauts autour — confiance élevée
- **Algèbre exacte.** AC PF du sous-réseau complet (bus frontières + bus éliminés) contre AC PF du réseau
  réduit, mêmes injections aux frontières : écart de tension **10⁻¹⁴ pu** (Montréal : 57 bus → 13 ; Saguenay :
  37 → 3). Les 47 + 10 branches frontière-frontière sautées sont bien toutes dans le modèle.
- **Défaut 1 (bug, effet faible).** `prepare_ac` appelle la préparation de Nathan, qui remplace le r de **tous**
  les transformateurs par x / (rapport X/R), y compris les 14 transformateurs équivalents de Kron (qui
  représentent des lignes) : leur r est divisé par 4 à 25. x est intact. Effet sur un cas test : 0,005 pu et
  0,15° aux frontières ; les pertes du réseau régional sont sous-estimées.
- **Défaut 2 (incohérence, effet à mesurer : voir B2).** `kron_equivalent` corrige les longueurs des lignes
  éliminées avec `T.factor` (géodésique / Mercator), c'est-à-dire **l'ancienne méthode, abandonnée pour le
  réseau principal**. Pour Montréal, 40 des 80 lignes éliminées ont un tracé non fiable ; avec la règle hybride,
  leurs impédances seraient **35 % plus grandes** en moyenne pondérée (Saguenay : +6 %). L'équivalent de
  Montréal est donc trop « fort ».
- **À savoir.** Le sous-réseau régional n'est pas connexe sans les branches du modèle : à Montréal, seuls
  13 + 3 des 30 bus frontières sont reliés entre eux par le réseau < 315 kV (Saguenay : 3 sur 9).
- **Preuve.** `a2_kron.py`, `a2b_kron_diag.py`.
- **Correctifs proposés** (non appliqués) :
```diff
--- reseau315.py  (kron_equivalent)
-        f = T.at[name, "factor"] if name in T.index else 1.0
+        if name in T.index:   # même règle que fix_lengths
+            L1 = T.at[name, "true_km"] if T.at[name, "reliable"] else sin_med * T.at[name, "gc_km"]
+            f = L1 / T.at[name, "length_stored"]
+        else:
+            f = 1.0
   (avec, en tête de fonction :  sin_med = float(T.loc[T.reliable, "sinuosity"].median()) )
--- reseau315.py  (prepare_ac)
+    kron_t = n.transformers.index[n.transformers.index.str.startswith("KRON")]
+    r_kron = n.transformers.loc[kron_t, "r"].copy()
     with quiet():
         run_pf.use_lopf_dispatch_as_setpoints(n)
         n.determine_network_topology()
         run_pf.prepare_for_ac_pf(n, pv_min_capacity=0.0, pv_max_load_ratio=1.0)
+    n.transformers.loc[kron_t, "r"] = r_kron   # la préparation de Nathan écrase le r des transformateurs
```

### A3. Cohérence PV/PQ : bonne. Mais le slack produit plus que sa centrale — confiance élevée
- **Cohérence.** Sur les 168 heures, les tensions et états enregistrés vérifient les équations à 10⁻³ MVAr
  près ; 0 régulateur PV hors limites Q ; **une seule** violation de complémentarité (> 0,001 pu) : REG 131 à
  la pointe (+0,0094 pu), déjà connue. Le Q du slack reste dans ses limites (−127 à +901 MVAr).
- **Problème nouveau : la puissance active du slack.** Le LOPF est sans pertes. En AC, toutes les pertes
  (608 à 1 801 MW, **5,2 % de la charge à la pointe**) sont produites par le bus slack, LG-2 (bus 339), au bout
  nord du réseau. Avec le dispatch du LOPF, le bus 339 produit jusqu'à **9 360 MW pour 7 722 MW installés** ;
  il dépasse sa puissance installée pendant **101 heures sur 168**.
- **Conséquence.** À la pointe, 1,8 GW de plus que prévu descend de la Baie-James vers le sud. Testé en B5 :
  répartir ces pertes sur les centrales change peu la marge (0 à +2 %). Ce n'est pas la cause principale.
- **Inductances.** 42,0 GVAr installés à k = 1 dans le modèle, 13,1 à 32,4 GVAr enclenchés ; l'enclenché
  dépasse les 25,5 GVAr installés d'HQ (chiffre de 1996) pendant 41 heures. Voir A1b.
- **Preuve.** `a3_coherence.py`.

### A4. Sources des équipements réactifs — confiance moyenne (une source relue)
- **Vérifié.** Le tableau 2 de la demande R-3712-2009 (HQT-1, doc. 1, p. 8) liste bien les compensateurs
  statiques : Albanel 2 × (−105/+300), Chamouchouane 2 × (−300/+300), Chibougamau 2 × (−105/+300),
  Laurentides 1 × (−100/+350), Némiscau 2 × (−105/+300), La Vérendrye 2 × (−105/+300), Lévis (déglaceur)
  1 × (−125/+250). Les 7 premières lignes de `REACTIVE_DEVICES` sont conformes.
- **Écart.** Chénier : le même tableau donne 2 × (−105/+300) ; le code met 2 × (−300/+300) en citant Hitachi.
  Seul Qmin change (−600 au lieu de −210 MVAr), sans effet à la pointe.
- **Supposé, et déjà dit comme tel dans le code.** Les compensateurs synchrones : seul Lévis a une puissance
  sourcée. Duvernay (3 CS), Abitibi (2 CS) et Manicouagan (1 CS ?) sont à « −200/+300, supposé = Lévis ». Les
  dossiers de la Régie confirment l'existence et le nombre à Duvernay (3) ; je n'ai pas trouvé les Mvar.
  **Duvernay pèse lourd** : +900 MVAr supposés au cœur de Montréal, là où la tension est la plus basse.
- **Non relu** (faute de temps) : D-2005-45 (Lévis), R-3684-2009, R-3810-2012.
- **Pas dans la liste.** Le même document mentionne un compensateur statique à Rimouski (absent du modèle).

### A5. pandas 3.x : ça tourne, résultats identiques — confiance élevée
- **Constat.** Avec pandas 3.1.1 (et numpy 2.5.5, PyPSA 1.3.0) : `build()` complet avec LOPF, 7 heures de
  `solve_week_reserve`, `validate_with_pypsa` et `to_ppc` passent sans erreur. Écart avec `resultats/`
  (pandas 2.3.3) : 5·10⁻¹⁵ pu sur les tensions. Non testé : pandapower sous pandas 3.
- **Preuve.** `a5_pandas3.py`, `a5_pandas3.log`.

---

## Partie B — la marge nulle à la pointe est-elle réelle ?

**Réponse courte : oui pour ce dispatch, et c'est même pire que « nulle ». Mais ce dispatch n'est pas
réaliste : deux défauts du LOPF l'expliquent (B4, B5), et aucun des deux n'est une question de réactif.**

### B1. Ce n'est pas un artefact de l'algorithme PV/PQ : la solution de la pointe est sur la branche instable — confiance élevée
- **Constat.** Critère de stabilité V-Q classique : injecter du réactif à un bus doit faire monter sa tension
  (dV/dQ > 0 ; valeurs propres du jacobien réduit J_R toutes positives). Sur les 168 solutions enregistrées
  dans `resultats/` : **167 heures sont stables ; la pointe (3 janv. 17 h) ne l'est pas** : une valeur propre
  négative de J_R (−149) et **59 bus à dV/dQ négatif**, dont 131 (Beauharnois), 136 (Chénier), 1, 42.
  La solution de la pointe est sur la **branche basse** de la courbe P-V : Newton converge, les limites Q
  sont respectées, les tensions sont dans 0,90-1,10, mais ce n'est pas un point de fonctionnement.
- **Ça explique la violation de complémentarité de REG 131** (à Qmax avec une tension au-dessus de sa
  consigne) : sur la branche basse, plus de Q fait baisser la tension.
- **Où est le nez ?** Courbe P-V autour de la pointe (shunts re-manœuvrés à chaque pas) : dV/dQ au bus 800
  vaut 0,004 pu/100 MVAr à 85 % de la charge, 0,011 à 98 %, 0,023 à 99,5 %, puis change de signe à 100 %.
  **Le nez est entre 99,5 % et 100 % de la charge de pointe (≈ 34,8 GW)** : la marge de ce cas est légèrement
  négative, pas nulle. Les deux autres heures basses (3 janv. 18 h, 4 janv. 10 h) sont stables.
- **Donc le bon énoncé du résultat est : 167/168 heures convergées au sens des trois critères ET stables ;
  la pointe n'a pas de solution stable avec ce dispatch.**
- **Méthodes essayées et abandonnées** (dites ici pour qu'on ne les refasse pas) : OPF AC de pandapower
  (ne converge sur aucun cas, même à 80 % de charge avec Q libre : problème de mise en données, pas de
  physique) ; continuation avec état PV/PQ conservé (trouve des « solutions » à +6 %, mais avec des
  régulateurs bloqués à Qmax à des tensions de 1,2 pu : fausses) ; régulateurs à statisme (Newton ne
  converge pas près du nez). Le critère dV/dQ est plus simple et tranche.
- **Preuve.** `b1d_stabilite_vq.py` → `b1d_stabilite_vq.csv` ; `b1c_nez.py` → `b1c_courbe_pv_pointe.csv` ;
  `b1_marge.py`.
- **Correctif proposé** (non appliqué) : dans `solve_week_reserve`, calculer ce critère à chaque heure et
  l'ajouter à la définition de « convergé » (10 lignes : jacobien à la solution, `np.linalg.eigvals` de J_R).

### B2. Les postes 735 kV sans transformation : Kron étendu au sud — confiance moyenne
- **Constat.** La réduction de Nathan a coupé tous les chemins < 315 kV ; les équivalents de Kron ne les
  rendent qu'au Saguenay et à Montréal. Étendu à tout le sud (lat 44,9-47,6 ; 460 bus éliminés, 104 bus
  frontières, 667 branches équivalentes), et à dispatch corrigé (B4) pour séparer les deux effets :

| Variante (toutes avec fil de l'eau prioritaire) | lim. Q | 0,95-1,05 | V min-max | V min pointe | marge stricte |
|---|---|---|---|---|---|
| Kron Saguenay + Montréal (= B4) | 168 | 168 | 0,965-1,037 | 0,965 | 1,070 |
| Kron Saguenay + tout le sud, boîte complète | 168 | 168 | 0,974-1,037 | 0,974 | 1,110 |
| Kron Saguenay + tout le sud, Québec seulement | 168 | 168 | 0,972-1,037 | 0,972 | 1,100 |
| Kron sur tout le Québec, Québec seulement | — | — | — | — | — |

  La dernière ligne a échoué : Newton ne converge à aucune heure (le script s'arrête faute de colonne
  `in_hard`). Non investigué faute de temps ; l'équivalent couvre alors aussi le nord, où les postes
  735 kV n'ont pas de réseau < 315 kV utile.

- **Lecture.** Rendre au modèle les chemins < 315 kV du sud ajoute **~3 % de marge** (1,07 → 1,10). C'est
  la réponse à la question « les 20 postes sans transformation coûtent-ils de la marge ? » : oui, un peu.
- **Piège trouvé.** La boîte rectangulaire du sud contient aussi ~500 bus de l'Ontario, de l'État de New York
  et du Maine (le réseau d'`elec_full` y est relié en AC). Kron en fait des chemins parallèles qui n'existent
  pas en exploitation. Les retirer (frontière approchée par une ligne brisée) coûte 1 % de marge. Les boîtes
  actuelles ne sont pas concernées : 0 bus hors Québec dans celles du Saguenay (50 bus) et de Montréal (84 bus).
- Kron avec la règle hybride des longueurs (défaut 2 de A2) : effet faible (tableau ci-dessous).
- **Vérifié avec tous les contrôles** (`b7_recommande_kron_sud.py` → `b7_*.csv`, `b7.log`) : fil de l'eau
  prioritaire + Kron Saguenay et sud (Québec seulement) donne **168/168 sur les cinq critères** (Newton,
  limites Q, 0,90-1,10, 0,95-1,05, stabilité V-Q), V 0,972-1,037, 0 violation de complémentarité, et une
  **marge stable de +10 % (38,5 GW)** ; +8 % (37,8 GW) avec le réactif des seules machines en marche.
- **Preuve.** `b_variantes3.py` → `b_variantes_part3*.csv`.

### B3. Capacité réactive des centrales : trop généreuse dans le sud — confiance élevée sur le constat
- **Constat.** À la pointe, le sud (latitude < 47,6°) dispose de 6 311 MVAr de Qmax « installés » et en
  fournit 5 689. Mais **1 734 MVAr (30 %) viennent de Beauharnois (REG 131, 975 MVAr) et de Carillon et ses
  voisines (REG 3373, 759 MVAr), deux groupes que le LOPF laisse à 0 MW à cette heure**. Des machines à
  l'arrêt ne fournissent pas de réactif (sauf marche en compensateur synchrone, non documentée ici).
- **Test.** Limites Q recalculées avec les seules machines en marche (P > 1 MW) : le Qmax du sud tombe à
  4 102 MVAr à la pointe. Newton converge encore, mais **18 régulateurs sortent de leurs limites** à la
  pointe (13 à 17 aux autres heures chargées), V min 0,892. Le cas final ne tient donc à la pointe que
  grâce au réactif de centrales arrêtées.
- **Mais** le vrai problème est en amont : ces centrales ne devraient pas être à l'arrêt (B4).
- **Preuve.** `b3_q_centrales.py` → `b3_q_centrales.csv`.

### B4. Le dispatch du LOPF arrête le fil de l'eau du sud à la pointe — confiance élevée
- **Constat.** À la pointe, le LOPF produit **1 602 MW dans le sud pour 6 644 MW installés**. Il laisse à
  0 MW **3 440 MW de fil de l'eau disponible** : Beauharnois (1 825 MW disponibles), Carillon (746), Paugan
  (216), Chelsea (152), Les Cèdres (111), Rapides-Farmer (104)… pendant que les réservoirs du nord tournent
  à 27 346 MW sur 28 467. Charge du sud : 31 210 MW sur 34 959.
- **Cause.** Toute l'hydraulique a un coût marginal nul (réservoirs et fil de l'eau) : le LOPF est
  indifférent entre Beauharnois et LG-2, et le solveur choisit au hasard parmi les optimums de même coût.
  C'est le « optimum dégénéré » déjà signalé, mais son effet est bien plus gros que 0,01-0,04 pu : il
  ajoute ~3,4 GW de transit nord → sud à la pointe.
- **Pourquoi c'est faux physiquement.** Une centrale au fil de l'eau ne stocke pas : l'eau non turbinée est
  perdue. En hiver, Beauharnois tourne. Le modèle donne déjà la priorité à l'éolien (coût −0,1) ; il faut la
  même chose pour le fil de l'eau. Ce n'est pas un réglage pour « caler » la marge : c'est la règle
  d'exploitation.
- **Preuve.** `b4_dispatch_sud.py` → `b4_dispatch_pointe.csv`. Effet mesuré : tableau ci-dessous.

### Effet de chaque correction sur la semaine et la marge (une hypothèse à la fois) — confiance élevée

Chaque ligne refait tout : LOPF, 168 heures d'AC avec limites Q, marge à la pointe. « Marge stricte » =
critère de `R.load_margin` (limites Q + complémentarité) ; « marge lim. Q » = sans le critère de
complémentarité. Ces marges ne vérifient pas encore la stabilité V-Q (voir B6 pour la marge stable).

| Variante | Newton | lim. Q | 0,90-1,10 | 0,95-1,05 | V min-max | P slack pointe (MW) | V min pointe | marge stricte | marge lim. Q |
|---|---|---|---|---|---|---|---|---|---|
| Cas final | 168 | 168 | 168 | 165 | 0,922-1,035 | 1 801 | 0,922 | 1,000 | 1,015 |
| **B4 fil de l'eau prioritaire** | 168 | 168 | 168 | **168** | 0,965-1,037 | 1 465 | 0,965 | **1,070** | 1,075 |
| B4 coûts perturbés de 0,1 % (graine 1) | 168 | 168 | 168 | 168 | 0,965-1,035 | 1 509 | 0,966 | 1,075 | 1,080 |
| B4 coûts perturbés de 0,1 % (graine 2) | 168 | 168 | 168 | 168 | 0,964-1,035 | 1 469 | 0,964 | 1,060 | 1,065 |
| B4 coûts perturbés de 0,1 % (graine 3) | 168 | 168 | 168 | 168 | 0,965-1,036 | 1 506 | 0,966 | 1,065 | 1,075 |
| B4 coûts perturbés de 0,1 % (graine 4) | 168 | 168 | 168 | 168 | 0,950-1,034 | 1 724 | 0,950 | 1,030 | 1,040 |
| B5 charges réduites des pertes | 168 | 168 | 168 | 168 | 0,963-1,036 | 1 670 | 0,963 | 1,035 | 1,040 |
| B5 pertes réparties sur les centrales | 168 | 168 | 168 | 165 | 0,915-1,034 | 85 | 0,915 | 1,000 | 1,020 |
| B5 les deux | 168 | 168 | 168 | 168 | 0,963-1,035 | 86 | 0,964 | 1,040 | 1,055 |
| B4 + B5 (pertes réparties) | 168 | 168 | 168 | 168 | 0,964-1,037 | 28 | 0,964 | 1,080 | 1,095 |
| B4 + B5 (les deux) | 168 | 168 | 168 | 168 | 0,971-1,035 | 13 | 0,972 | **1,115** | 1,160 |
| B4 + B5 (les deux), sans E3 | 168 | 168 | 168 | 168 | 0,972-1,035 | 15 | 0,972 | 1,105 | 1,110 |
| B2 Kron avec la règle hybride | 168 | 168 | 168 | 167 | 0,915-1,035 | 1 697 | 0,955 | 1,030 | 1,040 |
| A1b sans E3 (seul) | 168 | **166** | 166 | 165 | 0,932-1,033 | 1 885 | 0,956 | 1,015 | 1,015 |

(Fichiers : `b_variantes.py`, `b_variantes2.py` → `b_variantes_part*.csv`. Pour les variantes « charges
réduites », la charge de pointe du modèle est 33,2 GW au lieu de 35,0 ; λ multiplie la demande d'HQ dans les
deux cas, donc les marges se comparent directement.)

**Lecture.**
- **Le dispatch est le premier facteur.** Quatre perturbations aléatoires de 0,1 % des coûts (graines 1 à 4)
  donnent toutes 168/168 dans 0,95-1,05, une pointe à 0,950-0,966 pu et une marge de 1,03 à 1,075 : le cas
  final est tombé sur l'optimum dégénéré le plus défavorable des cinq. La marge dépend donc du hasard du
  solveur à ±2-3 %. Donner la priorité au fil de l'eau fixe ce choix de façon physique et reproductible, et
  se place dans le haut de la fourchette (1,07).
- **B5, les pertes.** Le LOPF est sans pertes, puis toutes les pertes AC (1,5-1,8 GW à la pointe) sont produites
  par LG-2. Deux corrections indépendantes :
  - *pertes réparties* sur les centrales en marche (au prorata de leur production, comme le réglage
    fréquence-puissance) : retire ~1,8 GW du transit nord → sud, mais change peu la marge (seul : 1,015 → 1,020 en
    « marge lim. Q » ; avec B4 : 1,075 → 1,095) ;
  - *charges réduites des pertes* : si la demande horaire d'HQ utilisée par Nathan contient déjà les pertes
    de transport (à vérifier dans la définition du jeu de données ; les chiffres de transport d'HQ sont publiés
    « avec pertes de transport »), le modèle compte les pertes deux fois (une fois dans la charge, une fois
    dans l'AC). Effet : +3,5 % de marge.
- **E3 (circuits de Nathan) n'est pas nécessaire une fois le dispatch corrigé** (B4 + B5 sans E3 : 168/168/168,
  marge 1,105 contre 1,115), et il ajoute 2 700 circuit-km que les données d'HQ ne justifient pas (A1b).
  Sans E3, les inductances enclenchées tombent à 26,8 GVAr au maximum (HQ : 25,5 installés) au lieu de 34.
- **Kron avec la règle hybride (B2)** : effet faible sur la semaine (167 au lieu de 165 dans 0,95-1,05) ; on
  devrait le corriger par cohérence, ce n'est pas lui qui fait la marge.

### B6. Configuration recommandée : cas final + fil de l'eau prioritaire — confiance élevée
Vérifiée avec tous les contrôles de la revue (`b6_recommande.py` → `b6_*.csv`, `b6.log`) :
- **168/168 Newton, 168/168 dans les limites Q, 168/168 dans 0,90-1,10, 168/168 dans 0,95-1,05, 168/168
  stables (V-Q)** ; V 0,965-1,037 pu ; 0 violation de complémentarité ; LOPF sans délestage.
- Le fil de l'eau du sud produit 4 844 MW à la pointe sur 4 929 installés (contre 1 405 dans le cas final).
- Bus slack (LG-2) : dépasse sa puissance installée 3 heures sur 168 (contre 101), de 640 MW au plus.
- Inductances : 16,5 à 33,8 GVAr enclenchés (toujours au-dessus des 25,5 d'HQ la nuit : A1b, E3).
- Réactif demandé à des machines arrêtées : il en reste (287 MVAr sur 7 régulateurs à la pointe). Mais avec
  les limites Q des **seules machines en marche**, la pointe converge aussi : 0 régulateur hors limites,
  V 0,960-1,035, stable.
- **Marge stable à la pointe** (Newton + limites Q + stabilité V-Q, shunts re-manœuvrés à chaque pas de 1 %) :
  **+7 % (37,4 GW)** avec les limites Q actuelles ; **+6 % (37,1 GW)** avec le réactif des seules machines en
  marche. Avec en plus V ≥ 0,90 : +6 % et +5 %. Le bus limitant reste le 800 (Montréal-Ouest).
- C'est encore loin des +23 % du réseau réel (43,1 GW en 2023), mais ce n'est plus un réseau au bord de
  l'effondrement. **Option plus forte** : y ajouter Kron étendu au sud (B2, `b7`) → marge stable +10 %. Avec aussi B5 (charges réduites des pertes, pertes réparties) : marge stricte 1,115 sans
  stabilité vérifiée.

**Correctif proposé** (non appliqué) :
```diff
--- reseau315.py  (lopf)
 def lopf(n, solver="highs", log=None):
     ...
+    # Fil de l'eau prioritaire : une centrale au fil de l'eau ne stocke pas d'eau ; à coût nul comme les
+    # réservoirs, le LOPF (optimum dégénéré) pouvait arrêter Beauharnois et Carillon à la pointe (B4).
+    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
     t0 = time.time()
```
et ajouter à la définition de « convergé » le critère de stabilité V-Q (`b1d_stabilite_vq.py`, ~10 lignes).

---

## Partie C — descendre sous 315 kV (exploratoire, dossier `reseau_complet/`)

**Réponse courte : ça converge à la pointe avec les limites Q, mais seulement après avoir laissé au 315 kV
7,2 GW de charge sur 35, et ~10 bus restent hors de 0,90-1,10. Ça ne remplace pas encore Kron, et la marge
ne remonte pas. Ce que ça apprend surtout : la répartition des charges de Nathan ne tient pas physiquement
sous 315 kV.** Confiance : moyenne (≈ 40 min de travail).

### Construction (`complet.py`, `build_complet`)
- Base : le modèle 315 kV de `R.build()` **sans** Kron (E1-E6, même LOPF, avec le fil de l'eau prioritaire
  de B4), plus le vrai réseau < 315 kV d'`elec_full.nc` : 597 bus, 716 lignes (12 474 km), 133 transformateurs.
  Mêmes corrections que dans `reseau315` : longueurs par la règle hybride, types 50 Hz (x et b × 1,2).
- Les charges reviennent à leur bus d'origine : `elec_full` et `elec_solved` ont les mêmes noms de charges,
  et la réduction de Nathan avait déplacé 539 charges sur 647 (30,3 GW sur 35 à la pointe) vers le bus
  ≥ 315 kV le plus proche. Les bancs de condensateurs sont replacés aux nouveaux bus de charge.
- Moteur : `SparseEngine` = `ACEngine` avec un Newton creux (796 bus).

### Constats
1. **Brut, ça ne converge pas** : même avec un réactif illimité, Newton diverge au-delà de **25 % de la
   charge de pointe** (`c4_rampe.py`, `c4b_ou.py`).
2. **Cause 1 : des charges physiquement impossibles** (`c3_diagnostic.py` → `c3_bus_charge.csv`). La
   répartition régionale de la demande (`rescale_demand_regional.py`) met de très grosses charges sur de
   petits bus : bus 3559 (69 kV) **851 MW au bout d'une seule ligne 69 kV dont la limite V²/x est 835 MW** ;
   bus 38 (69 kV, Île-Perrot) 1 046 MW ; bus 1045 (120 kV, antenne) 760 MW. 17 bus < 315 kV ont plus de
   300 MW. Ces charges ne passent pas, quel que soit le réactif.
3. **Cause 2 : des poches régionales faibles** : 120 kV des Hautes-Laurentides / Outaouais (bus 3204, 1228,
   1253 ; 46,4° N, −75,5°) et une antenne 230 kV au Témiscamingue (2734, 2712 ; −78 à −79°).
4. **Après élagage** (`c5_elagage.py`) : on laisse au 315 kV les charges des 8 bus « impossibles » (critère
   charge > 10 % de la somme des V²/x de leurs branches) et de 36 bus qui s'effondrent quand on monte la
   charge (liste : `c5_bus_exclus_ror.csv`, surtout du 120 kV autour de Montréal et de l'Outaouais) :
   **7,2 GW restent au 315 kV, 23,8 GW sont redescendus**. La pointe converge alors.
5. **Résultat à 6 heures de la semaine** (`c6_resultat.py` → `c6_heures_ror.csv`) :

| Heure | GW | Newton | régulateurs hors lim. Q | V min-max (tous bus) | bus hors 0,90-1,10 | V min-max 315-735 kV | stable V-Q | pertes (MW) |
|---|---|---|---|---|---|---|---|---|
| 3 janv. 17 h (pointe) | 35,0 | oui | 0 | 0,890-1,103 | 8 | 0,923-1,042 | oui | 2 394 |
| 3 janv. 18 h | 34,5 | oui | 0 | 0,904-1,104 | 10 | 0,940-1,042 | oui | 2 314 |
| 4 janv. 10 h | 34,2 | oui | 0 | 0,911-1,104 | 10 | 0,948-1,041 | oui | 2 276 |
| 2 janv. 16 h | 30,2 | oui | 0 | 0,947-1,108 | 10 | 0,977-1,046 | oui | 1 715 |
| 5 janv. 4 h | 27,4 | oui | 0 | 0,963-1,112 | 10 | 0,985-1,048 | oui | 1 422 |
| 1er janv. 3 h | 22,0 | oui | 0 | 0,979-1,125 | 11 | 0,987-1,053 | oui | 960 |

   - Au sens des trois critères : **Newton oui, limites Q oui, tensions non** (8 à 11 bus hors bande à chaque
     heure). Les surtensions (1,10-1,125 pu) sont au 161/230 kV de l'Abitibi-Témiscamingue : lignes longues peu
     chargées, sans inductance (le modèle n'en met qu'au 735 kV et sur les lignes 315 kV de plus de 100 km).
     La sous-tension est un bus 120 kV (3051).
   - Les tensions au 315-735 kV sont plus basses que dans le modèle 315 kV de même dispatch (pointe : 0,923
     contre 0,965) et les pertes plus fortes (2,4 GW contre 1,5 GW) : 24 GW transitent maintenant par le
     réseau < 315 kV, dont les transformateurs sont tous génériques (2 000 MVA, x = 0,1 pu, r = 0 dans
     `elec_full`, mis à x / 50 ici).
6. **Marge** : aucune. À la pointe, un bus est déjà sous 0,90 et +1 % de charge n'a plus de solution dans
   les limites Q (`c6_marge_ror.csv`).

### Ce qu'il faudrait pour que le réseau complet remplace Kron
1. Une vraie répartition des charges par poste (le problème principal, et il touche aussi le modèle 315 kV :
   ces charges-là y sont simplement agrégées au 315 kV, où elles ne posent pas de problème).
2. Des inductances sur les longues lignes 161/230 kV, et les vrais transformateurs (puissance, x par poste).
3. Les 45 bus < 315 kV non reliés au réseau principal (îlots d'`elec_full`) à examiner.
D'ici là, **les équivalents de Kron (corrigés, A2) restent le meilleur choix** pour le 315 kV : ils sont
exacts pour la partie passive du réseau et ne dépendent pas de la répartition des charges.

---

## Fichiers de ce dossier
| Fichier | Rôle |
|---|---|
| `_commun.py` | import de `reseau315` (dossier voisin), chargement du réseau final |
| `a1_longueurs.py`, `a1b_km735.py` | A1 : longueurs ; circuit-km à 735 kV |
| `a2_kron.py`, `a2b_kron_diag.py` | A2 : vérification numérique de Kron |
| `a3_coherence.py` | A3 : cohérence PV/PQ, slack, inductances |
| `a5_pandas3.py` (+ `.log`) | A5 : test sous pandas 3 (`PYTHONPATH=/tmp/pd3`) |
| `b1d_stabilite_vq.py`, `b1c_nez.py` | B1 : critère de stabilité V-Q, courbe P-V autour de la pointe |
| `b1_marge.py`, `b1b_statisme.py`, `b_outils.py` | B1 : essais abandonnés (gardés pour mémoire) |
| `b3_q_centrales.py`, `b4_dispatch_sud.py` | B3, B4 |
| `b_variantes.py`, `b_variantes2.py`, `b_variantes3.py` | tableau des variantes (B2, B4, B5, A1b) ; Kron étendu au sud |
| `b6_recommande.py` (+ `b6.log`) | B6 : configuration recommandée, tous les contrôles |
| `b7_recommande_kron_sud.py` (+ `b7.log`) | idem + Kron étendu au sud (Québec seulement) |
| `reseau_complet/complet.py` | C : construction du réseau complet + Newton creux |
| `reseau_complet/c0` à `c6` | C : inventaire, diagnostics, élagage, résultats |

Le fichier `reseau_complet/_complet_ror.pkl` (45 Mo, réseau complet sérialisé) n'est pas copié : il se refait
avec `c3b_dc.py`.
