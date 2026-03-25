#import "@preview/unequivocal-ams:0.1.2": ams-article, theorem, proof
#set math.equation(numbering: "(1)")
#set text(lang: "fr")

#show: ams-article.with(
  title: [Une méthode de différences finies pour la résolution numérique de l'équation de Fokker-Planck du modèle LIF],
  authors: (
    (
      name: "Jules Herrmann",
      department: [Master Mathématiques de la modélisation],
      organization: [Sorbonne Université],
      email: "jules.herrmann@etu.sorbonne-universite.fr",
    ),
  ),
  bibliography: bibliography("refs.bib"),
)

= Introduction

Le modèle Intergration-et-tir avec fuite, en anglais Leaky Integrate and Fire (LIF), est un modèle classique
en neuroscience. Introduit au début du XXeme siècle par Louis Lapicque, il décrit l'évolution du potentiel électrique
d'un neurone au cours du temps. Ce modèle a pu être étendu et généralisé, notament par l'ajout de 
comportement stochastique. Ainsi, son évolution non linéaire bruitée (NNLIF), basée sur
une équation différentielle stochastique, modélise les comportements de champ moyen d'un ensemble de neurones.

En l'absence de tir, l'évolution du potentiel $V$ est régie par l'équation stochastique suivante :

$
d V = -(V - V_L) d t + mu d t + sigma d B_t
$

Le terme $-(V - V_L) d t$ correspond à la relaxation du potentiel vers un potentiel de repos $V_L$.

Le terme $mu d t + sigma d B_t$ correspond à l'influence du courant synaptique généré par le reste des neurones, qui
comporte un terme de _drift_, controlé par $mu$, et un terme stochastique controlé par $sigma$ ($B_t$ est un mouvement Brownien).

En revanche, lorsqu'un neurone atteint le potentiel de tir $V_F$, son potentiel est réinitialisé à la valeur $V_R$.

Il est possible de déduire de ce modèle une équation différentielle de type Fokker Planck décrivant l'évolution de la 
densité de probabilité $p(v,t)$ de la variable aléatoire $V$.
$
cases(
  partial_t p + partial_v (h p) - a partial_(v v)p = 0  & forall v in \]-oo \, V_F \[\\{V_R},
  p(v,0) = p_0(v) & forall v in \]-oo \, V_F \[,
  p(- oo,t) = p(V_F,t) = 0 & forall t > 0,
  p(V_R^-,t) = p(V_R^+,t) & forall t > 0,
  partial_v p(V_R^-,t) = partial_v p(V_R^+,t) + N(t)/a & forall t > 0,
  N(t) = - a partial_v p(V_F,t) & forall t > 0
)
$ <eq:fokkerplanck>

- $p_0$ est la densité de probabilité initiale de $V$
- $h$ est le terme de drift
- $a$ est le coefficient de diffusion
- $N(t)$ est le taux moyen de tir au temps $t$

Le taux de tir $N$ du groupe de neurone a un impact sur le comportement de chaque neurone, ainsi, dans un réseau excitatoire,
un taux de tir élevé a pour conséquence un accroissement de la force du drift. En revanche, dans un réseau inhibiteur, la force
du drift diminue avec l'augmentation de $N$. Ainsi, les termes $a$ et $h$ sont
souvent définis par les relations suivantes #cite(<caceres_analysis_2011>):
$
h(v,N(t)) = -v + b N(t)
$
$
a(N(t)) = a_0 + a_1 N(t)
$

Avec ce choix, $b < 0$ correspond à un réseau inhibiteur, et $b>0$ à un réseau excitateur.

Dans leur article de 2011, #cite(<caceres_analysis_2011>,form:"prose"), ont posé les bases de l'étude de ce système. Ils ont caractérisé formellement le phénomène de _blow up_, dans lequel le taux de tir $N(t)$ diverge en temps fini. Si les blow up avaient été observés numériquement par le passé, ils ont prouvé leur existence dans le cas continu dans une famille de situations initiales. Ils ont également montré l'existence de situations présentant 0, 1 ou 2 solutions stationnaires.

De nombreuses propriétés de ce système sont toujours inconnues. La stabilité des solutions stationnaires et les comportements asymptotiques du système sont mal comprises. Dans le cas linéaire, il est possible de montrer que l'entropie relative de la solution est strictement décroissante, ce qui permet de montrer la convergence de la solution vers une solution stationnaire dans certains cas.

Ce manque de résultat théorique complique le développement de méthodes numériques capable d'approximer ce problème.

Une approche générale utilisée pour résoudre numériquement les équations de Fokker-Planck est l'approche dite de Scharfetter-Gummel. Cette approche reformule l'équation, comme une équation de balance des flux. #cite(<almeida_energy_2018>)

Dans leur article de 2020, #cite(<hu_structure_2020>,form:"prose") proposent une méthode numérique basée sur l'approche de Scharfetter-Gummel, et l'adapte pour prendre en compte le saut de flux. Cette approche satisfait plusieurs bonnes propriétés : elle est linéairement implicite et donc facile à résoudre numériquement, conserve la positivité de la densité et satisfait une version discrète de la propriété de l'entropie relative décroissante qui existe dans le cas continu.

Ce rapport présente premièrement la méthode de discrétisation étudiée par l'article.
Puis le détail de l'algorithme obtenu en utilisant la discrétisation proposée est exposé.
Une seconde partie se donne pour but de répliquer les expériences numériques de l'article et se concentre plus particulièrement sur le phénomène des _blow up_ et sur la décroissance de l'entropie relative.

Malheureusement, certaines tentatives de réplications des résultats de l'article n'ont pas donné des résultats identiques. Ces résultats ne sont pas inclus dans la @resultat, mais sont discutés dans la @conclusion.

#pagebreak()

= Méthode

== Formulation en équilibre de flux

La première ligne de @eq:fokkerplanck peut être reformulée en termes d'un flux.

$ partial_t p + partial_v F = 0 $ <eq:equiflux>

où 

$ F(v,t) = -a partial_v p + h p $

Cependant, le flux sortant en $V_F$, qui correspond à $N(t)$ doit sauter en $V_R$. Ce comportement peut être pris en compte en considérant le flux modifié :
$ tilde(F)(v,t) = F(v,t) - N(t) bb(1)_(v >= V_R) $ <eq:F_tilde>

== Reformulation de Scharfetter-Gummel

On pose 
$
M(v,t) = exp(- ((v - b N(t))^2)/(2 a(N(t))))
$

On remarque que 
$
F(v,t) = a(N(t)) M(v,t) partial_v (p(v,t)/M(v,t))
$ <eq:F_reform>

== Discrétisation

On choisit une valeur $V_("min")$ telle que $p(V_"min",t) approx 0$. On discrétise l'intervalle $\[V_"min" \, V_F \]$ en $n$ intervalles de longueur $h$, délimités par les points $v_0 \, v_1 dots, v_n$. Ce sont les points sur lesquels $p$ sera discrétisé. Le flux $F$, sera lui évalué sur les points intermédiaires $v_(1/2) \, v_(1+ 1/2), dots , v_(n- 1/2)$.

#figure(
  image("figs/discretisation.pdf", width: 100%),
  caption: [Schéma de discrétisation],
) <discretisation>

On discrétise le temps sur l'intervalle $\[0 \, T \]$ avec un pas de $tau$. On note alors $F_i^m$ et $p_i^m$ la valeur du flux et de la densité au point correspondant dans la discrétisation.

== Différence finie

@eq:equiflux peut être discrétisé avec une différence finie centrée $forall i in [|0,n|]$

$
(p_i^(m+1) - p_i^m)/ tau + (tilde(F)_(i+1/2)^m - tilde(F)_(i-1/2)^m)/h = 0
$

On discrétise également l'expression de $F$ dans @eq:F_reform avec une différence finie centrée. Avec @eq:F_tilde, on obtient, $forall i in [|1,n-2|]$,
$
tilde(F)_(i+1/2)^m = -a(N^m) M_(i+1/2)/h (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i))
- N^m bb(1)_(v_(i+1/2) >= V_R)
$<eq:F_discretisation>

$
tilde(F)_(1/2) = tilde(F)_(N-1/2) = 0
$

On remarque que, dans @eq:F_discretisation, les valeurs de $p$ sont prises au temps $m+1$, et les valeurs de $M$ au temps $i$. Cela permet d'obtenir un schéma implicite, qui nécessite la résolution d'une équation linéaire à chaque étape.

== Algorithme

Cette discrétisation peut être manipulée (voir @annexe) pour obtenir la procédure suivante pour calculer une itération de l'algorithme.

1. $N^m <- a_0 * p_(n-1)^m / (h - a_1 * p_(n-1)^m)$
2. $p_(n-1)^m <- p_(n-1)^m - tau/h N^m$
3. $p_(V_R)^m <- p_(V_R)^m + tau/h N^m$
4. Résoudre $A_m p^(m+1) = p^m$

Avec
$A_m = 
mat(
1 + alpha M_(1+1/2)/M^m_(1) ,
- alpha M_(1+1/2)/M^m_(2);

- alpha M_(1+1/2)/M^m_(1),
1 + alpha (M_(2+1/2) + M_(1+1/2))/M^m_(2),
- alpha M_(2+1/2)/M^m_(3);

,dots.down,dots.down,dots.down;

,,- alpha M_(N-3+1/2)/M^m_(N-3),
1 + alpha (M_(N-2+1/2) + M_(N-3+1/2))/M^m_(N-2),
- alpha M_(N-2+1/2)/M^m_(N-1);

,,,-alpha M_(N-2+1/2)/M^m_(N-2),
1 + alpha M_(N-2+1/2)/M^m_(N-1)
)$

== Implémentation

Cet algorithme a été implémenté en Python avec les librairies Numpy et Scipy. La matrice $A$ est stockée dans une structure de données aux lignes compressées, afin de profiter de son caractère creux.

== Distribution Stationnaire

Une solution stationnaire de @eq:fokkerplanck est une solution constante en temps. Une telle solution stationnaire satisfait donc l'équation différentielle 

$
cases(
  partial_v (h(v,N^oo) p^oo) - a(N^oo) partial_(v v)p^oo = 0  & forall v in \]-oo \, V_F \[\\{V_R},
  p^oo (- oo) = p^oo (V_F) = 0,
  p^oo (V_R^-) = p^oo (V_R^+),
  partial_v p^oo (V_R^-) = partial_v p^oo (V_R^+) + N^oo /a(N^oo),
  N^oo = - a(N^oo) partial_v p(V_F)
)
$ <eq:stationary>

D'après #cite(<caceres_analysis_2011>), si elle existe, une telle solution s'écrit sous la forme 
$
p^oo = N^oo /a(N^oo) exp(- h(v,N^oo )^2 /(2 a(N^oo))) integral_"max"(v,V_R)^V_F exp( h(omega ,N^oo )^2 /(2 a(N^oo))) d omega
$

Le même article à identifié qu'il peut exister 0,1 ou 2 valeurs $N^oo$ telles que $p^oo$ définit bien une densité de probabilité.

== Entropie relative

L'entropie relative de @eq:fokkerplanck est définie par
$
S(t) = sum_(i=1)^(n-1) h G(p_i^t/p_i^oo ) p_i^oo
$

où $G(x) = 1/2 (x-1)^2$ et $p^oo$ satisfait @eq:stationary.

Dans le cas $b=0$, il y a existence et unicité de l'équation stationnaire. Il n'y a donc pas de problème pour définir $S$.
Toujours dans ce cas, le théorème 3.4 de #cite(<hu_structure_2020>) montre que $S$ est décroissante au cours du temps.

#pagebreak()

= Résultats <resultat>

== Blow up

#cite(<caceres_analysis_2011>) a montré que la solution de l'équation peut exploser en temps fini. @fig:blowup_bottom et @fig:blowup_top montrent deux exemples de situations où un blow up a lieu. Ces blow up ont lieu dans le cas excitatoire, où le taux de tir est régi par une boucle de rétroaction positive. Une explication intuitive du phénomène est que lors d'un blow up, le taux de tir augmente, ce qui augmente la force du drift, et donc la vitesse à laquelle les potentiels se déplacent vers le potentiel de tir, ce qui contribue à augmenter le taux de tir.


#figure(
  image("figs/fig1_bottom.pdf", width: 100%),
  caption: [Blow up de la solution pour $a(N(t)) = 1$, $b=1.5$, avec condition initiale gaussienne avec $v_0 = 1.5$ et $sigma^2_0 = 5 times 10^(-3)$. \
  (Gauche) Évolution du taux de tir $N(t)$ au cours du temps. \
  (Droite) $p(v,t)$ à $t = 0.0325 \, 0.0365 \, 0.0405$. (#link("https://www.youtube.com/watch?v=jF3e-Q7j1z0")[vidéo])],
)<fig:blowup_bottom>

#figure(
  image("figs/fig1_top.pdf", width: 100%),
  caption: [Blow up de la solution pour $a(N(t)) = 1$, $b=3$, avec condition initiale gaussienne avec $v_0 = -1$ et $sigma^2_0 = 0.5$.\
  (Gauche) Évolution du taux de tir $N(t)$ au cours du temps.\
  (Droite) $p(v,t)$ à $t = 2.95\, 3.15\, 3.35$. (#link("https://www.youtube.com/watch?v=Sw56snxPUMg")[vidéo])],
)<fig:blowup_top>

== Entropie relative

La @fig:entro_lineaire montre un exemple de solution dans le cas $b=0$. Il existe une unique solution stationnaire du système $p^oo$, avec $N^oo = 0.1377$.
Comme prédite dans le théorème 3.4 de #cite(<hu_structure_2020>), l'entropie relative de la solution par rapport à $p^oo$ est décroissante.
Cela se traduit par une convergence de la solution vers $p^oo$.

#figure(
  image("figs/fig4.pdf", width: 140%),
  caption: [Décroissance de l'entropie pour $a(N(t)) = 1$, $b=0$, avec condition initiale gaussienne avec $v_0 = 0$ et $sigma^2_0 = 0.25$.\
  (Gauche) Évolution de la densité au cours du temps. L'unique distribution stationnaire avec $N^oo = 0.1377$ est indiquée en pointillés.\
  (Milieu) Évolution du taux de tir $N(t)$ au cours du temps.\
  (Droite) Évolution de l'entropie relative $S(t)$ au cours du temps.
]
)<fig:entro_lineaire>

Dans le cas $b=1.5$, il n'y a plus unicité de la solution stationnaire. Deux distributions stationnaires sont connues. La première, $p^oo_"stable"$, avec $N^oo = 0.1924$ est stable.
La seconde, $p^oo_"instable"$, avec $N^oo = 2.319$ est instable. On peut voir que seule l'entropie relative à $p^oo_"stable"$ est décroissante, ce qui témoigne du fait que la solution converge vers $p^oo_"stable"$.

#figure(
  image("figs/fig5.pdf", width: 140%),
  caption: [Décroissance de l'entropie pour $a(N(t)) = 1$, $b=1.5$, avec condition initiale gaussienne avec $v_0 = 0$ et $sigma^2_0 = 0.25$.\
  (Gauche) Évolution de la densité au cours du temps. Il existe deux distributions stationnaires, $p^oo_"stable"$ et $p^oo_"instable"$, avec respectivement $N^oo = 0.1924$ et $N^oo = 2.319$.\
  (Milieu) Évolution de l'entropie relative à la distribution stationnaire $p^oo_"stable"$.\
  (Droite) Évolution de l'entropie relative à la distribution stationnaire $p^oo_"instable"$],
)

#pagebreak()

= Conclusion <conclusion>

Ainsi, la méthode numérique présentée par #cite(<hu_structure_2020>) donne lieu à un algorithme assez simple, et son implémentation donne des résultats cohérents avec ceux de l'article original. Les phénomènes de blow up sont identiques, et on observe que la propriété de décroissance de l'entropie relative est vérifiée expérimentalement.

En revanche, l'estimation de l'ordre de convergence de la méthode ne donne pas de résultats cohérents. L'ordre de convergence en espace calculé expérimentalement est différent de celui exposé dans l'article. Pour ce qui est de l'ordre de temps, l'implémentation présentée dans ce rapport n'est même pas stable pour les paramètres utilisés dans l'article original.

Il y a fort à parier que mon implémentation est érronée. Malheureusement, le code source utilisé pour les expériences numériques de #cite(<hu_structure_2020>) n'a pas été publié, ce qui complique la comparaison des deux méthodes.

Je pense que les domaines des mathématiques appliquées dans lesquels les expériences numériques jouent un rôle important devraient prendre exemple des mesures qui ont été prises dans d'autres domaines scientifiques pour faire face à la crise de la réplicabilité, en normalisant le fait de publier tout code source ayant permis de produire des résultats expérimentaux.

#pagebreak()

= Annexe
<annexe>

Pour rappel, le schéma a été discrétisé de la manière suivante :

$(p_i^(m+1) - p_i^m)/ tau + (tilde(F)_(i+1/2)^m - tilde(F)_(i-1/2)^m)/h = 0$

$tilde(F)_(i+1/2)^m = -a(N^m) M_(i+1/2)/h (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i))
- N^m bb(1)_(v_(i+1/2) >= V_R)$

$tilde(F)_(1/2) = tilde(F)_(N-1/2) = 0$

On cherche à obtenir une expression de $p^(m+1)$ en fonction de $p^m$.

$p_i^m =  p_i^(m+1) + tau/h tilde(F)_(i+1/2)^m - tau/h tilde(F)_(i-1/2)^m$

En subsituant $tilde(F)$ avec son expression,

$p_i^m =  p_i^(m+1) \
+ tau/h (-a(N^m) M_(i+1/2)/h (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i)) 
- N^m bb(1)_(v_(i+1/2) >= V_R)) \
- tau/h (-a(N^m) M_(i-1/2)/h (p^(m+1)_(i)/M^m_(i) - p^(m+1)_(i-1)/M^m_(i-1))
- N^m bb(1)_(v_(i-1/2) >= V_R))$

En simplifiant et en faisant passer le terme de saut de flux du côté gauche de l'équation

$p_i^m + tau/h N^m (bb(1)_(v_(i+1/2) >= V_R) - bb(1)_(v_(i-1/2) >= V_R)) =  p_i^(m+1) \
- tau/h a(N^m) M_(i+1/2)/h (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i)) \
 + tau/h a(N^m) M_(i-1/2)/h (p^(m+1)_(i)/M^m_(i) - p^(m+1)_(i-1)/M^m_(i-1))$

On remarque que $bb(1)_(v_(i+1/2) >= V_R) - bb(1)_(v_(i-1/2) >= V_R) = bb(1)_(v_i = V_R)$, donc,

$p_i^m + tau/h N^m bb(1)_(v_i = V_R) =
p_i^(m+1) \
- tau/h^2 a(N^m) M_(i+1/2) (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i)) \
 + tau/h^2 a(N^m) M_(i-1/2) (p^(m+1)_(i)/M^m_(i) - p^(m+1)_(i-1)/M^m_(i-1))$

En factorisant les termes de $p^(m+1)$,

$p_i^m + tau/h N^m bb(1)_(v_i = V_R) = \
- tau/h^2 a(N^m) M_(i+1/2)/M^m_(i+1) p^(m+1)_(i+1) \
+ (1 + tau/h^2 a(N^m) (M_(i+1/2) + M_(i-1/2))/M^m_(i)) p_i^(m+1) \
 - tau/h^2 a(N^m) M_(i-1/2)/M^m_(i-1) p^(m+1)_(i-1)$

On s'intéresse maintenant au cas particuliers où $i=1$

$ p_1^m & =  p_1^(m+1) + tau/h tilde(F)_(1+1/2)^m - underbrace(tau/h tilde(F)_(1/2)^m ,=0) \

      & =  p_1^(m+1) - tau/h^2 a(N^m) M_(1+1/2) (p^(m+1)_(2)/M^m_(2) - p^(m+1)_(1)/M^m_(1))
- underbrace(N^m bb(1)_(v_(1+1/2) >= V_R),=0) \

& = - tau/h^2 a(N^m) M_(1+1/2)/M^m_(2) p^(m+1)_(2) + (1 + tau/h^2 a(N^m) M_(1+1/2)/M^m_(1))p^(m+1)_1$

On s'intéresse maintenant au cas particulier où $i=n-1$

$p_(N-1)^m & =  p_(N-1)^(m+1) + underbrace(tau/h tilde(F)_(N-1/2)^m,=0) - tau/h tilde(F)_(N-1-1/2)^m \
          & =  p_(N-1)^(m+1) - tau/h tilde(F)_(N-2+1/2)^m \
          & =  p_(N-1)^(m+1) + tau/h^2 a(N^m) M_(N-2+1/2) (p^(m+1)_(N-1)/M^m_(N-1) - p^(m+1)_(N-2)/M^m_(N-2))
+ tau/h N^m underbrace(bb(1)_(v_(N-2+1/2) >= V_R),=1)$

En faisant passer le terme de saut de flux du côté gauche de l'équation, on a bien un système linéaire.

$p_(N-1)^m - tau/h N^m 
 & = (1 + tau/h^2 a(N^m) M_(N-2+1/2)/M^m_(N-1)) p^(m+1)_(N-1) - tau/h^2 a(N^m) M_(N-2+1/2)/M^m_(N-2) p^(m+1)_(N-2)$


Toutes ces équations peuvent être résumée en un seul système linéaire.

$A_m p^(m+1) = tilde(p)^m$

Avec

$A_m = mat(
1 + alpha M_(1+1/2)/M^m_(1) ,
- alpha M_(1+1/2)/M^m_(2);

- alpha M_(1+1/2)/M^m_(1),
1 + alpha (M_(2+1/2) + M_(1+1/2))/M^m_(2),
- alpha M_(2+1/2)/M^m_(3);

,dots.down,dots.down,dots.down;

,,- alpha M_(N-3+1/2)/M^m_(N-3),
1 + alpha (M_(N-2+1/2) + M_(N-3+1/2))/M^m_(N-2),
- alpha M_(N-2+1/2)/M^m_(N-1);

,,,-alpha M_(N-2+1/2)/M^m_(N-2),
1 + alpha M_(N-2+1/2)/M^m_(N-1)
)$

Où on a noté $alpha = tau/h^2 a(N^m)$

Et avec 

$p^(m+1) = 
mat(
  p^(m+1)_1;
  p^(m+1)_2;
  dots.v;
  p^(m+1)_(N-2);
  p^(m+1)_(N-1);
)$
$
tilde(p)^m =
mat(
  p^(m)_1;
  p^(m)_2;
  dots.v;
  p^m_(V_R) - tau/h N^m;
  dots.v;
  p^(m)_(N-2);
  p^(m)_(N-1) + tau/h N^m;
)$

#pagebreak()
