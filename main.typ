#import "@preview/unequivocal-ams:0.1.2": ams-article, theorem, proof
#set math.equation(numbering: "(1)")
#show: ams-article.with(
  title: [Une méthode des différences finies pour la résolution numérique de l'équation de Fokker-Planck du modèle LIF],
  authors: (
    (
      name: "Jules Herrmann",
      department: [Master Math-Mod],
      organization: [Sorbonne Université],
      email: "jules.herrmann@etu.sorbonne-universite.fr",
    ),
  ),
  abstract : [test],
  bibliography: bibliography("refs.bib"),
)

= Introduction

Le modèle Intergration-et-tir avec fuite, en anglais Leaky Integrate and Fire (LIF), est un modèle classique
en neuroscience. Introduit au début du XXeme siècle par Louis Lapicque, il décrit l'évolution du potentiel electrique
d'un neurone au cours du temps. Ce modèle a pu être étendu et généralisé au cours du temps, notament par l'ajout de 
comportement stochastique. Ainsi, son évolution non linéaire bruité (NNLIF), basée sur
une équation différentielle stochastique, modèlise les comportement de champ moyen d'un ensemble de neurones.

En l'absence de tir, l'évolution du potentiel $V$ est régi par l'équation stochastique suivant :

$
d V = -(V - V_L) d t + mu d t + sigma d B_t
$

Le terme $-(V - V_L) d t$ correspond à la relaxation du potentiel vers un potentiel de repos $V_L$

Le terme $mu d t + sigma d B_t$ correspond à l'influence du courant synaptique généré par le reste des neurones, qui
comporte un terme de _drift_, controlé par $mu$, et un terme stochastique controlé par $sigma$ ($B_t$ est un mouvement Brownien).

En revanche, lorsqu'un neurone atteint le potentiel de tir $V_F$, son potentiel est réinitialisé à la valeur $V_R$.

Il est possible de déduire de ce modèle une équation différentielle de type Fokker Planck décrivant l'évolution de la 
densité de probabilité $p(v,t)$ de la variable aléatoire $V$.
$
cases(
  partial_t p + partial_v (h p) - a partial_(v v)p = 0  & forall v in \]-inf \, V_F \[\\{V_R},
  p(v,0) = p_0(v) & forall v in \]-inf \, V_F \[,
  p(- inf,t) = p(V_F,t) = 0 & forall t > 0,
  p(V_R^-,t) = p(V_R^+,t) & forall t > 0,
  partial_v p(V_R^-,t) = partial_v p(V_R^+,t) + N(t)/a & forall t > 0,
  N(t) = - a partial_v p(V_F,t) & forall t > 0
)
$ <eq:fokkerplanck>

- $p_0$ est la densité de probabilité initialte de $V$
- $h$ est le terme de drift
- $a$ est le coefficient de diffusion
- $N(t)$ est le taux moyen de tir au temps $t$

Le taux de tir $N$ du groupe de neurone a un impact sur le comportement de chaque neurone, ainsi, dans un réseaux excitatoire,
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

Dans leur article de 2011, #cite(<caceres_analysis_2011>,form:"prose"), ont posé les bases de l'étude de ce système. Ils ont caractérisé formelement le phénomène de "blow up", dans leque le taux de tir $N(t)$ diverge en temps fini. Si les blow up avaient été observés numériquement par le passé, ils ont prouvé leur existence dans le cas continue dans une famille de situation initiales. Ils ont également montré l'existence de situations présentant 0, 1 ou 2 solutions stationnaires.

De nombreuses propriété de ce système sont toujours inconnues. La stabilité des solutions stationnaires et les comportements asymptotiques du systèmes sont mal comprises. Dans le cas linéaire, il est possible de montrer que l'entropie relative de la solution est strictement décroissante, ce qui permet de montrer la convergence de la solution vers une solution stationnaire dans certains cas.

Ce manque de résultat théorique complique le développement de méthodes numériques capable d'approximer ce problème.

Une approche générale utilisée pour résoudre numériquement les équations de Fokker-Planck est l'approche dite de Scharfetter-Gummel. Cette approche reformule l'équation, comme une équation de balance des flux.

Dans leur article de 2020, #cite(<hu_structure_2020>,form:"prose") proposent une méthode numérique basée sur l'approche de Scharfetter-Gummel, et l'adapte pour prendre en compte le saut de flux. Cette approche satisfait plusieurs bonnes propriété : elle est linéairement implicite et donc facile à résoudre numériquement, conserve la positivité de la densité et satisfait une version discrete de la propriété de l'entropie relative décroissante qui existe dans le cas continue.

Ce rapport présente premièrement cette méthode numérique. Ainsi que les choix d'implémentations que j'ai effectué afin de pouvoir, dans une seconde partie, répliquer les expériences numériques de l'article.

= Méthode

== Formulation en équilibre de flux

La première ligne de @eq:fokkerplanck peut être reformulée en terme d'un flux.

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

On choisi une valeur $V_("min")$ telle que $p(V_"min",t) approx 0$. On discrétise l'intervalle $\[V_"min" \, V_F \]$ en $n$ intervalles de longueur $h$, délimités par les points $v_0 \, v_1 dots, v_n$. Ce sont les points sur lesquelles $p$ sera discrétisé. Le flux $F$, sera lui évalué sur les points intermédiaires $v_(1/2) \, v_(1+ 1/2), dots , v_(n- 1/2)$.

#figure(
  image("figs/discretisation.pdf", width: 100%),
  caption: [Schéma de discrétisation],
) <discretisation>

On discrétise le temps sur l'intervalle $\[0 \, T \]$ avec un pas de $tau$. On note lors $F_i^m$ et $p_i^m$ la valeur du flux et de la densité au point correspondant dans la discrétisation.

== Différence finie

@eq:equiflux peut être discrétisé avec une différence finie centrée $forall i in [|0,N|]$

$
(p_i^(m+1) - p_i^m)/ tau + (tilde(F)_(i+1/2)^m - tilde(F)_(i-1/2)^m)/h = 0
$

On discrétise également l'expression de $F$ dans @eq:F_reform avec une différence finie centrée. Avec @eq:F_tilde, on obtient, $forall i in [|1,N-2|]$,
$
tilde(F)_(i+1/2)^m = -a(N^m) M_(i+1/2)/h (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i))
- N^m bb(1)_(v_(i+1/2) >= V_R)
$

$
tilde(F)_(1/2) = tilde(F)_(N-1/2) = 0
$

= Résultat
https://www.youtube.com/watch?v=Sw56snxPUMg
https://www.youtube.com/watch?v=jF3e-Q7j1z0

= Conclusion

= Annexe

------
$
forall i in [|2,N-2|] \
p_i^m =  p_i^(m+1) + tau/h tilde(F)_(i+1/2)^m - tau/h tilde(F)_(i-1/2)^m
$

Substituting the $tilde(F)$ with their expression
$
p_i^m =  p_i^(m+1) + tau/h (-a(N^m) M_(i+1/2)/h (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i))
- N^m bb(1)_(v_(i+1/2) >= V_R)) - tau/h (-a(N^m) M_(i-1/2)/h (p^(m+1)_(i)/M^m_(i) - p^(m+1)_(i-1)/M^m_(i-1))
- N^m bb(1)_(v_(i-1/2) >= V_R))
$

Simplifying
$
p_i^m + tau/h N^m (bb(1)_(v_(i+1/2) >= V_R) - bb(1)_(v_(i-1/2) >= V_R))=  p_i^(m+1) - tau/h a(N^m) M_(i+1/2)/h (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i))
 + tau/h a(N^m) M_(i-1/2)/h (p^(m+1)_(i)/M^m_(i) - p^(m+1)_(i-1)/M^m_(i-1))
$

Since $bb(1)_(v_(i+1/2) >= V_R) - bb(1)_(v_(i-1/2) >= V_R) = bb(1)_(v_i = V_R)$, we get 
$
p_i^m + tau/h N^m bb(1)_(v_i = V_R) =
p_i^(m+1) - tau/h^2 a(N^m) M_(i+1/2) (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i))
 + tau/h^2 a(N^m) M_(i-1/2) (p^(m+1)_(i)/M^m_(i) - p^(m+1)_(i-1)/M^m_(i-1))
$
Distributing
$
p_i^m + tau/h N^m bb(1)_(v_i = V_R) =
- tau/h^2 a(N^m) M_(i+1/2)/M^m_(i+1) p^(m+1)_(i+1)
+ (1 + tau/h^2 a(N^m) (M_(i+1/2) + M_(i-1/2))/M^m_(i)) p_i^(m+1)
 - tau/h^2 a(N^m) M_(i-1/2)/M^m_(i-1) p^(m+1)_(i-1)
$

------

$
p_1^m & =  p_1^(m+1) + tau/h tilde(F)_(1+1/2)^m - underbrace(tau/h tilde(F)_(1/2)^m ,=0) \

      & =  p_1^(m+1) - tau/h^2 a(N^m) M_(1+1/2) (p^(m+1)_(2)/M^m_(2) - p^(m+1)_(1)/M^m_(1))
- underbrace(N^m bb(1)_(v_(1+1/2) >= V_R),=0) \

& = - tau/h^2 a(N^m) M_(1+1/2)/M^m_(2) p^(m+1)_(2) + (1 + tau/h^2 a(N^m) M_(1+1/2)/M^m_(1))p^(m+1)_1
$

------

$
p_(N-1)^m & =  p_(N-1)^(m+1) + underbrace(tau/h tilde(F)_(N-1/2)^m,=0) - tau/h tilde(F)_(N-1-1/2)^m \
          & =  p_(N-1)^(m+1) - tau/h tilde(F)_(N-2+1/2)^m \
          & =  p_(N-1)^(m+1) + tau/h^2 a(N^m) M_(N-2+1/2) (p^(m+1)_(N-1)/M^m_(N-1) - p^(m+1)_(N-2)/M^m_(N-2))
+ tau/h N^m underbrace(bb(1)_(v_(N-2+1/2) >= V_R),=1) \
$
$
p_(N-1)^m - tau/h N^m 
 & = (1 + tau/h^2 a(N^m) M_(N-2+1/2)/M^m_(N-1)) p^(m+1)_(N-1) - tau/h^2 a(N^m) M_(N-2+1/2)/M^m_(N-2) p^(m+1)_(N-2)
$

------

Avec $alpha = tau/h^2 a(N^m)$

$
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
)
mat(
  p^(m+1)_1;
  p^(m+1)_2;
  dots.v;
  p^(m+1)_(N-2);
  p^(m+1)_(N-1);
)
=
mat(
  tilde(p)^(m)_1;
  tilde(p)^(m)_2;
  dots.v;
  tilde(p)^(m)_(N-2);
  tilde(p)^(m)_(N-1);
)
$

Où $p^(m+1)$ est obtenu de $tilde(p)^(m+1)$ en soustrayant et en ajoutant $tau/h N^m$ aux posititons $N-1$ et $V_R$ respectivement.
