// #import "@preview/diatypst:0.9.1": *
// 
// #show: slides.with(
//   title: "Une méthode de différences finies pour la résolution numérique de l'équation de Fokker-Planck du modèle LIF", // Required
//   authors: ("Jules Herrmann"),
// 
//   // Optional (for more see docs at https://mdwm.org/diatypst/)
//   ratio: 16/9,
//   layout: "medium",
//   title-color: blue.darken(60%),
//   toc: true,
// )
#import "@preview/touying:0.6.3": *
#import themes.simple: *

#show: simple-theme.with(aspect-ratio: "16-9")
#set text(lang: "fr")

= Une méthode de différences finies pour la résolution numérique de l'équation de Fokker-Planck du modèle LIF

Jules Herrmann

UE de modélisation stochastique pour les neurosciences

== Nonlinear Noisy Leaky Integrate and Fire

#align(horizon,[

$
d V = -(V - V_L) d t + mu d t + sigma d B_t
$

#v(20pt)

Quand $V$ atteint $V_F$, il y a un saut à $V_R$
])

==
#scale(90%,reflow:true, [
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
])

==

$
h(v,N(t)) = -v + b N(t)
$
$
a(N(t)) = a_0 + a_1 N(t)
$

#v(20pt)

- $b < 0$ : réseau inhibiteur
- $b > 0$ : réseau excitateur.

==
#cite(<caceres_analysis_2011>,form:"full")

#align(left, [
- Existence des _blow up_
- Solutions stationnaires
- Décroissance de l'entropie relative
])

==
#cite(<hu_structure_2020>,form:"full")

#align(left, [
- Méthode numérique
- Décroissance de l'entropie relative dans le cas discret
])

== Reformulation en terme de flux

#align(horizon,[

$
partial_t p + partial_v (h p) - a partial_(v v)p = 0
$

#pause

$
arrow.b
$

$ partial_t p + partial_v F = 0 $ <eq:equiflux>

$ F(v,t) = -a partial_v p + h p $
])

== Reformulation de Scharfetter-Gummel

#align(horizon,[
$
M(v,t) = exp(- ((v - b N(t))^2)/(2 a(N(t))))
$

$
F(v,t) = a(N(t)) M(v,t) partial_v (p(v,t)/M(v,t))
$ 
])

==

Réinitialisation en $V_R$ lors d'un tir -> Le flux n'est pas continu

$
F(V_R^+,t) - F(V_R^-,t) = N(t)
$

#pause
#v(20pt)

$ tilde(F)(v,t) = F(v,t) - N(t) bb(1)_(v >= V_R) $ <eq:F_tilde>

$tilde(F)$ satisfait l'équation d'équilibre des flux et est continu

== Schéma de discrétisation

#align(horizon,[
#figure(
  image("figs/discretisation.pdf", width: 100%),
) <discretisation>
])

- $h$  : pas de discrétisation en espace
- $tau$ : pas de discrétisation en temps

== Différence finie

#align(horizon,[
$ partial_t p + partial_v F = 0 $ <eq:equiflux>

$
arrow.b
$

$
(p_i^(m+1) - p_i^m)/ tau + (tilde(F)_(i+1/2)^m - tilde(F)_(i-1/2)^m)/h = 0
$
])

== Différence finie

#align(horizon,[
$ tilde(F)(v,t) = F(v,t) - N(t) bb(1)_(v >= V_R) $ <eq:F_tilde>

$
F(v,t) = a(N(t)) M(v,t) partial_v (p(v,t)/M(v,t))
$ 

$
arrow.b
$

$
tilde(F)_(i+1/2)^m = -a(N^m) M_(i+1/2)/h (p^(#text(red,$m+1$))_(i+1)/M^m_(i+1) - p^(#text(red,$m+1$))_(i)/M^m_(i))
- N^m bb(1)_(v_(i+1/2) >= V_R)
$<eq:F_discretisation>
])

==

Après quelques manipulations

#align(center + horizon,[
1. $N^m <- a_0 * p_(n-1)^m / (h - a_1 * p_(n-1)^m)$
2. $p_(n-1)^m <- p_(n-1)^m - tau/h N^m$
3. $p_(V_R)^m <- p_(V_R)^m + tau/h N^m$
4. Résoudre $A_m p^(m+1) = p^m$
])

==

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


#figure(
  image("figs/fig1_bottom.pdf", width: 80%),
  caption: [$a(N(t)) = 1$, $b=1.5$, $p_0$ avec $v_0 = 1.5$ et $sigma^2_0 = 5 times 10^(-3)$. (#link("https://www.youtube.com/watch?v=jF3e-Q7j1z0")[vidéo])],
)<fig:blowup_bottom>

#figure(
  image("figs/fig1_top.pdf", width: 90%),
  caption: [$a(N(t)) = 1$, $b=3$, $p_0$ avec $v_0 = -1$ et $sigma^2_0 = 0.5$.
  (#link("https://www.youtube.com/watch?v=Sw56snxPUMg")[vidéo])],
)<fig:blowup_top>

== Distribution Stationnaire

$
cases(
  partial_v (h(v,N^oo) p^oo) - a(N^oo) partial_(v v)p^oo = 0  & forall v in \]-oo \, V_F \[\\{V_R},
  p^oo (- oo) = p^oo (V_F) = 0,
  p^oo (V_R^-) = p^oo (V_R^+),
  partial_v p^oo (V_R^-) = partial_v p^oo (V_R^+) + N^oo /a(N^oo),
  N^oo = - a(N^oo) partial_v p(V_F)
)
$ <eq:stationary>

#pause

- $b <= 0$ : Existence et unicité de $p^oo$
- $b > 0 $ : 0,1 ou 2 distribution stationnaires, ou plus (problème ouvert)

==

Si une distribution stationnaire existe, il existe $N^oo$ tel que

$
p^oo = N^oo /a(N^oo) exp(- h(v,N^oo )^2 /(2 a(N^oo))) integral_"max"(v,V_R)^V_F exp( h(omega ,N^oo )^2 /(2 a(N^oo))) d omega
$

== Entropie relative

$
S(t) = sum_(i=1)^(n-1) h G(p_i^t/p_i^oo ) p_i^oo
$

- $G(x) = 1/2 (x-1)^2$
- $p^oo$ une distribution stationnaire

#pause

Cas $a = 1$ et $b=0$ :

$ d/(d t) S(t) <= 0 $

#figure(
  image("figs/fig4.pdf", width: 100%),
  caption: [$a(N(t)) = 1$, $b=0$.\
$p_0$ avec $v_0 = 0$ et $sigma^2_0 = 0.25$.\
$p^oo$ avec $N^oo = 0.1377$
]
)<fig:entro_lineaire>

#figure(
  image("figs/fig5.pdf", width: 100%),
  caption: [
    $a(N(t)) = 1$, $b=1.5$.\
    $p_0$ avec $v_0 = 0$ et $sigma^2_0 = 0.25$.\
    $p^oo_"stable"$ et $p^oo_"instable"$, avec respectivement $N^oo = 0.1924$ et $N^oo = 2.319$.
  ]
)

==
#bibliography("refs.bib",title : none)
