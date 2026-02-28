#import "@preview/unequivocal-ams:0.1.2": ams-article, theorem, proof

#show: ams-article.with(
  title: [Title],
  authors: (
    (
      name: "Jules Herrmann",
    ),
  ),
  //bibliography: bibliography("refs.bib"),
)

$
forall i in [|0,N|] \
(p_i^(m+1) - p_i^m)/ tau + (tilde(F)_(i+1/2)^m - tilde(F)_(i-1/2)^m)/h = 0
$

$
forall i in [|1,N-2|] \
tilde(F)_(i+1/2)^m = -a(N^m) M_(i+1/2)/h (p^(m+1)_(i+1)/M^m_(i+1) - p^(m+1)_(i)/M^m_(i))
- N^m bb(1)_(v_(i+1/2) >= V_R)
$
$
tilde(F)_(1/2) = 0
$
$
tilde(F)_(N-1/2) = 0
$

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
  tilde(p)^(m+1)_1;
  tilde(p)^(m+1)_2;
  dots.v;
  tilde(p)^(m+1)_(N-2);
  tilde(p)^(m+1)_(N-1);
)
$

Où $p^(m+1)$ est obtenu de $tilde(p)^(m+1)$ en soustrayant et en ajoutant $tau/h N^m$ aux posititons $N-1$ et $V_R$ respectivement.
