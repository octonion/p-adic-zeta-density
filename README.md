# Multiweight arithmetic holonomy and density-one irrationality of p-adic zeta values

This repository contains a research draft on multiweight arithmetic holonomy for Kubota--Leopoldt $p$-adic zeta values. The manuscript develops genus-zero multiweight estimates, an intrinsic arbitrary-prime construction, exact mixed capacitary normalization and target-normal cancellation, an exact modular capacity-pair formula at every prime, explicit density-one bounds, and certified finite-packet consequences at several positive-genus prime levels.

> **Status.** This is a research draft. The manuscript records its internal audit status and notes that independent specialist review remains appropriate before publication.

## Files

- [`multiweight_arithmetic_holonomy_all_primes.tex`](multiweight_arithmetic_holonomy_all_primes.tex) — self-contained LaTeX source for the dependency-audited all-prime manuscript.
- [`multiweight_arithmetic_holonomy_all_primes.pdf`](multiweight_arithmetic_holonomy_all_primes.pdf) — compiled PDF generated from the LaTeX source above.
- [`allprime_capacity_finite_packet_certificate.py`](allprime_capacity_finite_packet_certificate.py) — outward-rounded integer fixed-point certificate for the positive-genus finite-packet bounds derived from the exact all-prime capacity formula.
- [`allprime_capacity_finite_packet_certificate_output.txt`](allprime_capacity_finite_packet_certificate_output.txt) — stored output of the capacity finite-packet certificate, including the certified margins and final `CERTIFICATE: PASS` result.

## Abstract

For each genus-zero prime

$$
\mathcal P_0=\{2,3,5,7,13\}
$$

and each finite set $\mathcal S$ of distinct odd integers at least $3$, we combine the Eisenstein--Eichler extensions attached to the values $\zeta_p(s)$, $s\in\mathcal S$, using one common homogenizing coordinate. Under the joint rationality hypothesis the resulting rank-$m$ horizontal leaf, where $m=1+\sum_{s\in\mathcal S}s$, has full functional rank. The blockwise genuine-source Hasse kernels and a common torsor-valued Cartier argument give the Hasse-improved multiweight type

$$
\tau^{\mathrm H}_{p,\mathcal S}(\xi)
=\frac{2}{m^2}\left[
\xi\sum_{s\in\mathcal S}s^2
-\left(\sum_{s\in\mathcal S}(s-1)\right)J_p(\xi)
+(\max\mathcal S)I^m_\xi(\xi)
\right].
$$

A thickness-sensitive Fitting refinement replaces the rectangular pole bound by a denominator polygon. The certified finite conclusions are unchanged: among the first $1200$ positive odd arguments, at least

$$
1194,\quad1191,\quad1185,\quad1180,\quad1061
$$

of the values $\zeta_p(s)$ are irrational for $p=2,3,5,7,13$, respectively.

We then develop an intrinsic-curve version that no longer assumes $X_0(p)\simeq\mathbf P^1$. On the fixed coarse modular curve, after removing the finite elliptic locus, a finite-map restriction-of-scalars argument gives a Shidlovsky zero estimate and exact $q\mapsto q^\ell$ divided Frobenius gives a general-curve $q$-jet Cartier window. A normalization-robust slopes argument applies finite-order capacitary estimates only after contraction to a scalar section. It yields

$$
\tau^{\mathrm{curve}}_{p,\mathcal S}(\xi)
\ge \Lambda_p^\circ-\frac{C_p}{m},
$$

where $\Lambda_p^\circ>0$ comes from the strict capacitary gain of a $p$-adic annular collar and $C_p$ is independent of the packet.

We also construct the canonical exact metric. Its finite components are the Bost--Chambert-Loir capacitary $\mathbf Q$-model metrics, while its Archimedean components are the Bost--Charles direct-image Green functions of the chosen continuation maps. After clearing the finitely many rational exponents, the capacitary tangent line is the normal line of a pseudoconcave formal-analytic surface. The Chen--Moriwaki determinant Hilbert--Samuel theorem and a uniform finite-order Green--Jensen estimate then give the exact source formula and the target-normal cancellation

$$
\widehat{\deg} M_D=\frac m2\mathscr H_\alpha D^2+o(D^2),
\qquad
\tau^{\mathrm{curve}}_{p,\mathcal S}(\xi)
\ge \Lambda_\alpha-\frac{\mathscr H_\alpha}{m}.
$$

The stronger exact theorem is not needed for the normalization-robust qualitative proof. We then compute the maximal modular capacity pair exactly at every prime. For $0<Y<1/p$,

$$
(\Lambda_p(Y),\mathscr H_p(Y))=
\left(
\frac{12\log p}{p-1}-2\pi Y,
\frac{12\log p}{p-1}-2\pi Y+C_p(Y)
\right).
$$

For $p\ge5$ the finite contribution is the effective resistance of the supersingular reduction graph, evaluated by the Eichler--Deuring mass formula; $p=2,3$ are treated directly with their integral Hauptmodul lemniscates. If $d_p$ is the degree of the finite map used in the general-curve zero estimate, then

$$
R_p(X):=\left\lvert{3\le s\le X:\ s\text{ odd and }\zeta_p(s)\in\mathbf Q}\right\rvert
\le\left(\frac{e\,d_p(p-1)}{12\log p}+o_p(1)\right)(\log X)^2.
$$

Thus the irrational values have natural density one with an explicit leading constant at every fixed prime. A fixed-point finite-packet certificate gives, among the first $1200$ odd arguments, at least

$$
1010,974,963,823,786,775,743,605,594
$$

irrational values for $p=11,17,19,23,29,31,37,41,43$, respectively.

## Reproducing the capacity finite-packet certificate

Run

```bash
python3 allprime_capacity_finite_packet_certificate.py > allprime_capacity_finite_packet_certificate_output.txt
```

The certificate uses outward-rounded integer fixed-point interval arithmetic and makes its certificate decisions without binary floating point. The stored output records the positive worst-case margins for the displayed prime levels and ends with `CERTIFICATE: PASS`.

## Building the manuscript

A standard LaTeX installation with the packages listed in the preamble can compile the manuscript. For example:

```bash
pdflatex multiweight_arithmetic_holonomy_all_primes.tex
pdflatex multiweight_arithmetic_holonomy_all_primes.tex
```

The second pass resolves the document's internal references.
