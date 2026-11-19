---
tags:
  - source-materials/lectures
  - bio/structural-biology
prerequisites: []
course date: Michaelmas 2025
lectuer: Prof. Peter Minary
share: true
---
# Lectures 13–15 Modelling DNA, RNA and protein structures

Readable notes from Minary, *Computational Biology*, Lectures 13–15. Slides: [13-15-modelling-dna-rna-protein.pdf](13-15-modelling-dna-rna-protein.pdf). Concept notes cite the headings below.

## Atomistic and coarse-grained models

The slides date the birth of computational structural biology to 1969: the first energy minimisation of an entire protein with an atomistic energy function, in Levitt and Lifson, *J. Mol. Biol.* 46, 269–279 (1969). That work proposed the functional form of the all-atom energy functions still used. Later improvements have mostly come from fitting the parameters to experiment.

An **atomistic** (all-atom) model represents the molecule with atomic detail. The conformation of an $N$-atom system is
$$
\mathbf{r}=(\mathbf{r}_1,\ldots,\mathbf{r}_N).
$$
A **coarse-grained** model replaces each residue by $k=1,2,3,\ldots$ centres. Proteins are built from amino-acid residues and nucleic acids from nucleotide residues. An $n$-residue chain then has $M=n\cdot k$ centres,
$$
\mathbf{c}=(\mathbf{c}_1,\ldots,\mathbf{c}_M).
$$
The example in the slides uses two centres per residue. The atomistic energy surface is drawn as the more rugged of the two.

![[./attachments/Pasted image 20260925152508.png|Pasted image 20260925152508.png]]

## Bonded and non-bonded energy

The atomistic energy splits into **bonded terms** and **non-bonded terms**. Bond stretch, angle bend, and dihedral torsion are the bonded part. Lennard-Jones and Coulomb interactions, summed over non-bonded pairs $i,j\in\mathrm{nb}$, are the non-bonded part. With $r_{ij}$ the distance between atoms $i$ and $j$,
$$
\begin{aligned}
U(\mathbf{r}_1,\ldots,\mathbf{r}_N)
&=\sum_{\mathrm{bonds}} K_r(r-r_0)^2
+\sum_{\mathrm{bends}} K_\theta(\theta-\theta_0)^2
+\sum_{\mathrm{dihedrals}} K_\phi\bigl[1+\cos(n\phi+\delta)\bigr]\\
&\quad+\sum_{i,j\in\mathrm{nb}}
\left(
4\epsilon_{ij}
\left[
\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12}
-\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}
\right]
+K_C\frac{q_i q_j}{r_{ij}}
\right).
\end{aligned}
$$

> ![[./attachments/slide-06.jpg|slide-06]]
> Small cartoons of bond stretch, angle bend, dihedral torsion, and a nonbonded pair, under the atomistic energy function.

Coarse-grained energies can keep the same pattern. The bonds are then pseudobonds between centres, and the nonbonded term $V_{\mathrm{nb}}$ may be physical, knowledge-based, or both. If $\ell$ and $\phi$ are the length and torsion of a pseudobond, $\theta$ is the bend angle at a centre, and $c_{ij}$ is the distance between centres $i$ and $j$,
$$
U(\mathbf{c})
=\sum_{\mathrm{bonds}} c_\ell(\ell-\ell_0)^2
+\sum_{\mathrm{bends}} c_\theta(\theta-\theta_0)^2
+\sum_{\mathrm{dihedrals}} c_\phi V(\phi)
+\sum_{i,j\in\mathrm{nb}} V_{\mathrm{nb}}(c_{ij}).
$$

## Metropolis Method
### Where the Metropolis method comes from

In a letter of 21 May 1947 from John von Neumann (Institute for Advanced Study) to Stanislav Ulam (Los Alamos), von Neumann described **acceptance-rejection sampling**. For a density $f$ on an interval $[a,b]$, draw $x_i$ uniformly from $[a,b]$ and $y_i$ uniformly from $[0,\max f]$. Keep $x_i$ when $y_i<f(x_i)$ and reject it when $y_i>f(x_i)$.

The same slides attribute the first pseudo-random generator to von Neumann’s **middle-square method**,
$$
R_{n+1}=\mathrm{mid}(R_n^2,m),
$$
where $m$ is the number of digits taken from the middle of $R_n^2$. The worked case is $m=3$ and $R_0=157$: $R_0^2=24649$, so $R_1=464$.

The paper that introduced Metropolis Monte Carlo studied $N$ interacting particles in two dimensions. Positions and momenta are $\mathbf{r}=(\mathbf{r}_1,\ldots,\mathbf{r}_N)$ and $\mathbf{p}=(\mathbf{p}_1,\ldots,\mathbf{p}_N)$, with $\mathbf{r}_i=(x_i,y_i)$. The particles sit in a square. To limit surface effects the whole substance is a periodic tiling of that square, every copy holding the same configuration. The pair distance $d_{ij}$ is the minimum distance from particle $i$ to any periodic image of particle $j$.

> ![[./attachments/slide-12.jpg|slide-12]]
> A square packed with $N$ particles, the repeating cell of the periodic system used to introduce Metropolis Monte Carlo.

The energy is a pair potential,
$$
E(\mathbf{r})=\frac12\sum_{i=1}^{N}\sum_{j\neq i} V(d_{ij}).
$$
Equilibrium averages use the canonical ensemble. For a quantity $F(\mathbf{r})$,
$$
\bar F=\frac{1}{Q}\int F(\mathbf{r})\,e^{-E(\mathbf{r})/kT}\,\mathrm{d}\mathbf{r}\,\mathrm{d}\mathbf{p},
\qquad
Q=\int e^{-E(\mathbf{r})/kT}\,\mathrm{d}\mathbf{r}\,\mathrm{d}\mathbf{p}.
$$
These are $4N$-dimensional integrals. Because the forces and $F(\mathbf{r})$ do not depend on velocity, the momentum integrals separate, and the average reduces to a $2N$-dimensional integral over configurations:
$$
\bar F
=\frac{\displaystyle\int F(\mathbf{r})\,e^{-E(\mathbf{r})/kT}\,\mathrm{d}\mathbf{r}}
{\displaystyle\int e^{-E(\mathbf{r})/kT}\,\mathrm{d}\mathbf{r}}.
$$
A naive Monte Carlo draw places every particle uniformly in the square. In a close-packed system most of those draws have a huge energy and a negligible Boltzmann weight. Metropolis Monte Carlo was the proposed remedy.

The journal page reproduced in the lecture states the idea: propose configurations with probability proportional to $e^{-E/kT}$ and then weight the retained samples evenly, rather than drawing configurations uniformly and reweighting afterwards.

> ![[./attachments/slide-20.jpg|slide-20]]
> The 1953 *Journal of Chemical Physics* page in which Metropolis, the Rosenbluths, and the Tellers state the modified Monte Carlo scheme.

> ![[./attachments/slide-21.jpg|slide-21]]
> Metropolis Monte Carlo walking through states of non-zero weight, contrasted with ordinary Monte Carlo scattered over the whole state space.

A later slide reproduces a *SIAM News* clipping that places this algorithm on an editors’ list of top algorithms of the twentieth century. The caption says the method, and the related Metropolis–Hastings algorithm, spread through operations research, economics, genomics, bioinformatics, and across physics, chemistry, biology, computer science, applied mathematics, neuroscience, sociology, and political science. It is also the basis of simulated annealing and other stochastic optimisers.

> ![[./attachments/slide-22.jpg|slide-22]]
> A *SIAM News* clipping listing the Metropolis algorithm among top algorithms of the twentieth century.

### Metropolis Monte Carlo

One iteration saves the current configuration, then tries a displacement of each particle. The maximum step is $\alpha$, and $\xi_1,\xi_2$ are uniform on $(-1,1)$. A trial is kept with probability
$$
P_{\mathrm{acc}}=\min\bigl\{1,\exp(-\Delta E/kT)\bigr\},
$$
where $\Delta E=E'-E$ is the energy change. Draw $\zeta\sim U(0,1)$. If $P_{\mathrm{acc}}>\zeta$, keep the new energy; otherwise restore the saved configuration.

```text
Metropolis Monte Carlo
  E ← E(r)
  for j = 1, …, M do
    r_s ← r
    for i = 1, …, N do
      x_i ← x_i + α ξ_1
      y_i ← y_i + α ξ_2
      E' ← E(r),  ΔE ← E' − E
      P_acc ← min{1, exp(−ΔE / kT)}
      if P_acc > ζ ∼ U(0,1) then E ← E'
      else r ← r_s
```

After $M$ iterations the equilibrium estimate of a quantity of interest is the plain average of its values on the visited configurations,
$$
\bar F=\frac1M\sum_{j=1}^{M} F_j.
$$

### Markov chain Monte Carlo

The same construction estimates a general weighted integral over a state space $\Omega$. For functions $\alpha$ and $\phi$,
$$
\bar\alpha
=\frac{\displaystyle\int_\Omega \alpha(x)\,\phi(x)\,\mathrm{d}x}
{\displaystyle\int_\Omega \phi(x)\,\mathrm{d}x}.
$$
The Boltzmann case is $x\equiv\mathbf{r}$, $\Omega$ the configuration space, and $\phi(x)=e^{-E(x)/kT}$.

**Markov chain Monte Carlo** samples $\Omega\subset\mathbb{R}^n$ against $f(x)=\phi(x)/N$, where $f$ is known only up to the normalising constant $N$. Proposals come from a symmetric density $q(x,y)$ in the second argument, so $\int q(x,y)\,\mathrm{d}y=1$ and $q(x,y)=q(y,x)$.

```text
MCMC
  x ← x_0
  for i = 1, …, M do
    x' ∼ q(x, x')
    P_acc(x ← x') ← min{1, φ(x')/φ(x)}
    if P_acc > ζ ∼ U(0,1) then x ← x'
    X_i ← x
```

In the physical applications the target is the Boltzmann distribution,
$$
f(x)=\frac{1}{N}\exp\bigl(-\beta E(x)\bigr),
\qquad
\beta=\frac{1}{kT},
\qquad
k=1.38\times 10^{-23}\,\mathrm{J/K},
$$
with $T$ the temperature.

### Hamiltonian Monte Carlo

Hamiltonian Monte Carlo still accepts or rejects a proposal, but the proposal is a short molecular-dynamics trajectory rather than a direct draw from $q$. At the start of an iteration, velocities are drawn from a Gaussian,
$$
\mathbf{v}\sim\mathcal{N}(\boldsymbol{\mu}_v,\Sigma_v),
\qquad
\boldsymbol{\mu}_v=\mathbf{0},
\qquad
\Sigma_v=\frac{kT}{m}I,
$$
with every mass $m_i=m$, and the momenta are $p_i=m_i v_i$ in three dimensions. The chain then takes $L$ molecular-dynamics steps. The slides write a half-step of the momenta with the current forces, a full step of the positions, a force update $\mathbf{F}\leftarrow\mathbf{F}(\mathbf{r})$, and a second half-step of the momenta. The proposed point is accepted using the Hamiltonian $H(\mathbf{r},\mathbf{p})$:
$$
P_{\mathrm{acc}}=\min\bigl\{1,\exp\bigl(-(H'-H)/kT\bigr)\bigr\}.
$$
On rejection the positions return to the saved configuration.

```text
Hamiltonian Monte Carlo
  H ← H(r, p)
  for j = 1, …, M do
    draw velocities v ∼ N(0, (kT/m) I) and set p_i ← m v_i
    r_s ← r
    for k = 1, …, L do
      for i = 1, …, N do
        p_i ← p_i + (Δt/2) F_i
        r_i ← r_i + Δt p_i / m_i
      F ← F(r)
      for i = 1, …, N do
        p_i ← p_i + (Δt/2) F_i
    H' ← H(r, p)
    P_acc ← min{1, exp(−(H' − H)/kT)}
    if P_acc > ζ ∼ U(0,1) then H ← H'
    else r ← r_s
```

### Energy barriers and rare events

A chain that only proposes local moves can fail to cross a high barrier in any practical number of steps. It then misses part of the region where $\phi$ is appreciable, and the Monte Carlo average of $\alpha$ is not the integral $\bar\alpha$. The slides draw that missed region in grey and assume $\phi(x)=\exp(-\beta E(x))$.

> ![[./attachments/slide-27.jpg|330]]
> A Markov chain trapped on one side of a high barrier, missing the grey part of the state space that still carries Boltzmann weight.

For sampling at the Boltzmann distribution $f(x)\propto\exp(-E(x)/kT)$, the rate of crossing a barrier is proportional to $\exp(-\Delta E/kT)$, where $\Delta E$ is the barrier height. On the one-dimensional examples, $\Delta E=E(0)-E(1)$ for each curve.

> ![[./attachments/slide-28.jpg|slide-28]]
> One-dimensional double wells of several barrier heights, next to a rugged surface in the plane.

## Advanced Sampling Methods
### Parallel tempering

Parallel tempering runs $K+1$ Markov chains at increasing temperatures,
$$
X^{(i)},\quad 0\le i\le K,
\qquad
f_i(x)\propto\exp\bigl(-h(x)/T_i\bigr),
\qquad
T_0<T_1<\cdots<T_K,
$$
with $h(x)=E(x)/k$ and $f(x)=N^{-1}\exp(-h(x)/T)$. The cold chain $X^{(0)}$ is usually at the temperature of interest. Hotter chains move more freely. *Adjacent chains may swap configurations, which feeds that exploration back to the cold chain*. The illustration is drawn for $x\in\mathbb{R}$; the same exchange extends to $\mathbb{R}^n$.

The chains are started from the hot end downward. *At each step one either proposes an ordinary move or an exchange of neighbours*, with
$$
p_{\mathrm{move}}+p_{\mathrm{exchange}}=1.
$$
The acceptance probabilities are
$$
P_{\mathrm{acc-exchange}}
=\min\left\{1,
\frac{f_i\bigl(X^{(i+1)}_n\bigr)\,f_{i+1}\bigl(X^{(i)}_n\bigr)}
{f_i\bigl(X^{(i)}_n\bigr)\,f_{i+1}\bigl(X^{(i+1)}_n\bigr)}
\right\},
$$
$$
P_{\mathrm{acc-move}}
=\min\left\{1,
\frac{f_i\bigl(X^{(i)}_{n+1}\bigr)}{f_i\bigl(X^{(i)}_n\bigr)}
\right\}.
$$

> ![[./attachments/slide-30.jpg|435]]
> A ladder of chains from $X^{(0)}$ to $X^{(K)}$, with a horizontal move inside one chain and a vertical exchange between neighbours.

> ![[./attachments/slide-31.jpg|slide-31]]
> The same ladder beside one-dimensional energy curves sampled at the different temperatures.

### Equi-energy Monte Carlo

Parellel tempering requires careful designs of the temperature ladders. Smaller temperature gaps between the chains gives better exchange probability, but is more computationally expensive. Bigger temperature gaps on the other hand could result in lower exchange rate due to less overlap between adjacent chains. 

Equi-energy Monte Carlo uses the same temperature ladder and the same family $f_i$. It adds energy levels
$$
H_0<H_1<\cdots<H_K<H_{K+1}=\infty,
\qquad
H_0<\inf_x h(x).
$$
An exchange is allowed only between adjacent chains and only when both states sit in the same energy window $[H_j,H_{j+1})$.

> ![[./attachments/slide-34.jpg|335]]
> Two states on a one-dimensional potential that occupy the same energy window and are therefore eligible to swap.

The hot chain $X^{(K)}$ is built first. Its visits are stored in energy disks
$$
D^{(K)}_j=\{x\in X^{(K)}: h(x)\in[H_j,H_{j+1})\},
\quad 0\le j\le K.
$$
While building $X^{(K-1)}$, with probability $p_{\mathrm{exch}}$ one tries to replace the current state by a state $y$ drawn from the matching disk of $X^{(K)}$; with probability $1-p_{\mathrm{exch}}$ one tries an ordinary move. The accept/reject rules are the usual Metropolis ratios of $f_{K-1}$ and $f_K$. Colder chains $X^{(K-2)},\ldots,X^{(0)}$ are generated in the same way, each using the disks of the chain above it. Building $X^{(i)}$ also fills that chain’s own disks.

On a one-dimensional periodic Fourier potential $h(x)$ (the slides set $k=1$) whose target $f$ is known at $T=1$, both equi-energy Monte Carlo and parallel tempering are compared with the analytical distribution. The histograms in the figure use statistics up to $40{,}000$ iterations, and a further panel shows convergence against the number of iterations.

> ![[./attachments/slide-37.jpg|slide-37]]
> A rugged one-dimensional potential, parallel-tempering and equi-energy histograms against the analytical density, and their convergence over Monte Carlo steps.

The same comparison is run on a model $\beta$-protein of $46$ pseudoatoms of three types: hydrophobic (B), hydrophilic (L), and neutral (N). Runs start $3$, $6$, and $9\,\text{Å}$ (RMSD) from the known native conformation. A vertical line on the probability plot marks structures within $1.5\,\text{Å}$ of native: about $0.25$ for parallel tempering and about $0.7$ for equi-energy Monte Carlo.



> ![[./attachments/slide-38.jpg|slide-38]]
> RMSD traces and native-recovery probabilities for parallel tempering and equi-energy Monte Carlo on the $46$-pseudoatom $\beta$-protein, with the model energy written beside the plots.

## Bypassing the energy barrier
### Reference-potential spatial warping

The optional REPSWA slides (reference potential spatial warping) change a torsion angle $\phi$ to a new coordinate $\phi_u$ so that a barrier is easier for Markov chain Monte Carlo to cross, and so that the map can be inverted. For a bond $a_1-a_2-a_3-a_4$, resolve the position of $a_4$ into spherical coordinates $(r,\theta,\phi)$ in a frame with $a_3$ at the origin, then
$$
\phi_u = c\int_0^{\phi}\exp\bigl(-\beta V_{\mathrm{ref}}(\phi')\bigr)\,\mathrm{d}\phi'.
$$
Rebuild the Cartesian position from $(r,\theta,\phi_u)$ and transform back. Repeated application covers the torsions of a chain. The numerical panels are a one-dimensional double well with a $10\,kT$ barrier, and end-to-end and central-dihedral fluctuations of a $400$-mer, each compared for plain MCMC and for MCMC with the torsional transformation.

> ![[./attachments/slide-39.jpg|slide-39]]
> The REPSWA sketch of the reference-potential warping.

> ![[./attachments/slide-41.jpg|slide-41]]
> Double-well and $400$-mer comparisons of MCMC in the original coordinate and MCMC after torsional potential smoothing.

### Extra dimensions, and three responses to a barrier

A related construction adds a coordinate without changing the marginal distribution of the original one. For a particle in a potential $V_x(x)$,
$$
P(x)=\frac{e^{-\beta V_x(x)}}{Q_x},
\qquad
Q_x=\int e^{-\beta V_x(x)}\,\mathrm{d}x.
$$
If $V(x,y)=V_x(x)+V_y(y)$, the joint partition function factors as $Q_x Q_y$, and integrating out $y$ returns the same $P(x)$. **The extra coordinate is there so that a path in the $(x,y)$ plane can get past a barrier that blocks the $x$ line**.

The slides collect three responses to a barrier that separates energy wells: overcome it, remove it, or go around it.

> ![[./attachments/slide-44.jpg|353]]
> Three sketches of a barrier: crossing over it, flattening it, and leaving the line to pass around it.

## Stochastic methods
### Simulated annealing

Stochastic optimisation looks for $x_{\min}\in\Omega\subset\mathbb{R}^n$ with $E(x_{\min})\le E(x)$ for every $x\in\Omega$. Simulated annealing is MCMC with a temperature that falls during the run. After each accepted or rejected step,
$$
T\leftarrow T\, r_a,\qquad r_a\in(0,1),
$$
and the loop stops at or below $T_{\mathrm{final}}$. The acceptance weight is $\phi(x)=\exp(-E(x)/kT)$ at the current temperature. For conformational optimisation of biomolecules the slides quote $T_{\mathrm{initial}}=500\,\mathrm{K}$, $T_{\mathrm{final}}=1\,\mathrm{K}$, and $r_a=0.995$.

```text
Simulated annealing
  x ← x_0
  T ← T_initial
  for i = 1, …, M do
    x' ∼ q(x, x')
    P_acc ← min{1, φ(x')/φ(x)}    with φ(x) = exp(−E(x)/kT)
    if P_acc > ζ ∼ U(0,1) then x ← x'
    X_i ← x
    T ← T r_a
    if T ≤ T_final then break
```

### Multi-canonical simulated annealing

Multi-canonical simulated annealing, as these slides define it, anneals every temperature in a parallel-tempering ladder at once. *Each chain is cooled independently*, and the figure uses one common annealing rate.

> ![[./attachments/slide-48.jpg|slide-48]]
> A parallel-tempering ladder in which every chain’s temperature is annealed, with a stop once a chain reaches the final temperature.

### Monte Carlo minimisation

Monte Carlo minimisation on $E$ is ordinary MCMC on a transformed energy $\tilde E$. The energy of each conformation is replaced by the energy of the nearest local minimum, so $\tilde E$ is flat across the basin of that minimum.

> ![[./attachments/slide-49.jpg|400]]
> A rugged energy $E$ (solid) and the staircase $\tilde E$ (dashed) that takes the energy of the nearest local minimum.
> The loci of the local minima may be lost depending how the new landscape is created

### Stochastic tunnelling

Stochastic tunnelling is MCMC on a transformation that keeps the order of the minima. If $E_0$ is the lowest energy found so far and $\gamma$ is the tunnelling parameter,
$$
\tilde E(x)=1-\exp\bigl[-\gamma\bigl(E(x)-E_0\bigr)\bigr].
$$
The map sends the interval $[E_0,E_{\max}]$ to $[0,1]$. Updating $E_0$ as deeper minima appear removes local minima that are no longer relevant. Wells that sit above the best minimum found so far are suppressed; wells that go deeper are kept and enhanced.

> ![[./attachments/slide-51.jpg|366]]
> A one-dimensional potential, then the tunnelled surface after the current best minimum, then the surface after a deeper minimum is found.

## Internal coordinates and torsion updates

Internal coordinates of a biomolecule are its bond lengths, bond angles, and torsion angles. A molecule of $N$ atoms has $3N-6$ of them.

To write atom $4$ of a chain $1-2-3-4$ in internal coordinates, build a local frame on $\mathbf{r}_1,\mathbf{r}_2,\mathbf{r}_3$:
$$
\begin{aligned}
\mathbf{e}_1&=(\mathbf{r}_3-\mathbf{r}_2)/|\mathbf{r}_3-\mathbf{r}_2|,\\
\mathbf{e}_2&=(\mathbf{r}_1-\mathbf{r}_2)/|\mathbf{r}_1-\mathbf{r}_2|,\\
\mathbf{e}_z&=\mathbf{e}_1,\\
\mathbf{e}_y&=(\mathbf{e}_1\times\mathbf{e}_2)/|\mathbf{e}_1\times\mathbf{e}_2|,\\
\mathbf{e}_x&=\mathbf{e}_y\times\mathbf{e}_z.
\end{aligned}
$$
Here $r=|\mathbf{r}_4-\mathbf{r}_3|$ is the length of bond $3-4$, $\theta$ is $\pi$ minus the bond angle at atom $3$, and $\phi$ is the torsion about bond $2-3$. Changing frame is a rotation: if $\mathbf{r}=x\mathbf{i}+y\mathbf{j}+z\mathbf{k}$ in the old orthonormal frame, the components along $(\mathbf{e}_x,\mathbf{e}_y,\mathbf{e}_z)$ are the dot products of $\mathbf{r}$ with those vectors, packed as a rotation matrix $R$.

> ![[./attachments/slide-56.jpg|slide-56]]
> The local frame on atoms $1,2,3$ in which atom $4$ is expressed by bond length, polar angle, and torsion.

The map from Cartesian coordinates to $(r,\theta,\phi)$ is invertible, in four conceptual steps and their inverses.

1. **Translate** so that atom $3$ is the origin: $\mathbf{r}_4'=\mathbf{r}_4-\mathbf{r}_3$.
2. **Rotate** into the local frame: $\mathbf{r}_4''=R\mathbf{r}_4'$.
3. **Resolve** $\mathbf{r}_4''=(r_{4,x},r_{4,y},r_{4,z})$ into spherical coordinates
$$
r=\sqrt{r_{4,x}^2+r_{4,y}^2+r_{4,z}^2},
\quad
\theta=\arccos(r_{4,z}/r),
\quad
\phi=\arccos\left(\frac{r_{4,x}}{\sqrt{r_{4,x}^2+r_{4,y}^2}}\right).
$$
If $r_{4,y}<0$, replace $\phi$ by $2\pi-\phi$.
4. **Resolve back** when returning to Cartesian coordinates:
$$
r_{4,x}=r\sin\theta\cos\phi,\qquad
r_{4,y}=r\sin\theta\sin\phi,\qquad
r_{4,z}=r\cos\theta,
$$
then apply $R^\mathsf{T}$ and translate by $\mathbf{r}_3$.

A torsion move uses the same pipeline, adds $\delta\phi$ to $\phi$, and inverts. The torsion is kept inside $[0,2\pi)$.

For a longer chain the same steps apply from the free end inward, then the inverse steps rebuild the Cartesian coordinates outward. They can change several torsions, and they show that the transformation is invertible. At atomistic resolution the atoms are moved in groups so that the stereochemistry stays correct.

## Reverse Monte Carlo

Reverse Monte Carlo fits a configuration to an experimental observable. For a monoatomic liquid the slides use the pair correlation $g(r)$. Let $\langle\alpha(r)\rangle$ be the expected number of atoms at a distance between $r$ and $r+\mathrm{d}r$ from a given atom; that count can be read off one configuration. Let $g_I(r)$ be the same count at density $\rho$ if the atoms did not interact,
$$
g(r)=\frac{\langle\alpha(r)\rangle}{g_I(r)},
\qquad
g_I(r)=\rho\,4\pi r^2\,\mathrm{d}r.
$$

> ![[./attachments/slide-62.jpg|slide-62]]
> A shell of radius $r$ and thickness $\mathrm{d}r$ around a central atom, the geometric picture behind $g(r)$.

A simulated configuration $x\in\mathbb{R}^n$ produces an estimate $(g_s(x))(r)$. In Cartesian coordinates $n=3N$. The experimental curve $g_e(r)$ and the simulated curve are compared by a $\chi^2$ statistic. With $n_r$ distance points and experimental error $\sigma_e$,
$$
\chi^2(x)
=\sum_{i=1}^{n_r}
\frac{\bigl[g_e(r_i)-(g_s(x))(r_i)\bigr]^2}{\sigma_e^2(r_i)}.
$$
A proposed move $x\to x'$ changes the fit by $\Delta\chi^2=\chi^2(x')-\chi^2(x)$.

> ![[./attachments/slide-63.jpg|slide-63]]
> Experimental and simulated pair-correlation curves plotted against distance in ångströms.

The reverse Monte Carlo chain targets $\phi(x)=\exp(-\chi^2(x))$ with the acceptance rule written on the slide,
$$
P_{\mathrm{acc}}(x\leftarrow x')
=\min\left\{1,\,\frac12\,\frac{\phi(x')}{\phi(x)}\right\}.
$$
Proposals still come from a symmetric $q$, and the visited states $X_i$ are the chain whose configurations are being driven toward $g_e$.

```text
Reverse Monte Carlo
  x ← x_0
  for i = 1, …, M do
    x' ∼ q(x, x')
    P_acc ← min{1, (1/2) φ(x')/φ(x)}    with φ(x) = exp(−χ²(x))
    if P_acc > ζ ∼ U(0,1) then x ← x'
    X_i ← x
```

# References

- Minary, Computational Biology, Lectures 13–15 (Michaelmas 2025)
- Levitt, M. and Lifson, S. *J. Mol. Biol.* 46, 269–279 (1969)
- Energy-term figure: *BMC Biology* 9, 71 (2011)
- Landscape figure: *Chem. Rev.* 116, 7898–7936 (2016)
- von Neumann to Ulam, 21 May 1947, on acceptance-rejection sampling
- Metropolis, Rosenbluth, Rosenbluth, Teller, and Teller, *J. Chem. Phys.* 21 (1953), as reproduced on the slides
- *The Annals of Statistics* 34, 1581–1619 (2006), https://doi.org/10.1214/009053606000000515
- *The Annals of Statistics* 34, 1636–1641 (2006), https://doi.org/10.1214/009053606000000470
- *SIAM J. Sci. Comput.* 30, 2055–2083 (2008), https://doi.org/10.1137/070686706
- *Physical Review Letters* 82, 3003 (1999), https://doi.org/10.1103/PhysRevLett.82.3003
- *Molecular Simulation* 1, 359–367 (1988), https://doi.org/10.1080/08927028808080958
