---
tags:
  - source-materials/lectures
  - bio/sequencing
prerequisites:
  - "[[NP-complete|NP-complete]]"
course date: Michaelmas 2025
lectuer: Prof. Peter Minary
share: true
---
# Lecture 2 Restriction Mapping

Readable notes from Minary, *Computational Biology*, Lecture 2. Slides: [02-restriction-mapping.pdf](02-restriction-mapping.pdf). Concept notes cite the headings below.

## 1 Restriction enzymes and restriction mapping

Restriction enzymes evolved in bacteria to protect against invading DNA. They restrict that DNA by cleaving it. A second class of enzymes methylates host DNA so the host is not cut; restriction enzymes cleave unmodified DNA.

A restriction enzyme cuts DNA at a short specific pattern, a **restriction site**. Most sites are palindromes in the molecular sense: the word equals its reverse complement, so both strands show the same sequence when each is read $5'\to 3'$.

A **restriction map** is the map of those sites along a DNA sequence. **Restriction mapping** infers the positions of the sites from experimental sets of fragment lengths.

Some enzymes leave blunt ends. HaeIII cuts `5' …GG↓CC… 3'` opposite `3' …CC↑GG… 5'`. AluI cuts `5' …AG↓CT… 3'` opposite `3' …TC↑GA… 5'`. Others leave overhangs. EcoRI cuts `5' …G↓AATTC… 3'` opposite `3' …CTTAA↑G… 5'`.
nh
Naming: the first three letters name the organism, an optional fourth letter names the strain, and Roman numerals index distinct enzymes from the same organism.

## 2 Gel electrophoresis

DNA length is measured in base pairs. Gel electrophoresis measures that length. The gel sits in an electric field. DNA is negatively charged and migrates toward the positive pole. In a fixed time, migration distance is a function of length.

A lane of DNA of known lengths calibrates the gel. Unknown fragments in other lanes are then read off that calibration.

![[./attachments/Pasted image 20260924200550.png|189]]

A practical model is
$$
D \approx a + b\log L,\qquad b<0,
$$
with intercept $a$ and slope $b$. $D$ is migration distance and $L$ is length. Smaller DNA threads through the gel more easily. Long DNA tangles in the matrix and does not migrate far.

## 3 Complete digest and partial digest

A **complete digest** cuts at every restriction site, so the fragments lie between consecutive sites. If the sites are at $0,2,4,7,10$, the fragment lengths are $2,2,3,3$.

![[./attachments/Pasted image 20260924200633.png|522]]

A **partial digest** does not cut at every site. Fragments appear between any pair of sites. The same positions give lengths $2,2,3,3,4,5,6,7,8,10$.

![[./attachments/Pasted image 20260924200743.png|526]]

## 4 Multiset notation

A **multiset** allows repeated elements, for example $\{2,2,3\}$. $C(D,x)$ is the number of times $x$ occurs in multiset $D$.

The difference of multisets $D$ and $L$ keeps the surplus counts:
$$
D\setminus L=\{x\mid x\in D \text{ and } C(D,x)-C(L,x)>0\}.
$$
The sum adds counts:
$$
D+L=\{x\mid C(D+L,x)=C(D,x)+C(L,x)>0\}.
$$

Differences between an element $y$ and a set $X$:
$$
\Delta(y,X)=\{|y-x_i|:x_i\in X\}.
$$
All pairwise differences of $X=\{x_1,\ldots,x_n\}$:
$$
\Delta X=\{|x_j-x_i|:1\le i<j\le n\}.
$$
A shift by $k$ is $X\oplus\{k\}=\{x+k:x\in X\}$. The reflection is $-X=\{-x:x\in X\}$.

## 5 Partial digest problem

**Partial Digest Problem (PDP).** Given a multiset of $N=\binom{n}{2}$ positive integers $D=\{d_1,\ldots,d_N\}$, find a set of $n$ integers $X=\{x_1,\ldots,x_n\}$ such that $\Delta X=D$.

The input is the multiset of all pairwise integer distances. The output is a set $X$ of $n$ integers with those distances. In the related Turnpike problem the same mathematics reconstructs highway-exit positions from distances in miles.

Two symmetries mean $X$ is not unique:

- $\Delta X=\Delta(X\oplus\{k\})$ for any shift $k$.
- $\Delta X=\Delta(-X)$. Reflection does not change the multiset of pairwise distances.

Example: $D=\{2,2,3,3,4,5,6,7,8,10\}$ and $X=\{0,2,4,7,10\}$. Then $X\oplus\{100\}=\{100,102,104,107,110\}$ and $-X=\{0,-2,-4,-7,-10\}$ give the same $D$.

Those distances can be written as a table with the elements of $X$ along the top and the side. The entries are the positive differences.

### 5.1 Brute force PDP

```text
PDP brute force(D, n)
  M ← max D
  for every set of n−2 integers 0 < x2 < … < x_{n−1} < M do
    X ← {0, x2, …, x_{n−1}, M}
    form ΔX from X
    if ΔX = D then return X
  output "No solution"
```

The search considers $\binom{M-1}{n-2}$ candidate sets. Endpoints are fixed at $0$ and $M=\max D$, which removes the shift freedom. Reflection can still remain.

### 5.2 The PLACE algorithm

A practical backtracking algorithm places the largest remaining distance at one end or the other.

```text
PDP(D)
  M ← max D
  delete M from D
  X ← {0, M}
  PLACE(D, X)

PLACE(D, X)
  if D is empty then output X and return
  y ← max D
  z ← M − y
  if Δ(y, X) ⊆ D then
    D ← D \ Δ(y, X)
    X ← X ∪ {y}
    PLACE(D, X)
    restore X and D
  if Δ(z, X) ⊆ D then
    D ← D \ Δ(z, X)
    X ← X ∪ {z}
    PLACE(D, X)
    restore X and D
  return
```

After each recursive call, $X$ and $D$ are restored, so both branches can be explored. The algorithm lists every $X$ on the search tree with $\Delta X=D$.

Tracing $D=\{2,2,3,3,4,5,6,7,8,10\}$ (after deleting $10$):

| Node | $D$ still to place | $X$ |
| --- | --- | --- |
| A | $\{2,2,3,3,4,5,6,7,8\}$ | $\{0,10\}$ |
| B | $\{2,3,3,4,5,6,7\}$ | $\{0,8,10\}$ |
| C | $\{2,3,4,6\}$ | $\{0,3,8,10\}$ |
| D | $\emptyset$ | $\{0,3,6,8,10\}$ |
| E | $\{2,3,3,4,5,6,7\}$ | $\{0,2,10\}$ |
| F | $\{2,3,4,6\}$ | $\{0,2,7,10\}$ |
| G | $\emptyset$ | $\{0,2,4,7,10\}$ |

A mark $x$ on a branch means $\Delta(y,X)\not\subseteq D$. A mark $x'$ means $\Delta(z,X)\not\subseteq D$. The two empty-$D$ leaves are the solutions $\{0,3,6,8,10\}$ and $\{0,2,4,7,10\}$, which are reflections of each other once the left endpoint is fixed at $0$.

![[./attachments/Pasted image 20260925075839.png|636]]

### 5.3 Cost of PDP

Let $T(N)$ be the maximum time for $|D|=N$. In the worst case both placements stay viable:
$$
T(N)=2T(N-1)+O(N)>2T(N-1),
$$
so $T(N)=2^k T(N-k)$. The running time is **exponential**.

If only one alternative is viable at every step,
$$
T(N)=T(N-1)+O(N)=c\frac{N(N+1)}{2}=O(N^2).
$$
The running time is then quadratic.

## 6 Double digest

A double digest records fragment lengths from three complete digests:

1. Restriction enzyme $A$ alone.
2. Restriction enzyme $B$ alone.
3. Enzymes $A$ and $B$ together.

*In each experiment the DNA is cut at every occurrence of the relevant sites(complete digest).*

For a molecule of length $L$, the experiments yield multisets (written in the slides as sorted lists)

$$
A=\{a_i:1\le i\le n\},\quad
B=\{b_i:1\le i\le m\},\quad
C=\{c_i:1\le i\le \ell\},
$$

with $\sum a_i=\sum b_i=\sum c_i=L$.

Let $\sigma$ be a permutation of $(1,\ldots,n)$ and $\mu$ a permutation of $(1,\ldots,m)$. The pair $(\sigma,\mu)$ is a **configuration**: it is a choice of order for the $A$ fragments and the $B$ fragments along the molecule. The cut-site locations are the partial sums
$$
S=\Bigl\{\sum_{j=1}^{r}a_{\sigma(j)}\Bigr\}\cup\Bigl\{\sum_{j=1}^{t}b_{\mu(j)}\Bigr\},
$$
including the empty sums $0$, for $0\le r\le n$ and $0\le t\le m$. Label the sorted sites $S=\{s_j:0\le j\le \ell\}$ with $s_i<s_j$ for $i<j$.

The double-digest lengths implied by that configuration are the consecutive gaps
$$
C(\sigma,\mu)=\{c_j(\sigma,\mu)=s_j-s_{j-1}:1\le j\le \ell\}.
$$

**Double Digest Problem (DDP).** Given experimental multisets $A$, $B$, and $C$, find a configuration $(\sigma,\mu)$ such that $C=C(\sigma,\mu)$ - Find the restriction sites of the two enzymes that can result in $C$ by ordering the complete digest segments.

### 6.1 Double digest solutions need not be unique

One example from the slides:
$$
A=\{1,3,3,12\},\quad B=\{1,2,3,3,4,6\},\quad C=\{1,1,1,1,2,2,2,3,6\}.
$$
Two different orders of the $B$ fragments are compatible with the same $A$ and the same $C$. Fixing $A$ in the order $1,3,12,3$ (positions $0,1,4,16,19$), one solution orders $B$ as $2,4,6,3,3,1$ and another orders $B$ as $3,3,6,1,2,4$. Both reproduce $C$ after the $A$ and $B$ sites are merged.

![[./attachments/Pasted image 20260925140519.png|548]]

### 6.2 DDP is NP-complete

NP-completeness needs two facts: a proposed solution can be **checked in polynomial time**, and some NP-complete problem reduces to DDP in polynomial time.

**Partition Problem.** Given a finite set $S$ with $|S|=n$ and a positive integer $f(a)$ for each $a\in S$, decide whether some subset $S'\subset S$ satisfies
$$
\sum_{a\in S'}f(a)=\sum_{a\in S\setminus S'}f(a).
$$
If $K=\sum_{a\in S}f(a)$ is odd, no such subset exists. Partition is NP-complete. A standard chain of reductions is $3\mathrm{SAT}\le_P\cdots\le_P\mathrm{VC}\le_P\mathrm{SUBSET\text{-}SUM}\le_P\mathrm{PP}$.

**DDP is in NP.** Given $(\sigma,\mu)$ and multisets of sizes $n$, $m$, and $\ell$, the merged site set $S$ can be built in time polynomial in those sizes. The implied lengths $C(\sigma,\mu)$ take $O(\ell)$ time. Sorting $C(\sigma,\mu)$ and $C$ takes less than $O(\ell^2)$, and comparing the ordered lists takes $O(\ell)$.

**Partition reduces to DDP.** Feed DDP the instance
$$
A=\{f(a):a\in S\},\qquad B=\{K/2,\,K/2\},\qquad C=A.
$$
A valid double-digest configuration exists exactly when a balanced partition exists. (The two $B$ fragments of length $K/2$ place a single interior cut that splits $A$ into two groups of equal total length, and the double-digest multiset stays equal to $A$ only for that kind of split.)

### 6.3 Interval graphs

Let $S$ be a set and $F=\{S_1,\ldots,S_n\}$ a family of distinct nonempty subsets of $S$. The **intersection graph** $I(F)$ has a vertex for each set in $F$, and an edge between $S_i$ and $S_j$ when $S_i\cap S_j\neq\emptyset$. An **interval graph** is an intersection graph in which every set in $F$ is an interval of the real line.

Each restriction enzyme cleaves DNA into intervals. Enzymes $A$ and $B$, used separately, produce two families of intervals. *If the experiment reports which $A$ interval overlaps which $B$ interval, the restriction map can be reconstructed.* The algorithmn **assumes that no interior $A$ site coincides with an interior $B$ site**. The ends of the molecule may still coincide.

### 6.4 The incidence matrix

Consider maps whose $A$ intervals are $A_2,A_3,A_1$ and whose $B$ intervals are $B_4,B_2,B_3,B_1$, written left to right. Let $S_A=\{A_1,A_2,A_3\}$ and $S_B=\{B_1,B_2,B_3,B_4\}$. The incidence matrix $I_{A,B}$ has $I_{A,B}(i,j)=1$ when interval $A_i$ overlaps interval $B_j$. For the labelling in the slides,
$$
I_{A,B}=\begin{pmatrix}1&0&1&0\\0&0&0&1\\0&1&1&1\end{pmatrix}.
$$

The bipartite graph $G_{A,B}$ has vertices $S_A\cup S_B$ and an edge $(A_i,B_j)$ whenever the intervals overlap. Highlight vertices of degree at least $2$. An end-interval vertex $a$ may still have degree at least $2$, but it always has $N_2(a)=1$, where $N_2(a)$ is the number of neighbours of $a$ that themselves have degree at least $2$.

### 6.5 Ordering the edges

```text
OrderEdges(G(V, E))
  L ← {v ∈ V : deg(v) ≥ 2}
  L0 ← L
  E0 ← E
  while E0 ≠ ∅ do
    if some v, u ∈ L0 have (u, v) ∈ E and N2(v) = 1 then
      v is an end
      v1 ← v, v2 ← u
      extend (v1, v2), …, (v_{k−1}, v_k) as far as possible
        with every vi ∈ L0
      append that path to E†
      delete those vertices from L0 and those edges from E0
    else
      for v ∈ L ∩ L0^c and u ∈ L^c with (u, v) ∈ E do
        insert (v, u) next to the earlier edge that touches v
        delete (v, u) from E0
  return E†
```

When several pendant edges must sit between the same pair of earlier edges, the algorithm uses one of the possible orders. If an edge cannot be inserted between previous edges, it is placed at the left or right end, beside the edge that already touches $v$.

On the example, $L_0=L=\{A_1,A_3,B_3,B_4\}$. Starting at the end $B_4$ and extending through high-degree vertices gives
$$
E^\dagger=((B_4,A_3),(B_3,A_3),(B_3,A_1))
$$
and then $L_0=\emptyset$.

The remaining edges have one end in $L$ (already removed from $L_0$) and one end outside $L$:

- $(A_3,B_2)$ is inserted between $(B_4,A_3)$ and $(B_3,A_3)$.
- $(A_1,B_1)$ is inserted to the right of $(B_3,A_1)$.
- $(B_4,A_2)$ is inserted to the left of $(B_4,A_3)$.

So
$$
E^\dagger=((B_4,A_2),(B_4,A_3),(A_3,B_2),(A_3,B_3),(B_3,A_1),(A_1,B_1)).
$$

That order is the order of the intervals on the two maps:
$$
A_2,\; A_3,\; A_1 \qquad\text{and}\qquad B_4,\; B_2,\; B_3,\; B_1.
$$
Drawing the edges between those intervals gives the relative positions of all restriction sites.

![[./attachments/Pasted image 20260925141318.png|561]]

# References

- Minary, Computational Biology, Lecture 2 (Michaelmas 2025)
- Jones & Pevzner, *An Introduction to Bioinformatics Algorithms* (complete-digest and partial-digest figures)
