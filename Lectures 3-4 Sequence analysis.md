---
tags:
  - source-materials/lectures
  - bio/sequencing
prerequisites: []
course date: Michaelmas 2025
lectuer: Prof. Peter Minary
share: true
---
# Lectures 3–4 Analysing nucleic acid and protein sequences

Readable notes from Minary, *Computational Biology*, Lectures 3–4. Slides: [03-04-nucleic-acid-protein-sequences.pdf](03-04-nucleic-acid-protein-sequences.pdf).


## 1 DNA sequencing of a few hundred bases

Nucleic-acid databases depend on being able to read DNA. The pioneering methods of 1976 were recognised by the Nobel Prize to Gilbert and Sanger in 1980. Those methods read strings of up to about $500$ base pairs. Sequencing by hybridization, an alternative, was developed in 1994.

For a read of about $350$–$1{,}000$ base pairs, identical copies of the same single-stranded DNA are **put through one reaction for each base**. Each reaction produces a collection of single-stranded molecules that **end at one of A, G, C, or T**. For `5'-ATCCACGGCCATAGG-3'`, the data are all $5'$ prefixes that end in a given base. Gel electrophoresis separates those fragments. Reading across the four lanes, from the shortest band at the $5'$ end toward longer bands, recovers the original strand from $5'$ to $3'$.

![[./attachments/Pasted image 20260925152719.png|360]]

## 2 Shotgun sequencing

Shotgun sequencing uses reads of length $\ell\approx 500$–$700$, which can be acquire from the method described in the last section.

1. Obtain a **large sample of** genomic DNA.
2. Randomly break each piece into inserts of that length.
3. Sequence the inserts.
4. **Assemble the reads computationally into one long sequence.**

This is a simplified description. The computational problem is: a molecule $a=a_1\ldots a_L$ is to be sequenced; the data are strings $\{f_1,\ldots,f_N\}$ from unknown locations in $a$; recover $a$. 

*We assume all strings were read from the same orientation, and reads contain no error.*

> We need to consider reverse complement if we want to consider which strand(orientation) the fragment came from

### 2.1 Shortest superstring problem


> [!def] **Shortest Superstring Problem (SSP).**
> Given $F=\{f_1,\ldots,f_N\}$, no string in $F$ a substring of another, find a shortest string $S$ such that every $f_i$ is a substring of $S$.

**Prefix graph.** For $f_i=uv$ and $f_j=vw$, take the longest such $v$. The empty prefix $u$ is not allowed(*no substring*); the empty overlap $v$ is allowed. Then $\mathrm{overlap}(f_i,f_j)=|v|$ and the prefix length $\mathrm{pf}(f_i,f_j)=|u|$. The prefix graph has a vertex for each fragment and a directed edge $(f_i,f_j)$ of weight $\mathrm{pf}(f_i,f_j)$. Superstrings correspond to cycles in this graph.

> [!example] 
> $F=\{\mathtt{ATAT},\mathtt{TATT},\mathtt{TTAT},\mathtt{TATA}\}$ has a cycle whose superstring is $\mathtt{ATATTATA}$, of length $8$. Another set $\{\mathtt{ACT},\mathtt{CTA},\mathtt{AGT}\}$ gives the cycle superstring $\mathtt{ACTAGT}$ of length $6$, with a zero overlap on the closing edge. Prefix length $0$ is forbidden; overlap $0$ is allowed.
> 
> Prefix graph of $\{\mathtt{ATAT},\mathtt{TATT},\mathtt{TTAT},\mathtt{TATA}\}$. One vertex on each corner of the square. Bold edges are the cycle $\mathtt{ATAT}\to\mathtt{TATT}\to\mathtt{TTAT}\to\mathtt{TATA}$, prefix weights $1,2,1,1$, superstring $\mathtt{ATATTATA}$. Grey edges are the other pairs.
> 
> ```tikz
> \begin{document}
> \begin{tikzpicture}[x=2.6cm,y=2.6cm,>=latex,
>   v/.style={circle,draw,minimum size=10mm,inner sep=0pt,font=\scriptsize\ttfamily},
>   cyc/.style={->,very thick},
>   oth/.style={->,draw=black!45},
>   w/.style={font=\scriptsize,fill=white,inner sep=1pt}]
>   \node[v] (ATAT) at (0,1) {ATAT};
>   \node[v] (TATT) at (1,1) {TATT};
>   \node[v] (TTAT) at (1,0) {TTAT};
>   \node[v] (TATA) at (0,0) {TATA};
>   \draw[cyc] (ATAT) -- node[w,below] {$1$} (TATT);
>   \draw[cyc] (TATT) -- node[w,left] {$2$} (TTAT);
>   \draw[cyc] (TTAT) -- node[w,above] {$1$} (TATA);
>   \draw[cyc] (TATA) -- node[w,right] {$1$} (ATAT);
>   \draw[oth] (TATT) to[bend right=30] node[w,above] {$4$} (ATAT);
>   \draw[oth] (TTAT) to[bend right=30] node[w,right] {$1$} (TATT);
>   \draw[oth] (TATA) to[bend right=30] node[w,below] {$4$} (TTAT);
>   \draw[oth] (ATAT) to[bend right=30] node[w,left] {$1$} (TATA);
>   \draw[oth] (ATAT) to[bend left=14] node[w,pos=0.32,auto=left] {$3$} (TTAT);
>   \draw[oth] (TTAT) to[bend left=14] node[w,pos=0.32,auto=left] {$2$} (ATAT);
>   \draw[oth] (TATT) to[bend left=14] node[w,pos=0.32,auto=left] {$3$} (TATA);
>   \draw[oth] (TATA) to[bend left=14] node[w,pos=0.32,auto=left] {$2$} (TATT);
> \end{tikzpicture}
> \end{document}
> ```
> 
> Prefix graph of $\{\mathtt{ACT},\mathtt{CTA},\mathtt{AGT}\}$. Bold edges spell $\mathtt{ACTAGT}$. The closing edge $\mathtt{AGT}\to\mathtt{ACT}$ has prefix weight $3$ and overlap $0$.
> 
> ```tikz
> \begin{document}
> \begin{tikzpicture}[x=2.4cm,y=2.4cm,>=latex,
>   v/.style={circle,draw,minimum size=9mm,inner sep=0pt,font=\scriptsize\ttfamily},
>   cyc/.style={->,very thick},
>   oth/.style={->,draw=black!45},
>   w/.style={font=\scriptsize,fill=white,inner sep=1pt}]
>   \node[v] (ACT) at (0,1.5) {ACT};
>   \node[v] (CTA) at (1.3,0) {CTA};
>   \node[v] (AGT) at (-1.3,0) {AGT};
>   \draw[cyc] (ACT) -- node[w,auto=right] {$1$} (CTA);
>   \draw[cyc] (CTA) -- node[w,auto=right] {$2$} (AGT);
>   \draw[cyc] (AGT) -- node[w,auto=right] {$3$} (ACT);
>   \draw[oth] (CTA) to[bend right=28] node[w,auto=right] {$2$} (ACT);
>   \draw[oth] (AGT) to[bend right=28] node[w,auto=right] {$3$} (CTA);
>   \draw[oth] (ACT) to[bend right=28] node[w,auto=right] {$3$} (AGT);
> \end{tikzpicture}
> \end{document}
> ```
> 
If $S_{\mathrm{MWHC}}(G)$ is the total prefix weight of[[ | ]]a minimum-weight Hamiltonian cycle in $G$, then $S_{\mathrm{MWHC}}(G)\le |S_{\mathrm{opt}}|$. A Hamiltonian cycle visits each vertex once, so the minimum-weight one is a travelling-salesman problem. The optional proof compares that cycle with the cycle of an **optimal superstring**: the cycle weight omits the final overlap, so it is at most the superstring length.

SSP can therefore be cast as TSP on the prefix graph. TSP is NP-hard, so the reduction is not an efficient algorithm. A greedy algorithm repeatedly merges the pair with the largest overlap until one string remains. Early DNA assemblers used this greedy strategy. It is at most $2.75$ times optimal.

## 3 Sequencing by hybridization

Shotgun sequencing did not have a good computational solution, to improve upon this, DNA chips were invented.

SBH was suggested in 1988. Build a two-dimensional array of all DNA $k$-tuples. A distinct probe sits at each position. A labelled single-stranded target is applied. Each $k$-tuple present in the target hybridizes to its reverse complement on the array. For example the probe `ACCGTGGA` ($k=8$) hybridizes to `CCCTGGCACCTA` because it is complementary to a substring. Unhybridized DNA is washed off, and the hybridized probes are detected from the label. The reverse complements of those probes are the $k$-mer composition of the target, from which the target is reconstructed.

### 3.1 Spectrum and the SBH problem

The **spectrum** $S(a,k)$ of a string $a$ of length $n$ is the multiset of its $n-k+1$ consecutive $k$-mers. If $k=3$ and $a=\mathtt{TATGGTGC}$, then
$$
S(a,3)=\{\mathtt{TAT},\mathtt{ATG},\mathtt{TGG},\mathtt{GGT},\mathtt{GTG},\mathtt{TGC}\}.
$$


> [!def] **SBH problem.**
> Given the set $F$ of all $k$-mers of an unknown string $a$, reconstruct a string with $S(a,k)=F$.

### 3.2 SBH as a Hamiltonian path

$k$-mers $p$ and $q$ overlap when $\mathrm{overlap}(p,q)=k-1$: the last $k-1$ letters of $p$ are the first $k-1$ letters of $q$. The graph $H$ has a vertex for each $k$-mer and a directed edge $(p,q)$ when they overlap. A Hamiltonian path visits every vertex and spells a reconstruction.

For $S(a,3)=\{\mathtt{ATG},\mathtt{AGG},\mathtt{TGC},\mathtt{TCC},\mathtt{GTC},\mathtt{GGT},\mathtt{GCA},\mathtt{CAG}\}$ there is one Hamiltonian path, giving $a=\mathtt{ATGCAGGTCC}$. For $\{\mathtt{ATG},\mathtt{TGG},\mathtt{TGC},\mathtt{GTG},\mathtt{GGC},\mathtt{GCA},\mathtt{GCG},\mathtt{CGT}\}$ there are two, giving $\mathtt{ATGCGTGGCA}$ and $\mathtt{ATGGCGTGCA}$. Hamiltonian path is NP-complete, so this reduction is not a practical algorithm.

One Hamiltonian path, $a=\mathtt{ATGCAGGTCC}$. Vertices are in spectrum order.

```tikz
\begin{document}
\begin{tikzpicture}[x=1.15cm,y=1cm,>=latex,
  v/.style={circle,draw,minimum size=7.5mm,inner sep=0pt,font=\scriptsize}]
  \node[v] (ATG) at (0,0) {ATG};
  \node[v] (AGG) at (1,0) {AGG};
  \node[v] (TGC) at (2,0) {TGC};
  \node[v] (TCC) at (3,0) {TCC};
  \node[v] (GTC) at (4,0) {GTC};
  \node[v] (GGT) at (5,0) {GGT};
  \node[v] (GCA) at (6,0) {GCA};
  \node[v] (CAG) at (7,0) {CAG};
  \draw[->] (GGT) -- (GTC);
  \draw[->] (GTC) -- (TCC);
  \draw[->] (GCA) -- (CAG);
  \draw[->] (ATG) .. controls (0.45,0.5) and (1.55,0.5) .. (TGC);
  \draw[->] (TGC) .. controls (2.7,1.05) and (5.3,1.05) .. (GCA);
  \draw[->] (AGG) .. controls (2.1,-0.7) and (4.2,-0.7) .. (GGT);
  \draw[->] (CAG) .. controls (5.6,-1.15) and (2.4,-1.15) .. (AGG);
\end{tikzpicture}
\end{document}
```

Two Hamiltonian paths on the second spectrum. Red edges belong only to $a=\mathtt{ATGCGTGGCA}$, green edges only to $a=\mathtt{ATGGCGTGCA}$, and black edges to both.

```tikz
\begin{document}
\begin{tikzpicture}[x=1.15cm,y=1cm,>=latex,
  v/.style={circle,draw,minimum size=7.5mm,inner sep=0pt,font=\scriptsize},
  onlya/.style={->,red},
  onlyb/.style={->,green!55!black},
  both/.style={->,black}]
  \node[v] (ATG) at (0,0) {ATG};
  \node[v] (TGG) at (1,0) {TGG};
  \node[v] (TGC) at (2,0) {TGC};
  \node[v] (GTG) at (3,0) {GTG};
  \node[v] (GGC) at (4,0) {GGC};
  \node[v] (GCA) at (5,0) {GCA};
  \node[v] (GCG) at (6,0) {GCG};
  \node[v] (CGT) at (7,0) {CGT};
  \draw[onlyb] (ATG) -- (TGG);
  \draw[onlyb] (GTG) -- (TGC);
  \draw[onlya] (GGC) -- (GCA);
  \draw[both] (GCG) -- (CGT);
  \draw[onlya] (ATG) .. controls (0.4,0.48) and (1.55,0.48) .. (TGC);
  \draw[onlyb] (GGC) .. controls (4.4,0.5) and (5.6,0.5) .. (GCG);
  \draw[both] (TGG) .. controls (1.7,0.95) and (3.4,0.95) .. (GGC);
  \draw[onlya] (TGC) .. controls (2.8,1.3) and (5.2,1.3) .. (GCG);
  \draw[onlya] (GTG) .. controls (2.55,-0.45) and (1.45,-0.45) .. (TGG);
  \draw[onlyb] (TGC) .. controls (2.6,-0.8) and (4.4,-0.8) .. (GCA);
  \draw[both] (CGT) .. controls (6.2,-1.25) and (3.8,-1.25) .. (GTG);
\end{tikzpicture}
\end{document}
```
This construction is NP-Complete, still not very useful.
### 3.3 SBH as an Eulerian path

Build a different graph $G$. Vertices are the $(k-1)$-mers that appear inside the spectrum (not a multiset, although in reality the spectrum is a multiset). Direct an edge from $u$ to $v$ when some $k$-mer in the spectrum has prefix $u$ and suffix $v$. Each $k$-mer is an edge. An Eulerian path visits every edge and spells a reconstruction. On the second spectrum above, the two Eulerian paths are the same two strings.

Vertices are the $(k-1)$-mers. Each spectrum $k$-mer is one directed edge, from its prefix to its suffix.

```tikz
\begin{document}
\begin{tikzpicture}[x=1.05cm,y=1cm,>=latex,
  v/.style={circle,draw,minimum size=7.5mm,inner sep=0pt,font=\scriptsize}]
  \node[v] (GT) at (1.7,1.55) {GT};
  \node[v] (CG) at (3.4,1.55) {CG};
  \node[v] (AT) at (0,0) {AT};
  \node[v] (TG) at (1.7,0) {TG};
  \node[v] (GC) at (3.4,0) {GC};
  \node[v] (CA) at (5.1,0) {CA};
  \node[v] (GG) at (2.55,-1.45) {GG};
  \draw[->] (AT) -- (TG);
  \draw[->] (TG) -- (GC);
  \draw[->] (GC) -- (CA);
  \draw[->] (CG) -- (GT);
  \draw[->] (GT) -- (TG);
  \draw[->] (GC) -- (CG);
  \draw[->] (TG) -- (GG);
  \draw[->] (GG) -- (GC);
\end{tikzpicture}
\end{document}
```

Both Eulerian paths use every edge, so each edge is drawn twice. Red, read in order, spells $a=\mathtt{ATGCGTGGCA}$ ($\mathtt{ATG},\mathtt{TGC},\mathtt{GCG},\mathtt{CGT},\mathtt{GTG},\mathtt{TGG},\mathtt{GGC},\mathtt{GCA}$). Green spells $a=\mathtt{ATGGCGTGCA}$ ($\mathtt{ATG},\mathtt{TGG},\mathtt{GGC},\mathtt{GCG},\mathtt{CGT},\mathtt{GTG},\mathtt{TGC},\mathtt{GCA}$).

```tikz
\begin{document}
\begin{tikzpicture}[x=1.05cm,y=1cm,>=latex,
  v/.style={circle,draw,fill=white,minimum size=7.5mm,inner sep=0pt,font=\scriptsize},
  reda/.style={->,red},
  grnb/.style={->,green!55!black}]
  \node[v] (GT) at (1.7,1.55) {GT};
  \node[v] (CG) at (3.4,1.55) {CG};
  \node[v] (AT) at (0,0) {AT};
  \node[v] (TG) at (1.7,0) {TG};
  \node[v] (GC) at (3.4,0) {GC};
  \node[v] (CA) at (5.1,0) {CA};
  \node[v] (GG) at (2.55,-1.45) {GG};
  \draw[reda] ([yshift=1.5pt]AT.east) -- ([yshift=1.5pt]TG.west);
  \draw[grnb] ([yshift=-1.5pt]AT.east) -- ([yshift=-1.5pt]TG.west);
  \draw[reda] ([yshift=1.5pt]TG.east) -- ([yshift=1.5pt]GC.west);
  \draw[grnb] ([yshift=-1.5pt]TG.east) -- ([yshift=-1.5pt]GC.west);
  \draw[reda] ([yshift=1.5pt]GC.east) -- ([yshift=1.5pt]CA.west);
  \draw[grnb] ([yshift=-1.5pt]GC.east) -- ([yshift=-1.5pt]CA.west);
  \draw[reda] ([yshift=1.5pt]CG.west) -- ([yshift=1.5pt]GT.east);
  \draw[grnb] ([yshift=-1.5pt]CG.west) -- ([yshift=-1.5pt]GT.east);
  \draw[reda] ([xshift=-1.5pt]GT.south) -- ([xshift=-1.5pt]TG.north);
  \draw[grnb] ([xshift=1.5pt]GT.south) -- ([xshift=1.5pt]TG.north);
  \draw[reda] ([xshift=-1.5pt]GC.north) -- ([xshift=-1.5pt]CG.south);
  \draw[grnb] ([xshift=1.5pt]GC.north) -- ([xshift=1.5pt]CG.south);
  \draw[reda] (TG) to[bend right=16] (GG);
  \draw[grnb] (TG) to[bend left=16] (GG);
  \draw[reda] (GG) to[bend right=16] (GC);
  \draw[grnb] (GG) to[bend left=16] (GC);
\end{tikzpicture}
\end{document}
```

For a directed graph whose nonzero-degree vertices lie in one connected component of the underlying undirected graph, there is an Eulerian path from $s$ to $t$ if and only if every other vertex is balanced ($\mathrm{in}(v)=\mathrm{out}(v)$), $\mathrm{in}(s)-\mathrm{out}(s)=-1$, and $\mathrm{in}(t)-\mathrm{out}(t)=1$. Adding an edge from $t$ to $s$ balances every vertex.


> [!theorem] 
> A connected graph has an Eulerian cycle if and only if every vertex is balanced.

`\begin{proof}`
If a cycle uses every edge, every arrival is matched by a departure, so every vertex is balanced. Conversely, start at any vertex and walk unused edges. Balance forces the walk to stop only at its start, so the walk is a cycle. If unused edges remain, they form a balanced subgraph, and some vertex $w$ on the first cycle still has an unused edge. Walk a second cycle from $w$ and splice it in: follow the first cycle from $v$ to $w$, the second cycle around $w$, then the rest of the first cycle. Repeat until no edges remain. This runs in time linear in the number of edges.

SBH therefore has a linear-time reconstruction algorithm once it is posed as an Eulerian path.
`\end{proof}`

**Fragment assembly (optional).** Assembling a genome from reads of length $500$–$1{,}000$ has three steps: find overlaps, lay the reads out along the DNA, and derive a consensus sequence. SSP is the most relevant problem covered here, but it is too simple to capture fragment assembly.
## 4 Motif finding

Short DNA strings upstream of genes are regulatory motifs. Transcription factors bind them and can recruit RNA polymerase, which transcribes the downstream gene. In eukaryotes, Pol II makes mRNA precursors, Pol I makes ribosomal RNA, and Pol III makes tRNA, 5S rRNA, and other small RNAs. The computational problem is that one sees only the longer sequences that contain the motif, not the motif itself.

A pattern $P$ of length $l=8$ planted in seven sequences of length about $32$–$40$ is unlikely to arise by chance: fewer than $280/4^8\approx 0.004$ random $8$-mers of that total length are expected to match a fixed pattern, while $P$ appears at least seven times. Real motifs are *mutated*, so one string cannot represent every copy.

Let $H$ be $m$ DNA sequences of length $n$. A vector of start positions $s=(s_1,\ldots,s_m)$, with $1\le s_i\le n-l+1$, picks one $l$-mer from each sequence. Those $l$-mers form an $m\times l$ alignment matrix $A$, where $A(i,j)$ is nucleotide $s_i+j-1$ of sequence $i$. The worked example uses $l=8$, $m=7$, $n=40$, and $s=(8,19,3,5,31,27,15)$.

![[./attachments/Pasted image 20261020091713.png|Pasted image 20261020091713.png]]
> Worked example, 8-mers with mutations in sequences of length 40

The $4\times l$ **profile** $P$ counts how often each nucleotide appears in each column of $A$. The **consensus string** takes the most frequent nucleotide in each column. The **consensus score** is
$$
S_c=\sum_{j=1}^{l}\max_i P(i,j).
$$
In the example the column maxima are $5,5,6,4,5,5,6,6$, so $S_c=42$.

![[./attachments/Pasted image 20261020091858.png|332]]
> Alignment matrix and consensus score


> [!def] **Motif Finding Problem** 
> Given an $m\times n$ sequence matrix $M_{\mathrm{DNA}}$ and a length $l$, find start positions $s$ maximising $S_c(s,M_{\mathrm{DNA}})$.

**Brute force solution.** tries every $s$ in $\{1,\ldots,n-l+1\}^m$. There are $(n-l+1)^m$ such vectors, and each score costs $O(l)$ after the profile is built, so the complexity is $O(l n^{m})$, which is how the slide writes $O(ln^{m})$. The exact count of start vectors is $(n-l+1)^{m}$.

## 5 Median string problem

> [!def] **Hamming distance**
> The Hamming distance $d_H(v,w)$ between two $l$-mers is the number of positions at which they differ. 
> For example $d_H(\mathtt{GCATTA},\mathtt{GCAATT})=2$. For a start vector $s$ and an $l$-mer $v$,
> $$
> d_H(v,s)=\sum_{i=1}^{m} d_H(v,s_i),
> $$
> where $s_i$ here means the $l$-mer starting at $s_i$ in sequence $i$. 

Then
$$
d(v,M_{\mathrm{DNA}})=\min_s d_H(v,s),
$$
and the minimum can be taken *separately* in each of the $m$ sequences. 

The **median string** is
$$
v_m=\arg\min_v d(v,M_{\mathrm{DNA}}),
$$
the minimum running over all $4^l$ strings of length $l$.

> [!def] **Median String Problem.** 
> Given $M_{\mathrm{DNA}}$ and $l$, return a string $v$ of length $l$ minimising $d(v,M_{\mathrm{DNA}})$.

Motif finding and the median-string problem are equivalent. If $w$ is the consensus of the alignment at starts $s$, then
$$
d_H(w,s)=lm-S_c(s,M_{\mathrm{DNA}}).
$$
In the example each of the seven rows differs from `ATGCAACT` in two positions, so $d_H(w,s)=14$, and $14=8\cdot 7-42$. The consensus minimises $d_H(v,s)$ over $v$, and because $l$ and $m$ are fixed,
$$
\min_s\min_v d_H(v,s)=lm-\max_s S_c(s,M_{\mathrm{DNA}}).
$$
The median string solves motif finding: in each sequence, take the substring closest in Hamming distance to the median string. The motif-finding solution solves the median-string problem: its consensus is the median string.

**Brute force** for the median string scans all $4^l$ candidates. Each distance costs $O(nm)$, so the complexity is $O(4^l nm)$. For $l$ about $8$–$15$ and sequence regions of length about $1{,}000$, *this is the more efficient of the two brute-force searches when $l$ and $m$ are comparable*.

## 6 Sequence comparison

Sequence comparison measures similarity. A newly sequenced gene can be assigned a function by similarity to a gene of known function. A protein sequence can be assigned structural information by similarity to a sequence of known structure.

Hamming distance compares strings of equal length. Strings of different lengths need edit distance: the minimum number of insertions, deletions, and substitutions that turn one string into the other.

### 6.1 Alignment as a path

An alignment of $v$ ($n$ characters) and $w$ ($m$ characters) is a two-row matrix. The characters of each string stay in order, spaces may be inserted, and a column may not be a space in both rows. The matrix therefore has at most $n+m$ columns.

```text
A T – G T T A T –
A T C G T – A – C
```

Columns with the same letter are matches. Columns with different letters are mismatches. A space on the top row is an insertion; a space on the bottom row is a deletion. Insertions and deletions together are indels.

Writing how many symbols have been consumed gives coordinates. `AT-GTTAT-` becomes $1,2,2,3,4,5,6,7,7$ and `ATCGT-A-C` becomes $1,2,3,4,5,5,6,6,7$. The columns are points on an $n\times m$ grid, and the alignment is a path from $(0,0)$ to $(n,m)$. Every path is an alignment and every alignment is a path. Diagonal, horizontal, and vertical edges are the columns of the alignment matrix.

The edit graph for $v=\mathtt{ATGTTAT}$ and $w=\mathtt{ATCGTAC}$. Horizontal edges consume a letter of $w$, vertical edges a letter of $v$, and diagonal edges one letter of each. The bold path is one alignment,
```text
A T G T T - A T
A T - C G T A C
```

```tikz
\begin{document}
\begin{tikzpicture}[x=0.9cm,y=0.88cm,>=latex,
  dot/.style={circle,draw,fill=white,inner sep=0pt,minimum size=1.55mm}]
  \foreach \r in {0,...,7}
    \foreach \c in {0,...,7}
      \coordinate (p\c\r) at (\c,-\r);
  \foreach \r in {0,...,7}
    \foreach \c [evaluate=\c as \n using int(\c+1)] in {0,...,6}
      \draw[->,gray,thin,shorten >=0.9mm,shorten <=0.9mm] (p\c\r) -- (p\n\r);
  \foreach \c in {0,...,7}
    \foreach \r [evaluate=\r as \n using int(\r+1)] in {0,...,6}
      \draw[->,gray,thin,shorten >=0.9mm,shorten <=0.9mm] (p\c\r) -- (p\c\n);
  \foreach \r [evaluate=\r as \nr using int(\r+1)] in {0,...,6}
    \foreach \c [evaluate=\c as \nc using int(\c+1)] in {0,...,6}
      \draw[->,gray,thin,shorten >=0.9mm,shorten <=0.9mm] (p\c\r) -- (p\nc\nr);
  \draw[line width=1.05pt]
    (p00) -- (p11) -- (p22) -- (p32) -- (p43) -- (p54) -- (p55) -- (p66) -- (p67) -- (p77);
  \foreach \r in {0,...,7}
    \foreach \c in {0,...,7}
      \node[dot] at (p\c\r) {};
  \foreach \c/\r in {0/0,1/1,2/2,3/2,4/3,5/4,5/5,6/6,7/7}
    \node[dot,fill] at (p\c\r) {};
  \foreach \c/\t in {1/A,2/T,3/C,4/G,5/T,6/A,7/C}
    \node[font=\scriptsize] at (\c,0.85) {\t};
  \foreach \r/\t in {1/A,2/T,3/G,4/T,5/T,6/A,7/T}
    \node[font=\scriptsize] at (-0.85,-\r) {\t};
  \node[font=\scriptsize] at (3.5,1.45) {$w$};
  \node[font=\scriptsize] at (-1.45,-3.5) {$v$};
  \foreach \k in {0,...,7} {
    \node[font=\tiny] at (\k,0.38) {\k};
    \node[font=\tiny] at (-0.38,-\k) {\k};
  }
\end{tikzpicture}
\end{document}
```

An alignment score is the score of that path. The score “$+1$ for a match and $0$ otherwise” is the longest common subsequence problem.

### 6.2 Longest common subsequence

> [!def] **LCS problem** 
A subsequence keeps characters in order but need **not** keep them consecutive. `ACTGA` is a subsequence of `AGGCTCGA`; `TCGG` is not, because that order of letters does not occur.

Given strings $v$ and $w$, find a longest common subsequence.

Let $s_{i,j}$ be the length of an LCS of the $i$-prefix of $v$ and the $j$-prefix of $w$. Then $s_{i,0}=s_{0,j}=0$ and
$$
s_{i,j}=\max\begin{cases}s_{i-1,j}\\ s_{i,j-1}\\ s_{i-1,j-1}+1 & \text{if }v_i=w_j.\end{cases}
$$
The length of the LCS is $s_{n,m}$. Backtracking pointers $\uparrow$, $\leftarrow$, and $\nwarrow$ record which term achieved the maximum. `PrintLCS` follows $\nwarrow$ and emits $v_i$, or follows $\uparrow$ or $\leftarrow$ without emitting. The call `PrintLCS(b,v,n,m)` prints the subsequence. For $v=\mathtt{ATCTGAT}$ and $w=\mathtt{TGCATA}$ the LCS is `TCTA`.

One backtracking path for $v=\mathtt{ATCTGAT}$ and $w=\mathtt{TGCATA}$. Diagonal steps emit the matched letter; the path spells $\mathtt{TCTA}$, and $s_{7,6}=4$.

```tikz
\begin{document}
\begin{tikzpicture}[x=0.58cm,y=0.48cm,>=latex,font=\scriptsize]
  \foreach \x/\t in {1/T,2/G,3/C,4/A,5/T,6/A}
    \node at (\x,1.5) {\t};
  \foreach \x in {0,...,6}
    \node at (\x,0.78) {\x};
  \foreach \y/\lab/\idx/\a/\b/\c/\d/\e/\f/\g in {
    0//0/0/0/0/0/0/0/0,
    -1/A/1/0/0/0/0/1/1/1,
    -2/T/2/0/1/1/1/1/2/2,
    -3/C/3/0/1/1/2/2/2/2,
    -4/T/4/0/1/1/2/2/3/3,
    -5/G/5/0/1/2/2/2/3/3,
    -6/A/6/0/1/2/2/3/3/4,
    -7/T/7/0/1/2/2/3/4/4} {
    \node at (-1.15,\y) {\lab};
    \node at (-0.55,\y) {\idx};
    \node at (0,\y) {\a};
    \node at (1,\y) {\b};
    \node at (2,\y) {\c};
    \node at (3,\y) {\d};
    \node at (4,\y) {\e};
    \node at (5,\y) {\f};
    \node at (6,\y) {\g};
  }
  \node at (6,-7) {\textbf{4}};
  \draw (-1.55,1.88) -- (6.38,1.88);
  \draw (-1.55,0.38) -- (6.38,0.38);
  \draw (-1.55,-7.34) -- (6.38,-7.34);
  \draw (-0.28,1.88) -- (-0.28,-7.34);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (6,-7) -- (6,-6);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (6,-6) -- (5,-5);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (5,-5) -- (5,-4);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (5,-4) -- (4,-3);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (4,-3) -- (3,-3);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (3,-3) -- (2,-2);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (2,-2) -- (1,-2);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (1,-2) -- (0,-1);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (0,-1) -- (0,0);
\end{tikzpicture}
\end{document}
```

### 6.3 Edit distance by dynamic programming

The lecture’s edit-distance recurrence allows only insertions and deletions. Let $d_{i,j}$ be that distance between the $i$-prefix of $v$ and the $j$-prefix of $w$. Then $d_{i,0}=i$, $d_{0,j}=j$, and
$$
d_{i,j}=\min\begin{cases}d_{i-1,j}+1\\ d_{i,j-1}+1\\ d_{i-1,j-1} & \text{if }v_i=w_j.\end{cases}
$$
The distance is $d_{n,m}$. The same three backtracking arrows reconstruct the operations. For $v=\mathtt{ATCTGAT}$ and $w=\mathtt{TGCATA}$, one transformation deletes A, G, and T and inserts G and A.

One backtracking path for the same pair. It deletes $\mathtt{A}$, $\mathtt{G}$, and $\mathtt{T}$, and inserts $\mathtt{G}$ and $\mathtt{A}$. The distance is $d_{7,6}=5$.

```tikz
\begin{document}
\begin{tikzpicture}[x=0.58cm,y=0.48cm,>=latex,font=\scriptsize]
  \foreach \x/\t in {1/T,2/G,3/C,4/A,5/T,6/A}
    \node at (\x,1.5) {\t};
  \foreach \x in {0,...,6}
    \node at (\x,0.78) {\x};
  \foreach \y/\lab/\idx/\a/\b/\c/\d/\e/\f/\g in {
    0//0/0/1/2/3/4/5/6,
    -1/A/1/1/2/3/4/3/4/5,
    -2/T/2/2/1/2/3/4/3/4,
    -3/C/3/3/2/3/2/3/4/5,
    -4/T/4/4/3/4/3/4/3/4,
    -5/G/5/5/4/3/4/5/4/5,
    -6/A/6/6/5/4/5/4/5/4,
    -7/T/7/7/6/5/6/5/4/5} {
    \node at (-1.15,\y) {\lab};
    \node at (-0.55,\y) {\idx};
    \node at (0,\y) {\a};
    \node at (1,\y) {\b};
    \node at (2,\y) {\c};
    \node at (3,\y) {\d};
    \node at (4,\y) {\e};
    \node at (5,\y) {\f};
    \node at (6,\y) {\g};
  }
  \node at (6,-7) {\textbf{5}};
  \draw (-1.55,1.88) -- (6.38,1.88);
  \draw (-1.55,0.38) -- (6.38,0.38);
  \draw (-1.55,-7.34) -- (6.38,-7.34);
  \draw (-0.28,1.88) -- (-0.28,-7.34);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (6,-7) -- (6,-6);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (6,-6) -- (5,-5);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (5,-5) -- (5,-4);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (5,-4) -- (4,-3);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (4,-3) -- (3,-3);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (3,-3) -- (2,-2);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (2,-2) -- (1,-2);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (1,-2) -- (0,-1);
  \draw[->,shorten >=1.5mm,shorten <=1.5mm] (0,-1) -- (0,0);
\end{tikzpicture}
\end{document}
```

### 6.4 Global alignment

Extend the alphabet by the gap character `-`. For an alphabet of $k$ letters this gives a $(k+1)\times(k+1)$ scoring matrix $\delta$. Here $k=4$ for DNA and $k=20$ for proteins. The score of a column $\binom{x}{y}$ is $\delta(x,y)$, and the alignment score is the sum of its columns.

> [!def] **Global Alignment Problem** 
Given strings $v$, $w$, and scoring matrix $\delta$, find an alignment of maximum score.

The score $s_{i,j}$ of an optimal alignment of the prefixes satisfies
$$
s_{i,j}=\max\begin{cases}
s_{i-1,j}+\delta(v_i,-)\\
s_{i,j-1}+\delta(-,w_j)\\
s_{i-1,j-1}+\delta(v_i,w_j),
\end{cases}
$$
with $s_{0,0}=0$ and the borders equal to the total gap score of the corresponding prefix. The global score is $s_{n,m}$.

If a match scores $+1$, a mismatch costs $\mu$, and an indel costs $\sigma$,
$$
s_{i,j}=\max\begin{cases}
s_{i-1,j}-\sigma\\
s_{i,j-1}-\sigma\\
s_{i-1,j-1}-\mu & \text{if }v_i\neq w_j\\
s_{i-1,j-1}+1 & \text{if }v_i=w_j.
\end{cases}
$$

### 6.5 Local alignment

During the evolution of species, a large part of their genes would diverge. Still, it would be interesting to find segments of their genes that are can matchup, this can be formulated as the local alignment problem.

A substring alignment can score higher than the best alignment of the full strings, and that local similarity is often the biologically relevant part. Genes that regulate the same function in different species may differ everywhere except in one conserved region.

On the same pair, the global alignment is a path from one corner of the edit graph to the opposite corner. A local alignment is a path that can start and stop at any vertices, so it need not pay for the dissimilar flanks.

![[./attachments/Pasted image 20260928202520.png|521]]

> [!def] **Local Alignment Problem** 
Given $v$, $w$, and $\delta$, find substrings from $v$ and $w$ whose global alignment score is the best among all pairs of substrings.

Global alignment is a longest path from $(0,0)$ to $(n,m)$ in the edit graph. Local alignment is a longest path between some pair of vertices $(i,j)$ and $(i',j')$. Adding a zero-weight edge from the source $(0,0)$ to every vertex reduces this to longest paths from the source.

Zero-weight edges from the source $(0,0)$ to the other vertices, drawn dashed. Solid edges are the ordinary edit-graph steps. The recurrence is the global one, **with an extra option of $0$ that restarts the alignment**:
$$
s_{i,j}=\max\begin{cases}
0\\
s_{i-1,j}+\delta(v_i,-)\\
s_{i,j-1}+\delta(-,w_j)\\
s_{i-1,j-1}+\delta(v_i,w_j),
\end{cases}
$$
and $s_{i,0}=s_{0,j}=0$. The best local score is the largest $s_{i,j}$ in the table.

### 6.6 Matrices $\delta$ for sequence comparison
- For DNA, $\delta$ is fixed by the mismatch penalty $\mu$ and the indel penalty $\sigma$. 
- For proteins, $\delta(i,j)$ records how often amino acid $i$ replaces amino acid $j$ in alignments of related sequences. Building those alignments seems to require $\delta$ already. One either starts from very similar sequences or from an initial binary matrix. 
	- PAM matrices score closely related proteins. BLOSUM matrices score more divergent proteins.

## 7 Multiple sequence alignment

Comparing many sequences can show similarities that no pair shows on its own. An alignment of $\ell$ strings $v_1,\ldots,v_\ell$ of lengths $n_1,\ldots,n_\ell$ is an $\ell\times n$ matrix with $n\ge\max n_i$. Entries lie in the alphabet plus `-`. Row $i$ contains the characters of $v_i$ in order, with $n-n_i$ spaces. Every column contains at least one real symbol.

The score is the sum of column scores. The scoring function is an $\ell$-dimensional matrix $\delta$. The optimal alignment maximises that score.

**Multiple Sequence Alignment(MSA)** Given $\ell$ strings and $\delta$, find a maximum-score alignment.

Dynamic programming uses an $\ell$-dimensional edit graph. For three strings the cell $s_{i,j,k}$ is the best score of the corresponding prefixes, and the recurrence maximises over all $2^3-1=7$ ways of consuming a symbol or a gap from each string.

The cell $(i,j,k)$ has an incoming edge from each of the $2^3-1=7$ predecessors. Dashed edges are the rest of the cube.

- [ ] Add scoring matrix

```tikz
\begin{document}
\begin{tikzpicture}[>=latex,font=\scriptsize,
  dot/.style={circle,draw,fill=white,inner sep=0pt,minimum size=1.8mm}]
  \node[dot,label=above:{$(i{-}1,j{-}1,k{-}1)$}] (A) at (0.15,3.15) {};
  \node[dot,label=above:{$(i{-}1,j,k{-}1)$}] (B) at (3.15,3.15) {};
  \node[dot,label=left:{$(i{-}1,j{-}1,k)$}] (C) at (0,1.85) {};
  \node[dot,label=right:{$(i{-}1,j,k)$}] (D) at (3.0,1.85) {};
  \node[dot,label=left:{$(i,j{-}1,k{-}1)$}] (E) at (1.15,2.15) {};
  \node[dot,label={[xshift=2pt]above right:{$(i,j,k{-}1)$}}] (F) at (4.15,2.15) {};
  \node[dot,label=left:{$(i,j{-}1,k)$}] (G) at (1.0,0.7) {};
  \node[dot,label=below:{$(i,j,k)$}] (H) at (4.0,0.55) {};
  \draw[dashed] (A) -- (B);
  \draw[dashed] (A) -- (C);
  \draw[dashed] (A) -- (E);
  \draw[dashed] (B) -- (D);
  \draw[dashed] (B) -- (F);
  \draw[dashed] (C) -- (D);
  \draw[dashed] (C) -- (G);
  \draw[dashed] (E) -- (F);
  \draw[dashed] (E) -- (G);
  \draw[->] (C) -- (H);
  \draw[->] (D) -- (H);
  \draw[->] (G) -- (H);
  \draw[->] (A) -- (H);
  \draw[->] (B) -- (H);
  \draw[->] (E) -- (H);
  \draw[->] (F) -- (H);
\end{tikzpicture}
\end{document}
```
> Example of one step for $l=3$

With $\ell$ sequences of length about $n$, the running time is $O((2n)^\ell)$. Practical methods are *heuristics*. Some combine a compatible subset of optimal pairwise alignments into one multiple alignment.

# References

- Minary, Computational Biology, Lectures 3–4 (Michaelmas 2025)
- Jones & Pevzner, *An Introduction to Bioinformatics Algorithms*
