---
tags:
  - source-materials/lectures
  - physics/statistical-mechanics
  - physics/thermodynamics
prerequisites:
  - "[[limit|limit]]"
  - "[[derivative|derivative]]"
  - "[[taylor's series|taylor's series]]"
  - "[[partial derivative|partial derivative]]"
  - "[[multivariable chain rule|multivariable chain rule]]"
  - "[[nabla operator|nabla operator]]"
  - "[[integral|integral]]"
course date: Michaelmas 2025
lectuer: Prof. Peter Minary
share: true
---
# Lectures 8–10 Towards the statistical mechanics of biomolecules

Readable notes from Minary, *Computational Biology*, Lectures 8–10. Slides: [08-10-statistical-mechanics-biomolecules.pdf](08-10-statistical-mechanics-biomolecules.pdf).

These lectures are one pipeline. A biomolecule is many atoms, and an experiment reports a few numbers: a folded structure, a free energy, an average distance. Newton's laws move the atoms. Phase space and hte Hamiltonian turn that motion into one trajectory of known energy. 

Thermodynamics names the quantities an experiment actually controls, and an ensemble is how those quantities are averages over atomic arrangements. 

Molecular dynamics computes the average by a trajectory. Verlet and the propagator are how a finite time step can be taken without abandoning the mechanics, and the canonical ensemble is the version of the average that matches a laboratory held at fixed temperature. Chain models and folding landscapes are what the pipeline is for.

## 1 Newton's laws and the N-particle system

A simulation has to move the atoms. Newton's second law is that rule, written once for every atom. Newton's *Principia* (1687) states three laws for interacting bodies.

1. With no external force, a body stays at rest or moves in a straight line at constant velocity $v$.
2. An external force $F$ produces acceleration $a = F/m$, where $m$ is the mass.
3. If body $A$ exerts a force on body $B$, then $B$ exerts an equal and opposite force on $A$: $F_{BA} = -F_{AB}$.

The classical motion of $N$ particles in three dimensions is Newton's second law written for every particle. With no external forces, the force on particle $i$ is the sum of contributions from the other particles. If those forces are pairwise additive,
$$
F_i(r_1,\ldots,r_N) = \sum_{j\neq i} f_{ij}(r_i-r_j),
$$
and
$$
m_i\ddot r_i = F_i(r_1,\ldots,r_N), \qquad i = 1,\ldots,N.
$$
The momentum is $p_i = m_i v_i = m_i \dot r_i$, so the same law reads $F_i = \dot p_i$. Dots are time derivatives: $\dot r_i = dr_i/dt$ and $\ddot r_i = d^2 r_i/dt^2$.

## 2 Phase space

A protein has thousands of coordinates. Keeping them as separate functions is unusable. Phase space packs every position and every momentum into one point, so the molecule is one state and its motion is one curve.

Positions and momenta(plural, momentum) at one instant determine $\{r_1(t),\ldots,r_N(t)\}$ and $\{p_1(t),\ldots,p_N(t)\}$ at every later time. That is $2dN$ functions. In three dimensions $d = 3$, so there are $6N$ of them.

Those $2dN$ numbers are the microscopic state at time $t$. As an ordered tuple they are one point in a $2dN$-dimensional space, the **phase space**. In the usual Cartesian picture the space has $6N$ orthogonal axes, and a point is the phase-space vector
$$
x = (r_1,\ldots,r_N,p_1,\ldots,p_N).
$$
Solving Newton's equations draws a trajectory
$$
x_t = (r_1(t),\ldots,r_N(t),p_1(t),\ldots,p_N(t)).
$$
The motion of the whole system is the motion of that single point.
	 
A scalar indexing of the same vector is often easier:
$$
x = (q_1,\ldots,q_{3N},p_1,\ldots,p_{3N}),
$$
with $q_1 = r_{1,x}$, $q_2 = r_{1,y}$, and so on through $p_{3N} = p_{N,z}$. Along a trajectory each component is a function of time.

### 2.1 Lagrangian formulation

The forces in a molecular model come from a potential energy. The Lagrangian is the bookkeeping that turns that potential back into Newton's second law, and it is available once the force is conservative. A geometric description draws state as a point and motion as a curve. 

To have a description of the whole energy function at every point in the space, we need an algebraic description.

An algebraic description with position $r$ and velocity $\dot r$ (sitting in between the configuration space of position alone and position + momenta) is available for **conservative** forces. A conservative force is the gradient of a scalar potential energy $U(r_1,\ldots,r_N)$ (Also written as $V$),
$$
F_i(r_1,\ldots,r_N) = -\nabla_{r_i} U(r_1,\ldots,r_N).
$$
The work done by $F_i$ between two points is then **independent of the path**. This formulation sits closer to the phase-space picture than Newton's laws do.

The kinetic energy $K$ (Also written as $T$) of the $N$-particle velocities is
$$
K(\dot r_1,\ldots,\dot r_N) = \frac12\sum_{i=1}^N m_i \dot r_i^2.
$$
The **Lagrangian** is
$$
L = K(\dot r_1,\ldots,\dot r_N) - U(r_1,\ldots,r_N).
$$

It generates the equations of motion through the Euler–Lagrange equation
$$
\frac{d}{dt}\left(\frac{\partial L}{\partial \dot r_i}\right) - \frac{\partial L}{\partial r_i} = 0.
$$
Substituting $L$ recovers Newton's second law. We have $\partial L/\partial\dot r_i = m_i\dot r_i$, so the time derivative $\frac{d}{dt}\left(\frac{\partial L}{\partial \dot r_i}\right)$ is $m_i\ddot r_i$, while $\partial L/\partial r_i = -\partial U/\partial r_i = F_i$. The Euler–Lagrange equation is therefore $m_i\ddot r_i - F_i = 0$.
### 2.2 Hamiltonian and the Poisson bracket

With the Lagrangian, we now have everything to derive a formulation in the phase space.

Molecular dynamics and the propagator step the energy written in positions and momenta. That function is the Hamiltonian. The Poisson bracket is the same idea for any observable: it is the rate at which that observable changes along the trajectory.

The **Hamiltonian** is the total energy as a function of positions and momenta. It is related to the Lagrangian by
$$
\begin{align}
H(r_1,\ldots,r_N,p_1,\ldots,p_N)
&= \sum_{i=1}^N p_i\cdot\dot r_i(p_i)
- L(r_1,\ldots,r_N,\dot r_1(p_1),\ldots,\dot r_N(p_N))\\
&= \sum_{i=1}^N \frac{p_i^2}{2m_i} + U(r_1,\ldots,r_N)\\ 
&= K(p_{1}, \ldots, p_{N}) + U(r_1,\ldots,r_N)
\end{align}
$$
> $\dot{r}_{i}(p_{i})$ means $\dot{r_{i}}$ rewritten as a function of $p_{i}$ using $p_{i}=\frac{ \partial L }{ \partial \dot{r_{i}} }=m_{i}\dot{r_{i}}$, which gives $\frac{p_{i}}{m_{i}}$

A standard example is the **harmonic oscillator**: a mass $m$ pulled back toward $x=0$ by a spring of stiffness $k$. The force is $F=-kx$, so the pull grows in proportion to the displacement and always points back toward the rest point. The stored potential is $\tfrac12 kx^2$, and
$$
H = \frac{p^2}{2m} + \frac12 kx^2.
$$
Released from rest, the mass swings through $0$ and back at one angular frequency $\omega=\sqrt{k/m}$, which is the same statement as $k=m\omega^2$. The motion is the single cosine $x(t)=A\cos(\omega t+\phi)$, and the period is $2\pi/\omega$. "Harmonic" names that pure single-frequency swing.

Hamilton's equations are
$$
\dot r_i = \frac{\partial H}{\partial p_i} = \frac{p_i}{m_i},
\qquad
\dot p_i = -\frac{\partial H}{\partial r_i} = -\frac{\partial U}{\partial r_i} = F_i(r_1,\ldots,r_N).
$$
Newton's second law follows again: $m_i\ddot r_i = \dot p_i = F_i$.

These equations conserve $H$. In the scalar coordinates $q_i$,
$$
\frac{dH}{dt}
= \sum_{i=1}^{3N}\left(
\frac{\partial H}{\partial q_i}\dot q_i + \frac{\partial H}{\partial p_i}\dot p_i
\right)
= \sum_{i=1}^{3N}\left(
\frac{\partial H}{\partial q_i}\frac{\partial H}{\partial p_i}
- \frac{\partial H}{\partial p_i}\frac{\partial H}{\partial q_i}
\right)
= 0,
$$
after substituting $\dot q_i = \partial H/\partial p_i$ and $\dot p_i = -\partial H/\partial q_i$. The two products in each term cancel.

Let $a(x)$ be any function of the phase-space vector. Along the trajectory $x_t$,
$$
\frac{da}{dt}
= \sum_{i=1}^{3N}\left(
\frac{\partial a}{\partial q_i}\frac{\partial H}{\partial p_i}
- \frac{\partial a}{\partial p_i}\frac{\partial H}{\partial q_i}
\right)
\equiv \{a,H\}.
$$
The right-hand side is the **Poisson bracket** of $a$ with $H$. For two phase-space functions $a$ and $b$,
$$
\{a,b\}
= \sum_{i=1}^{3N}\left(
\frac{\partial a}{\partial q_i}\frac{\partial b}{\partial p_i}
- \frac{\partial a}{\partial p_i}\frac{\partial b}{\partial q_i}
\right).
$$
This measures that rate of changes of $a$ when the motion (velocity) is generated by $b$
## 3 The three laws of thermodynamics

An experiment does not set the coordinates of every atom. It sets pressure, volume, temperature, and the amount of material. Thermodynamics is the description in those variables. Entropy is introduced so that heat, which depends on the path taken, still determines the change in a state function. The laws do not refer to the microscopic constituents. The system is in **thermodynamic equilibrium** when its thermodynamic state does not change in time.

The parameters that define a thermodynamic state are **pressure $P$, volume $V$, temperature $T$**, and **the number of moles $n$**. One mole contains $6.02214076\times 10^{23}$ particles. 

An **equation of state** is a relation $g(n,P,V,T) = 0$ among those parameters, describing how they change from one equilibrium state to another. 

> For an ideal gas the relation is $PV - nRT = 0$. $R$ is the ideal gas constant, about $8.314\,\mathrm{J\,mol^{-1}\,K^{-1}}$.

Counted per particle the same law is $PV=NkT$. One mole is $N_A$ particles, so the constant per mole and the constant per particle are the same number in different units: $R=N_A k$. The Boltzmann constant $k$ is given in section 4.2.

A **thermodynamic transformation** changes the thermodynamic state. It may be carried out reversibly or irreversibly. In a *reversible* transformation the outside conditions change more slowly than the system can catch up, so the system stays in equilibrium at every step and can be led back through the same states. An *irreversible* transformation does not allow that adjustment.

A piston makes the difference concrete. Raise the external pressure in tiny steps and wait after each one: the gas pressure equals the external pressure throughout, and lowering the pressure through the same steps returns the gas along the same path. Jump the external pressure up in one step and the gas next to the piston is squeezed before the rest of the gas has moved. During that push the cylinder has no single pressure, so there is no equilibrium path to retrace.

A **state function** has the form $f(n,P,V,T)$. Its change in any transformation depends *only on the initial and final states*. Internal energy $E$ is one example. Entropy, written $S$ from here on, is another. The letter $S$ is the entropy; "state function" is the property that its change ignores the path.

In a reversible transformation, if the volume changes by $dV$ and the number of moles by $dn$, the work performed on the system is
$$
dW_{\mathrm{rev}} = -P\,dV + \mu\,dn,
$$
where $\mu$ is the **chemical potential**. 
 
*Changing the volume or the number of moles therefore requires work*. 

The heat required for a *reversible* temperature change $dT$ is
$$
dQ_{\mathrm{rev}} = C\,dT,
$$
where $C$ is the **heat capacity**: the heat required to change the temperature of $1.0$ mole of a substance by $1\,\mathrm{K}$. A rise of $1\,\mathrm{K}$ is the same temperature interval as a rise of $1^\circ\mathrm{C}$.

### 3.1 **First law**
In any thermodynamic transformation, if the system absorbs heat $\Delta Q$ and work $\Delta W$ is performed on it, the internal energy changes by
$$
\Delta E = \Delta Q + \Delta W.
$$
$E$ is a state function, so a reversible path and an irreversible path between the same end states give the same $\Delta E$.

### 3.2 **Second law** 
The entropy of the universe, meaning the system and its surroundings, must increase or stay the same, so $dS_{\mathrm{universe}} \ge 0$. Then $dS_{\mathrm{universe}} > 0$ for an irreversible process and $dS_{\mathrm{universe}} = 0$ for a reversible process. This $dS_{\mathrm{universe}}$ is not the entropy change of the system alone.

The entropy of the system is defined along a reversible path by
$$
dS = \frac{dQ_{\mathrm{rev}}}{T}.
$$
It is an exact differential, so the entropy $S$ is a state function of the system. The change between two states is this integral along any reversible path connecting them. If the process that actually occurred was irreversible, $\Delta S$ is still that number: the heat absorbed on the irreversible path is smaller than $T\,dS$ and is not what goes into the formula. When the system does absorb heat $dQ_{\mathrm{rev}}$, its entropy rises by $dQ_{\mathrm{rev}}/T$ and the surroundings lose that same heat, so their entropy falls by $dQ_{\mathrm{rev}}/T$. The two changes cancel and $dS_{\mathrm{universe}} = 0$. A reversible heating therefore changes the system's entropy while leaving the entropy of the universe fixed.


> [!faq] Exact differential and state function
> An exact differential is the small change in a function of the state alone. Integrated from state $A$ to state $B$ it equals the function at $B$ minus the function at $A$, whichever path is taken. Altitude works this way: every route from the foot of a hill to the summit changes height by the same amount. 
> 
> The heat absorbed does not. For one mole of ideal gas, $dQ_{\mathrm{rev}} = C\,dT + RT\,dV/V$. Heating from $T_1$ to $T_2$ and then expanding from $V_1$ to $V_2$ absorbs $C(T_2-T_1) + RT_2\ln(V_2/V_1)$. Expanding first and then heating absorbs $RT_1\ln(V_2/V_1) + C(T_2-T_1)$. The expansion contributes heat proportional to the temperature at which it is carried out, so the two totals differ. Dividing each contribution by that temperature cancels the extra factor of $T$. Both routes then give
> $$
> \int\frac{dQ_{\mathrm{rev}}}{T} = C\ln\frac{T_2}{T_1} + R\ln\frac{V_2}{V_1},
> $$
> which depends only on the end states. That path independence is what makes $dQ_{\mathrm{rev}}/T$ the differential of the entropy $S$.

### 3.3 **Third law** 
The entropy of the system at absolute zero is a universal constant, which can be taken to be zero.

Absolute zero is $T = 0$ on the Kelvin scale, a temperature that can never be physically reached. 

## 4 The ensemble concept

> One temperature, or one energy, is shared by enormously many atomic arrangements. A laboratory number is an average over them. An ensemble is that collection, and it is what makes "the probability of a conformation" a definite object. *Temperature inferred from the average kinetic energy is an example: many assignments of the particle velocities share the same total kinetic energy.*

An **ensemble** is a collection of systems with the **same microscopic *interactions* and the same macroscopic *properties***, for example the same total energy $E$, volume $V$, and particle number $N$.

### 4.1 The microcanonical ensemble

The simplest ensemble is an isolated molecule, with the energy fixed. The trajectory then stays on one surface, which is the setting in which the Dirac delta and Liouville's theorem are stated. A laboratory sample is usually held at fixed temperature instead; that case is section 9. The isolated case is the starting point because energy conservation is already in Hamilton's equations. Consider $N$ identical particles in a container of volume $V$ with fixed internal energy $E$. The triple $(E,V,N)$ is a set of macroscopic control variables, and it fixes the thermodynamic properties. The state function of an ensemble controlled by $E$, $V$, and $N$ is a function of those three variables.

Take a reversible transformation that adds heat $dQ$ and does work $dW$ on the system. The first law gives $dE = dQ_{\mathrm{rev}} + dW_{\mathrm{rev}}$. With $dS = dQ_{\mathrm{rev}}/T$, the heat is $dQ_{\mathrm{rev}} = T\,dS$.

A small volume change contributes *mechanical* work $-P\,dV$, and a small change $dN$ in the particle number contributes *chemical* work $\mu\,dN$. 

Particle number changes by integers, but the changes considered here are tiny compared with $N \approx 10^{23}$, so *$N$ is treated as continuous*. Therefore
$$
\color{red}{dE = T\,dS - P\,dV + \mu\,dN}
$$
Rearrangement gives
$$
dS = \frac{1}{T}\,dE + \frac{P}{T}\,dV - \frac{\mu}{T}\,dN.
$$
The state function this identifies is the entropy $S = S(N,V,E)$.

The same differential, written from the partial derivatives of $S$, is
$$
dS
= \left(\frac{\partial S}{\partial E}\right)_{N,V} dE
+ \left(\frac{\partial S}{\partial V}\right)_{N,E} dV
+ \left(\frac{\partial S}{\partial N}\right)_{E,V} dN.
$$
The subscript on a partial derivative names the variables held fixed. Since $S=S(N,V,E)$, each derivative changes one argument and keeps the other two constant. $\bigl(\partial S/\partial E\bigr)_{N,V}$ is the slope along $E$ at constant particle number and volume. The same derivative taken at constant temperature would be a different number, which is why the subscript is written out.

Matching coefficients,
$$
\frac{1}{T} = \left(\frac{\partial S}{\partial E}\right)_{N,V},
\quad
\frac{P}{T} = \left(\frac{\partial S}{\partial V}\right)_{N,E},
\quad
\frac{\mu}{T} = -\left(\frac{\partial S}{\partial N}\right)_{E,V}.
$$

### 4.2 Boltzmann's entropy

Boltzmann (1877) connected the entropy to the number of microscopic states $\Omega(N,V,E)$. Planck put the relation in its present form around 1900:
$$
S(N,V,E) = k\ln\Omega(N,V,E),
$$
where $k$ is the Boltzmann constant,
$$
k = 1.380650524\times 10^{-23}\,\mathrm{J\cdot K^{-1}}.
$$
If $\Omega$ can be computed from the microscopic description, the macroscopic observables follow.

The microcanonical ensemble is a collection of isolated systems obeying Hamilton's equations. **Those equations conserve the Hamiltonian, so $H(x) = E$**. The constraint **removes one degree of freedom**: an $N$-particle system in three dimensions has $6N-1$ independent phase-space degrees of freedom, or $2dN-1$ in $d$ dimensions.

Each allowed microstate is a phase-space point with $H(x) = E$. The ensemble distribution function that enforces this is
$$
F\bigl(H(x)\bigr) = \mathcal{N}\,\delta\bigl(H(x)-E\bigr),
$$
> where $\mathcal{N}$ is a normalisation constant (the slides denote it by $N$) and $\delta$ is the Dirac delta. The delta vanishes for a nonzero argument, integrates to $1$, and is the limit of the rectangular spikes $\delta_n$ that equal $n/2$ on $[-1/n,1/n]$ and equal $0$ outside. **In $\delta\bigl(H(x)-E\bigr)$ the argument is zero only on the surface $H(x)=E$, so a phase-space point of any other energy contributes nothing**. The average written below therefore adds $a(x)$ only on that surface; the denominator is the surface's area, up to the constant in $\Omega$.

Every phase-space point of the ensemble lies on the *constant-energy hypersurface*, and every point on that surface is equally probable. Points off the surface have probability zero. The quantity $1/\Omega(N,V,E)$ measures the probability of drawing a microstate at random on the constant-energy hypersurface.

$\Omega(N,V,E)$ measures the phase space available to the system. It is proportional to the $(6N-1)$-dimensional area of the constant-energy hypersurface. Precisely,
$$
\Omega(N,V,E) = M\int dx\,\delta\bigl(H(x)-E\bigr),
$$
with a constant factor $M$. 
> $M$ come from the consideration that  if all particles are equivalent, we will be overcounting the states by $n!$. 

The same integral, written over momenta and positions, is
$$
\Omega(N,V,E)
= M\int d^N p\int_{D(V)} d^N r\,\delta\bigl(H(r,p)-E\bigr).
$$
Each momentum is integrated over the whole real line. $D(V)$ is the container. For a cube of side $L$ in the positive octant with a corner at the origin, each Cartesian coordinate of each position lies in $[0,L]$. The notation $d^N r$ and $d^N p$ stands for the product of the $N$ vector measures. When an observable is averaged, $M$ cancels.

Thermodynamic quantities at given $N$, $V$, and $E$ are then derivatives of $\ln\Omega$:
$$
\frac{1}{kT} = \left(\frac{\partial\ln\Omega}{\partial E}\right)_{N,V},
\quad
\frac{P}{kT} = \left(\frac{\partial\ln\Omega}{\partial V}\right)_{N,E},
\quad
\frac{\mu}{kT} = -\left(\frac{\partial\ln\Omega}{\partial N}\right)_{E,V}.
$$

The ensemble average of an observable $A$ is the average of its phase-space function $a(x)$ on the energy surface:
$$
\begin{aligned}
\langle a\rangle
&= \frac{M}{\Omega(N,V,E)}\int dx\,a(x)\,\delta\bigl(H(x)-E\bigr)  \texttt{substitute }\Omega\\
&= \frac{\int dx\,a(x)\,\delta\bigl(H(x)-E\bigr)}{\int dx\,\delta\bigl(H(x)-E\bigr)}\\
\end{aligned}
$$
However, this is not realistic to solve this integral in practice.

The lecture uses this formula to recover $S = k\ln\Omega$. The observable $1/T$ must itself be such an average, and it must also equal $(\partial S/\partial E)_{N,V}$. If one assumes only that $S = CG(\Omega)$ for some constant $C$ and some function $G$, matching the two expressions for $(\partial S/\partial E)_{N,V}$ forces $CG(\Omega) = k\ln\Omega$. The steps identify the phase-space function whose average is $1/T$ as $(1/k)\,\partial[\ln\delta(H(x)-E)]/\partial E$. The conclusion is the Boltzmann formula above. The integrals that define $\langle a\rangle$ cannot be done analytically in most cases, so the averages are computed numerically.

- [ ] Update this section with up to date slides 🔺 

## 5 Molecular dynamics

The phase-space integrals that define an ensemble average cannot be evaluated by hand for a protein. Molecular dynamics replaces the integral by a trajectory: **step Hamilton's equations, and average the observable along the steps**. The output is a trajectory: an ordered set of microscopic configurations, the evolution of the system in time. Many numerical experiments can be run by changing the control parameters. The results are only as good as the numerical model, and a run of fixed length may fail to visit enough microstates.

In the microcanonical ensemble the system follows
$$
\dot q_i = \frac{\partial H}{\partial p_i},
\qquad
\dot p_i = -\frac{\partial H}{\partial q_i}.
$$
These equations conserve energy, $dH/dt = 0$, so $H$ stays equal to the constant $E$. The configurations on the trajectory belong to the microcanonical ensemble of energy $E$.

### 5.1 Ergodicity and time averages

A trajectory is one path. It can stand in for the ensemble average only if a long enough run visits the energy surface in the same proportions the ensemble would. That assumption is ergodicity, and it is the reason a time average is used as the measurement. The system is then called **ergodic**, and the trajectory is an ergodic trajectory. ***Microcanonical phase-space averages may be replaced by time averages along that trajectory***:
$$
\frac{\int dx\,a(x)\,\delta\bigl(H(x)-E\bigr)}{\int dx\,\delta\bigl(H(x)-E\bigr)}
= \lim_{\tau\to\infty}\frac{1}{\tau}\int_0^\tau dt\,a(x_t)
\equiv \overline{a}.
$$

A numerical trajectory starts from an initial phase-space vector $x_0$ and a time step $\Delta t$, and it produces vectors $x_{n\Delta t}$ for $n = 1,\ldots,M$. The ensemble average is estimated by the discrete time average
$$
\langle a\rangle = \frac{1}{M}\sum_{n=1}^M a(x_{n\Delta t}).
$$
This equality holds when the molecular-dynamics trajectory is ergodic.

### 5.2 Numerical Integrators
#### 5.2.1 Verlet and velocity Verlet

A computer cannot take the limit of a vanishing time step. Verlet is a finite step that remains time-reversible, as Hamilton's equations are, so the numerical trajectory stays close to a surface of constant energy. The simplest integrators come from the Taylor expansion about $t$ 
$$
r_i(t+\Delta t)
= r_i(t) + \Delta t\,v_i(t) + \frac{\Delta t^2}{2}\frac{F_i(t)}{m_i} + O(\Delta t^3),
$$
where Newton's law has replaced the acceleration. The expansion about $t-\Delta t$ is the same series with the signs of the odd powers of $\Delta t$ reversed. Adding the two expansions cancels the velocities:
$$
r_i(t+\Delta t) + r_i(t-\Delta t)
= 2r_i(t) + \Delta t^2\frac{F_i(t)}{m_i} + O(\Delta t^3).
$$
The lecture then drops the explicit $O(\Delta t^3)$ term and notes that the computed positions still carry an error of order $\Delta t^3$.

Solving for the forward position gives the **Verlet algorithm** (1967),
$$
r_i(t+\Delta t) = 2r_i(t) - r_i(t-\Delta t) + \Delta t^2\frac{F_i(t)}{m_i}.
$$
Verlet advances positions only. Velocities along the trajectory are recovered by
$$
v_i(t) = \frac{r_i(t+\Delta t) - r_i(t-\Delta t)}{2\Delta t}.
$$

**Velocity Verlet** advances positions and velocities together:
$$
r_i(t+\Delta t)
= r_i(t) + \Delta t\,v_i(t) + \frac{\Delta t^2}{2}\frac{F_i(t)}{m_i},
$$
$$
v_i(t+\Delta t)
= v_i(t) + \frac{\Delta t}{2m_i}\bigl[F_i(t) + F_i(t+\Delta t)\bigr].
$$
The velocity update is obtained by writing the backward Taylor expansion about $t+\Delta t$, substituting the forward expansion of $r_i(t+\Delta t)$, and solving for $v_i(t+\Delta t)$.

Both algorithms are time-reversible, which is a symmetry of Hamilton's equations. If the positions and velocities at $t+\Delta t$ are taken as a new initial state and the integrator is run with step $-\Delta t$, it returns the state at time $t$.

#### 5.2.2 The classical propagator and Trotter factorization

> The velocity Verlet formulas of §5.2 are one step of the exact evolution, written so that a computer can carry it out. Hamilton's equations move a phase-space point by the classical propagator $e^{iLt}$. The Liouville operator inside that exponential splits into a position piece and a momentum piece, and the two pieces do not commute, so the exponential of the sum is not the product of the exponentials. 
> 
> The symmetric Trotter factorization replaces one exact step by three shifts that are exact on their own: half a step of the momentum, a full step of the position, and another half-step of the momentum. Applying those three shifts recovers the velocity Verlet updates, and the remainder $O(\Delta t^3)$ is the error in the factorization(second-order accurate). Placing the position piece on the outside of the same pattern gives the drift–kick–drift step, and nesting the pattern for fast and slow forces gives the multiple-time-step integrator.

The construction starts from the time derivative of a phase-space function. With Hamilton's equations, $da/dt = \{a,H\}$.

Define an operator $iL$ by $iLa = \{a,H\}$, so $da/dt = iLa$. Abstractly, $iL = \{\ldots,H\}$, with the operand written in place of the dots. As a differential operator,
$$
iL = \sum_{i=1}^{3N}\left(
\frac{\partial H}{\partial p_i}\frac{\partial}{\partial q_i}
- \frac{\partial H}{\partial q_i}\frac{\partial}{\partial p_i}
\right).
$$
The symbol $iL$ is the conventional name of the Liouville operator. The factor $i$ is there by analogy with quantum mechanics and plays no mathematical role in the classical formula.

Solving $da/dt = iLa$ gives
$$
a(x_t) = e^{iLt}a(x_0).
$$
The operator $e^{iLt}$ is the **classical propagator**. Setting $a(x) = x$ solves Hamilton's equations in the form $x_t = e^{iLt}x_0$. The factor $i$ makes this expression resemble the quantum propagator $e^{-i\hat H t/\hbar}$. The action of $e^{iLt}$ on $x_0$ cannot be evaluated exactly.

Split $iL = iL_1 + iL_2$, where $iL_1$ contains the derivatives with respect to the coordinates and $iL_2$ contains the derivatives with respect to the momenta. The two pieces do not commute, $[iL_1,iL_2] \neq 0$, so $e^{iLt}$ is not the product $e^{iL_1 t}e^{iL_2 t}$. That product would be useful, because each factor separately ca n often be applied exactly.

> $[A, B] = AB - BA$, $[A,B]\neq0$ means order of the two operator matters. 

The **Trotter theorem**, also called the symmetric Trotter theorem or **Strang splitting**, says that if $[A,B] \neq 0$ and $P$ is an integer, then
$$
e^{A+B} = \lim_{P\to\infty}\Bigl(e^{B/(2P)}e^{A/P}e^{B/(2P)}\Bigr)^P
$$
> Note that the addition still commutes, so changing the order of splitting by moving $A$ to the outside also works


The failure of $e^{A+B}=e^{A}e^{B}$ is visible in the Taylor series of the exponential.

`\begin{proof}`
For a number $x$,
$$
e^{x}=\sum_{n=0}^{\infty}\frac{x^{n}}{n!}=1+x+\frac{x^{2}}{2!}+\frac{x^{3}}{3!}+\cdots,
$$
and the infinite sum equals $e^{x}$ exactly. A matrix or an operator has sums and products, so the same series defines the exponential, with $1$ replaced by the identity:
$$
e^{X}=I+X+\frac{X^{2}}{2!}+\frac{X^{3}}{3!}+\cdots.
$$
Through second order the two sides of the disputed identity are
$$
\begin{align}
e^{A+B}
&= I+(A+B)+\frac{(A+B)^{2}}{2}+\cdots \\
&= I+A+B+\frac{A^{2}+AB+BA+B^{2}}{2}+\cdots
\end{align}
$$
and
$$
\begin{align}
e^{A}e^{B}
&=\Bigl(I+A+\frac{A^{2}}{2}+\cdots\Bigr)\Bigl(I+B+\frac{B^{2}}{2}+\cdots\Bigr) \\
&= I\cdot I + I\cdot B + I\cdot\frac{B^{2}}{2} \\
&\quad + A\cdot I + A\cdot B + A\cdot\frac{B^{2}}{2} \\
&\quad + \frac{A^{2}}{2}\cdot I + \frac{A^{2}}{2}\cdot B + \frac{A^{2}}{2}\cdot\frac{B^{2}}{2}+\cdots.
\end{align}
$$
$A\cdot B^{2}/2$, $A^{2}B/2$, and $A^{2}B^{2}/4$ each contain three or four factors, so they sit in the omitted higher-order terms. What remains through second order is
$$
e^{A}e^{B}= I+A+B+\frac{A^{2}}{2}+AB+\frac{B^{2}}{2}+\cdots.
$$
$A^{2}/2$ comes only from the left factor and $B^{2}/2$ only from the right factor. The cross term is one $A$ from the left and one $B$ from the right, so its coefficient is $1$ and its order is $AB$. The other order $BA$ does not appear: the left factor supplies no $B$, and the right factor supplies no $A$.

The linear terms agree. In the quadratic terms, $(A+B)^{2}$ produces both orders, so $e^{A+B}$ contains the average $(AB+BA)/2$. The product $e^{A}e^{B}$ contains only $AB$, because every factor taken from $e^{A}$ stands to the left of every factor taken from $e^{B}$. Subtracting the two series leaves
$$
e^{A}e^{B}-e^{A+B}=\frac{AB-BA}{2}+\cdots=\frac{[A,B]}{2}+\cdots.
$$
The two sides differ whenever $[A,B]\neq 0$. If $[A,B]=0$, then $AB=BA$, the extra terms vanish, and the series agree. Ordinary numbers commute, which is why $e^{a+b}=e^{a}e^{b}$ holds for them.

A matrix pair with $[A,B]\neq 0$ makes the same gap explicit. Take
$$
A=\begin{pmatrix}0&1\\0&0\end{pmatrix},\qquad
B=\begin{pmatrix}0&0\\1&0\end{pmatrix}.
$$
Then $A^{2}=B^{2}=0$, so $e^{A}=I+A$ and $e^{B}=I+B$, while $(A+B)^{2}=I$, so $e^{A+B}=\cosh(1)\,I+\sinh(1)\,(A+B)$. The two results are
$$
e^{A}e^{B}=\begin{pmatrix}2&1\\1&1\end{pmatrix},\qquad
e^{A+B}\approx\begin{pmatrix}1.543&1.175\\1.175&1.543\end{pmatrix}.
$$
`\end{proof}`

Applied to the propagator, with time step $\Delta t = t/P$,
$$
e^{iLt}
= \lim_{P\to\infty,\,\Delta t\to 0}
\Bigl(e^{iL_2\Delta t/2}e^{iL_1\Delta t}e^{iL_2\Delta t/2}\Bigr)^P.
$$
At finite $P$ this is an approximation,
$$
e^{iLt}
= \Bigl(e^{iL_2\Delta t/2}e^{iL_1\Delta t}e^{iL_2\Delta t/2}\Bigr)^P + O(P\Delta t^3),
$$
so one step of the factorized propagator is
$$
e^{iL\Delta t} = e^{iL_2\Delta t/2}e^{iL_1\Delta t}e^{iL_2\Delta t/2} + O(\Delta t^3).
$$

For one particle in one dimension, $H = p^2/(2m) + U(x)$, the pieces are
$$
iL_1 = \frac{p}{m}\frac{\partial}{\partial x},
\qquad
iL_2 = F(x)\frac{\partial}{\partial p}.
$$
The operator identity $e^{c\partial/\partial q}f(q) = f(q+c)$ turns each factor into a pure shift of one coordinate, holding the other fixed. The factors act from the right, so the rightmost exponential is applied first. In a momentum shift the force is computed from whatever position the preceding factors have already produced. Applying the three factors in that order, and writing $v = p/m$, *recovers velocity Verlet*:
$$
x(\Delta t) = x(0) + \Delta t\,v(0) + \frac{\Delta t^2}{2m}F\bigl(x(0)\bigr),
$$
$$
v(\Delta t) = v(0) + \frac{\Delta t}{2m}\Bigl[F\bigl(x(0)\bigr) + F\bigl(x(\Delta t)\bigr)\Bigr].
$$

We have seemingly achieved nothing, but this derivation gives us a very powerful ***three-step procedure***:
- A half-step of momentum;
- A full step of position;
- and another half-step of momentum:
$$
p(\Delta t/2) = p(0) + \frac{\Delta t}{2}F\bigl(x(0)\bigr),
$$
$$
x(\Delta t) = x(0) + \frac{\Delta t}{m}p(\Delta t/2),
$$
$$
p(\Delta t) = p(\Delta t/2) + \frac{\Delta t}{2}F\bigl(x(\Delta t)\bigr).
$$
Each operator is one instruction. The first shifts the momentum. The second shifts the position, after which the force must be recomputed. The third shifts the momentum with the new force. One integration step is therefore

```text
p ← p + (Δt/2) F
x ← x + Δt (p/m)
F ← F(x)
p ← p + (Δt/2) F
```

The translation avoids deriving a finite-difference formula by hand, and the time-reversible step preserves symmetries of classical mechanics.

For $N$ particles in three dimensions,
$$
H = \sum_{i=1}^N \frac{p_i^2}{2m_i} + U(r_1,\ldots,r_N),
$$
and
$$
iL = \sum_{i=1}^N \frac{p_i}{m_i}\cdot\frac{\partial}{\partial r_i} + \sum_{i=1}^N F_i\cdot\frac{\partial}{\partial p_i}.
$$
All terms inside the first sum commute with each other, and all terms inside the second sum commute with each other, so the same Trotter factorization applies. The code updates every Cartesian component of every momentum by a half-step, then every component of every position, recomputes the force from the new positions, and finishes the momentum half-step.

#### 5.2.3 Multiple time scales

A timestep small enough to follow a bond vibration cannot reach a folding event: the run would be almost all vibration. Motions in a biomolecule are separated by orders of magnitude, so the slow forces can be recomputed much less often than the fast ones. That separation is what makes a long trajectory affordable.

![[./attachments/Pasted image 20260924201442.png|450]]

From fast to slow the slide shows bond vibrations (fs), side-chain rotations (ps), backbone fluctuations (ns), loop motion (ns–ms), ligand binding and unbinding (μs), and domain rearrangements (μs–ms). The units are $1\,\mathrm{fs} = 10^{-15}\,\mathrm{s}$, $1\,\mathrm{ps} = 10^{-12}\,\mathrm{s}$, $1\,\mathrm{ns} = 10^{-9}\,\mathrm{s}$, $1\,\mu\mathrm{s} = 10^{-6}\,\mathrm{s}$, and $1\,\mathrm{ms} = 10^{-3}\,\mathrm{s}$.

A standard biomolecular potential separates bonded terms from nonbonded terms:
$$
\begin{aligned}
U(r_1,\ldots,r_N)
&= \sum_{\mathrm{bonds}} K_r(r-r_0)^2
+ \sum_{\mathrm{bends}} K_\theta(\theta-\theta_0)^2
+ \sum_{\mathrm{dihedrals}} K_\phi\bigl[1+\cos(n\phi+\delta)\bigr] \\
&\quad + \sum_{i,j\in\mathrm{nb}}\left(
4\varepsilon_{ij}\left[
\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12}
- \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}
\right]
+ k_e\frac{q_i q_j}{r_{ij}}
\right),
\end{aligned}
$$
with $\sigma_{ij} = (\sigma_i+\sigma_j)/2$.

![[./attachments/Pasted image 20260924201525.png|Pasted image 20260924201525.png]]

The first sum is over covalent pairs and penalises departure from the equilibrium bond length $r_0$. The second is over bond angles of neighbouring covalent bonds and penalises departure from the equilibrium angle $\theta_0$. The third is the conformational energy of each dihedral (torsion) angle $\phi$; the parameters $n$ and $\delta$ are constants. Those three sums are the bonded energy.

The nonbonded sum contains *van der Waals* interactions, written as a Lennard-Jones potential between spheres of radii $\sigma_i$ and $\sigma_j$ at distance $r_{ij}$ with well depth $\varepsilon_{ij}$, and Coulomb interactions between charges $q_i$ and $q_j$. The factor $k_e$ is a constant.

The force on particle $i$ is $F_i = -\partial U/\partial r_i$. Bonded forces are large and vary rapidly, and they are cheap to evaluate. *Nonbonded forces vary slowly, but a direct sum over pairs costs $O(N^2)$ work*. On the timescale of the fast forces, the slow forces barely change. A single velocity-Verlet step $\Delta t$ is limited by the fast forces. If bond vibrations require $\Delta t = 1\,\mathrm{fs}$, a trajectory of $1\,\mathrm{ms}$ takes $10^{12}$ steps, and both the fast and the slow forces are evaluated $10^{12}$ times. An integrator that uses the separation of timescales avoids that cost.

The construction is shown for one particle in one dimension whose potential is a slow part plus a fast part, so
$$
F(x) = F_{\mathrm{slow}}(x) + F_{\mathrm{fast}}(x),
$$
with $F_{\mathrm{slow}} = -dU_{\mathrm{slow}}/dx$ and $F_{\mathrm{fast}} = -dU_{\mathrm{fast}}/dx$. The Liouville operator is streaming plus the total force. Regrouping gives
$$
iL_{\mathrm{fast}} = \frac{p}{m}\frac{\partial}{\partial x} + F_{\mathrm{fast}}(x)\frac{\partial}{\partial p},
\qquad
iL_{\mathrm{slow}} = F_{\mathrm{slow}}(x)\frac{\partial}{\partial p}.
$$
A step $\Delta t$ appropriate to the slow force is factorized by the Trotter theorem as a slow half-step, a full fast interval, and a slow half-step. The fast interval is itself split into $n$ steps of size $\delta t = \Delta t/n$. The resulting operator is a slow momentum half-step, then $n$ copies of a fast velocity-Verlet step, then a slow momentum half-step. In code:

```text
p ← p + (Δt/2) F_slow
for i = 1, …, n do
  p ← p + (δt/2) F_fast
  x ← x + δt (p/m)
  F_fast ← F_fast(x)
  p ← p + (δt/2) F_fast
F_slow ← F_slow(x)
p ← p + (Δt/2) F_slow
```

The extension to $N$ particles in three dimensions follows the single-step integrator. **More than two scales are allowed**. A split $F = F_{\mathrm{fast}} + F_{\mathrm{intermed}} + F_{\mathrm{slow}}$ uses three steps $\delta t$, $\Delta t = n\,\delta t$, and $\Delta T = k\,\Delta t$. Each such scheme can be made time-reversible.

## 6 The canonical ensemble and Helmholtz free energy

A folding experiment is held at fixed temperature, not sealed at one energy. The canonical ensemble is the probability distribution for that situation. Once the controls are $N$, $V$, and $T$, the state function that decides whether a change is favorable is the Helmholtz free energy $A=E-TS$: a decrease in $A$ is the thermodynamic sign that the change, folding included, is favored. The **canonical ensemble** uses the control variables $N$, $V$, and $T$. It describes a system in thermal contact with an infinite heat source. It is also the starting point for ensembles with other controls, including the isothermal–isobaric ensemble ($NPT$) and the grand-canonical ensemble ($\mu PT$).

For the microcanonical ensemble, $dE = T\,dS - P\,dV + \mu\,dN$, so the energy can be read as a function $E = E(N,V,S)$ by inverting $S = S(N,V,E)$. Then
$$
T = \left(\frac{\partial E}{\partial S}\right)_{N,V},
\quad
P = -\left(\frac{\partial E}{\partial V}\right)_{N,S},
\quad
\mu = \left(\frac{\partial E}{\partial N}\right)_{S,V}.
$$

The passage from $E(N,V,S)$ to a function of $N$, $V$, and $T$ is a Legendre transform. The Legendre transform of a function $f(s)$ is $\tilde f(t) = f(s(t)) - t\,s(t)$, where $t = f'(s)$. Applied here,
$$
A(N,V,T) = E - TS.
$$
$A$ is the **Helmholtz free energy**. It has an energetic part and an entropic part. A process that decreases the free energy is thermodynamically favorable, and a process that increases it is unfavorable. Thermodynamics can decide which of those two outcomes occurs. It does not determine the timescale of the process.

Differentiating $A = E - TS$ and substituting $dE = T\,dS - P\,dV + \mu\,dN$ gives
$$
dA = -S\,dT - P\,dV + \mu\,dN.
$$
Since $A = A(N,V,T)$, the coefficients are the partial derivatives
$$
\mu = \left(\frac{\partial A}{\partial N}\right)_{V,T},
\quad
P = -\left(\frac{\partial A}{\partial V}\right)_{N,T},
\quad
S = -\left(\frac{\partial A}{\partial T}\right)_{V,N}.
$$

### 6.1 The Boltzmann distribution

The canonical ensemble is usable only once each conformation has a weight. The Boltzmann factor $e^{-\beta H}$ is that weight, and the partition function $Q$ is the integral that normalises it. Averages, pressures, and free energies are then derivatives or ratios of $Q$, which is why the free-particle integral on the problem sheet is computed at all. Up to normalisation, the canonical phase-space distribution is the Boltzmann distribution
$$
f(x) \propto e^{-\beta H(x)},
\qquad
\beta = \frac{1}{kT}.
$$
The **canonical partition function** counts the accessible microscopic states,
$$
Q(N,V,T) = C_c\int dx\,e^{-\beta H(x)},
$$
where $C_c$ is a constant that cancels in ensemble averages. The average of a phase-space function $a(x)$ is
$$
\langle a\rangle
= \int dx\,a(x)\,f(x)
= \frac{\int dx\,a(x)\,e^{-\beta H(x)}}{\int dx\,e^{-\beta H(x)}}.
$$
This is the microcanonical average with the Dirac delta replaced by the Boltzmann factor. There, $\delta(H-E)$ kept a single energy. Here every phase-space point contributes, weighted by $e^{-\beta H}$, and the energy fluctuates. The symbol $dx$ stands for the product of every $dq_i\,dp_i$. When the limits of one variable do not depend on the others, the integrals are done from the inside outward: momenta first, then positions. For one free particle in a box of length $L$, with $H=p^2/(2m)$, the Boltzmann factor does not depend on $x$. It cancels between numerator and denominator, and the average position is the midpoint of the box,
$$
\langle x\rangle
= \frac{\int_0^L x\,dx}{\int_0^L dx}
= \frac{L}{2}.
$$
The link to thermodynamics is
$$
A(N,V,T) = -kT\ln Q(N,V,T).
$$
The macroscopic quantities are therefore
$$
\mu = -kT\left(\frac{\partial\ln Q}{\partial N}\right)_{V,T},
\quad
P = kT\left(\frac{\partial\ln Q}{\partial V}\right)_{N,T},
$$
$$
S = k\ln Q + kT\left(\frac{\partial\ln Q}{\partial T}\right)_{V,N}.
$$

### 6.2 Canonical sampling

A plain molecular-dynamics trajectory conserves energy, so it samples the microcanonical surface. 

The Boltzmann weight describes a system that *can* exchange energy with a heat bath. Canonical sampling is how a simulation produces that exchange, so the saved configurations can be averaged as a laboratory sample at temperature $T$ would be. The energy is not a conserved quantity. Molecular dynamics in this ensemble has to produce those fluctuations, so that a macromolecule can be simulated at the control variables of an ordinary experiment rather than as an isolated system at fixed $N$, $V$, and $E$.

The methods sought are canonical dynamics: trajectories whose microstates are distributed according to the canonical distribution. Generating a realistic trajectory is not required. Equilibrium averages computed from the sample can be compared with experiment. The lecture gives one method of each of two types.

#### 6.2.1 Extended phase space method

**Extended phase space.** Extra variables are added to the physical phase space so that they mimic a heat bath. One picture is an agent that compares the instantaneous kinetic energy with the value expected at the target temperature and rescales the velocities. For $N$ particles in $d$ dimensions the extended Hamiltonian uses a new coordinate $s$, its momentum $p_s$, and a parameter $Q$ that sets the timescale of the agent:
$$
H_e
= \sum_{i=1}^N \frac{\mathbf{p_{i}}^2}{2m_i s^2}
+ U(\mathbf{r_{1}},\ldots, \mathbf{r_{n}})
+ \frac{p_s^2}{2Q}
+ gkT\ln s.
$$
The term $gkT\ln s$ is chosen so that a microcanonical distribution on the extended space, which has $2dN+2$ dimensions, induces the canonical distribution on the physical space of $2dN$ dimensions. Carrying out that reduction fixes
$$
g = dN + 1.
$$
With that choice, integrating out $s$ and $p_s$ leaves the canonical partition function of the physical system, up to a constant prefactor. The equations of motion are Hamilton's equations for the extended Hamiltonian,
$$
\dot r_i = \frac{\partial H_e}{\partial p_i},
\quad
\dot p_i = -\frac{\partial H_e}{\partial r_i},
\quad
\dot s = \frac{\partial H_e}{\partial p_s},
\quad
\dot p_s = -\frac{\partial H_e}{\partial s}.
$$

#### 6.2.2 Isokinetic Ensemble
**Isokinetic ensemble.** An alternative constrains the momenta to a surface of constant kinetic energy and keeps a canonical weight for the positions. Its partition function is
$$
Q(N,V,T)
= C_I\int d^N p\,d^N r\,
\delta\!\left(\sum_{i=1}^N \frac{p_i^2}{2m_i} - K\right)
e^{-U(r_1,\ldots,r_N)/kT},
$$
where $C_I$ is a constant and $K$ is the kinetic energy. The position-dependent factor is the part that matters for equilibrium properties. A canonical distribution of the momenta is of little importance for those properties. The lecture states that a target temperature $T$ requires
$$
\sum_{i=1}^N \frac{p_i^2}{m_i} = (dN-1)kT,
$$
where $dN-1$ is the number of kinetic degrees of freedom left by the one constraint.

The equations of motion that enforce $\sum_i p_i^2/(2m_i) = K$ are
$$
\dot r_i = \frac{p_i}{m_i},
\qquad
\dot p_i = F_i - \alpha p_i,
$$
with $\alpha$ chosen so that the constraint holds at all times. Differentiating the constraint and substituting $\dot p_i$ gives
$$
\alpha
= \frac{\sum_{i=1}^N p_i\cdot F_i/m_i}{\sum_{i=1}^N p_i^2/m_i}.
$$
Hence
$$
\dot p_i
= F_i - \frac{\sum_{j=1}^N p_j\cdot F_j/m_j}{\sum_{j=1}^N p_j^2/m_j}\,p_i.
$$
A numerical step comes from the same Liouville splitting used for velocity Verlet: half a step of the constrained momentum update, a full streaming step of the positions, a force update, and another momentum half-step. While the positions, and therefore the forces, are held fixed, the momentum update is available in closed form through an auxiliary function $s(t)$ obeying $\ddot s = a + bs$ with $a$ and $b$ constant. The lecture's integrator uses
$$
p_i \leftarrow p_i + F_i\,\frac{s(\Delta t/2)}{\dot s(\Delta t/2)},
\qquad
r_i \leftarrow r_i + \Delta t\,\frac{p_i}{m_i},
$$
then recomputes the forces and repeats the momentum update.

## 7 Chain models

Questions about the size and stiffness of a long molecule do not need every atom. A chain of beads and links keeps the feature that matters, connectivity, and gives the end-to-end distance and the persistence length in closed form. A bead may be an atom or a monomer, such as a residue. A link may be a chemical bond or a pseudobond, for example a link between $\mathrm{C}_\alpha$ atoms that are not joined by a real bond. Bead $i$ sits at $r_i = (x_i,y_i,z_i)$, and the link from bead $i$ to bead $i+1$ is $l_i = r_{i+1}-r_i$. The average link length is
$$
l = \frac{1}{n}\sum_{i=1}^n \sqrt{l_i\cdot l_i},
$$
and the models below take $\langle l_i\cdot l_i\rangle = l^2$ for every $i$.

If bead $1$ is placed at $r_1$, then $r_i = r_1 + \sum_{j=1}^{i-1} l_j$. The **end-to-end vector** is
$$
h \equiv r_N - r_1 = \sum_{i=1}^{N-1} l_i.
$$
The **end-to-end distance** is the root-mean-square length $h \equiv \langle |h|^2\rangle^{1/2}$. Expanding the square separates a universal piece from a model-dependent piece:
$$
h^2
= nl^2 + \sum_{i=1}^{N-1}\sum_{j\neq i}\langle l_i\cdot l_j\rangle.
$$
The first term depends only on the number of links and their average length.

The **persistence length** is
$$
\xi_p = \frac{1}{l}\sum_{j=i}^{n}\langle l_i\cdot l_j\rangle.
$$
Flory (*Statistical Mechanics of Chain Molecules*, 1989) describes $\xi_p$ as the average sum of the projections of all bonds $j\ge i$ onto an arbitrary bond $i$ in an infinitely long chain. It is the average distance the chain travels before it has turned $90^\circ$. The sum from $j = i$ to $n$ equals the sum from $j = 1$ to $i$. Averaging the two forms over the choice of link $i$ relates the persistence length to the end-to-end distance:
$$
\xi_p = \frac{h^2}{2nl} + \frac{l}{2}.
$$

**Freely jointed chain.** The simplest model puts no restriction on bond angles and allows free rotation about every bond. The only restriction is the link length, taken to be the same value $l$ for all $n$ links. Then $\langle l_i\cdot l_j\rangle = l^2\delta_{ij}$, where $\delta_{ij}$ is $1$ if $i = j$ and $0$ otherwise, and the cross terms vanish. The model predicts $h^2 = nl^2$, and therefore $\xi_p = l$. The chain has no persistence beyond a single link.

**Freely rotating chain.** In a biomolecule the bond angles stay near their equilibrium values, apart from small high-frequency deviations, which a potential usually represents by stiff harmonic angle terms. The freely rotating chain keeps every bond angle at one fixed value $\theta'$ and every link at length $l$. A chain of $n$ links has $n-1$ bond angles. Adjacent links then satisfy
$$
l_i\cdot l_{i+1} = l^2\cos\theta,
$$
where $\theta$ is the supplement of the bond angle $\theta'$.

![[./attachments/Pasted image 20260924201853.png|231]]

A third link can be split into a part parallel to the previous link and a part perpendicular to it. The perpendicular part points along the dihedral angle and averages to zero because that angle is free. Only the parallel part survives, and in general
$$
\langle l_i\cdot l_{i+k}\rangle = l^2\cos^k\theta.
$$

![[./attachments/Pasted image 20260924201915.png|361]]

Substituting the cosine rule into $h^2$ and summing the resulting geometric series, for large $n$, gives
$$
h^2 \approx nl^2\frac{1+\cos\theta}{1-\cos\theta}.
$$
The approximation uses $\sum_{k=1}^{n} x^k \approx x/(1-x)$ and $\sum_{k=1}^{n} kx^k \approx x/(1-x)^2$. At $\theta = \pi/2$ one recovers $h^2 = nl^2$, the freely jointed result. If $\theta < \pi/2$ then $h^2 > nl^2$, and if $\theta > \pi/2$ then $h^2 < nl^2$. At $\theta = 0$ the chain is straight, so the end-to-end length should be $h = nl$, but the large-$n$ formula predicts an infinite length. Evaluating the sum before that approximation, at $\theta = 0$, restores the expected result $h^2 = n^2 l^2$.

**Rotational isomeric states.** A chain of carbon atoms is better described by tetrahedral bond angles, about $109.5^\circ$. Rotation about each rotatable link is taken to begin and end in one of three states: $t$, $g^+$, and $g^-$, with torsion angles $180^\circ$, $60^\circ$, and $-60^\circ$. These are the trans, gauche$^+$, and gauche$^-$ states. The first and last bonds are not rotatable, so the state of bond $i$ is $\phi_i\in\Gamma$ for $1 < i < n$, where $\Gamma=\{t,g^+,g^-\}$.

A rotational-state model assigns an energy $U(\phi)$ to the state vector $\phi = (\phi_2,\ldots,\phi_{n-1})$. The partition function is the sum of Boltzmann factors over every state,
$$
Z = \sum_{\phi_2\in\Gamma}\cdots\sum_{\phi_{n-1}\in\Gamma} \exp\bigl(-U(\phi)/kT\bigr).
$$
The free energy and entropy per chain follow at once:
$$
A_c = -kT\ln Z,
\qquad
S_c = -\frac{\partial A_c}{\partial T} = k\ln Z + kT\frac{\partial\ln Z}{\partial T}.
$$

If the bonds are independent, $U(\phi) = \sum_{i=2}^{n-1} u(\phi_i)$, the multiple sum factors. With energies $E_t$, $E_{g^+}$, and $E_{g^-}$ for the three states,
$$
Z = \Bigl(
e^{-E_t/kT} + e^{-E_{g^+}/kT} + e^{-E_{g^-}/kT}
\Bigr)^{n-2}.
$$
Any one of the bonds $i\in\{2,\ldots,n-1\}$ may be used to write the single-bond sum.

If neighbouring rotatable bonds interact, $U(\phi) = \sum_{i=3}^{n-1} u(\phi_{i-1},\phi_i)$, and the partition function is the sum, over all states, of the product of the pair factors $\exp(-u(\phi_{i-1},\phi_i)/kT)$.

## 8 Protein folding and Levinthal's paradox

The counting argument is why folding cannot be a blind walk through conformations, and why a picture with direction, the energy landscape of the next section, is introduced at all. **Protein folding** is the process in which a protein molecule organises itself into a three-dimensional native structure, also called the native conformation or the native state. That folded structure is the tertiary structure of the protein domain. Several folded domains may organise further into a larger complex; that step is the formation of quaternary structure.

Secondary structure and tertiary structure form together, and they remain distinct. Secondary structure is often fixed by local geometry. In an $\alpha$-helix, hydrogen bonds join monomers $k$ and $k+4$ along the main chain. Tertiary structure is shaped by the whole molecule.

Levinthal's counting argument runs as follows. A polypeptide of $100$ amino acids has $198$ torsion angles $(\phi,\psi)$, and each angle is supposed to have only $3$ stable rotational states. The chain then has at most
$$
3^{198}
$$
distinct conformations. If the protein reached its native structure by sampling those conformations one after another, at a rate of one per nanosecond ($10^{-9}\,\mathrm{s}$) or even one per picosecond ($10^{-12}\,\mathrm{s}$), the time required would exceed the age of the universe. Many proteins fold spontaneously on a microsecond timescale. That contradiction is **Levinthal's paradox**.

One proposed resolution is a pathway from denatured conformations to the native conformation, possibly through a folding intermediate,
$$
U \rightleftharpoons I \rightleftharpoons N,
$$
where $U$, $I$, and $N$ are the unfolded, intermediate, and native states. A denatured conformation is fully unfolded, close to an extended chain. The alternative resolution is the energy-landscape, or folding-funnel, picture of folding.

### 8.1 Energy landscapes

Levinthal's count says a random search cannot finish. The landscape is the replacement picture: conformations are arranged so the chain is guided downhill, which is what makes a fold on a microsecond scale possible. Simulation and the folding mechanisms below are read off this surface. An energy landscape plots the internal free energy of one chain configuration against the degrees of freedom, such as the backbone torsion angles. That internal free energy excludes the conformational entropy. It includes everything else in the conformation: hydrogen bonds, ion pairs, torsion-angle energies, and hydrophobic and solvation free energies. Later slides call this quantity the energy. The lateral area of the landscape at a given depth represents the number of conformations, the conformational entropy, that share that internal free energy.

Folding, as a kinetic process, is a ball rolling on this surface. Hills, valleys, and plateaus guide the ball. The ball can move downhill or uphill, and uphill steps are less frequent.

The landscape version of Levinthal's paradox is a golf course: the native conformation is found only by a random search.

![[./attachments/Pasted image 20260924201959.png|400]]

The pathway resolution is a tunnel from a denatured conformation to the native state. The search is directed, and folding is faster than a random search.

![[./attachments/Pasted image 20260924202017.png|400]]

A pathway is tied to a few specific intermediate structures, which corresponds to a narrow tunnel. Ensembles of intermediate structures are represented by funnels. Once the protein is caught by a folding funnel, it changes conformation so that its energy falls. The fall in energy compensates the loss of conformational entropy that comes from moving deeper in the landscape.

![[./attachments/Pasted image 20260924202551.png|627]]

A. In an ideal funnel every conformation moves rapidly toward the native state, and every denatured conformation has the same smooth access to it.

B. A rugged funnel has energy barriers and kinetic traps, and it is linked to folding through more than one state.

C. A moat can account for proteins that fold in either of two ways: a fast process, or a slow process that enters kinetic traps.

D. Conformational entropy can also slow folding by creating free-energy barriers. On a champagne-glass landscape the search across the flat plateau is undirected, and that search is the rate-limiting step.
### 8.2 Folding mechanisms

A funnel does not say which contacts form first in a given protein. These three mechanisms are the concrete sequences of events the landscape picture is broad enough to contain. The landscape picture is not meant as a full account of any one protein's folding mechanism. Most mechanisms are compatible with it. The lecture names three.

1. **Hydrophobic collapse.** The unfolded chain first makes a fast transition to a globule, with hydrophobic residues inside and hydrophilic residues outside. A slower second stage builds secondary structure until the tertiary native structure is reached. *On a funnel, the first stage is the more downhill of the two*.
2. **Diffusion–collision.** Secondary structure, such as $\alpha$-helices and $\beta$-strands, forms rapidly. The chain, with that secondary structure in place, then assembles into the native tertiary structure. On a funnel, the second stage is the more downhill move.
3. **Nucleation–condensation.** A small part of the chain collapses first and forms a nucleus. The nucleus may contain non-local contacts and some partial secondary structure. From that nucleus the tertiary structure grows together with the rest of the secondary structure. *Both stages are downhill moves on the landscape*. 
	- Most successful of the three, constant compensation for loss of conformational entropy by losing energy

![[./attachments/Pasted image 20260924202302.png|423]]
> Illustration of the three protein folding mechanism:
> (Top) Hydrophobic-collapse, (Middle) Nucleation-condensation,
> (Bottom) Diffusion-collision.
### 8.3 Limits of folding simulation

The earlier sections are the checklist a folding simulation has to meet: a potential, a canonical sample, and a reversible multiple-timescale integrator, run for long enough. This section is why meeting the checklist can still fail to show a folded structure.

1. An accurate biomolecular potential.
2. Equations of motion that sample microstates with their canonical probabilities, or more generally with the probabilities of the chosen control variables, such as $N$, $V$, and $T$.
3. Time-reversible multiple-timescale integrators. The time step is still limited by resonance, but resonance-free time-reversible schemes exist and allow outer steps on the order of $100\,\mathrm{fs}$ (*Phys. Rev. Lett.* 93: 150201, 2004).
4. Enough computing to run long trajectories.

These conditions can still fail to produce a folded structure. The chosen sequence may not fold spontaneously, or the simulated trajectory may be shorter than the time on which that sequence folds. Unusually powerful computing resources can change that outcome.

Folding@home, started by Stanford professor V. J. Pande and co-workers, is a distributed-computing project. It uses GPUs, CPUs, and ARM processors on machines volunteered by users, and it is a platform for protein-folding trajectories on millisecond scales. Volunteers contribute at `https://foldingathome.org/`.

Anton is a special-purpose computer for molecular-dynamics simulations of biological macromolecules, built by a company headquartered in New York City.

![[./attachments/Pasted image 20260924202233.png|626]]

# References

- Minary, Computational Biology, Lectures 8–10 (Michaelmas 2025)
- Flory, *Statistical Mechanics of Chain Molecules* (1989), cited for the persistence length
- *J. Chem. Phys.* 97, 1990 (1992), https://doi.org/10.1063/1.463137
- *J. Chem. Phys.* 81, 511 (1984), https://doi.org/10.1063/1.447334
- *J. Chem. Phys.* 118, 2510 (2003), https://doi.org/10.1063/1.1534582
- *Nature Structural Biology* 4, 10–19 (1997), https://www.nature.com/articles/nsb0197-10
- *Phys. Rev. Lett.* 93, 150201 (2004), cited for resonance-free outer steps
- *Int. J. Mol. Sci.* 21, 2713 (2020), source of the timescale figure
- *BMC Biology* 9, 71 (2011), source of the bonded-term figure
