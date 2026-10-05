# One-and-Done Hyperbolic / Euclidean LC Board
## Revised architecture for a single fabrication + assembly

## Goal

Collapse the staged program into **one PCB, one stencil, one reflow assembly**, while preserving as much mathematical range as possible.

The board should provide:

1. a full 85-node Euclidean \(\{3,6\}\) network;
2. a full 85-node hyperbolic \(\{3,7\}\) network;
3. identical electrical scaling on both halves;
4. node-by-node probing;
5. a small number of built-in boundary/defect controls;
6. enough measurement access that the same hardware supports:
   - eigenvalue/eigenvector comparison;
   - heat kernel / heat trace;
   - spectral zeta function;
   - Green's functions / resolvent measurements;
   - time-domain propagation;
   - geodesic-spreading comparison;
   - disorder/localization analysis.

The guiding principle is:

> Do not make the topology itself electronically programmable. Put both reference geometries on the same board and make the analysis programmable.

This avoids hundreds of analog switches, their parasitics, and another debugging project.

---

# 1. Single-board architecture

Use two independent networks on one PCB:

- **Left:** 85-node Euclidean \(\{3,6\}\)
- **Right:** 85-node hyperbolic \(\{3,7\}\)

They share:

- ground;
- connector style;
- component families;
- measurement conventions;
- silkscreen conventions;
- excitation circuitry/interface;
- board environment.

They do **not** share graph nodes.

This gives a cleaner controlled comparison than using separate fabrication lots.

## Board dimensions

Target:

\[
\boxed{150\text{ mm} \times 130\text{ mm maximum}}
\]

Prefer smaller if routing permits.

Orient 130 mm along the T-962's 180 mm drawer dimension, leaving about 25 mm nominal clearance on each side.

Do not exceed 150 x 130 mm without explicit justification.

Use 2 layers unless routing proves impossible. Both graph families are planar, so crossings should not fundamentally require extra layers.

---

# 2. Component choice

For a one-and-done build, prioritize robust spectral resolution over saving the last ~$10.

Starting electrical scale:

\[
L \approx 1\,\mu\text H,\qquad C\approx330\,\text{pF}
\]

giving approximately a 1.4–9.6 MHz translated spectral window relative to the published 10 µH / 1 nF implementation.

The agent must compare at least:

### Robust option
A ~1 µH, ~5% inductor with Q around 50 near 8–10 MHz.

### Compact ultra-cheap option
Abracon AIML-0805-1R0K-T or similar:
- 0805
- 1 µH
- Q min ~45 at 10 MHz
- compact enough to make a dual-network PCB easier
- tolerance is ~10%, which is the main drawback.

The simulation must determine whether 10% L tolerance destroys mode identification. If it does, use the tighter-tolerance part.

Because this is intended as a one-time build, spending an extra ~$10–30 total for tighter L tolerance is preferable to needing another PCB revision.

Use C0G/NP0 capacitors. Prefer 1–2% if the cost delta is small.

---

# 3. Do not use a giant switch matrix

Do **not** attempt to switch every capacitor edge electronically.

An 85-node triangulated graph contains hundreds of edges. Making the union of Euclidean and hyperbolic edge sets electronically selectable would require hundreds of switching channels.

At ~1–10 MHz those switches introduce:

- on resistance;
- off capacitance;
- channel capacitance;
- routing complexity;
- control wiring;
- power rails;
- firmware;
- a new failure mode comparable in complexity to the experiment itself.

The cost advantage of sharing 85 cheap inductors is not worth it.

Two complete passive networks on one PCB are the simpler one-and-done solution.

---

# 4. Built-in reconfiguration: only where high-value

Reconfiguration should be sparse and intentional.

## 4.1 Excitation selector

Provide a simple header or jumper selection arrangement allowing a signal source to be routed manually to several preselected representative nodes on either geometry.

No active RF switching is required unless simulation shows a clear advantage.

Select sources corresponding to:

- center;
- intermediate shell;
- boundary;
- one or two symmetry-related points.

All nodes must still remain individually probeable.

## 4.2 Controlled defect bank

On each geometry reserve approximately 4–8 scientifically chosen edges for alternate coupling.

For each chosen edge provide footprints permitting:

- nominal \(C\);
- weaker coupling;
- stronger coupling;
- open edge.

Prefer solder-jumper / 0-ohm configuration unless a tiny DIP-switch implementation is clearly cleaner.

The board need not make every edge mutable.

## 4.3 Boundary-condition bank

Do not add one switch per boundary node.

Instead, reserve alternate termination footprints and group a *small number* of representative boundary-condition experiments where mathematically meaningful.

Default assembly should implement the same Dirichlet-like boundary condition used by the reference experiment.

The PCB should be modifiable later without cutting traces, but V1 measurement should require no modification.

---

# 5. What one board gives us

The following previously staged items are now **one experiment**.

## A. Euclidean vs. hyperbolic spectra

Measure resonances and recover the characteristic spectral reordering.

## B. Spatial eigenmodes

Measure node amplitudes/transfer functions at each resonance.

## C. Heat kernel and heat trace

Derived from the measured eigenpairs:

\[
K(t)=e^{-tQ}
\]

and

\[
\operatorname{Tr}e^{-tQ}
=\sum_\beta e^{-tq_\beta}.
\]

No additional hardware.

## D. Spectral zeta

\[
\zeta_Q(s)=\sum_{q_\beta\ne0}q_\beta^{-s}.
\]

No additional hardware.

## E. Green's function / resolvent

Driven response already gives access to the inverse circuit Laplacian / resolvent-like response.

Measure transfer functions between selected node pairs.

## F. Time-domain propagation

Inject the same localized pulse into both halves and scan the wavefront.

## G. Geodesic spreading

Compare wavefront spreading / arrival structure from nearby starting points in flat and negatively curved networks.

## H. Disorder / localization

Use measured component disorder plus the defect bank.

All of these should be treated as outputs of the same physical board, not future PCB revisions.

---

# 6. What should NOT be forced onto this PCB

Do not try to include:

- full Hodge \(L_1\) hardware;
- topolectrical nodal-knot circuitry;
- event-horizon / nonlinear transmission-line circuitry;
- active hyperbolic topological-band hardware;
- a sphere/torus/genus-2 network merely to claim extensibility.

Those require sufficiently different operators/couplings that including them would compromise the simplicity of the main experiment.

However, the software pipeline should generalize to arbitrary graph Laplacians so that a later board, if ever justified, is mostly a topology change.

---

# 7. Simulation-before-layout requirements

Before PCB freeze, answer:

1. Can two 85-node networks fit cleanly in <=150 x 130 mm?
2. Can they fit with the preferred tighter-tolerance inductors?
3. Does an 0805 Q~45 / 10% inductor remain usable after Monte Carlo disorder?
4. If not, what is the cheapest 5% part whose Q is sufficient?
5. Is 330 pF still the best C after probe and trace capacitance are included?
6. What resonances are expected across the accessible frequency range?
7. Which 4–8 defect edges are most mathematically informative?
8. Which excitation nodes best expose the Euclidean/hyperbolic contrast?
9. What measurement sequence reconstructs enough eigenvectors for the derived heat-kernel/zeta/Green-function analysis?

Do not freeze PCB dimensions until graph placement has been automatically attempted.

---

# 8. Layout requirements

Use scripted graph placement as the starting point.

For each half:

- preserve recognizable graph shells;
- label every node;
- place the inductor at/near its node;
- place coupling capacitors along graph edges;
- expose a probe pad for every node;
- include nearby ground access;
- keep high-parasitic test structures away from sensitive nodes;
- avoid large copper beneath long high-impedance graph edges unless simulation shows it is harmless.

Add clear silkscreen:

- `EUCLIDEAN`
- `HYPERBOLIC`
- node IDs
- shell/radius index
- designated source nodes
- designated defect edges.

---

# 9. Agent deliverables

Produce:

- one KiCad project for the combined board;
- one schematic containing both networks;
- one PCB <=150 x 130 mm;
- BOM;
- stencil/paste outputs;
- manufacturing outputs;
- graph adjacency files;
- automatic graph-to-KiCad generation where practical;
- automatic PCB-netlist-vs-target-graph verification;
- Monte Carlo/Q/probe-loading simulations;
- time-domain simulations;
- a measurement script/procedure;
- post-processing for spectrum, eigenvectors, heat trace, spectral zeta, Green's functions, and propagation.

---

# 10. Decision criterion

The desired endpoint is:

\[
\boxed{\text{One PCB + one reflow session + many mathematical experiments}}
\]

Do not optimize for the absolute smallest BOM if doing so creates a significant chance of another hardware revision.

The design should favor:
- passive components;
- generous test access;
- robust Q;
- controlled tolerances;
- simple measurement;
- analysis flexibility.

The main distinction from the prior plan is that **Euclidean and hyperbolic hardware are fabricated together**, while nearly all later “extensions” are extracted from the same measured eigenpairs and transfer functions rather than implemented as separate boards.
