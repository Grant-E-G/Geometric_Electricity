# Geometric Electricity: design and measurement plan

Design baseline: 2026-10-05. The original idea remains in
[the project outline](one_and_done_hyperbolic_euclidean_board.md).

## Question and scope

Can one passive PCB distinguish the low-order spatial spectra and propagation
of Euclidean and hyperbolic disks using the equipment already available, and
support reproducible finite-graph spectral analysis?

Keep one board, one assembly, two independent 85-node networks, common ground,
and no active switch matrix. The intended comparison is between induced regular
triangular disks with coordination six and seven. Electrical traces do not set
the simulated metric; connectivity does. Physical placement coordinates and
mathematical embedding coordinates must remain separate.

## Equipment and measurement contract

The user confirmed a Rigol DHO924S (four channels, 250 MHz), JDS6600-60M generator,
PVP3150 probes, DE-5000 with TL-21/TL-22 fixtures, and T-962 oven. There is no VNA.
Use 10× probes, modeled as 10 MΩ in parallel with 10 pF at **each** connected
network node. Also simulate 15 pF uncertainty and 50 pF loading. CH1 measures
actual voltage before the injection resistor; CH2 remains at the source node;
CH3 scans the response. CH1 capacitance is on the generator side of the resistor,
not another network-node load. Ground every probe to the common ground.

The JDS6600's sine-frequency limit does not establish pulse rise time or arbitrary
waveform bandwidth. Confirm its delivered pulse on CH1 and use the measured
waveform for propagation interpretation. The DE-5000 is useful for loose-part
screening, but its 100 kHz ceiling cannot certify MHz Q. Do not treat measured
100 kHz Q or DC resistance as the loss at a circuit resonance.

## Operator and exact boundary convention

Let A be the adjacency matrix and d the retained degree. For each missing neighbor,
install one extra coupling capacitor to ground. Then

\[
Q=\operatorname{diag}(d)-A+\operatorname{diag}(q-d)=qI-A.
\]

This is a positive Dirichlet operator, not the zero-mode-bearing free-boundary
Laplacian. The default hardware has 330 pF per retained edge and **330 pF per
missing neighbor**, including multiple capacitors in parallel where necessary.
No special 990 pF or 1320 pF purchase is required. Removing all capacitive
boundary terminations would produce the free graph Laplacian; removing only a
few makes a local boundary defect and must not be called a global Neumann change.

For uniform L and C,

\[
J(\omega)=i\omega C Q+\frac{1}{i\omega L}I,
\qquad f_k=\frac{1}{2\pi\sqrt{LCq_k}}.
\]

Ascending q means **descending electrical frequency**. Low-order drum modes
are the highest resonances. With losses, the simulator stamps the full series
R-L branch rather than adding an arbitrary damping term to a graph Laplacian.
The Q sweeps use constant R = 2π(8 MHz)L/Q; consequently Q varies with frequency.
These are engineering sensitivity scenarios, not manufacturer broadband models.

With inductance disorder, the measured voltage eigenvectors are generalized
modes of C_matrix v = (1/ω²) diag(1/L_i) v. Ordinary orthonormal eigenvectors are
obtained after dividing voltage coordinates by sqrt(L_i). The self-adjoint
operator is S = diag(sqrt(L_i)) C_matrix diag(sqrt(L_i))/(L_nom C_nom).
Do not apply uniform-L orthogonality or heat-kernel formulas blindly to raw
voltage maps from a disordered circuit.

## Graph construction and mathematical validation

The Euclidean disk consists of axial triangular-lattice coordinates satisfying
 a²+ab+b² ≤ 21. It has 85 nodes, 222 retained edges and 66 missing-neighbor
capacitors. Graph-distance shells contain 1, 6, 12, 18, 24, 24 nodes.

The hyperbolic disk is the radius-three graph ball of the regular {3,7}
tessellation. Construct neighbors with Poincaré disk isometries and sevenfold
rotations. Its edge distance satisfies cosh(ell)=cos(2π/7)/(1-cos(2π/7)); the
embedding uses curvature -1. It has shells 1, 7, 21, 56, 196 retained edges and
203 missing-neighbor capacitors. Distances for the curvature -4 convention in
the reference are half these distances. This graph is generated independently;
we have not claimed equality to the reference's original node numbering or its
archived Euclidean adjacency file.

Checks cover node count, connectivity, degree, cyclic symmetry, positivity,
row sums, eigen-equation residuals, independently exported schematic pin nets,
and PCB capacitor incidence plus boundary multiplicities. Pin/net parity alone
does not prove that copper is routed: KiCad DRC must separately pass.

## Electrical scale and simulation findings

At 1 µH / 330 pF the nominal full spectrum is:

| Geometry | Minimum electrical frequency | Maximum electrical frequency |
| --- | ---: | ---: |
| Euclidean | 2.9455 MHz | 15.6483 MHz |
| Hyperbolic | 2.8070 MHz | 7.6601 MHz |

Thus the outline's estimated 1.4–9.6 MHz window is not the spectrum of these
specific disks. Both computed windows are accessible to the confirmed generator
and scope. Start with a 2–17 MHz sweep and refine to ≤5 kHz spacing near selected
peaks. The continuum interpretation is restricted to low-order modes.

The first six distinct nominal resonance clusters, in increasing q order, are:

| Geometry | Electrical frequencies, MHz |
| --- | --- |
| Euclidean | 15.6483, 9.9262, 7.5118, 7.0204, 6.1951, 6.1108 |
| Hyperbolic | 7.6601, 6.1478, 5.4316, 4.9893, 4.6011, 4.4550 |

These frequency lists by themselves do not label radial/angular modes. Determine
spatial profiles and angular content before identifying a reordered radial mode.

Seed 20261005, 200 trials per combination, independent uniform L tolerance,
independent edge C tolerance (including the budget ±5% case), conservative
lumped boundary-bank tolerance,
and 2 pF uniform stray capacitance:

| L tolerance | C tolerance | Scenario Q at 8 MHz | 5th percentile worst first-six-cluster subspace overlap, E / H |
| --- | --- | ---: | ---: |
| ±10% | ±5% | 36 (budget baseline) | 0.827 / 0.852 |
| ±10% | ±2% | 45 | 0.850 / 0.870 |
| ±5% | ±2% | 50 | 0.962 / 0.968 |
| ±2% | ±2% | 50 | 0.994 / 0.991 |

All scenarios use the same random draws within each geometry for a fair
comparison. Q affects the linewidth screen, not the lossless eigenspace overlap.
For the budget case, the 95th percentile maximum first-12 frequency shift is
about 1.6% for E and 1.7% for H. This is a simulation of bounded independent
uniform variation, not a manufacturing yield guarantee. Correlated lot errors
can move the common frequency scale and must be checked experimentally.

Degenerate eigenspaces are compared by principal angles after eigenvector-based
assignment; an arbitrary rotated basis inside a degenerate space is not counted
as a failure. Disorder and spectral resolution are separate issues. The cheap
10% option does not destroy all low-order identification, but gives worse spatial
fidelity. Higher-order modes are not covered by these low-order overlap numbers.
The gap-versus-linewidth screen is deliberately conservative: Q≈50 does not
separate all first-six clusters; the hyperbolic 4.601/4.455 MHz pair is especially
close. Q≈100 improves this screen substantially. It does not follow that 85
individual peak maps will be recoverable with either part.

With 2 pF stray at every node and two 10 pF probes at the designated intermediate
source and boundary node, the largest first-12 frequency shifts are about 1.08%
for E and 0.26% for H. Uniform stray is a simplification; copper extraction,
uneven node capacitance, inductor SRF, mutual inductance, capacitor ESR/ESL,
connector/cable capacitance and channel calibration still need assessment.
The simulations sweep C = 220, 330, 470, 680 and 1000 pF. Retain 330 pF as the
baseline because it preserves useful MHz bandwidth with modest modeled probe
loading. Changing C alone does not improve normalized mode spacing or Q.

Independent ngspice AC checks agree with the matrix solution to <10⁻⁶ relative
error. The ideal current-pulse transient agrees within 0.1% of peak voltage.
A second ngspice transient includes a 50 Ω generator, 1 kΩ injection resistance
and source/scan probe capacitance. Synthetic waveforms are predictions, not
measurements of the JDS6600's pulse or a fabricated board.

## Parts decision and assembly

The budget is a hard **under $200 USD all-in** limit for one populated board,
including spare parts, bare PCB order, stencil, paste/consumables, shipping, tax
and import charges. Existing measurement equipment and the owned oven are not
new purchases. Do not spend the allowance on precision components before obtaining
landed fabrication prices.

The cost-controlled CAD baseline is **Abracon AIML-0805-1R0K-T**, 1 µH,
±10%, Q≥45 at 10 MHz, 0805, and **Yageo CC0603JRNPO9BN331**, 330 pF,
±5%, C0G/NP0, 50 V, 0603. The current Abracon revision F specifies minimum
SRF 95 MHz for this value. A Q specification at 10 MHz does not establish Q
throughout the measurement band: the simulations also include Q=36 at 8 MHz,
corresponding to the same constant series resistance as Q=45 at 10 MHz.
Use standard KiCad 0805/0603 footprints, front-side assembly and one stencil.
Headers are installed afterward by hand. No custom inductor footprint is needed.

The following is a spending allocation, **not a fabrication quote**. Supplier
component listings were checked 2026-10-05; stock/prices can change. The board is
150×130 mm, outside common 100×100 promotional pricing. Obtain a quote at the
actual outline, 2 layers, ordinary FR-4, 1 oz copper and standard finish. Populate
one board; budget any mandatory minimum bare-board order in the fabrication row.

| Item | Buy quantity / allowance | USD |
|---|---:|---:|
| Abracon inductors, cut tape (170 used + 30 spare) | 200 × $0.09380 | 18.76 |
| Yageo capacitors (687 used; spare/defect stock included) | 800 × $0.013 | 10.40 |
| 1 kΩ 1% 0603 resistors and two 1×02 headers, with spares | allowance | 5.00 |
| Bare PCBs and one front stencil | maximum allocation | 70.00 |
| All supplier shipping combined | maximum allocation | 35.00 |
| Paste, hookup wire, probe ground adapters, consumables | allowance | 15.00 |
| Tax/import charges | allowance | 20.00 |
| Contingency | reserve | 20.00 |
| **Planned ceiling** | **$5.84 below the hard limit** | **194.16** |

Component sources: [DigiKey US cut-tape listing](https://www.digikey.com/en/products/detail/abracon-llc/AIML-0805-1R0K-T/2662996),
[LCSC capacitor listing](https://www.lcsc.com/product-detail/Multilayer-Ceramic-Capacitors-MLCC-SMD-SMT_YAGEO-CC0603JRNPO9BN331_C62784.html),
and [Yageo part specification](https://yageogroup.com/download/specsheet/CC0603JRNPO9BN331).
Buy cut tape rather than a full reel or paid Digi-Reel. Do not count a promotional
PCB price as the landed cost. If actual invoices exceed an allocation, use the
contingency only while keeping the total strictly below $200. If the total still
exceeds the ceiling, the design fails its cost acceptance gate; revise suppliers,
stencil choice or board size before ordering. No order has been placed.

Bourns CM453232-1R0JL (5%, Q≥50 at 7.96 MHz, 1812) remains a simulation
comparison, not an alternate that fits the present footprint. Its current
recommended land pattern differs from the standard 1812 footprint, so any future
upgrade requires a new verified footprint and layout. The expensive 1% capacitor
option is rejected. DE-5000 measurements can record component variation, but do
not assume its 100 kHz inductance measurement establishes RF Q or exact RF L.
Use a noninductive 1% injection resistor and calibrate its effective complex
impedance if phase precision warrants it.

T-962 drawer clearance must be checked against the actual usable tray, not only
nominal oven specifications. Keep 130 mm along the 180 mm drawer direction.
Use a thermocouple on the populated board and a profile compatible with the
actual solder paste and the inductor's specified reflow limits. Oven presets do
not establish the component temperature profile.

## Placement, routing and sparse controls

The KiCad project uses a 150 × 130 mm outline and two copper layers. Start with
recognizable elliptical graph shells, then relax actual component courtyards.
A successful placement is not evidence of a completed RF layout; record all DRC
results, copper continuity, label readability and return paths before release.
Node signal pads and nearby ground access support manual probing. Some labels
are displaced up to 12 mm to avoid pads and other silkscreen. Use the generated
`build/probe_map.pdf` and `build/probe_map.csv` to identify the **actual signal
pad**, and find its explicit testpoint reference in KiCad. Never assume the pad
closest to displaced text is the corresponding node. `graphs.json` records actual
placed signal coordinates separately from the mathematical metric and includes
nearby ground coordinates; some ground contacts are exposed return vias. A removable
patch selects one source; unused nodes are not connected to a switch matrix.
Source choices are center, shell two, boundary and a rotated shell-two partner.
Permanent long traces to a bank of duplicate node pads are avoided. The removable
source patch still has capacitance and inductance: use the shortest practical
connection, keep its geometry fixed, and calibrate it with the injection fixture.
At these frequencies, current inferred across the 1 kΩ resistor includes current
into the patch as well as the source probe. Subtract a measured fixture admittance
when extracting network impedance; do not attribute that current to the LC disk.
Record any extra source-node load in the model.

Six chosen retained edges on each geometry have a nominal capacitor plus an
unpopulated parallel footprint. Replace the nominal capacitor with ~165 pF for
weaker coupling, populate the parallel 330 pF for stronger coupling, or remove
the nominal capacitor for an open edge. Changes leave traces intact. Rank
candidate edges by low-mode voltage-difference sensitivity and distribute them
across graph shells; record the exact selected edges in graphs.json and BOM.
`simulation.json` contains predictions for opening, halving and doubling each
selected edge, plus removing one boundary capacitor. It reports frequency shifts
and low-mode inverse participation ratios; changed sorted indices are not assumed
to preserve mode identity. The assembly PDF marks the parallel DNP footprints.
For local boundary experiments, remove/repopulate selected individual boundary
capacitors. Changing one of several parallel capacitors gives a known step in
termination strength. It does not change the global boundary condition.

Short ground-return vias and ground pours on both sides connect the reference
around routed signal paths. Pours use full pad connections, 0.25 mm clearance
and 0.15 mm minimum copper width; the many return vias stitch the two sides.
Full ground connections increase heat sinking at assembly, reinforcing the need
for a measured reflow profile. The routing exchange excludes pours, preserves
existing return copper and refills the native KiCad pours after import. The initial
signal pass may skip the ground class; the final connectivity check includes it.
The 0–10 pF stray sweep is a first capacitance budget; it does not validate
RF ground impedance or layout capacitance extraction. Inspect actual routed
lengths and ground-plane continuity before treating this as the final placement.

The saved board is routed. Native KiCad 9.0.9 checks on 2026-10-05 reported
zero DRC violations, zero unconnected items and zero ERC violations; independently
reconstructed graph edges, boundary banks and schematic/PCB pins agree. Routing
used Freerouting 1.9.0 followed by local copper repairs and native pour refill.
A back-layer detour permits the last ground stitch, using a 0.50 mm via with
0.30 mm drill; other return vias use 0.60/0.30 mm. Selected 0.60 mm ground vias
in the probe map have explicit front solder-mask openings; the small stitch is
kept tented. The saved board preserves these routing edits. Regenerating
placement does not reproduce them automatically.
The routed layout contains 31 signal vias and 590 ground vias, with 53 selected
ground vias exposed on the front. The longest individual signal segment is
7.10 mm; this is not an extracted pad-to-pad path length or inductance bound.
Mapped ground-contact distances are at most 6.86 mm for E and 8.17 mm for H.
Check the actual probe spring or short ground adapter against these coordinates;
long clip leads require separate fixture characterization.

## Measurement sequence with existing equipment

1. Calibrate probe compensation, 10× scope scaling, common grounding, relative
   channel gain and phase using the same sinusoid over 2–17 MHz. Record calibration
   data and connected-probe count. Select acquisition bandwidth/sample rate to
   capture ≥10 points per cycle; use ≥20 cycles and averaging when needed.
2. Check the injection interface against known R/C loads. Use CH1 at the actual
   resistor input and CH2 at its output. Start at 200 mVpp generator amplitude,
   zero DC offset, and check linearity by repeating at half amplitude.
3. Measure representative loose L parts in a repeatable parallel-LC fixture at
   the relevant MHz frequencies using the generator and scope. Account for fixture
   and probe C and source-resistor loss. This qualifies RF Q without buying a VNA.
   Loose capacitors can be screened by the compensated DE-5000. Measuring network
   components in circuit will include parallel paths and is not a part-value test.
4. Sweep each designated source and repeat for both geometries. Compute complex
   source impedance rather than relying on visible voltage peaks alone. The
   generator's front-panel voltage is not the calibrated input current.
5. At each selected resonance, scan all 85 node pads with CH3 while CH2 remains
   at the source. Revisit reference nodes periodically to detect drift. Include
   narrow frequency windows around overlapping modes, not just one frequency.
   Scan multiple sources so an eigenvector zero at one source is not missed.
6. Fit common poles and residues across complex response maps. For separated
   poles, a residue vector estimates a spatial eigenvector up to normalization.
   Overlapping/degenerate poles require multiple frequencies/sources and subspace
   fitting. Full 85-mode reconstruction is an additional identifiability/noise
   problem, not a promised consequence of 85 accessible pads. The `fit` command
   estimates common conjugate pole pairs and residues by variable projection.
   It has recovered a noisy synthetic overlapping pair; it has not been validated
   against this equipment or a fabricated board. Moving-probe deembedding and
   channel calibration are prerequisites. Inspect fit residuals, conditioning,
   multiple seed sets and held-out frequencies before using fitted eigenpairs.
7. Repeat identical pulse injections at a boundary node and capture every node
   with the same trigger/reference. Characterize the delivered pulse first,
   choose a repeat period longer than the observed decay, and compare envelopes
   and instantaneous phases before boundary reflections dominate. Use metric
   coordinates, not PCB millimeters, to compare geodesics.
8. Change one defect or one local boundary bank at a time, record the assembly
   state and repeat the relevant measurements. Analyze localization through
   participation ratio and eigenvector/subspace changes, accounting for probe load.

Let Vin, Vs and Vi denote simultaneously measured complex voltages. With a real
sense resistor R, the current delivered past it is (Vin−Vs)/R. Remove source-probe
current Ys Vs, where Ys=1/10 MΩ+iωCprobe. Then Zss=Vs/Igraph and Zis=Vi/Igraph.
This is the calculation implemented by the CSV phasor command. It removes the
source probe's shunt current; it does **not** remove the moving scan probe from
the rest of the network.

To deembed the moving scan probe exactly in the lumped linear model, first
measure each Zii without CH3, correcting the CH2 source probe as above. With CH3
at i, the measured transfer is Zis_loaded = Zis/(1+Yi Zii). Thus
Zis = Zis_loaded (1+Yi Zii). Do not leave two probes on the same diagonal node
without counting both. Close to resonance this correction can amplify noise;
propagate calibration uncertainty and compare probe-loaded simulations.

## Derived mathematics and acceptance criteria

For a complete orthonormal eigenbasis of the specified self-adjoint graph operator,
compute K(t)=U diag(exp(-tq)) Uᵀ, heat trace, ζ(s)=Σq⁻ˢ, Q⁻¹ and participation
ratios. These are finite-disk observables, not an infinite-plane spectral zeta or
continuum heat kernel. Positive q avoids an arbitrary zero-mode cutoff.

If only m low modes are measured, report partial sums/kernels. When all lowest m
modes (including partners) are known to have been recovered, the unmeasured heat
trace tail is bounded by (85−m) exp(-t q_m); for s>0, the zeta tail is bounded by
(85−m) q_m⁻ˢ. For incomplete mode selection these bounds are unavailable. The
analysis command requires an explicit flag before producing these bounds.
A resolvent measurement is a frequency-dependent inverse admittance and includes
loss/source/probe effects; Q⁻¹ from eigenpairs is a separate mathematical object.

Fabrication gates are: complete routing and native DRC/ERC; independent graph/pin
parity; readable node/source/defect maps and usable ground access; actual part
MPNs and stencil verification; and a landed purchasing total strictly below
$200 including contingency. Keep manufacturing output marked as previews until
these gates are met. Screen representative loose inductors in the RF fixture
before committing to a full order where practical.

Routing and electrical connectivity gates pass for the current saved board.
Landed supplier quotes, physical ground-spring/patch fit, stencil inspection and
the measured assembly profile remain outstanding. The native checks do not
establish RF parasitics or experimental recovery of the complete spectrum.

After assembly, acceptance requires a measured/calibrated RF loss model,
adequate loaded signal/noise and repeatability, and demonstrated low-mode
separation or validated multi-pole extraction. These are physical experiment
steps, not results already obtained by simulation. Full-spectrum trace/zeta
claims additionally require complete eigenpair recovery or a stated truncation
bound. A fabricated prototype does not by itself validate the scientific claims.

## Primary references

- [Lenggenhager et al., paper and supplementary methods](https://arxiv.org/pdf/2109.01148).
- [Published article and archived data/code link](https://www.nature.com/articles/s41467-022-32042-4).
- [Abracon AIML-0805 specifications, revision F](https://abracon.com/Magnetics/chips/AIML-0805.pdf).
- [Bourns CM453232 datasheet](https://www.bourns.com/pdfs/cm.pdf) and
  [2025 change notice](https://www.bourns.com/docs/technical-documents/product-change-notifications/Bourns_IC25125_Chip_Inductors_PCN.pdf).
- [TDK NLV32 specifications](https://product.tdk.com/en/search/inductor/inductor/smd/info?part_no=NLV32T-1R0J-EF).
- [Obsolete TDK NL453232 status](https://product.tdk.com/en/search/inductor/inductor/smd/info?part_no=NL453232T-1R0J-PF).
- [DE-5000 manual](https://www.ietlabs.com/pdf/Manuals/DE_5000_im.pdf).
- [Rigol programming guide](https://download.rigol.com/en/Manual/Digital%20Oscilloscope/DHO900/DHO800900_ProgrammingGuide_EN.pdf).
- [JDS6600 supplier specifications](https://www.joy-it.net/en/products/JT-JDS6600).
