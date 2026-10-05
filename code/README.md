# Reproducible LC design

`design.py` keeps graph generation, electrical models, KiCad generation and
measurement post-processing together. Python is the backup language allowed by
the repository guidelines: SciPy and KiCad's native `pcbnew` interface make it
appropriate here. No extra simulator framework or application is introduced.

Use KiCad 9, ngspice, system Python 3.10+ with `pcbnew`, and a virtual environment
that can see system packages:

```bash
uv venv --python /usr/bin/python3 --system-site-packages .venv
uv pip install --python .venv/bin/python numpy==2.2.6 scipy==1.15.3 matplotlib==3.10.9 black==26.10.0 mypy==2.4.0 ruff==0.16.10
.venv/bin/python code/design.py build --trials 200
.venv/bin/python code/design.py check
.venv/bin/black --check code/design.py
.venv/bin/ruff check code/design.py
.venv/bin/mypy --ignore-missing-imports code/design.py
```

`build` regenerates the schematic, **placed PCB**,
graphs and BOM. It replaces existing PCB edits/routing; preserve manual changes
before using it. `check` verifies the current board and schematic against the
independently reconstructed graphs without overwriting the copper.
`simulate` reruns the mathematics, tolerances, AC response and two transient
models without regenerating the KiCad files. All intermediate files, plots,
ngspice decks and manufacturing previews are in ignored `build/`.

Generated outputs include:

- `simulation.json` and `simulation.png`: seed, Monte Carlo scenarios, frequency
  clusters, probe/stray sweeps, C choices and illustrative responses.
- `E/H_simulation.npz`: q, voltage eigenvectors, frequencies and loaded AC maps.
- `E/H_pulse.npz`: ideal Gaussian current pulse with series inductor loss.
- `E/H_loaded_pulse.npz`: Gaussian voltage pulse, 50 Ω generator, 1 kΩ injection,
  and two 10 pF probes, simulated independently in ngspice.
- `spice_validation.json`: independent matrix/modal-versus-ngspice errors.
- `verification.json`, `erc.json`, `drc.json`: graph/pin and native rule checks.
- `budget.json`: component quantities, spending allocations and unconfirmed quote gates.
- `probe_map.pdf` / `probe_map.csv`: actual signal-pad coordinates, KiCad reference
  and nearest ground contact. Ground contacts include return vias and ground test pads.
- `manufacturing-preview/`: Gerbers including F.Paste, Excellon drills and an
  explicit NOT_FOR_FABRICATION marker; `positions.csv` and `schematic.pdf` / `assembly.pdf` (DNP footprints crossed out).

Record the actual assembly state, L/C/Q calibration, source, scanned node,
probe count, frequency and channel calibration with acquired data. The documented
measurement procedure is in `notes/research_plan.md`.

For one stopped scope acquisition, the optional LAN capture command sets the
waveform readout to RAW ASCII and reads CH1/CH2/CH3 from the same frozen record:

```bash
.venv/bin/python code/design.py capture --host SCOPE_IP --csv build/capture.csv
```

It assumes a SCPI socket at TCP 5555. The command has **not been tested against
the user's scope**; validate this interface and scaling with a known sine first,
or export CSV manually. The front panel must already have 10× probe scaling,
channels enabled and suitable acquisition settings. It does not configure or
control the JDS6600. The first 10,000 RAW samples are read; ensure that interval
contains the desired signal. SCPI failures are reported rather than hidden.

For manual or scripted captures, use exactly these columns (volts already scaled
to the probe tip): `time_s,vin_v,vsource_v,vscan_v`. Analyze a settled sine window:

```bash
.venv/bin/python code/design.py phasor --csv build/capture.csv --frequency 6000000 --sense 1000 --probe-pf 10
```

This fits sine/cosine/DC terms and subtracts source-probe current. If fixture
calibration supports a lumped extra source capacitance, pass `--fixture-pf VALUE`
to subtract it as well (default zero; do not guess it). Distributed fixture
impedance requires a fuller calibration. It assumes a
real sense resistance and calibrated relative channel phase/gain. It leaves the
moving scan probe's network loading in the transfer response; see the explicit
deembedding procedure in the design plan. A synthetic capture checks current
orientation and source-probe subtraction analytically during build/check.

For calibrated, current-normalized complex frequency responses, supply an NPZ
with `frequency_hz[n]` and `impedance_ohm[n,channels]`. Correct moving-probe loading
first; a moving load does not have exact common poles. Seed distinct resonance
clusters from a coarse sweep (one pole pair per degenerate cluster):

```bash
.venv/bin/python code/design.py fit --responses build/responses.npz --initial-mhz 4.42 4.63 --output build/poles.npz
```

Variable projection fits real-impulse-response conjugate pole pairs with complex
residues and a smooth background. Output includes damping, undamped frequencies,
residues, residual RMS and Jacobian condition. A synthetic two-channel overlap
case checks pole/damping recovery independently. Test several seed sets and
withhold frequencies to assess model error. Convergence is not completeness or
identifiability; residues still need calibrated spatial factorization before
becoming the real orthonormal vectors required below.

For **fitted** real orthonormal eigenvectors, provide an NPZ containing `q[m]`
and `vectors[85,m]`. Uniform-L voltage vectors can be used directly; for disordered
L, convert to coordinates of the self-adjoint operator described in the plan.
Raw voltage maps at unresolved peaks are not eigenvectors.

```bash
.venv/bin/python code/design.py derive --eigenpairs build/E_simulation.npz --output build/E_analysis.npz
```

This produces the finite-graph heat kernel/trace, zeta at s=1,2,3, Q inverse and
participation ratios. Fewer than 85 modes are explicitly marked partial. Only
when recovery of **all lowest m modes** has been established may
`--lowest-modes-complete` be used to add rigorous heat-trace/zeta tail bounds.
This command is not a common-pole fitting implementation.

For uniform-time pulse scans, provide `time_s[time]`, `voltage_v[time,85]` and
`source` in an NPZ:

```bash
.venv/bin/python code/design.py propagate --waveforms build/H_loaded_pulse.npz --geometry H --output build/H_propagation.npz
```

This extracts Hilbert envelopes, masked instantaneous phases and distances in the
true metric and embedding. Compare only reproducibly triggered scans before
reflections dominate. Hilbert phase alone does not establish causal arrival time.

The reproducible assertions and native checks are meaningful electrical/topology
checks, not substitutes for physical calibration, noise characterization or
fabrication readiness. Read the current design status before using any Gerbers.

Export/import routing through a Specctra-compatible router without regenerating
the schematic or placement:

```bash
.venv/bin/python code/design.py route-export
# Route build/routing.dsn in a compatible router, saving build/routing.ses.
# Preserve the ground returns; in Freerouting use -inc GND and -mt 0.
.venv/bin/python code/design.py route-import --session build/routing.ses
.venv/bin/python code/design.py check
```

The exported design assigns ground to its own class and protects its existing
return copper and excludes pours from the router's copy. The initial signal pass
may skip that class; a final repair pass must also connect remaining ground items.
Native refill/connectivity checks remain mandatory.
The saved board was routed with Freerouting 1.9.0, then repaired locally in KiCad
copper geometry. The initial router command was `java -jar freerouting-1.9.0.jar
-da -de build/routing.dsn -do build/routing.ses -mp 20 -mt 0 -inc GND`; the ground
repair pass omitted `-inc GND`. This version required a local graphical display.
Preserve the saved routed PCB: `build` produces placement, and a new router
run can produce different copper. `check` is the authority for the current native
connectivity result. Refresh exports after copper edits using
`python -c 'import sys; sys.path.insert(0, "code"); import design; design.export_preview()'`
with the same virtual environment, then run `check` again. Export refresh also
updates the ground-contact map and opens selected standard ground vias on the
front solder mask.
Import refills the ground plane. A successful router import is not a passing
DRC or permission to fabricate. Current manufacturing output is explicitly a
preview, and the budget's fabrication/shipping figures remain quote allowances.
