# Geometric Electricity

A passive LC circuit experiment comparing 85-node Euclidean and hyperbolic disks.

The [project outline](notes/one_and_done_hyperbolic_euclidean_board.md) describes the intended board. The [design and measurement plan](notes/research_plan.md) records calculations, equipment, decisions and remaining validation. Open [the KiCad project](hardware/geometric_electricity.kicad_pro) for the current CAD draft; [reproduction commands](code/README.md) generate simulations and verify connectivity.

The hardware is a design draft. Generated manufacturing files remain previews until routing, procurement and assembly checks pass. Bench validation follows assembly.

## Repository Layout

- `notes/`: research plans, reading notes, project logs, and informal planning.
- `math/`: definitions, theorem targets, proof sketches, examples, and open questions.
- `latex/`: manuscript source, macros, bibliography, and TeX inputs.
- `latex/render/`: rendered PDFs and TeX build products. This directory is ignored by git.
- `code/`: one consolidated simulation, CAD-generation and measurement-analysis script.
- `hardware/`: one KiCad project, sparse adjacency/metric data, BOM.
- `build/`: ignored simulation reports, figures, native checks and manufacturing previews.
- `sources/`: source material used during research.
- `sources/pdfs/`: local PDF references. PDF files are ignored by git.

## Working Sequence

1. Read the design/measurement plan and the original outline.
2. Reproduce simulations, probe loading, tolerance and defect sweeps.
3. Generate the KiCad schematic and placement; route and run native checks.
4. Confirm a landed total **below $200**, including PCB/stencil, shipping, tax,
   consumables and spares. Component precision must not break this constraint.
5. Verify the stencil/assembly map, assemble with a measured reflow profile and
   calibrate the RF fixture and probes.
6. Measure low-order responses and pulses, fit overlapping modes where justified,
   and state partial-spectrum limits explicitly.

All generated reports/previews stay in ignored `build/`. The refreshed spending
allocation is $188.72 with contingency. A recommended TDK inductor alternate
reduces the proposed allocation to $177.58 and improves tolerance; its purchasing
and assembly qualification are recorded in the design plan. Fabrication and
shipping are allowances, not confirmed quotes. No purchasing or physical
measurements have occurred.

## Coding Guidelines

- Prefer Rust for implementation work. Use Python as a backup, with type checking.
- Prefer functional code when possible. Relax this guideline when it causes performance issues.
- Long functions are fine.
- Avoid inheritance. Object-oriented design is generally discouraged.
- Avoid excessive code fragmentation. Longer files are acceptable when they keep related logic together.
- Lint code aggressively. Use Black for Python formatting.
- Shorter code is almost always better code.
