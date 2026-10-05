# Theory Project Template

Template repository for pure math and theory-heavy research projects.

## Repository Layout

- `notes/`: research plans, reading notes, project logs, and informal planning.
- `math/`: definitions, theorem targets, proof sketches, examples, and open questions.
- `latex/`: manuscript source, macros, bibliography, and TeX inputs.
- `latex/render/`: rendered PDFs and TeX build products. This directory is ignored by git.
- `code/`: optional computational checks, symbolic experiments, and figure-generation scripts.
- `sources/`: source material used during research.
- `sources/pdfs/`: local PDF references. PDF files are ignored by git.

## Suggested Workflow

1. Start by editing `notes/research_plan.md`.
2. Extract stable definitions, claims, and proof obligations into `math/`.
3. Promote mature material from `math/` into `latex/`.
4. Put references and reading notes under `sources/` and `notes/`.
5. Use `code/` only when computation helps check an example, verify algebra, generate figures, or support reproducibility.

## Theory-First Convention

The central path is:

```text
notes/research_plan.md -> math/ -> latex/
```

Code is optional and should remain subordinate to the mathematical argument unless the project explicitly becomes computational.

## Coding Guidelines

- Prefer Rust for implementation work. Use Python as a backup, with type checking.
- Prefer functional code when possible. Relax this guideline when it causes performance issues.
- Long functions are fine.
- Avoid inheritance. Object-oriented design is generally discouraged.
- Avoid excessive code fragmentation. Longer files are acceptable when they keep related logic together.
- Lint code aggressively. Use Black for Python formatting.
- Shorter code is almost always better code.
