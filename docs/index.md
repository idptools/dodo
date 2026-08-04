# DODO

**re**Design **A**lphaF**o**ld **D**isordered regi**o**ns

DODO takes a predicted protein structure, works out which parts are folded domains and which are
intrinsically disordered regions, and rebuilds the disordered parts so they adopt realistic
polymer dimensions instead of AlphaFold's characteristic extended spaghetti.

AlphaFold isn't trying to represent a disordered region as an ensemble — an IDR has no single
structure to predict. But for figures and presentations it is useful to have those regions *look*
like what they are.

```bash
pip install git+https://github.com/idptools/dodo.git
pip install git+https://github.com/idptools/sparrow.git   # strongly recommended, see below
```

```bash
dodo fetch P04637 -o p53_dodo.pdb
```

```python
import dodo

report = dodo.rebuild("AF-P04637-F1-model_v6.pdb", mode="expanded", n_models=10, seed=0)
print(report.summary())
dodo.write_pdb(report.models, "p53_dodo.pdb")
```

```{toctree}
:maxdepth: 2
:caption: Contents

algorithm
guide
validation
api/index
```

## Install sparrow

The base install is numpy, scipy, getSequence and tqdm -- no torch. Everything works without
sparrow, but the *dimensions*
DODO builds to are then estimated from an analytical polymer scaling law that is blind to sequence
composition. With sparrow, they come from ALBATROSS, which predicts them from the sequence itself.

The difference is not subtle. For a 100-residue poly-glutamate region ALBATROSS predicts 122.2 Å
and the fallback estimates 68.8 Å; for poly-glycine the fallback *over*estimates, 68.8 Å against
54.9 Å. DODO warns when it has fallen back, but install sparrow.

sparrow is a separate step rather than an extra because it is not published on PyPI, and a PyPI
extra is not permitted to reference a git URL.

## What DODO is not

It is a geometric sampler, not a force field, and a multi-model file is not a trajectory. The
conformations are physically plausible and correctly sized; they are not sampled from a Boltzmann
distribution. It is a geometric sampler, not a substitute for a simulation.

## Where things are

- {doc}`algorithm` — what DODO actually does, in seven steps, and why the order matters
- {doc}`guide` — the command line and the Python API
- {doc}`validation` — how DODO checks its own output, and what the reference numbers are
- {doc}`api/index` — the generated reference

Source, issues and the changelog live at
[github.com/idptools/dodo](https://github.com/idptools/dodo).
