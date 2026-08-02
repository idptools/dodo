# Changelog

This changelog starts at 2.0. Earlier releases predate it.

It records what changes **for you as a user** — the commands and functions you call, the
structures that come out, and the things you have to do differently. Internal restructuring is
not listed, however large, unless you can see it from outside.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] — unreleased

DODO 2.0 is a rewrite. **It breaks the 1.x API deliberately**, so read the migration notes below
before upgrading a script.

The scientific behaviour changed too, in ways you can see in your output. The two that matter
most: build modes now mean something different, and a multi-model run is now a real ensemble.

### The install name changed

```bash
pip install "idptools-dodo[starling]"     # was: dodo[starling]
```

`import dodo` is unchanged and always will be. Only the distribution name moved, because PyPI's
`dodo` belongs to an unrelated project from 2014.

There is no longer an `[albatross]` extra. It pointed at a PyPI package called `sparrow` that is
**not** the sparrow DODO needs — PyPI's `sparrow` is an RDF/SPARQL library, and installing that
extra quietly put the wrong package on your import path. sparrow is not published on PyPI, and an
extra is not allowed to reference a git URL, so it is now an explicit step:

```bash
pip install git+https://github.com/idptools/sparrow.git
```

Install it. Without it DODO falls back to a sequence-blind approximation, and says so.

The `[viz]` extra is also gone. It installed matplotlib for a plotter that did not exist.

### Build modes mean something different — check your figures

In 1.x a mode was **Ångströms per residue**: `normal` was 0.8 Å/residue, so a 100-residue IDR got
an 80 Å target and a 500-residue IDR got 400 Å. Real IDR end-to-end distance scales as roughly
N^0.55, not linearly, so a fixed per-residue value can only be right at one length.

In 2.0 a mode is a **multiplier on the ALBATROSS-predicted end-to-end distance**, so it is
length-independent. `expanded` means 1.3× whatever this particular sequence is predicted to be.

**The same mode name now gives a different structure, and the difference runs in both
directions.** Measured on a generic disordered composition, across all six modes:

| Region length | What 2.0 gives, relative to 1.x |
|---|---|
| ~11 residues | **2.0–2.4× larger** |
| ~20 residues | 1.7–1.9× larger |
| ~40 residues | 1.3–1.4× larger |
| ~75–95 residues | about the same (the crossover) |
| ~150 residues | 0.74–0.81× |
| ~280 residues | 0.56–0.62× |
| ~1200 residues | 0.31–0.34× |

Short regions — most loops, and many linkers — come out **bigger** than 1.x made them. Long tails
come out considerably smaller. If you are regenerating a figure from a 1.x script, expect it to
look different, and treat 2.0's version as the correct one.

The ratio also depends on composition now, which is the point: at 20 residues the same
`super_compact` request gives 1.65× for poly-glycine and 2.85× for poly-proline, because 2.0 asks
ALBATROSS what *this* sequence does instead of applying one number to every sequence.

Also: `mode='normal'` and `mode='predicted'` are now exact synonyms. In 1.x they differed, by
1.70× on a 282-residue IDR.

### Multi-model output is a real ensemble now

`-n 10` writes ten conformers. In 1.x every model shared one arrangement *and* essentially one
end-to-end distance, so it was an ensemble at fixed dimensions — the one thing an IDR ensemble
should not be. Each model now draws its own dimensions from the physical distribution.

Folded domains are positioned once and held fixed across every model, so flicking between frames
in a viewer shows only the disordered regions moving.

### Added

- **A single `dodo` command** with subcommands, replacing three separate scripts. `dodo rebuild`,
  `dodo fetch`, `dodo sequence`, `dodo regions`, and `dodo validate`.
- **`dodo validate`** checks bond lengths, steric clashes and CONECT records on any PDB or mmCIF
  file, not just DODO's output. Its reference values were measured from 105,299,848 bonds across
  23,587 AlphaFold structures. Findings on atoms DODO did not build are labelled as inherited, so
  a defect that arrived in your input is not reported as DODO's.
- **`dodo regions`** prints what DODO thinks your structure looks like, with the score profile and
  threshold behind every boundary, without rebuilding anything.
- **Reports instead of print statements.** `dodo.rebuild()` returns a `RebuildReport` with
  `.ok`, `.models`, `.failures`, `.assignments` and `.outcomes`. In 1.x these functions returned
  `None` and printed, so a script could not tell whether a region had actually been rebuilt.
- **Region identification strategies**, selectable with `-s`: `density` (DODO's original all-atom
  metric, and the default), `contact` (a CA-only alternative), `plddt`, and `metapredict`.
  `auto` picks between the first two based on whether your input has side chains, and tells you.
- **Granular control from Python.** Assign regions however you like, then pass
  `strategy='preset'` and DODO builds exactly those:

  ```python
  import dodo
  structure = dodo.read_structure("model.pdb")
  dodo.assign_regions_from_spec(structure, {"A": [("idr", 1, 60), ("folded", 61, 912)]})
  report = dodo.rebuild(structure, strategy="preset")
  ```

  This replaces 1.x's `regions_dict=`, which took a separate description of the structure that
  could disagree with the real one.
- **mmCIF reading**, alongside PDB. Both accept `.gz`.
- **Reproducibility.** `seed=` / `--seed` makes output bit-identical.
- **`--ca-only`** and **`-b`/`--annotate-regions`**, which were 1.x's `-f/--no_FD_atoms` and
  `-b/--beta_for_FD_IDR`.

### Changed

- **Downloads no longer break.** 1.x built AlphaFold URLs with a hardcoded `model_v4`, which now
  returns HTTP 404 for every accession. The model version is resolved from the AFDB API instead.
- **Rebuilt regions land where the sequence says.** Folded domains are now translated and rotated
  so that the gap between them matches the predicted end-to-end distance of the linker joining
  them. AlphaFold packs domains 2–3.6× closer than the linker predicts, so without this the
  linker is built into a gap unrelated to its sequence.
- **Structures render correctly in VMD and PyMOL.** CONECT records are written by default, atom
  names occupy the correct columns, and the element column is populated. 1.x's CA-only output was
  read by MDTraj as *calcium* atoms.
- **Failures are explicit.** A region that cannot be built is reported with a reason and keeps its
  input coordinates. 1.x could return coordinate arrays of exact `(0, 0, 0)` rows, or NaN.
- **Requires Python 3.10+, numpy 2.0+ and scipy 1.13+.**
- Exceptions all derive from `DodoError`. Bad arguments raise `InvalidParameterError`, which is
  both a `DodoError` and a `ValueError`, so either `except` clause works.

### Removed

- `build.pdb_from_name()`, `build.pdb_from_pdb()`, `build.pdb_from_sequence()`, and the
  `pdb-from-name` / `pdb-from-pdb` / `pdb-from-sequence` commands. See the migration table.
- `regions_dict=`. Use `strategy='preset'` as shown above.
- `attempts_per_region` / `attempts_per_coord` / `attempts_per_residue` / `attempts_per_IDR`.
  Retry budgets are internal and reported when a region fails.
- `graph=True` and `plot_structure()`. No plotting is included.
- `linear_placement`. Folded domains are placed by facing the next domain and perturbing, then
  rejecting arrangements that clash.
- `just_fds=True`. Not reimplemented.
- `--all-atom` and `--sidechains`. Building backbone or side-chain atoms **for rebuilt regions**
  is not in 2.0. Note this is *not* 1.x's "all atom" behaviour, which meant keeping the atoms of
  the regions DODO does **not** rebuild — 2.0 does that unconditionally, so folded domains keep
  full atomic detail and only rebuilt regions are alpha-carbon only.

### Migrating from 1.x

| 1.x | 2.0 |
|---|---|
| `build.pdb_from_pdb(path, out_path=...)` | `dodo.rebuild(path)` then `dodo.write_pdb(...)` |
| `build.pdb_from_name(name, out_path=...)` | `dodo.fetch_alphafold(accession)` then `dodo.rebuild(...)` |
| `build.pdb_from_sequence(seq, out_path=...)` | `dodo.build_from_sequence(seq)` |
| `pdb-from-pdb in.pdb -o out.pdb` | `dodo rebuild in.pdb -o out.pdb` |
| `pdb-from-name "human p53" -o out.pdb` | `dodo fetch "human p53" -o out.pdb` |
| `pdb-from-sequence SEQ -o out.pdb` | `dodo sequence SEQ -o out.pdb` |
| `-o`, `--out_path` | `-o`, `--out` |
| `-c`, `--no_CONECT_lines` | `--no-conect` (no short form) |
| `-s`, `--silent` | `-q`, `--quiet` (`-s` now means `--strategy`) |
| `-u`, `--use_metapredict` | `-s metapredict` |
| `-f`, `--no_FD_atoms` | `--ca-only` |
| `-b`, `--beta_for_FD_IDR` | `-b`, `--annotate-regions` |
| `num_models=N` | `n_models=N`, or `-n N` |
| `mode='normal'` (0.8 Å/residue) | `mode='predicted'` — **different result**, see above |
| `use_metapredict=True` | `strategy='metapredict'` |
| `regions_dict={...}` | `assign_regions_from_spec(...)` + `strategy='preset'` |
| `beta_for_FD_IDR=True` | `write_pdb(..., annotate_regions=True)` |
| functions returned `None` and printed | functions return a `RebuildReport` |
| `except dodoException` | `except DodoError` |

Multi-word protein names now need quoting: `dodo fetch "human p53"`.

Exit status: `0` on success, `2` if some regions could not be built, `1` on error.
