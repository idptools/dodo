# Known failures on the human proteome

Every structure DODO does not handle cleanly, kept so the next person does not have to
rediscover them. Generated from a full run over the AlphaFold human proteome —
**23,587 structures**, seed 0, backbone on (the default) — on 2026-08-19.

**Everything here is rare and none of it is silent.** The run completed on 100% of structures
with no crash, hang or out-of-memory, rebuilt 63,082 of 63,149 regions (99.894%), and produced
zero bond defects in DODO's own geometry on every single structure. The cases below are the
tail, and they were judged shippable.

## Reproducing any of these

```bash
PROTEOME=your_path_to_human_proteome
dodo rebuild $PROTEOME/AF-<ACCESSION>-F1-model_v6.pdb -o /tmp/out.pdb --seed 0
dodo validate /tmp/out.pdb
```

Seed 0 throughout. A different seed will usually move these — most are the walk running out of
room in a crowded region, not a deterministic geometric impossibility — so reproduce with the
seed before concluding a fix worked. **Section 1 in particular churns freely between seeds**;
see the note at the end of it.

Fragment numbering matters for the long proteins. AlphaFold splits anything over 1,400 residues,
and several entries below are a fragment rather than `F1` — the accession alone will not
reproduce them. The exact file is named in each table.

## 1. Introduced impossible contacts

**12 structures (0.051%).** An atom pair closer than 1.00 Å, the floor below which no bond in
any protein exists, that was *not* already in the input. This is the one hard invariant DODO
otherwise holds everywhere, so these are the highest-value cases to chase.

Measure this as *introduced*, never as a raw count on the output. Several AlphaFold inputs ship
impossible pairs of their own — `AF-Q9BTC0-F1` arrives with 92 — and DODO rebuilds regions
containing some of them, so a number of structures come out **cleaner than they went in**.
Counting output pairs conflates DODO's defects with AlphaFold's and roughly triples the figure.

| File | Residues | Regions | Contact |
|---|---:|---:|---|
| AF-A2VEC9-F12 | 1,400 | 7 | `A:ASP610 CA` / `A:PHE820 CD1` at **0.422 Å** |
| AF-A2VEC9-F12 | 1,400 | 7 | `A:SER607 O` / `A:ARG782 NH2` at **0.561 Å** |
| AF-A2VEC9-F12 | 1,400 | 7 | `A:GLN619 CD` / `A:THR667 O` at **0.730 Å** |
| AF-A2VEC9-F12 | 1,400 | 7 | `A:SER611 CB` / `A:ARG782 CB` at **0.751 Å** |
| AF-A2VEC9-F12 | 1,400 | 7 | `A:SER611 CA` / `A:ARG782 CB` at **0.894 Å** |
| AF-Q86XZ4-F1 | 545 | 3 | `A:THR540 N` / `A:SER545 O` at **0.628 Å** |
| AF-Q6NXP2-F1 | 309 | 3 | `A:LEU241 N` / `A:ALA309 O` at **0.704 Å** |
| AF-Q12830-F9 | 1,400 | 8 | `A:GLU341 O` / `A:CYS369 N` at **0.819 Å** |
| AF-Q9NQS1-F1 | 362 | 1 | `A:MET337 N` / `A:SER362 O` at **0.841 Å** |
| AF-P84550-F1 | 965 | 4 | `A:HIS952 N` / `A:PRO965 O` at **0.852 Å** |
| AF-Q9BWT6-F1 | 205 | 3 | `A:ALA165 O` / `A:PHE201 N` at **0.859 Å** |
| AF-Q96M86-F2 | 1,400 | 8 | `A:PRO981 N` / `A:GLN995 O` at **0.899 Å** |
| AF-Q5T7V8-F1 | 369 | 1 | `A:GLY360 CA` / `A:THR369 O` at **0.919 Å** |
| AF-Q8WXT5-F1 | 416 | 2 | `A:PRO109 O` / `A:ARG206 N` at **0.921 Å** |
| AF-P58107-F16 | 1,400 | 6 | `A:LEU16 N` / `A:GLN167 O` at **0.972 Å** |
| AF-Q9P0J7-F1 | 381 | 1 | `A:ASN364 N` / `A:LEU381 O` at **0.982 Å** |

`AF-A2VEC9-F12` is worth starting with: it is the only structure that is simultaneously an
impossible contact, a blocked region, and the clash outlier of the entire proteome (101
introduced clashes). If one structure is pathological, it is that one.

**Do not read a change in this table as a regression without a seed control.** These 16 contacts
come from 23,587 structures and the set is not stable: any change that alters the random stream
reshuffles which structures land here. Measured over two independent build seeds, the count runs
11 and 15 for one steering policy and 16 and 15 for another — the between-policy gap is smaller
than the between-seed gap for either policy alone. Compare rates over at least two seeds, never
one A/B.

The contacts split into two mechanisms. Most are **long-range** — a rebuilt region landing on a
distant part of the structure, tens to hundreds of residues away in sequence. A minority are
**short-range**, the chain doubling back onto itself within 5–25 residues; those became relatively
more common when the steering width was loosened to fix the flat-zig-zag artifact (median
sequence separation moved from ~92 to ~48 residues at an unchanged total). Neither is diagnosed.
Both trace to the same exposure: `CA_CLASH_DISTANCE` is 3.2 Å, and two alpha carbons at 3.2 Å can
carry backbone N and O atoms far closer than any bond.

## 2. Blocked regions

**49 structures (0.21%), one blocked region each.** These are *honest*
failures, not silent ones: `report.ok` is `False`, the CLI exits 2, and the region keeps its
input coordinates rather than being replaced with anything invented. They are concentrated in
large proteins.

This set is unchanged from the previous run — the same 49 structures failing the same regions
for the same reasons.

| Cause | Count | What it means |
|---|---:|---|
| `ExhaustedAttemptsError` | 38 | The walk could not thread the region through the space available after 40 attempts. Usually a crowded region; often seed-dependent. |
| `GeometryError` | 11 | The request is geometrically impossible as stated — either the anchors are further apart than the residue count can span, or the input chain is already broken at the anchor. Not DODO's to fix. |

| Accession | Residues | Failed region | Length | Cause |
|---|---:|---|---:|---|
| A2VEC9 | 1,400 | 605–648 (connecting IDR) | 44 | `ExhaustedAttemptsError` |
| A6NIE9 | 313 | 1–67 (terminal IDR) | 67 | `ExhaustedAttemptsError` |
| B2RTY4 | 2,548 | 1169–1495 (connecting IDR) | 327 | `ExhaustedAttemptsError` |
| O15020 | 2,390 | 2107–2221 (connecting IDR) | 115 | `ExhaustedAttemptsError` |
| O15050 | 1,325 | 1293–1309 (loop in FD 1174-1310) | 17 | `GeometryError` |
| O15496 | 165 | 1–42 (terminal IDR) | 42 | `ExhaustedAttemptsError` |
| O60673 | 1,330 | 1212–1241 (connecting IDR) | 30 | `ExhaustedAttemptsError` |
| O75445 | 1,400 | 290–323 (connecting IDR) | 34 | `ExhaustedAttemptsError` |
| P24158 | 256 | 1–28 (terminal IDR) | 28 | `ExhaustedAttemptsError` |
| P24821 | 2,201 | 1111–1177 (connecting IDR) | 67 | `ExhaustedAttemptsError` |
| P39877 | 138 | 1–21 (terminal IDR) | 21 | `ExhaustedAttemptsError` |
| P78509 | 1,400 | 52–70 (loop in FD 39-422) | 19 | `GeometryError` |
| Q15027 | 740 | 359–407 (connecting IDR) | 49 | `ExhaustedAttemptsError` |
| Q15149 | 1,400 | 1173–1216 (connecting IDR) | 44 | `ExhaustedAttemptsError` |
| Q15149 | 1,284 | 973–1016 (connecting IDR) | 44 | `ExhaustedAttemptsError` |
| Q2V2M9 | 1,422 | 938–957 (loop in FD 895-1265) | 20 | `ExhaustedAttemptsError` |
| Q58FG1 | 418 | 198–214 (loop in FD 79-238) | 17 | `ExhaustedAttemptsError` |
| Q5CZC0 | 1,400 | 1297–1312 (loop in FD 1173-1400) | 16 | `ExhaustedAttemptsError` |
| Q5CZC0 | 1,400 | 1097–1112 (loop in FD 973-1296) | 16 | `ExhaustedAttemptsError` |
| Q5CZC0 | 1,400 | 897–912 (loop in FD 773-1096) | 16 | `ExhaustedAttemptsError` |
| Q5CZC0 | 1,307 | 697–712 (loop in FD 312-896) | 16 | `ExhaustedAttemptsError` |
| Q5H8C1 | 2,179 | 1850–1920 (connecting IDR) | 71 | `GeometryError` |
| Q5SZK8 | 1,400 | 705–715 (loop in FD 7-1400) | 11 | `ExhaustedAttemptsError` |
| Q5VST9 | 1,400 | 820–845 (connecting IDR) | 26 | `ExhaustedAttemptsError` |
| Q5VST9 | 1,400 | 712–729 (loop in FD 646-1248) | 18 | `GeometryError` |
| Q5VST9 | 1,400 | 512–527 (loop in FD 67-690) | 16 | `GeometryError` |
| Q6W4X9 | 2,439 | 755–764 (loop in FD 43-1188) | 10 | `ExhaustedAttemptsError` |
| Q702N8 | 1,843 | 315–484 (connecting IDR) | 170 | `ExhaustedAttemptsError` |
| Q7Z3E5 | 818 | 625–818 (terminal IDR) | 194 | `ExhaustedAttemptsError` |
| Q7Z442 | 2,459 | 544–564 (loop in FD 35-861) | 21 | `ExhaustedAttemptsError` |
| Q8IZH2 | 1,706 | 1177–1442 (connecting IDR) | 266 | `ExhaustedAttemptsError` |
| Q8NF91 | 1,400 | 838–851 (loop in FD 721-1371) | 14 | `GeometryError` |
| Q8NF91 | 1,400 | 321–351 (connecting IDR) | 31 | `ExhaustedAttemptsError` |
| Q8TDW7 | 1,400 | 56–117 (connecting IDR) | 62 | `ExhaustedAttemptsError` |
| Q8WXG9 | 1,400 | 1025–1039 (loop in FD 72-1393) | 15 | `GeometryError` |
| Q8WZ42 | 1,400 | 815–824 (loop in FD 277-1343) | 10 | `ExhaustedAttemptsError` |
| Q8WZ42 | 1,400 | 564–581 (loop in FD 1-672) | 18 | `ExhaustedAttemptsError` |
| Q8WZ42 | 1,400 | 873–883 (loop in FD 1-951) | 11 | `GeometryError` |
| Q8WZ42 | 1,400 | 19–28 (loop in FD 1-1385) | 10 | `ExhaustedAttemptsError` |
| Q8WZ42 | 1,400 | 1128–1138 (loop in FD 1-1382) | 11 | `ExhaustedAttemptsError` |
| Q92824 | 1,860 | 110–121 (loop in FD 33-602) | 12 | `ExhaustedAttemptsError` |
| Q96RW7 | 1,400 | 540–558 (loop in FD 46-1380) | 19 | `ExhaustedAttemptsError` |
| Q96S53 | 571 | 385–399 (loop in FD 384-417) | 15 | `GeometryError` |
| Q9BTC0 | 2,240 | 1393–1408 (loop in FD 1388-1413) | 16 | `GeometryError` |
| Q9H251 | 1,400 | 734–747 (loop in FD 733-1354) | 14 | `GeometryError` |
| Q9NQ94 | 594 | 301–447 (connecting IDR) | 147 | `ExhaustedAttemptsError` |
| Q9P2H5 | 1,018 | 600–760 (connecting IDR) | 161 | `ExhaustedAttemptsError` |
| Q9P2N4 | 1,935 | 1294–1310 (loop in FD 1183-1582) | 17 | `ExhaustedAttemptsError` |
| Q9Y613 | 1,164 | 673–692 (loop in FD 628-1018) | 20 | `ExhaustedAttemptsError` |

Every `GeometryError` case is the input's problem, and no rebuild strategy fixes them:

- **O15050**, residues 1293–1309: the chain is broken next to this region's anchor — the anchor and the fixed residue beyond it are 128.82 Å apart, where consecutive alpha carbons are 3.81 Å.
- **P78509**, residues 52–70: cannot bridge anchors 78.04 Å apart with 19 residues; 20 CA-CA bonds at pseudo-angles of at most 161° span 75.15 Å.
- **Q5H8C1**, residues 1850–1920: chain broken at the anchor — 5.26 Å where 3.81 Å is required.
- **Q5VST9**, residues 712–729: cannot bridge anchors 72.93 Å apart with 18 residues (maximum span 71.40 Å).
- **Q5VST9**, residues 512–527: cannot bridge anchors 72.01 Å apart with 16 residues (maximum span 63.88 Å).
- **Q8NF91**, residues 838–851: chain broken at the anchor — 5.02 Å where 3.81 Å is required.
- **Q8WXG9**, residues 1025–1039: cannot bridge anchors 60.90 Å apart with 15 residues (maximum span 60.12 Å).
- **Q8WZ42**, residues 873–883: cannot bridge anchors 46.71 Å apart with 11 residues (maximum span 45.09 Å).
- **Q96S53**, residues 385–399: chain broken at the anchor — 89.29 Å where 3.81 Å is required.
- **Q9BTC0**, residues 1393–1408: the two fixed residues flanking the region are 3.04 Å apart, closer than the 3.20 Å two non-bonded alpha carbons can be, and they are 18 residues apart along the chain.
- **Q9H251**, residues 734–747: chain broken at the anchor — 88.48 Å where 3.81 Å is required.

## 3. High introduced-clash outliers

**3 structures with 20 or more introduced steric clashes** that are otherwise clean.
For scale: 67.8% of the proteome introduces zero, the 90th percentile is 3 and the 99th is 8.
These are soft van der Waals overlaps, not impossible contacts — visible to `dodo validate`,
but nothing is on top of anything.

| File | Residues | Regions | Introduced clashes |
|---|---:|---:|---:|
| AF-A2VEC9-F12 | 1,400 | 7 | 101 |
| AF-Q7Z460-F1 | 1,538 | 5 | 22 |
| AF-Q8IZQ1-F8 | 1,400 | 9 | 22 |

This section is where the steering-width change did the most good. Under the previous flat width
this list held **12 structures** at seed 0 and **20** at seed 1; it now holds 3 and 4. Total
introduced clashes across the proteome fell from 26,823 to 22,595, the 99th percentile from 10 to
8, and the share of structures introducing none rose from 65.6% to 67.8%. `AF-A2VEC9-F12` is
unmoved at 101 and remains the single worst structure in the proteome.

## Previously fixed, kept for the lesson

An earlier pass found **86 structures (0.365%)** with introduced impossible contacts. 72 of them
were the same defect: adjacent carbonyl oxygens, `O(i)` against `O(i+1)`, at 0.42–0.99 Å. The
coupled-clash polish excluded pairs by *residue distance*, which hides a pair that is five
covalent bonds apart and free to collide, so it rotated peptide units into collisions it could
not see. The polish now counts real covalent separation, matching what
`backbone_refine._bond_separation` already did — **86 → 7** at the time.

Worth knowing when adding any future clash-scoring path: reuse `_bond_separation` rather than
re-deriving an exclusion rule. That is exactly how this bug was introduced.
