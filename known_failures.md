# Known failures on the human proteome

Every structure DODO does not handle cleanly, kept so the next person does not have to
rediscover them. Generated from a full run over the AlphaFold human proteome —
**23,587 structures**, seed 0, backbone on (the default) — on 2026-08-18.

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
seed before concluding a fix worked.

## 1. Introduced impossible contacts

**7 structures (0.030%).** An atom pair closer than 1.00 Å, the floor
below which no bond in any protein exists, that was *not* already in the input. This is the one
hard invariant DODO otherwise holds everywhere, so these are the highest-value cases to chase.

All of them are **long-range**: a rebuilt region landing on a distant part of the structure,
tens to hundreds of residues away in sequence. That is a different mechanism from the
adjacent-backbone class fixed in the same run (see the note at the bottom), and it is **not
diagnosed**.

| Accession | Residues | Regions | Contact |
|---|---:|---:|---|
| A2VEC9 | 1,400 | 7 | `CA` / `A:PHE820 CD1` at **0.422 Å** |
| A2VEC9 | 1,400 | 7 | `O` / `A:ARG782 NH2` at **0.561 Å** |
| A2VEC9 | 1,400 | 7 | `CD` / `A:THR667 O` at **0.730 Å** |
| A2VEC9 | 1,400 | 7 | `CB` / `A:ARG782 CB` at **0.751 Å** |
| A2VEC9 | 1,400 | 7 | `CA` / `A:ARG782 CB` at **0.894 Å** |
| Q12830 | 1,400 | 8 | `O` / `A:CYS369 N` at **0.819 Å** |
| Q6ISU1 | 281 | 2 | `O` / `A:SER117 N` at **0.739 Å** |
| Q92625 | 1,134 | 4 | `CA` / `A:ASN1134 O` at **0.843 Å** |
| Q96M86 | 1,400 | 8 | `N` / `A:GLN995 O` at **0.899 Å** |
| Q9BUR5 | 198 | 1 | `CA` / `A:LYS198 O` at **0.921 Å** |
| Q9BWT6 | 205 | 3 | `O` / `A:PHE201 N` at **0.859 Å** |

`A2VEC9` is worth starting with: it is the only structure that is simultaneously an impossible
contact, a blocked region, and the clash outlier of the entire proteome (101 introduced
clashes). If one structure is pathological, it is that one.

## 2. Blocked regions

**49 structures (0.21%), one blocked region each.** These are *honest*
failures, not silent ones: `report.ok` is `False`, the CLI exits 2, and the region keeps its
input coordinates rather than being replaced with anything invented. They are concentrated in
large proteins.

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

Two `GeometryError` cases are worth reading in full, because they are the input's problem and
no rebuild strategy fixes them:

- **O15050**, residues 1293–1309: The chain is broken next to this region's anchor: the anchor and the fixed residue beyond it are 128.82 A apart, where consecutive alpha carbons are 3.81 A. That residue exists only to pin the pseudo-
- **P78509**, residues 52–70: Cannot bridge anchors 78.04 A apart with 19 residue(s): 20 CA-CA bonds constrained to pseudo-angles of at most 161 degrees span at most 75.15 A. The region needs more residues, or the flanking domains
- **Q5H8C1**, residues 1850–1920: The chain is broken next to this region's anchor: the anchor and the fixed residue beyond it are 5.26 A apart, where consecutive alpha carbons are 3.81 A. That residue exists only to pin the pseudo-an
- **Q5VST9**, residues 712–729: Cannot bridge anchors 72.93 A apart with 18 residue(s): 19 CA-CA bonds constrained to pseudo-angles of at most 161 degrees span at most 71.40 A. The region needs more residues, or the flanking domains
- **Q5VST9**, residues 512–527: Cannot bridge anchors 72.01 A apart with 16 residue(s): 17 CA-CA bonds constrained to pseudo-angles of at most 161 degrees span at most 63.88 A. The region needs more residues, or the flanking domains
- **Q8NF91**, residues 838–851: The chain is broken next to this region's anchor: the anchor and the fixed residue beyond it are 5.02 A apart, where consecutive alpha carbons are 3.81 A. That residue exists only to pin the pseudo-an
- **Q8WXG9**, residues 1025–1039: Cannot bridge anchors 60.90 A apart with 15 residue(s): 16 CA-CA bonds constrained to pseudo-angles of at most 161 degrees span at most 60.12 A. The region needs more residues, or the flanking domains
- **Q8WZ42**, residues 873–883: Cannot bridge anchors 46.71 A apart with 11 residue(s): 12 CA-CA bonds constrained to pseudo-angles of at most 161 degrees span at most 45.09 A. The region needs more residues, or the flanking domains
- **Q96S53**, residues 385–399: The chain is broken next to this region's anchor: the anchor and the fixed residue beyond it are 89.29 A apart, where consecutive alpha carbons are 3.81 A. That residue exists only to pin the pseudo-a
- **Q9H251**, residues 734–747: The chain is broken next to this region's anchor: the anchor and the fixed residue beyond it are 88.48 A apart, where consecutive alpha carbons are 3.81 A. That residue exists only to pin the pseudo-a

## 3. High introduced-clash outliers

**11 structures with 20 or more introduced steric clashes** that are otherwise clean.
For scale: 65.6% of the proteome introduces zero, the 90th percentile is 4 and the 99th is 10.
These are soft van der Waals overlaps, not impossible contacts — visible to `dodo validate`,
but nothing is on top of anything.

| Accession | Residues | Regions | Introduced clashes |
|---|---:|---:|---:|
| Q6P3W6 | 1,400 | 1 | 29 |
| Q6P3W6 | 1,400 | 1 | 29 |
| P51826 | 1,226 | 3 | 26 |
| P49790 | 1,475 | 1 | 25 |
| Q7Z6E9 | 1,792 | 3 | 24 |
| A0A087WUL8 | 1,400 | 1 | 22 |
| O95996 | 2,303 | 5 | 22 |
| Q7Z460 | 1,538 | 5 | 22 |
| Q8IZQ1 | 1,400 | 9 | 22 |
| P08572 | 1,712 | 1 | 20 |
| P20930 | 1,400 | 1 | 20 |

## What was already fixed by this run

The first pass over the proteome found **86 structures (0.365%)** with introduced impossible
contacts. 72 of them were the same defect: adjacent carbonyl oxygens, `O(i)` against `O(i+1)`,
at 0.42–0.99 Å. The coupled-clash polish excluded pairs by *residue distance*, which hides a
pair that is five covalent bonds apart and free to collide, so it rotated peptide units into
collisions it could not see. The polish now counts real covalent separation, matching what
`backbone_refine._bond_separation` already did — **86 → 7**, which is the 7 in section 1.

Worth knowing when adding any future clash-scoring path: reuse `_bond_separation` rather than
re-deriving an exclusion rule. That is exactly how this bug was introduced.

