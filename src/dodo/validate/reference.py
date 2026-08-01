"""Reference bond geometry, measured across the human proteome.

Every number here was MEASURED from 23,587 AlphaFold DB structures (UniProt reference
proteome UP000005640, Homo sapiens, model v6), covering 105,299,823 intra-residue bond
measurements, 15,065,706 peptide bonds and 15,065,401 CA-CA virtual bonds.

Nothing here is transcribed from a table. v1 shipped a hand-written side-chain template library
whose own comment read "these might not be accurate TBH"; it was wrong by up to 2.8 A, and one
residue's coordinates had been copy-pasted from another's. Measuring is cheap -- the full scan
takes 17 seconds -- so there is no reason to guess.

Why AlphaFold structures are the right reference
-----------------------------------------------
DODO's input is predicted structures, and folded-domain atoms pass through untouched. So the
question a validator answers is "did DODO preserve the geometry it was given, and is what it
built consistent with it" -- which makes AlphaFold-derived statistics the correct baseline.

These are NOT experimental bond-length distributions and must not be cited as such. The
difference is measurable. An earlier version of this table was derived from four structures, of
which the largest by far was 6kn7, an experimental EM structure, so it was effectively
experiment-weighted. Replacing it with the proteome moved several bonds systematically:

    TYR CE2-CZ   1.3781 -> 1.4098      GLU CD-OE1  1.2217 -> 1.2526
    TRP CH2-CZ3  1.3757 -> 1.4061      ASP CG-OD2  1.2225 -> 1.2512

All aromatic C-C and carboxylate C-O, all shifted the same way, which is a systematic difference
between AlphaFold's restraint targets and experimental refinement rather than sampling noise.
Standard deviations also narrowed by a median of 0.007 A, because predicted geometry is more
uniform than measured geometry.

Two practical consequences. Tolerances derived from these standard deviations are STRICTER than
experimental spread would justify, which is what you want for catching a defect DODO introduced.
But validating an experimental structure against them will flag geometry that is perfectly real,
so a caller working with crystal or EM input should loosen the tolerance and expect to.

How bonded pairs are identified
------------------------------
By distance, which is unambiguous. Measured over 3,066,489 intra-residue pairs closer than 3.0 A
in 300 random structures, covalent bonds end at 1.85 A and 1-3 contacts begin at 2.15 A, with
just 23 pairs in the whole 1.85-2.10 A span. Any cutoff in that gap works; 1.90 is used.

The low-confidence caveat, which matters for tolerances
-----------------------------------------------------
AlphaFold geometry degrades where AlphaFold is unsure. 2.48% of CA-CA virtual bonds fall below
3.3 A, and those have a mean pLDDT of 38.8 against 72.2 for the rest. Part of that population is
genuine cis-proline -- proline is 25.6% of short pairs against a ~5% baseline, a 5x enrichment --
and the remainder is the model producing compressed geometry where it is guessing.

Two consequences for validation. A strict validator run against an unmodified AlphaFold structure
will flag a few percent of its bonds, and those belong to AlphaFold rather than to DODO, so a
validator has to distinguish a defect introduced from one faithfully preserved. And because DODO
rebuilds precisely the low-confidence regions, and builds them to exact bond lengths, rebuilt
regions should measure *cleaner* than the input they replaced.
"""

from __future__ import annotations

from typing import Final

#: Distance below which two atoms in one residue are taken to be covalently bonded, in
#: Angstroms. See the module docstring for the measurement that justifies it.
BOND_CUTOFF: Final[float] = 1.90

#: Number of structures the reference was measured from.
N_REFERENCE_STRUCTURES: Final[int] = 23587

#: Inter-residue peptide C(i)-N(i+1) bond, in Angstroms.
#: Measured n=15,065,706: mean 1.3376, sd 0.0086, 0.1-99.9 percentile 1.320-1.399.
#:
#: This is +0.0086 A from :data:`dodo.constants.C_N_PEPTIDE_BOND_LENGTH`
#: (1.329), which sits near the low edge of the observed distribution.
PEPTIDE_BOND_LENGTH: Final[float] = 1.3376
PEPTIDE_BOND_SD: Final[float] = 0.0086

#: CA(i)-CA(i+1) virtual bond for a TRANS peptide, in Angstroms.
#: Measured n=14,697,875 (bonds >= 3.3 A): mean 3.8317, sd 0.0842.
#:
#: +0.0217 A from :data:`dodo.constants.CA_CA_BOND_LENGTH` (3.81), well inside
#: one standard deviation. The 3.81 DODO builds to is the package author's specified value and
#: stands; this is recorded so a validator can set its tolerance from data.
CA_CA_TRANS_LENGTH: Final[float] = 3.8317
CA_CA_TRANS_SD: Final[float] = 0.0842

#: CA(i)-CA(i+1) for the short population: cis-peptides plus low-confidence geometry.
#: Measured n=367,526 (2.440% of all): mean 3.1912, sd 0.0776, min 1.339.
#: A validator must not treat these as errors when they came from the input.
CA_CA_SHORT_LENGTH: Final[float] = 3.1912
CA_CA_SHORT_FRACTION: Final[float] = 0.02440
#: Threshold separating the two CA-CA populations, in Angstroms.
CA_CA_CIS_THRESHOLD: Final[float] = 3.30

#: Mean pLDDT of residues whose CA-CA bond is short, versus the rest. Evidence that the short
#: population is partly AlphaFold uncertainty and not only real cis-peptides.
PLDDT_AT_SHORT_CA_CA: Final[float] = 38.8
PLDDT_AT_NORMAL_CA_CA: Final[float] = 72.2

#: Intra-residue bonds: ``{residue: {(atom_a, atom_b): (mean, sd, n, p0.1, p99.9)}}``.
#: Atom-name pairs are sorted, so look up with ``tuple(sorted((a, b)))``.
#: Only pairs observed at least 100 times are included, so a mis-parsed record cannot
#: introduce a phantom bond. Percentiles are given because several distributions are asymmetric,
#: and a symmetric sd-based tolerance misrepresents those.
RESIDUE_BONDS: Final[dict[str, dict[tuple[str, str], tuple[float, float, int, float, float]]]] = {
    "ALA": {
        ("C", "CA"): (1.5422, 0.0095, 1030991, 1.5189, 1.5994),
        ("C", "O"): (1.2298, 0.0040, 1031001, 1.2141, 1.2501),
        ("C", "OXT"): (1.2505, 0.0206, 1236, 1.2410, 1.6406),
        ("CA", "CB"): (1.5337, 0.0037, 1031001, 1.5155, 1.5640),
        ("CA", "N"): (1.4642, 0.0096, 1030988, 1.4398, 1.5214),
    },
    "ARG": {
        ("C", "CA"): (1.5318, 0.0103, 817707, 1.5079, 1.5901),
        ("C", "O"): (1.2268, 0.0043, 817710, 1.2093, 1.2482),
        ("C", "OXT"): (1.2476, 0.0205, 1256, 1.2368, 1.6298),
        ("CA", "CB"): (1.5386, 0.0062, 817710, 1.5205, 1.5628),
        ("CA", "N"): (1.4747, 0.0099, 817706, 1.4494, 1.5306),
        ("CB", "CG"): (1.5355, 0.0048, 817710, 1.5166, 1.5599),
        ("CD", "CG"): (1.5327, 0.0052, 817704, 1.5120, 1.5588),
        ("CD", "NE"): (1.4751, 0.0048, 817708, 1.4534, 1.4950),
        ("CZ", "NE"): (1.3208, 0.0045, 817707, 1.3013, 1.3438),
        ("CZ", "NH1"): (1.3037, 0.0044, 817708, 1.2830, 1.3704),
        ("CZ", "NH2"): (1.3015, 0.0049, 817704, 1.2790, 1.3665),
    },
    "ASN": {
        ("C", "CA"): (1.5411, 0.0090, 538574, 1.5185, 1.5961),
        ("C", "O"): (1.2297, 0.0042, 538583, 1.2134, 1.2513),
        ("C", "OXT"): (1.2481, 0.0030, 914, 1.2376, 1.2591),
        ("CA", "CB"): (1.5405, 0.0054, 538580, 1.5222, 1.5614),
        ("CA", "N"): (1.4665, 0.0092, 538581, 1.4420, 1.5214),
        ("CB", "CG"): (1.5180, 0.0043, 538584, 1.4992, 1.5340),
        ("CG", "ND2"): (1.3110, 0.0032, 538583, 1.2973, 1.3681),
        ("CG", "OD1"): (1.2190, 0.0039, 538581, 1.2026, 1.2341),
    },
    "ASP": {
        ("C", "CA"): (1.5419, 0.0100, 725887, 1.5176, 1.6053),
        ("C", "O"): (1.2319, 0.0042, 725895, 1.2164, 1.2554),
        ("C", "OXT"): (1.2527, 0.0031, 1139, 1.2405, 1.2651),
        ("CA", "CB"): (1.5404, 0.0061, 725891, 1.5198, 1.5643),
        ("CA", "N"): (1.4689, 0.0097, 725887, 1.4430, 1.5315),
        ("CB", "CG"): (1.5335, 0.0064, 725890, 1.5089, 1.5568),
        ("CG", "OD1"): (1.2524, 0.0040, 725895, 1.2368, 1.2677),
        ("CG", "OD2"): (1.2512, 0.0039, 725893, 1.2355, 1.2665),
    },
    "CYS": {
        ("C", "CA"): (1.5373, 0.0083, 339450, 1.5154, 1.5918),
        ("C", "O"): (1.2303, 0.0039, 339450, 1.2158, 1.2522),
        ("C", "OXT"): (1.2512, 0.0270, 718, 1.2411, 1.7554),
        ("CA", "CB"): (1.5344, 0.0049, 339449, 1.5133, 1.5614),
        ("CA", "N"): (1.4639, 0.0085, 339448, 1.4407, 1.5192),
        ("CB", "SG"): (1.8137, 0.0061, 339437, 1.7798, 1.8517),
    },
    "GLN": {
        ("C", "CA"): (1.5382, 0.0101, 718460, 1.5148, 1.5978),
        ("C", "O"): (1.2292, 0.0041, 718467, 1.2109, 1.2509),
        ("C", "OXT"): (1.2499, 0.0213, 1149, 1.2353, 1.6807),
        ("CA", "CB"): (1.5385, 0.0053, 718462, 1.5212, 1.5611),
        ("CA", "N"): (1.4684, 0.0100, 718461, 1.4428, 1.5267),
        ("CB", "CG"): (1.5369, 0.0041, 718465, 1.5199, 1.5609),
        ("CD", "CG"): (1.5261, 0.0037, 718465, 1.5058, 1.5378),
        ("CD", "NE2"): (1.3140, 0.0032, 718466, 1.3000, 1.3697),
        ("CD", "OE1"): (1.2195, 0.0036, 718468, 1.2070, 1.2353),
    },
    "GLU": {
        ("C", "CA"): (1.5402, 0.0105, 1086309, 1.5161, 1.6026),
        ("C", "O"): (1.2307, 0.0042, 1086320, 1.2126, 1.2509),
        ("C", "OXT"): (1.2518, 0.0031, 1241, 1.2413, 1.2635),
        ("CA", "CB"): (1.5398, 0.0063, 1086324, 1.5174, 1.5648),
        ("CA", "N"): (1.4705, 0.0102, 1086316, 1.4446, 1.5308),
        ("CB", "CG"): (1.5382, 0.0050, 1086321, 1.5191, 1.5620),
        ("CD", "CG"): (1.5378, 0.0055, 1086322, 1.5135, 1.5589),
        ("CD", "OE1"): (1.2526, 0.0039, 1086326, 1.2402, 1.2678),
        ("CD", "OE2"): (1.2525, 0.0037, 1086325, 1.2395, 1.2672),
    },
    "GLY": {
        ("C", "CA"): (1.5336, 0.0087, 976065, 1.5086, 1.5881),
        ("C", "O"): (1.2286, 0.0043, 976076, 1.2147, 1.2465),
        ("C", "OXT"): (1.2484, 0.0032, 885, 1.2381, 1.2624),
        ("CA", "N"): (1.4670, 0.0089, 976072, 1.4425, 1.5192),
    },
    "HIS": {
        ("C", "CA"): (1.5380, 0.0102, 388135, 1.5149, 1.5955),
        ("C", "O"): (1.2294, 0.0042, 388140, 1.2140, 1.2522),
        ("C", "OXT"): (1.2506, 0.0183, 781, 1.2409, 1.3705),
        ("CA", "CB"): (1.5363, 0.0059, 388137, 1.5161, 1.5619),
        ("CA", "N"): (1.4667, 0.0100, 388133, 1.4422, 1.5242),
        ("CB", "CG"): (1.5071, 0.0047, 388137, 1.4880, 1.5248),
        ("CD2", "CG"): (1.3684, 0.0046, 388125, 1.3560, 1.3749),
        ("CD2", "NE2"): (1.3959, 0.0082, 388118, 1.3673, 1.4105),
        ("CE1", "ND1"): (1.3423, 0.0047, 388125, 1.3296, 1.3530),
        ("CE1", "NE2"): (1.3406, 0.0052, 388125, 1.3276, 1.3527),
        ("CG", "ND1"): (1.3779, 0.0071, 388119, 1.3623, 1.4042),
    },
    "ILE": {
        ("C", "CA"): (1.5345, 0.0096, 663012, 1.5127, 1.5970),
        ("C", "O"): (1.2298, 0.0038, 663020, 1.2128, 1.2508),
        ("C", "OXT"): (1.2491, 0.0033, 1054, 1.2359, 1.2638),
        ("CA", "CB"): (1.5482, 0.0056, 663018, 1.5315, 1.5748),
        ("CA", "N"): (1.4711, 0.0090, 663012, 1.4486, 1.5302),
        ("CB", "CG1"): (1.5468, 0.0037, 663018, 1.5297, 1.5625),
        ("CB", "CG2"): (1.5409, 0.0036, 663019, 1.5204, 1.5650),
        ("CD1", "CG1"): (1.5329, 0.0040, 663017, 1.5059, 1.5638),
    },
    "LEU": {
        ("C", "CA"): (1.5371, 0.0092, 1474796, 1.5158, 1.5947),
        ("C", "O"): (1.2292, 0.0039, 1474802, 1.2126, 1.2512),
        ("C", "OXT"): (1.2491, 0.0031, 2562, 1.2390, 1.2622),
        ("CA", "CB"): (1.5410, 0.0047, 1474802, 1.5254, 1.5622),
        ("CA", "N"): (1.4670, 0.0090, 1474796, 1.4448, 1.5223),
        ("CB", "CG"): (1.5397, 0.0037, 1474799, 1.5236, 1.5584),
        ("CD1", "CG"): (1.5284, 0.0033, 1474797, 1.5094, 1.5616),
        ("CD2", "CG"): (1.5291, 0.0036, 1474794, 1.5100, 1.5613),
    },
    "LYS": {
        ("C", "CA"): (1.5335, 0.0105, 858603, 1.5099, 1.5967),
        ("C", "O"): (1.2271, 0.0043, 858624, 1.2094, 1.2495),
        ("C", "OXT"): (1.2474, 0.0128, 1686, 1.2391, 1.2623),
        ("CA", "CB"): (1.5389, 0.0063, 858618, 1.5198, 1.5638),
        ("CA", "N"): (1.4747, 0.0103, 858608, 1.4497, 1.5369),
        ("CB", "CG"): (1.5361, 0.0047, 858621, 1.5201, 1.5597),
        ("CD", "CE"): (1.5321, 0.0051, 858617, 1.5112, 1.5605),
        ("CD", "CG"): (1.5340, 0.0052, 858618, 1.5162, 1.5606),
        ("CE", "NZ"): (1.4837, 0.0055, 858618, 1.4643, 1.5560),
    },
    "MET": {
        ("C", "CA"): (1.5349, 0.0098, 310671, 1.5070, 1.5934),
        ("C", "O"): (1.2288, 0.0044, 310675, 1.2118, 1.2506),
        ("C", "OXT"): (1.2486, 0.0029, 496, 1.2401, 1.2586),
        ("CA", "CB"): (1.5372, 0.0050, 310674, 1.5190, 1.5602),
        ("CA", "N"): (1.4692, 0.0094, 310675, 1.4456, 1.5248),
        ("CB", "CG"): (1.5326, 0.0040, 310676, 1.5137, 1.5584),
        ("CE", "SD"): (1.8113, 0.0053, 310674, 1.7774, 1.8337),
        ("CG", "SD"): (1.8137, 0.0052, 310675, 1.7866, 1.8428),
    },
    "PHE": {
        ("C", "CA"): (1.5350, 0.0087, 530596, 1.5140, 1.5905),
        ("C", "O"): (1.2300, 0.0038, 530598, 1.2152, 1.2528),
        ("C", "OXT"): (1.2505, 0.0220, 1070, 1.2384, 1.7191),
        ("CA", "CB"): (1.5354, 0.0053, 530597, 1.5167, 1.5626),
        ("CA", "N"): (1.4652, 0.0086, 530599, 1.4436, 1.5203),
        ("CB", "CG"): (1.5160, 0.0044, 530601, 1.4970, 1.5358),
        ("CD1", "CE1"): (1.4055, 0.0048, 530565, 1.3941, 1.4138),
        ("CD1", "CG"): (1.4076, 0.0050, 530567, 1.3967, 1.4150),
        ("CD2", "CE2"): (1.4053, 0.0049, 530564, 1.3934, 1.4129),
        ("CD2", "CG"): (1.4070, 0.0050, 530567, 1.3959, 1.4139),
        ("CE1", "CZ"): (1.4049, 0.0048, 530565, 1.3952, 1.4099),
        ("CE2", "CZ"): (1.4049, 0.0049, 530566, 1.3949, 1.4097),
    },
    "PRO": {
        ("C", "CA"): (1.5417, 0.0136, 953619, 1.5161, 1.6151),
        ("C", "O"): (1.2297, 0.0041, 953640, 1.2166, 1.2521),
        ("C", "OXT"): (1.2504, 0.0134, 1413, 1.2407, 1.2644),
        ("CA", "CB"): (1.5328, 0.0039, 953439, 1.5156, 1.5607),
        ("CA", "N"): (1.4656, 0.0112, 953448, 1.4420, 1.5319),
        ("CB", "CG"): (1.5288, 0.0040, 953442, 1.5126, 1.5435),
        ("CD", "CG"): (1.5276, 0.0045, 953438, 1.5089, 1.5437),
        ("CD", "N"): (1.4576, 0.0086, 953453, 1.4385, 1.6452),
    },
    "SER": {
        ("C", "CA"): (1.5360, 0.0106, 1308986, 1.5113, 1.5965),
        ("C", "O"): (1.2292, 0.0043, 1309000, 1.2143, 1.2509),
        ("C", "OXT"): (1.2483, 0.0110, 2347, 1.2372, 1.2622),
        ("CA", "CB"): (1.5339, 0.0048, 1308995, 1.5138, 1.5599),
        ("CA", "N"): (1.4720, 0.0110, 1308991, 1.4446, 1.5322),
        ("CB", "OG"): (1.4144, 0.0045, 1308998, 1.3961, 1.4477),
    },
    "THR": {
        ("C", "CA"): (1.5326, 0.0112, 893870, 1.5083, 1.5981),
        ("C", "O"): (1.2290, 0.0042, 893887, 1.2126, 1.2486),
        ("C", "OXT"): (1.2487, 0.0152, 1162, 1.2362, 1.2656),
        ("CA", "CB"): (1.5393, 0.0061, 893881, 1.5219, 1.5640),
        ("CA", "N"): (1.4750, 0.0114, 893875, 1.4475, 1.5396),
        ("CB", "CG2"): (1.5278, 0.0034, 893882, 1.5080, 1.5610),
        ("CB", "OG1"): (1.4085, 0.0048, 893884, 1.3880, 1.4514),
    },
    "TRP": {
        ("C", "CA"): (1.5333, 0.0083, 175784, 1.5119, 1.5848),
        ("C", "O"): (1.2299, 0.0037, 175784, 1.2146, 1.2516),
        ("C", "OXT"): (1.2492, 0.0033, 274, 1.2383, 1.2612),
        ("CA", "CB"): (1.5332, 0.0055, 175783, 1.5144, 1.5612),
        ("CA", "N"): (1.4654, 0.0082, 175784, 1.4440, 1.5181),
        ("CB", "CG"): (1.4971, 0.0052, 175783, 1.4775, 1.5196),
        ("CD1", "CG"): (1.3523, 0.0053, 175782, 1.3393, 1.3821),
        ("CD1", "NE1"): (1.3804, 0.0039, 175769, 1.3677, 1.4043),
        ("CD2", "CE2"): (1.4264, 0.0053, 175756, 1.4124, 1.4320),
        ("CD2", "CE3"): (1.4114, 0.0050, 175769, 1.3958, 1.4176),
        ("CD2", "CG"): (1.4596, 0.0045, 175764, 1.4309, 1.4745),
        ("CE2", "CZ2"): (1.4011, 0.0055, 175771, 1.3848, 1.4098),
        ("CE2", "NE1"): (1.3772, 0.0041, 175772, 1.3608, 1.3874),
        ("CE3", "CZ3"): (1.4076, 0.0051, 175768, 1.3922, 1.4139),
        ("CH2", "CZ2"): (1.4046, 0.0049, 175769, 1.3907, 1.4099),
        ("CH2", "CZ3"): (1.4061, 0.0048, 175767, 1.3939, 1.4266),
    },
    "TYR": {
        ("C", "CA"): (1.5347, 0.0086, 388235, 1.5138, 1.5874),
        ("C", "O"): (1.2301, 0.0038, 388234, 1.2141, 1.2519),
        ("C", "OXT"): (1.2494, 0.0035, 701, 1.2397, 1.2667),
        ("CA", "CB"): (1.5353, 0.0054, 388234, 1.5174, 1.5627),
        ("CA", "N"): (1.4651, 0.0084, 388233, 1.4437, 1.5172),
        ("CB", "CG"): (1.5162, 0.0045, 388236, 1.4978, 1.5358),
        ("CD1", "CE1"): (1.4074, 0.0053, 388209, 1.3945, 1.4148),
        ("CD1", "CG"): (1.4085, 0.0050, 388208, 1.3974, 1.4157),
        ("CD2", "CE2"): (1.4072, 0.0054, 388210, 1.3940, 1.4149),
        ("CD2", "CG"): (1.4079, 0.0051, 388207, 1.3966, 1.4149),
        ("CE1", "CZ"): (1.4101, 0.0052, 388206, 1.3940, 1.4156),
        ("CE2", "CZ"): (1.4098, 0.0051, 388206, 1.3941, 1.4155),
        ("CZ", "OH"): (1.3574, 0.0032, 388233, 1.3384, 1.3907),
    },
    "VAL": {
        ("C", "CA"): (1.5315, 0.0106, 909353, 1.5094, 1.6016),
        ("C", "O"): (1.2292, 0.0039, 909360, 1.2136, 1.2508),
        ("C", "OXT"): (1.2484, 0.0029, 1503, 1.2380, 1.2615),
        ("CA", "CB"): (1.5436, 0.0057, 909359, 1.5280, 1.5692),
        ("CA", "N"): (1.4734, 0.0101, 909355, 1.4493, 1.5405),
        ("CB", "CG1"): (1.5337, 0.0037, 909357, 1.5123, 1.5621),
        ("CB", "CG2"): (1.5346, 0.0037, 909357, 1.5150, 1.5640),
    },
}
