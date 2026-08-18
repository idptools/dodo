"""Place backbone N, C and O on a CA-only trace, from the neighbouring alpha carbons.

The advantage DODO has over general backbone prediction is that it already knows *every* alpha
carbon. That turns backbone placement from a modelling problem into a lookup, because of one fact
about protein geometry:

**The peptide unit is rigid.** Measured over 19,302 peptide units from 100 frames of all-atom IDR
simulation, its internal bonds do not move:

===============  =================  =============================
bond             measured           :mod:`dodo.constants` value
===============  =================  =============================
CA(i)-C(i)       1.525 +/- 0.004 A  ``CA_C_BOND_LENGTH`` 1.525
C(i)-N(i+1)      1.329 +/- 0.004 A  ``C_N_PEPTIDE_BOND_LENGTH`` 1.329
N(i+1)-CA(i+1)   1.458 +/- 0.004 A  ``N_CA_BOND_LENGTH`` 1.458
CA(i)-CA(i+1)    3.813 +/- 0.033 A  ``CA_CA_BOND_LENGTH`` 3.81
===============  =================  =============================

So given both flanking alpha carbons, the unit between them has exactly **one** degree of freedom:
a rotation about the CA-CA axis. Everything else is fixed. The whole problem is predicting that one
angle, and the neighbouring alpha carbons are what predict it.

Why four and not three, and why five beats four
------------------------------------------------
Three consecutive alpha carbons give the CA-CA-CA pseudo-*angle*. The peptide plane's rotation is
set by the pseudo-*dihedral*, which needs a fourth. Measured, held out on alternating frames:

==========  =========  =========  =========
predictor   C(i)       N(i+1)     O(i)
==========  =========  =========  =========
3 CAs       0.471 A    0.311 A    1.495 A
4 CAs       0.274 A    0.192 A    0.843 A
==========  =========  =========  =========

The fourth alpha carbon removes 53% of the peptide plane's angular uncertainty, from 64.4 to
30.5 degrees. It roughly halves the error on C and N.

A *fifth* adds a second, smaller increment. Conditioning a unit on both flanking pseudo-dihedrals
rather than only its own (:data:`_C_BY_BIN_2D`) is worth a further 5.1% on C and 3.8% on N in
held-out placed-atom error, each with a paired 95% CI excluding zero. That is the shipped
predictor for every interior unit with a fifth alpha carbon; the last interior unit of a region has
none and uses the 4-CA table above. Diminishing returns are real, though -- see the module note on
what was measured and rejected, and ``backbone_plan.md`` Phase BB-3.

Three things this module does *not* do, each for a measured reason
------------------------------------------------------------------
**It does not average positions.** Binning the local-frame coordinates and using the bin mean --
which is the obvious implementation -- produces C-N bonds of 1.254 A against an ideal 1.329, because
averaging points scattered on a sphere lands inside it. Every placed atom is instead rescaled onto
its exact ideal bond length, which costs nothing in accuracy (C 0.279 -> 0.281 A) and makes the
output valid by construction.

**It does not predict the carbonyl O.** O is fully determined by CA(i), C(i) and N(i+1): placed by
ideal geometry from the *true* three, it lands within **0.013 A**. Its 0.84 A error above is
entirely inherited from C and N, amplified because O sits 1.769 A from the CA-CA axis where C sits
0.556 A -- so the same angular uncertainty displaces it three times as far. Predicting it
separately is both redundant and worse.

**It does not invent alpha carbons at the termini.** Synthesizing a missing CA by reflection, or by
extrapolating at the mean pseudo-angle, was measured at 3.9-5.3 A error -- larger than the 3.81 A
bond itself, so a synthesized alpha carbon carries no information about the real one. Nor does the
adjacent peptide unit help: copying its plane rotation leaves a 99.0 degree residual, *worse* than
the 64.4 degree marginal, because successive peptide planes alternate rather than track. So a
terminal unit falls back to the marginal for its frame -- the honest 3-CA answer -- and the two
atoms no peptide unit covers are placed by ideal internal geometry.

Provenance
----------
The table is MEASURED from 100 frames of all-atom IDR simulation supplied by the author
(``subset_frames``), 19,302 peptide units, binned at 20 degrees. Every bin holds at least 146
observations. It is IDR data on purpose: this module exists to rebuild disordered regions.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Final

import numpy as np

from ..constants import (
    C_N_PEPTIDE_BOND_LENGTH,
    C_O_BOND_LENGTH,
    CA_C_BOND_LENGTH,
    CA_C_O_ANGLE,
    N_CA_BOND_LENGTH,
    N_CA_C_ANGLE,
    N_CA_C_WINDOW_MAX,
    N_CA_C_WINDOW_MIN,
)
from ..exceptions import GeometryError

if TYPE_CHECKING:  # pragma: no cover
    from ..structure import Structure

__all__ = [
    "BackboneResult",
    "CABackboneResult",
    "SeamStrain",
    "add_backbone_to_rebuilt",
    "backbone_from_ca",
]

#: Width of a table bin, in degrees of CA pseudo-dihedral.
#:
#: MEASURED choice. Narrowing below 20 degrees stops helping -- held-out C error is 0.283 A at 30
#: degrees, 0.281 at 20, 0.281 at 10 and 0.281 at 5 -- because the residual is the peptide plane's
#: intrinsic spread, not binning coarseness. 20 keeps every bin well populated.
_BIN_WIDTH_DEGREES: Final[int] = 20

_C_BY_BIN: Final[tuple[tuple[float, float, float], ...]] = (
    (+1.4223, -0.2170, +0.4056),  # n=423
    (+1.4230, -0.2540, +0.3861),  # n=613
    (+1.4226, -0.2919, +0.3709),  # n=787
    (+1.4217, -0.3135, +0.2853),  # n=855
    (+1.4207, -0.2905, +0.0443),  # n=656
    (+1.4192, -0.1793, -0.2858),  # n=558
    (+1.4169, -0.0321, -0.4767),  # n=948
    (+1.4173, +0.0406, -0.5059),  # n=1478
    (+1.4172, +0.1045, -0.5091),  # n=2144
    (+1.4183, +0.1581, -0.4886),  # n=2614
    (+1.4189, +0.1972, -0.4665),  # n=2824
    (+1.4195, +0.2447, -0.4339),  # n=2277
    (+1.4209, +0.3013, -0.3824),  # n=1487
    (+1.4200, +0.3691, -0.2863),  # n=676
    (+1.4216, +0.3166, -0.0690),  # n=291
    (+1.4210, +0.0592, +0.3254),  # n=146
    (+1.4180, -0.1227, +0.4088),  # n=201
    (+1.4208, -0.2036, +0.4170),  # n=324
)
_N_BY_BIN: Final[tuple[tuple[float, float, float], ...]] = (
    (+2.4047, +0.1532, -0.2535),  # n=423
    (+2.4065, +0.1816, -0.2439),  # n=613
    (+2.4035, +0.2039, -0.2350),  # n=787
    (+2.4051, +0.2202, -0.1704),  # n=855
    (+2.4055, +0.2061, -0.0169),  # n=656
    (+2.4024, +0.1234, +0.1837),  # n=558
    (+2.3988, +0.0380, +0.2837),  # n=948
    (+2.3980, -0.0154, +0.3049),  # n=1478
    (+2.3965, -0.0550, +0.3108),  # n=2144
    (+2.3979, -0.0995, +0.2998),  # n=2614
    (+2.3979, -0.1298, +0.2879),  # n=2824
    (+2.3982, -0.1737, +0.2651),  # n=2277
    (+2.4003, -0.2111, +0.2318),  # n=1487
    (+2.4006, -0.2538, +0.1724),  # n=676
    (+2.4019, -0.2071, +0.0242),  # n=291
    (+2.4023, -0.0430, -0.2006),  # n=146
    (+2.4043, +0.0830, -0.2624),  # n=201
    (+2.4062, +0.1408, -0.2570),  # n=324
)

#: Used when the fourth alpha carbon is missing but the preceding one is present -- the last
#: peptide unit of a chain. This is the 3-CA answer: C 0.471 A, N 0.311 A.
_C_MARGINAL: Final[tuple[float, float, float]] = (+1.4193, +0.0752, -0.2886)
_N_MARGINAL: Final[tuple[float, float, float]] = (+2.3998, -0.0459, +0.1756)

#: Used for the FIRST peptide unit, whose frame must be built forward from CA(i), CA(i+1), CA(i+2)
#: because there is no CA(i-1). Measured in that convention over 19,402 units: C 0.343 A, N 0.225 A.
_C_FORWARD_MARGINAL: Final[tuple[float, float, float]] = (+1.4193, -0.0140, +0.4396)
_N_FORWARD_MARGINAL: Final[tuple[float, float, float]] = (+2.3997, +0.0130, -0.2806)


#: 2D peptide-plane table, indexed ``[tau_i bin][tau_{i+1} bin]`` -- effectively a 5-CA
#: predictor, conditioning each interior unit on BOTH adjacent CA pseudo-dihedrals rather than
#: only its own. MEASURED win over the 1D table above (5-fold CV on the committed 19,302 units,
#: placed-atom error, i.e. after the exact bond projections): C 0.3139->0.2978 A (-5.1%),
#: N 0.1941->0.1868 A (-3.8%); the paired per-unit |error| improvement's 95% CI excludes 0 for
#: both. Used for every interior unit that has a fifth alpha carbon (CA[i+3]); the last interior
#: unit, which does not, falls back to the 1D ``_C_BY_BIN`` / ``_N_BY_BIN`` row for ``tau_i``.
#:
#: Cells no training unit reached are backfilled from that same 1D row AT DERIVATION TIME, so
#: every cell is defined and the placement path never branches on an empty bin. A min-count
#: threshold on sparse cells was measured and REJECTED: held-out error rose monotonically with
#: it (C 0.2978 at no threshold, 0.2988 at 10, 0.3015 at 20, 0.3105 at 75), so the sparse cells
#: carry signal rather than noise -- do not add one. Finer bins, bin interpolation and a
#: separate O predictor were all measured to give nothing; see backbone_plan.md Phase BB-3.
#:
#: Re-derive with ``scripts/derive_peptide_table.py``; ``tests/unit/test_backbone_table.py``
#: pins the regeneration and the placement accuracy.
_C_BY_BIN_2D: Final[tuple[tuple[tuple[float, float, float], ...], ...]] = (
    (  # tau_i bin 0 [-180, -160) deg, n=419
        (+1.4365, +0.0219, -0.3061), (+1.4619, -0.4479, +0.1134), (+1.4246, -0.2640, +0.4784),
        (+1.4224, -0.1853, +0.3316), (+1.4177, -0.2258, +0.4742), (+1.4222, -0.1998, +0.3897),
        (+1.4212, -0.2820, +0.4463), (+1.4199, -0.2070, +0.3624), (+1.4234, -0.1597, +0.3365),
        (+1.4216, -0.2124, +0.4673), (+1.4230, -0.2341, +0.4183), (+1.4242, -0.2319, +0.4243),
        (+1.4207, -0.2068, +0.4020), (+1.4241, -0.2558, +0.4638), (+1.4046, -0.1704, -0.0554),
        (+1.4262, -0.3627, +0.4063), (+1.4201, +0.1735, -0.5271), (+1.4223, -0.2170, +0.4056),
    ),
    (  # tau_i bin 1 [-160, -140) deg, n=608
        (+1.4196, -0.1293, +0.0434), (+1.4320, -0.3719, +0.1911), (+1.4220, -0.2856, +0.1025),
        (+1.4238, -0.2055, +0.1489), (+1.4299, -0.2768, +0.3745), (+1.4157, -0.1927, +0.5202),
        (+1.4256, -0.2686, +0.3927), (+1.4201, -0.2378, +0.3696), (+1.4252, -0.2667, +0.4479),
        (+1.4215, -0.2418, +0.3843), (+1.4242, -0.2540, +0.4235), (+1.4243, -0.2770, +0.4192),
        (+1.4205, -0.2848, +0.3997), (+1.4174, -0.2334, +0.3860), (+1.4212, -0.0783, +0.1495),
        (+1.4230, -0.2540, +0.3861), (+1.4075, -0.0996, -0.2903), (+1.4054, -0.2826, +0.5145),
    ),
    (  # tau_i bin 2 [-140, -120) deg, n=783
        (+1.4173, -0.3819, +0.2225), (+1.4162, -0.2475, +0.1431), (+1.4207, -0.2706, +0.1085),
        (+1.4288, -0.3612, +0.2062), (+1.4235, -0.3163, +0.3531), (+1.4234, -0.2930, +0.3380),
        (+1.4178, -0.2880, +0.3407), (+1.4216, -0.2990, +0.3780), (+1.4240, -0.2902, +0.4273),
        (+1.4233, -0.3153, +0.4010), (+1.4236, -0.2819, +0.4280), (+1.4234, -0.3023, +0.3925),
        (+1.4221, -0.2783, +0.4024), (+1.4199, -0.2517, +0.2739), (+1.4292, +0.0985, -0.2411),
        (+1.4095, -0.3979, +0.4202), (+1.4221, +0.0341, -0.5162), (+1.4164, -0.1598, +0.2268),
    ),
    (  # tau_i bin 3 [-120, -100) deg, n=851
        (+1.4282, -0.2346, -0.0438), (+1.4208, -0.3080, +0.0361), (+1.4219, -0.3568, -0.0031),
        (+1.4150, -0.2031, -0.0042), (+1.4223, -0.1762, +0.0478), (+1.4202, -0.2884, +0.2813),
        (+1.4276, -0.3480, +0.2611), (+1.4232, -0.2964, +0.3474), (+1.4198, -0.3452, +0.3989),
        (+1.4223, -0.3488, +0.3495), (+1.4217, -0.3338, +0.3539), (+1.4225, -0.3453, +0.3568),
        (+1.4209, -0.2485, +0.2614), (+1.4167, -0.2587, +0.1153), (+1.4178, -0.4972, +0.2039),
        (+1.4115, +0.1152, -0.1324), (+1.4332, -0.2305, -0.1829), (+1.4154, -0.1359, -0.1616),
    ),
    (  # tau_i bin 4 [-100, -80) deg, n=655
        (+1.4171, -0.1250, -0.4602), (+1.4120, -0.1226, -0.3936), (+1.4192, -0.1839, -0.4299),
        (+1.4167, -0.1503, -0.3489), (+1.4165, -0.2165, -0.1262), (+1.4249, -0.3160, +0.0391),
        (+1.4218, -0.3460, +0.0450), (+1.4248, -0.3731, +0.2185), (+1.4214, -0.3720, +0.1890),
        (+1.4232, -0.3886, +0.2404), (+1.4187, -0.3615, +0.2590), (+1.4233, -0.3392, +0.2184),
        (+1.4239, -0.3179, +0.1933), (+1.4261, -0.0644, -0.2747), (+1.4191, -0.0689, -0.2024),
        (+1.4276, +0.0020, -0.4285), (+1.4131, +0.0134, -0.4391), (+1.4240, +0.0588, -0.4417),
    ),
    (  # tau_i bin 5 [-80, -60) deg, n=555
        (+1.4230, -0.0373, -0.4007), (+1.4148, -0.0877, -0.4944), (+1.4212, -0.0400, -0.5038),
        (+1.4157, -0.0931, -0.4731), (+1.4162, -0.1643, -0.4022), (+1.4108, +0.0286, -0.4930),
        (+1.4160, -0.2074, -0.1926), (+1.4176, -0.2866, -0.1180), (+1.4204, -0.3610, -0.0593),
        (+1.4200, -0.3389, -0.0873), (+1.4199, -0.3017, -0.1091), (+1.4182, -0.3166, -0.0443),
        (+1.4208, -0.2529, -0.2463), (+1.4183, -0.1073, -0.3975), (+1.4211, -0.0246, -0.4556),
        (+1.4362, -0.1887, -0.3792), (+1.4205, +0.1386, -0.4816), (+1.4248, +0.0171, -0.4775),
    ),
    (  # tau_i bin 6 [-60, -40) deg, n=942
        (+1.4231, +0.0183, -0.4853), (+1.4192, +0.0804, -0.5026), (+1.4185, +0.0047, -0.5129),
        (+1.4163, +0.0874, -0.5178), (+1.4159, -0.0321, -0.5102), (+1.4136, -0.0994, -0.5170),
        (+1.4180, -0.0797, -0.4414), (+1.4164, -0.0884, -0.4451), (+1.4144, -0.1410, -0.4452),
        (+1.4125, -0.1389, -0.4371), (+1.4144, -0.1113, -0.4389), (+1.4175, -0.1494, -0.4140),
        (+1.4129, -0.0047, -0.4956), (+1.4220, +0.0059, -0.4551), (+1.4094, -0.0184, -0.5331),
        (+1.4157, +0.0527, -0.4460), (+1.4247, +0.0392, -0.5166), (+1.4211, +0.0732, -0.5195),
    ),
    (  # tau_i bin 7 [-40, -20) deg, n=1467
        (+1.4195, +0.1352, -0.4626), (+1.4182, +0.1047, -0.5014), (+1.4160, +0.0930, -0.5058),
        (+1.4174, +0.0879, -0.5209), (+1.4160, +0.0391, -0.5219), (+1.4184, +0.0382, -0.4974),
        (+1.4174, +0.0083, -0.5177), (+1.4173, -0.0109, -0.5036), (+1.4166, +0.0043, -0.5154),
        (+1.4166, -0.0250, -0.5024), (+1.4161, +0.0190, -0.5304), (+1.4162, -0.0113, -0.5176),
        (+1.4181, +0.0153, -0.5234), (+1.4149, +0.0527, -0.5241), (+1.4178, +0.0691, -0.3948),
        (+1.4228, +0.0474, -0.4823), (+1.4219, +0.0684, -0.3794), (+1.4207, +0.1167, -0.4846),
    ),
    (  # tau_i bin 8 [-20, +0) deg, n=2130
        (+1.4207, +0.1935, -0.4522), (+1.4179, +0.1564, -0.4761), (+1.4178, +0.1373, -0.4944),
        (+1.4185, +0.1610, -0.5033), (+1.4203, +0.1317, -0.5029), (+1.4169, +0.1146, -0.5199),
        (+1.4162, +0.1033, -0.5281), (+1.4181, +0.0743, -0.5086), (+1.4136, +0.0787, -0.5238),
        (+1.4171, +0.0687, -0.5205), (+1.4175, +0.0788, -0.5192), (+1.4164, +0.0721, -0.5191),
        (+1.4151, +0.1062, -0.5229), (+1.4190, +0.0840, -0.5112), (+1.4196, +0.1548, -0.4661),
        (+1.4195, +0.1128, -0.4432), (+1.4210, +0.1405, -0.4610), (+1.4209, +0.1652, -0.4753),
    ),
    (  # tau_i bin 9 [+0, +20) deg, n=2598
        (+1.4182, +0.2440, -0.4562), (+1.4155, +0.1893, -0.4734), (+1.4196, +0.2249, -0.4583),
        (+1.4200, +0.2278, -0.4569), (+1.4185, +0.1924, -0.4777), (+1.4209, +0.2157, -0.4713),
        (+1.4178, +0.1402, -0.4952), (+1.4205, +0.1320, -0.5017), (+1.4203, +0.1208, -0.4957),
        (+1.4176, +0.1249, -0.5078), (+1.4180, +0.1312, -0.5104), (+1.4166, +0.1533, -0.5011),
        (+1.4167, +0.1653, -0.4875), (+1.4151, +0.1334, -0.4761), (+1.4188, +0.1910, -0.4479),
        (+1.4231, +0.1810, -0.3394), (+1.4200, +0.2132, -0.4210), (+1.4182, +0.2287, -0.4239),
    ),
    (  # tau_i bin 10 [+20, +40) deg, n=2813
        (+1.4222, +0.3288, -0.3440), (+1.4191, +0.2743, -0.4322), (+1.4202, +0.2733, -0.4388),
        (+1.4173, +0.2771, -0.4357), (+1.4185, +0.2301, -0.4187), (+1.4171, +0.1969, -0.4072),
        (+1.4181, +0.1955, -0.4898), (+1.4187, +0.1752, -0.4850), (+1.4191, +0.1894, -0.4778),
        (+1.4189, +0.1816, -0.4845), (+1.4185, +0.1787, -0.4840), (+1.4191, +0.1765, -0.4741),
        (+1.4188, +0.1795, -0.4786), (+1.4197, +0.1927, -0.4622), (+1.4196, +0.2510, -0.4311),
        (+1.4288, +0.2485, -0.2573), (+1.4171, +0.2642, -0.3578), (+1.4192, +0.2392, -0.3969),
    ),
    (  # tau_i bin 11 [+40, +60) deg, n=2266
        (+1.4173, +0.3450, -0.3999), (+1.4213, +0.2952, -0.3249), (+1.4238, +0.3517, -0.3267),
        (+1.4192, +0.3588, -0.3515), (+1.4187, +0.3210, -0.3830), (+1.4148, +0.2991, -0.4359),
        (+1.4181, +0.2243, -0.4621), (+1.4187, +0.2265, -0.4443), (+1.4198, +0.2362, -0.4576),
        (+1.4192, +0.2371, -0.4517), (+1.4199, +0.2337, -0.4470), (+1.4207, +0.2255, -0.4296),
        (+1.4189, +0.2123, -0.4591), (+1.4211, +0.2696, -0.4050), (+1.4198, +0.2486, -0.3243),
        (+1.4203, +0.3942, -0.3769), (+1.4199, +0.3140, -0.2031), (+1.4127, +0.3144, -0.3039),
    ),
    (  # tau_i bin 12 [+60, +80) deg, n=1479
        (+1.4168, +0.4558, -0.1766), (+1.4227, +0.4152, -0.0722), (+1.4174, +0.4324, -0.1981),
        (+1.4213, +0.4613, -0.1965), (+1.4215, +0.3291, -0.2106), (+1.4238, +0.2948, -0.3829),
        (+1.4176, +0.2898, -0.3776), (+1.4201, +0.2807, -0.4006), (+1.4209, +0.3024, -0.4028),
        (+1.4201, +0.3139, -0.3978), (+1.4205, +0.2803, -0.4146), (+1.4221, +0.2815, -0.4209),
        (+1.4234, +0.3198, -0.3844), (+1.4213, +0.2716, -0.3958), (+1.4182, +0.2581, -0.2721),
        (+1.4299, +0.4364, -0.2651), (+1.4245, +0.1188, -0.1241), (+1.4247, +0.3449, -0.0460),
    ),
    (  # tau_i bin 13 [+80, +100) deg, n=675
        (+1.4342, +0.4589, -0.1962), (+1.4255, +0.1075, +0.0035), (+1.4263, +0.4370, -0.0155),
        (+1.4249, +0.3344, -0.0474), (+1.4129, -0.0392, +0.2784), (+1.4194, +0.3547, -0.1334),
        (+1.4137, +0.4164, -0.2724), (+1.4210, +0.3752, -0.3191), (+1.4196, +0.3698, -0.3111),
        (+1.4234, +0.3740, -0.3094), (+1.4225, +0.3856, -0.2850), (+1.4209, +0.3674, -0.3263),
        (+1.4186, +0.3689, -0.3487), (+1.4244, +0.3705, -0.2587), (+1.4193, +0.4247, -0.2255),
        (+1.4308, +0.0650, +0.0613), (+1.4211, -0.0306, +0.3399), (+1.3218, +0.2293, -0.0756),
    ),
    (  # tau_i bin 14 [+100, +120) deg, n=291
        (+1.4216, +0.3166, -0.0690), (+1.4324, +0.1485, +0.2407), (+1.4294, -0.1580, +0.4251),
        (+1.4266, +0.0375, +0.3330), (+1.4288, -0.0211, -0.0115), (+1.4390, +0.3526, -0.0712),
        (+1.4173, +0.2814, -0.0910), (+1.4200, +0.3783, -0.1589), (+1.4186, +0.3407, -0.1135),
        (+1.4210, +0.3909, -0.0945), (+1.4213, +0.4283, -0.1454), (+1.4239, +0.4176, -0.1364),
        (+1.4161, +0.3737, -0.2557), (+1.4273, +0.2367, -0.0234), (+1.4229, +0.0931, +0.0579),
        (+1.4219, -0.1772, +0.2474), (+1.4287, -0.0478, +0.5244), (+1.4257, +0.1182, +0.4407),
    ),
    (  # tau_i bin 15 [+120, +140) deg, n=146
        (+1.4184, -0.2343, +0.4827), (+1.4314, -0.2931, +0.4411), (+1.4299, -0.2671, +0.4484),
        (+1.4250, -0.1467, +0.3772), (+1.4274, +0.0561, +0.3976), (+1.4032, +0.0382, +0.4343),
        (+1.4306, +0.1938, +0.3564), (+1.4175, +0.1429, +0.1825), (+1.4158, +0.2619, +0.1526),
        (+1.4278, +0.2414, +0.1546), (+1.4136, +0.1241, +0.3638), (+1.4160, +0.1931, +0.3572),
        (+1.4234, +0.0376, +0.3534), (+1.4201, -0.0966, +0.3039), (+1.4279, -0.1248, +0.4820),
        (+1.4221, -0.1333, +0.5274), (+1.4158, -0.2436, +0.5021), (+1.4272, -0.1857, +0.4979),
    ),
    (  # tau_i bin 16 [+140, +160) deg, n=201
        (+1.4227, -0.2908, +0.3384), (+1.4280, -0.1538, +0.4963), (+1.4143, -0.3074, +0.4488),
        (+1.4237, -0.2312, +0.4481), (+1.4283, -0.1791, +0.3530), (+1.4289, +0.0521, +0.3890),
        (+1.4186, +0.0432, +0.2555), (+1.4232, -0.1397, +0.4713), (+1.4172, -0.1136, +0.5278),
        (+1.4140, -0.0719, +0.4556), (+1.4204, -0.1166, +0.4234), (+1.4229, -0.0962, +0.4043),
        (+1.4232, -0.1107, +0.3915), (+1.4171, -0.1691, +0.4627), (+1.4184, -0.3030, +0.3068),
        (+1.2839, -0.2581, +0.0466), (+1.4337, -0.1430, +0.4678), (+1.4258, -0.1822, +0.4959),
    ),
    (  # tau_i bin 17 [+160, +180) deg, n=323
        (+1.4042, -0.3965, +0.4096), (+1.4218, -0.2215, +0.4790), (+1.4132, -0.1364, +0.5419),
        (+1.4244, -0.2739, +0.2710), (+1.4295, -0.2117, +0.4754), (+1.4272, -0.2249, +0.3180),
        (+1.4242, -0.1231, +0.3391), (+1.4226, -0.1597, +0.3077), (+1.4202, -0.1784, +0.3843),
        (+1.4193, -0.1817, +0.4267), (+1.4217, -0.1955, +0.4471), (+1.4209, -0.2114, +0.4588),
        (+1.4221, -0.2141, +0.4527), (+1.4150, -0.2326, +0.4801), (+1.4202, -0.2616, +0.4657),
        (+1.4099, -0.3243, +0.0211), (+1.4317, -0.1417, +0.1351), (+1.4175, -0.2175, +0.5154),
    ),
)
_N_BY_BIN_2D: Final[tuple[tuple[tuple[float, float, float], ...], ...]] = (
    (  # tau_i bin 0 [-180, -160) deg, n=419
        (+2.4115, -0.0014, +0.1694), (+2.4683, +0.3988, -0.0784), (+2.4012, +0.2641, -0.2404),
        (+2.4098, +0.1724, -0.1611), (+2.3912, +0.2183, -0.2431), (+2.3946, +0.1428, -0.2568),
        (+2.4079, +0.1838, -0.2725), (+2.4045, +0.1239, -0.2300), (+2.4108, +0.1549, -0.2027),
        (+2.4045, +0.1283, -0.3022), (+2.3992, +0.1706, -0.2635), (+2.4099, +0.1362, -0.2852),
        (+2.4030, +0.1424, -0.2483), (+2.4175, +0.2086, -0.2683), (+2.3301, +0.1401, +0.0260),
        (+2.3519, +0.2209, -0.3536), (+2.4008, -0.1571, +0.3131), (+2.4047, +0.1532, -0.2535),
    ),
    (  # tau_i bin 1 [-160, -140) deg, n=608
        (+2.3971, +0.1719, -0.0247), (+2.4222, +0.3249, -0.0564), (+2.4074, +0.2120, -0.0312),
        (+2.4037, +0.1353, -0.0884), (+2.4187, +0.2310, -0.2132), (+2.4009, +0.2274, -0.2661),
        (+2.4104, +0.2002, -0.2677), (+2.4081, +0.1514, -0.2361), (+2.4086, +0.1975, -0.2869),
        (+2.4045, +0.1672, -0.2419), (+2.4050, +0.1833, -0.2727), (+2.4088, +0.1849, -0.2722),
        (+2.4044, +0.1869, -0.2586), (+2.4024, +0.1809, -0.2154), (+2.4021, +0.0884, -0.0585),
        (+2.4065, +0.1816, -0.2439), (+2.3983, +0.0390, +0.1533), (+2.3905, +0.3308, -0.1418),
    ),
    (  # tau_i bin 2 [-140, -120) deg, n=783
        (+2.3976, +0.2600, -0.1049), (+2.3866, +0.2770, -0.0525), (+2.3979, +0.2018, -0.0634),
        (+2.4139, +0.2593, -0.1401), (+2.4059, +0.2533, -0.2173), (+2.3993, +0.2136, -0.1903),
        (+2.3996, +0.2046, -0.2141), (+2.4061, +0.2022, -0.2355), (+2.4028, +0.2036, -0.2754),
        (+2.4039, +0.2136, -0.2560), (+2.4065, +0.1941, -0.2775), (+2.4020, +0.2012, -0.2587),
        (+2.3999, +0.1960, -0.2551), (+2.4081, +0.1829, -0.1720), (+2.4295, -0.0536, +0.1012),
        (+2.4083, +0.2345, -0.1632), (+2.3897, +0.0018, +0.3515), (+2.4054, +0.1621, -0.0500),
    ),
    (  # tau_i bin 3 [-120, -100) deg, n=851
        (+2.4100, +0.1799, +0.0264), (+2.4059, +0.2323, +0.0139), (+2.4031, +0.2639, +0.0226),
        (+2.4028, +0.1358, +0.0332), (+2.4128, +0.1359, -0.0493), (+2.4044, +0.1937, -0.1762),
        (+2.4044, +0.2204, -0.1944), (+2.4089, +0.2130, -0.1898), (+2.4004, +0.2325, -0.2369),
        (+2.4065, +0.2473, -0.2168), (+2.4020, +0.2368, -0.2154), (+2.4101, +0.2352, -0.2190),
        (+2.4006, +0.1959, -0.1553), (+2.4051, +0.1769, -0.0762), (+2.3834, +0.2960, -0.0735),
        (+2.3990, +0.0042, +0.1345), (+2.4251, +0.1193, +0.1114), (+2.3969, +0.1319, +0.1760),
    ),
    (  # tau_i bin 4 [-100, -80) deg, n=655
        (+2.3895, +0.0918, +0.2794), (+2.3976, +0.0808, +0.2166), (+2.4081, +0.1467, +0.2557),
        (+2.3967, +0.1241, +0.2120), (+2.4036, +0.1381, +0.0845), (+2.4020, +0.2363, +0.0362),
        (+2.4092, +0.2403, -0.0150), (+2.4111, +0.2721, -0.1267), (+2.4056, +0.2525, -0.1107),
        (+2.4078, +0.2750, -0.1254), (+2.4029, +0.2581, -0.1382), (+2.4143, +0.2417, -0.1247),
        (+2.4070, +0.2057, -0.1234), (+2.3922, -0.0414, +0.1858), (+2.4065, +0.0880, +0.1499),
        (+2.4025, -0.0560, +0.3022), (+2.3975, -0.0481, +0.3065), (+2.4155, +0.0120, +0.2789),
    ),
    (  # tau_i bin 5 [-80, -60) deg, n=555
        (+2.3997, -0.0053, +0.2664), (+2.3934, +0.0535, +0.2914), (+2.3998, +0.0405, +0.3105),
        (+2.3945, +0.0471, +0.2876), (+2.4018, +0.1022, +0.2336), (+2.3921, +0.0337, +0.2610),
        (+2.4041, +0.1257, +0.1052), (+2.4032, +0.2211, +0.0754), (+2.4105, +0.2198, +0.0583),
        (+2.4050, +0.2185, +0.0940), (+2.4058, +0.2249, +0.0883), (+2.3991, +0.2068, +0.0498),
        (+2.4086, +0.1681, +0.1459), (+2.3978, +0.0922, +0.2546), (+2.4105, +0.0511, +0.2541),
        (+2.4112, +0.1128, +0.2608), (+2.3918, -0.0619, +0.3192), (+2.4150, +0.0183, +0.3067),
    ),
    (  # tau_i bin 6 [-60, -40) deg, n=942
        (+2.4041, -0.0110, +0.2972), (+2.3998, -0.0316, +0.3089), (+2.3973, -0.0090, +0.3094),
        (+2.3981, -0.0174, +0.3133), (+2.4012, +0.0296, +0.3001), (+2.3967, +0.0597, +0.2787),
        (+2.4001, +0.1000, +0.2747), (+2.4017, +0.0745, +0.2623), (+2.3963, +0.0981, +0.2598),
        (+2.3960, +0.1004, +0.2352), (+2.4005, +0.1100, +0.2460), (+2.3982, +0.1324, +0.2630),
        (+2.3926, +0.0051, +0.2841), (+2.3976, +0.0123, +0.2842), (+2.3859, -0.0558, +0.3144),
        (+2.3958, -0.0158, +0.2646), (+2.4140, +0.0337, +0.3294), (+2.4021, -0.0572, +0.3296),
    ),
    (  # tau_i bin 7 [-40, -20) deg, n=1467
        (+2.3939, -0.0820, +0.3024), (+2.3946, -0.0562, +0.3148), (+2.3944, -0.0833, +0.2996),
        (+2.3974, -0.0386, +0.3121), (+2.3974, +0.0056, +0.3141), (+2.3955, -0.0149, +0.3108),
        (+2.3998, -0.0052, +0.3063), (+2.3971, +0.0443, +0.3086), (+2.4006, +0.0199, +0.3025),
        (+2.3974, +0.0240, +0.2984), (+2.3979, +0.0043, +0.3155), (+2.4000, +0.0183, +0.3122),
        (+2.4015, -0.0039, +0.3074), (+2.3980, -0.0261, +0.3058), (+2.4098, -0.0170, +0.2102),
        (+2.4061, -0.0636, +0.3075), (+2.4008, -0.0438, +0.2334), (+2.3966, -0.0622, +0.3068),
    ),
    (  # tau_i bin 8 [-20, +0) deg, n=2130
        (+2.3980, -0.1139, +0.2875), (+2.3980, -0.1058, +0.2943), (+2.3967, -0.1018, +0.3034),
        (+2.3967, -0.1071, +0.3088), (+2.4006, -0.0965, +0.3161), (+2.3986, -0.0487, +0.3010),
        (+2.3948, -0.0407, +0.3155), (+2.3970, -0.0319, +0.3207), (+2.3929, -0.0299, +0.3159),
        (+2.3962, -0.0314, +0.3185), (+2.3980, -0.0325, +0.3154), (+2.3970, -0.0134, +0.3156),
        (+2.3940, -0.0512, +0.3126), (+2.3982, -0.0324, +0.3128), (+2.4031, -0.0970, +0.2862),
        (+2.3961, -0.0881, +0.2683), (+2.3999, -0.0748, +0.2619), (+2.3927, -0.1169, +0.3110),
    ),
    (  # tau_i bin 9 [+0, +20) deg, n=2598
        (+2.4006, -0.1865, +0.2574), (+2.3936, -0.1710, +0.2592), (+2.3977, -0.1799, +0.2702),
        (+2.4021, -0.1697, +0.2692), (+2.4007, -0.1333, +0.2799), (+2.4001, -0.1315, +0.3042),
        (+2.3967, -0.0882, +0.3010), (+2.3982, -0.0916, +0.3159), (+2.4008, -0.0729, +0.3136),
        (+2.3984, -0.0660, +0.3117), (+2.3958, -0.0714, +0.3203), (+2.3965, -0.0777, +0.3097),
        (+2.3953, -0.0944, +0.3007), (+2.3962, -0.0774, +0.2787), (+2.4040, -0.0939, +0.2658),
        (+2.4016, -0.1610, +0.2257), (+2.3952, -0.0951, +0.2947), (+2.3940, -0.1613, +0.2750),
    ),
    (  # tau_i bin 10 [+20, +40) deg, n=2813
        (+2.4016, -0.2274, +0.2125), (+2.3970, -0.2266, +0.2317), (+2.3970, -0.2200, +0.2427),
        (+2.3960, -0.1952, +0.2599), (+2.3959, -0.1935, +0.2389), (+2.3943, -0.1383, +0.2480),
        (+2.3999, -0.1268, +0.2964), (+2.3972, -0.1147, +0.3032), (+2.3968, -0.1204, +0.2983),
        (+2.3988, -0.1083, +0.3047), (+2.3983, -0.1114, +0.2981), (+2.4001, -0.1090, +0.3009),
        (+2.3969, -0.1189, +0.3000), (+2.3996, -0.1090, +0.2849), (+2.3959, -0.1288, +0.2824),
        (+2.4033, -0.1778, +0.2013), (+2.3895, -0.1810, +0.2075), (+2.3946, -0.2118, +0.2244),
    ),
    (  # tau_i bin 11 [+40, +60) deg, n=2266
        (+2.3914, -0.2556, +0.2015), (+2.3982, -0.2425, +0.1845), (+2.4000, -0.2415, +0.2143),
        (+2.3891, -0.2604, +0.1948), (+2.4048, -0.1922, +0.2397), (+2.3923, -0.1927, +0.2637),
        (+2.3988, -0.1760, +0.2638), (+2.3992, -0.1628, +0.2757), (+2.3973, -0.1685, +0.2825),
        (+2.3977, -0.1749, +0.2720), (+2.3977, -0.1564, +0.2791), (+2.4007, -0.1603, +0.2673),
        (+2.3988, -0.1527, +0.2760), (+2.4012, -0.1868, +0.2536), (+2.4011, -0.1727, +0.2360),
        (+2.3976, -0.2215, +0.2489), (+2.4105, -0.2225, +0.1468), (+2.3854, -0.2469, +0.1338),
    ),
    (  # tau_i bin 12 [+60, +80) deg, n=1479
        (+2.4051, -0.2836, +0.1215), (+2.4064, -0.2774, +0.0161), (+2.3989, -0.2691, +0.1039),
        (+2.4061, -0.2708, +0.1390), (+2.4095, -0.2512, +0.1208), (+2.4070, -0.2153, +0.2318),
        (+2.3962, -0.2063, +0.2074), (+2.4025, -0.2029, +0.2376), (+2.3981, -0.2155, +0.2430),
        (+2.3994, -0.2056, +0.2482), (+2.3988, -0.2078, +0.2506), (+2.3984, -0.2030, +0.2651),
        (+2.4034, -0.2154, +0.2491), (+2.4034, -0.1940, +0.2270), (+2.4046, -0.1702, +0.1271),
        (+2.3915, -0.3816, +0.0558), (+2.4130, -0.1260, +0.0788), (+2.3997, -0.2728, -0.0102),
    ),
    (  # tau_i bin 13 [+80, +100) deg, n=675
        (+2.4055, -0.3306, +0.1470), (+2.4056, -0.1282, -0.0390), (+2.4060, -0.3020, -0.0116),
        (+2.4139, -0.2090, -0.0430), (+2.4047, +0.0091, -0.1732), (+2.3978, -0.2053, +0.1085),
        (+2.3906, -0.2819, +0.1379), (+2.3996, -0.2706, +0.1761), (+2.3966, -0.2641, +0.1772),
        (+2.4053, -0.2546, +0.1979), (+2.4047, -0.2673, +0.1903), (+2.4021, -0.2516, +0.2021),
        (+2.3995, -0.2339, +0.2084), (+2.4048, -0.2656, +0.1986), (+2.3894, -0.2571, +0.1773),
        (+2.3903, -0.0355, -0.0584), (+2.3829, -0.0464, -0.2689), (+2.3605, -0.2030, -0.1923),
    ),
    (  # tau_i bin 14 [+100, +120) deg, n=291
        (+2.4019, -0.2071, +0.0242), (+2.4238, -0.1269, -0.2241), (+2.3981, +0.1495, -0.2667),
        (+2.4105, +0.0046, -0.1985), (+2.3983, -0.0920, -0.0254), (+2.4051, -0.2928, +0.0038),
        (+2.4051, -0.1688, -0.0020), (+2.4025, -0.2263, +0.0842), (+2.4008, -0.2443, +0.0436),
        (+2.4003, -0.2565, +0.0357), (+2.4023, -0.2813, +0.0859), (+2.3973, -0.2838, +0.1012),
        (+2.3916, -0.2498, +0.1302), (+2.4161, -0.1678, -0.0108), (+2.4058, -0.0189, -0.0725),
        (+2.3997, +0.1843, -0.1709), (+2.3982, +0.1140, -0.3481), (+2.4141, -0.0412, -0.3274),
    ),
    (  # tau_i bin 15 [+120, +140) deg, n=146
        (+2.4230, +0.1566, -0.2933), (+2.4027, +0.2324, -0.2821), (+2.4085, +0.1344, -0.3443),
        (+2.4114, +0.1241, -0.2418), (+2.4263, -0.0326, -0.2485), (+2.3962, -0.0419, -0.2377),
        (+2.4259, -0.1087, -0.2352), (+2.3973, -0.1316, -0.1085), (+2.4050, -0.1492, -0.0594),
        (+2.4160, -0.1832, -0.0983), (+2.3896, -0.0843, -0.1901), (+2.3894, -0.1013, -0.2261),
        (+2.3957, -0.0307, -0.2417), (+2.3946, +0.0373, -0.1789), (+2.3952, +0.0530, -0.3523),
        (+2.4099, -0.0113, -0.3261), (+2.3930, +0.1832, -0.2844), (+2.3969, +0.1494, -0.3190),
    ),
    (  # tau_i bin 16 [+140, +160) deg, n=201
        (+2.4128, +0.2214, -0.2404), (+2.4191, +0.1978, -0.2714), (+2.3838, +0.1773, -0.2839),
        (+2.4124, +0.0389, -0.3285), (+2.4023, +0.1973, -0.2081), (+2.4001, -0.0056, -0.2827),
        (+2.4049, -0.0577, -0.2023), (+2.4085, +0.0804, -0.3140), (+2.4099, +0.1326, -0.2744),
        (+2.3976, +0.0645, -0.2605), (+2.4041, +0.0851, -0.2672), (+2.4105, +0.0542, -0.2612),
        (+2.4075, +0.0737, -0.2315), (+2.4021, +0.0916, -0.2671), (+2.4000, +0.1613, -0.1895),
        (+2.3627, -0.0046, -0.5066), (+2.4225, +0.1956, -0.3303), (+2.4165, +0.2116, -0.2994),
    ),
    (  # tau_i bin 17 [+160, +180) deg, n=323
        (+2.3891, +0.1917, -0.2564), (+2.3864, +0.2142, -0.3311), (+2.3855, +0.2148, -0.2855),
        (+2.4125, +0.2316, -0.1838), (+2.4186, +0.2237, -0.2814), (+2.4076, +0.1958, -0.2120),
        (+2.3999, +0.0554, -0.2310), (+2.4072, +0.1033, -0.1917), (+2.4029, +0.0913, -0.2597),
        (+2.4053, +0.1264, -0.2530), (+2.4070, +0.1246, -0.2859), (+2.4098, +0.1612, -0.2711),
        (+2.4080, +0.1429, -0.2834), (+2.4018, +0.1507, -0.2769), (+2.3984, +0.1938, -0.2564),
        (+2.3978, +0.1885, -0.0414), (+2.4316, +0.2186, -0.0303), (+2.3599, +0.0921, -0.3643),
    ),
)

#: The per-bin tables as arrays, for the vectorized placement path. Identical values to the tuples
#: above; the NumPy form lets a whole chain's bins be gathered in one indexing operation.
_C_TABLE: Final[np.ndarray] = np.array(_C_BY_BIN, dtype=np.float64)
_N_TABLE: Final[np.ndarray] = np.array(_N_BY_BIN, dtype=np.float64)

#: The 2D table as an array, for the vectorized placement path (see _C_TABLE / _N_TABLE).
_C_TABLE_2D: Final[np.ndarray] = np.array(_C_BY_BIN_2D, dtype=np.float64)
_N_TABLE_2D: Final[np.ndarray] = np.array(_N_BY_BIN_2D, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class CABackboneResult:
    """Backbone atoms placed for a CA-only trace.

    Attributes
    ----------
    n_xyz, c_xyz, o_xyz
        ``(n_residues, 3)`` coordinates. Every row is finite: the two atoms no peptide unit covers
        (N of the first residue, C and O of the last) are placed by ideal internal geometry.
    source
        Per peptide unit, which predictor was used: ``"table"``, ``"marginal"`` or ``"forward"``.
        Length ``n_residues - 1``.
    notes
        What happened, for a caller that wants to report it.
    """

    n_xyz: np.ndarray
    c_xyz: np.ndarray
    o_xyz: np.ndarray
    source: tuple[str, ...]
    notes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SeamStrain:
    """A seam peptide bond that could not be satisfied, so was left long and reported.

    Where a rebuilt region meets a residue DODO did not rebuild, the peptide bond crosses from an
    atom this pass placed to one that was already there. When the fixed neighbour -- part of a
    folded domain moved rigidly in step 3 -- sits farther than the ``2.854 A`` a peptide unit can
    bridge, no placement satisfies both bonds, so the atom is aimed as close as the residue's own
    N-CA-C angle allows and the bond is left long. This records that it happened, on the residue it
    happened to, so a caller can surface it rather than leaving the validator to rediscover it.

    Attributes
    ----------
    residue
        0-based index of the rebuilt boundary residue whose placed atom was strained.
    side
        ``"C"`` -- this residue's carbon reaching the next residue's nitrogen -- or ``"N"`` -- this
        residue's nitrogen reaching the previous residue's carbon.
    bond_length
        The resulting peptide-bond length, in Angstroms; ideal is
        :data:`~dodo.constants.C_N_PEPTIDE_BOND_LENGTH` (1.329).
    """

    residue: int
    side: str
    bond_length: float


@dataclass(frozen=True, slots=True)
class BackboneResult:
    """What :func:`add_backbone_to_rebuilt` produced: the structure and what it could not close.

    Attributes
    ----------
    structure
        The rebuilt structure with N, C and O placed on every rebuilt region; folded domains are
        returned exactly as they arrived.
    strained_seams
        Every seam left with a long peptide bond, in residue order. Empty when every seam closed --
        always the case for a from-sequence build, which has no folded domains to seam against.
    """

    structure: Structure
    strained_seams: tuple[SeamStrain, ...] = ()


def _unit(v: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(v))
    if norm < 1e-9:
        raise GeometryError("Cannot normalize a zero-length vector while placing a backbone.")
    return v / norm


def _frame(along: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, bool]:
    """Right-handed frame: x along ``along``, z perpendicular to the along/reference plane.

    Returns the frame and whether it had to be invented. The second element matters: three
    collinear alpha carbons define no plane, and a caller that cannot tell the difference between
    a measured peptide plane and an arbitrary one will report the arbitrary one as though it were
    determined by the trace.
    """
    x = _unit(along)
    normal = np.cross(reference, along)
    degenerate = float(np.linalg.norm(normal)) < 1e-6
    if degenerate:
        # Any perpendicular will do; the peptide plane is genuinely undetermined here, and the
        # caller passes that on through `notes`.
        fallback = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(fallback, x))) > 0.9:
            fallback = np.array([0.0, 1.0, 0.0])
        normal = np.cross(fallback, along)
    z = _unit(normal)
    return np.stack([x, np.cross(z, x), z]), degenerate


def _pseudo_dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    axis = _unit(b1)
    v = b0 - float(np.dot(b0, axis)) * axis
    w = b2 - float(np.dot(b2, axis)) * axis
    return float(np.degrees(np.arctan2(float(np.dot(np.cross(axis, v), w)), float(np.dot(v, w)))))


def _rotate(v: np.ndarray, axis: np.ndarray, radians: float) -> np.ndarray:
    axis = _unit(axis)
    cos, sin = float(np.cos(radians)), float(np.sin(radians))
    rotated: np.ndarray = (
        v * cos + np.cross(axis, v) * sin + axis * float(np.dot(axis, v)) * (1.0 - cos)
    )
    return rotated


def _on_two_spheres(
    centre_a: np.ndarray,
    radius_a: float,
    centre_b: np.ndarray,
    radius_b: float,
    toward: np.ndarray,
) -> np.ndarray | None:
    """Point at exact distances from two centres, on that circle, nearest ``toward``.

    Two spheres meet in a circle, so demanding both bonds leaves one free parameter -- which the
    table's prediction then resolves by picking the nearest point on it. Returns ``None`` when the
    spheres do not intersect, which cannot happen for a physically plausible CA-CA separation but
    is checked rather than assumed.
    """
    separation = float(np.linalg.norm(centre_b - centre_a))
    if (
        separation < 1e-9
        or separation > radius_a + radius_b
        or separation < abs(radius_a - radius_b)
    ):
        return None
    axis = (centre_b - centre_a) / separation
    along = (radius_a**2 - radius_b**2 + separation**2) / (2.0 * separation)
    height_squared = radius_a**2 - along**2
    if height_squared < 0.0:
        return None
    centre = centre_a + along * axis
    spoke = toward - centre
    spoke = spoke - float(np.dot(spoke, axis)) * axis
    if float(np.linalg.norm(spoke)) < 1e-9:
        arbitrary = np.array([1.0, 0.0, 0.0])
        spoke = arbitrary - float(np.dot(arbitrary, axis)) * axis
    resolved: np.ndarray = centre + np.sqrt(height_squared) * _unit(spoke)
    return resolved


def _place_carbonyl_oxygen(ca: np.ndarray, c: np.ndarray, n_next: np.ndarray) -> np.ndarray:
    """O is fully determined by CA(i), C(i) and N(i+1); from the true three it lands to 0.013 A."""
    to_ca = _unit(ca - c)
    to_n = _unit(n_next - c)
    normal = np.cross(to_ca, to_n)
    if float(np.linalg.norm(normal)) < 1e-6:
        normal = np.array([0.0, 0.0, 1.0])
    placed: np.ndarray = c + C_O_BOND_LENGTH * _rotate(to_ca, normal, -np.radians(CA_C_O_ANGLE))
    return placed


_Z_AXIS: Final[np.ndarray] = np.array([0.0, 0.0, 1.0])
_Y_AXIS: Final[np.ndarray] = np.array([0.0, 1.0, 0.0])
_X_AXIS: Final[np.ndarray] = np.array([1.0, 0.0, 0.0])


def _unit_rows(v: np.ndarray) -> np.ndarray:
    """Normalize along the last axis, flooring the divisor so no row can divide by zero.

    The scalar :func:`_unit` raises on a zero-length vector; here the only genuinely zero row would
    be a zero virtual bond, which :func:`_backbone_atoms` rejects up front, so the floor never bites
    on valid input and simply keeps the discarded branch of a :func:`numpy.where` finite.
    """
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    unit: np.ndarray = v / np.maximum(norm, 1e-12)
    return unit


def _rotate_rows(v: np.ndarray, axis: np.ndarray, radians: float) -> np.ndarray:
    """Rodrigues rotation of every row of ``v`` about the matching row of ``axis``."""
    axis = _unit_rows(axis)
    cos, sin = float(np.cos(radians)), float(np.sin(radians))
    dot = np.sum(axis * v, axis=-1, keepdims=True)
    rotated: np.ndarray = v * cos + np.cross(axis, v) * sin + axis * dot * (1.0 - cos)
    return rotated


def _frames_rows(
    along: np.ndarray, reference: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Batched :func:`_frame`: x along ``along``, z perpendicular to the along/reference plane.

    Returns the three basis vectors and a boolean mask of the units whose alpha carbons were
    collinear -- the same ``degenerate`` signal the scalar path reports.
    """
    x = _unit_rows(along)
    normal = np.cross(reference, along)
    degenerate = np.linalg.norm(normal, axis=-1) < 1e-6
    use_y = np.abs(x[..., 2]) > 0.9
    fallback = np.where(use_y[..., None], _Y_AXIS, _Z_AXIS)
    normal = np.where(degenerate[..., None], np.cross(fallback, along), normal)
    z = _unit_rows(normal)
    y = np.cross(z, x)
    return x, y, z, degenerate


def _dihedral_rows(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_pseudo_dihedral`, in degrees."""
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    axis = _unit_rows(b1)
    v = b0 - np.sum(b0 * axis, axis=-1, keepdims=True) * axis
    w = b2 - np.sum(b2 * axis, axis=-1, keepdims=True) * axis
    y = np.sum(np.cross(axis, v) * w, axis=-1)
    x = np.sum(v * w, axis=-1)
    degrees: np.ndarray = np.degrees(np.arctan2(y, x))
    return degrees


def _two_spheres_rows(
    centre_a: np.ndarray,
    radius_a: float,
    centre_b: np.ndarray,
    radius_b: float,
    toward: np.ndarray,
) -> np.ndarray:
    """Batched form of :func:`_on_two_spheres`.

    Where the spheres do not intersect -- the ``None`` case in the scalar version -- this returns a
    point one ``radius_a`` bond from ``centre_a`` toward ``toward``, exactly the fallback the scalar
    caller applies.
    """
    diff = centre_b - centre_a
    sep = np.linalg.norm(diff, axis=-1)
    sep_safe = np.maximum(sep, 1e-12)
    axis = diff / sep_safe[..., None]
    along = (radius_a**2 - radius_b**2 + sep**2) / (2.0 * sep_safe)
    height_squared = radius_a**2 - along**2
    valid = (
        (sep >= 1e-9)
        & (sep <= radius_a + radius_b)
        & (sep >= abs(radius_a - radius_b))
        & (height_squared >= 0.0)
    )
    centre = centre_a + along[..., None] * axis
    spoke = toward - centre
    spoke = spoke - np.sum(spoke * axis, axis=-1, keepdims=True) * axis
    arbitrary = _X_AXIS - np.sum(_X_AXIS * axis, axis=-1, keepdims=True) * axis
    spoke = np.where(np.linalg.norm(spoke, axis=-1, keepdims=True) < 1e-9, arbitrary, spoke)
    resolved = centre + np.sqrt(np.maximum(height_squared, 0.0))[..., None] * _unit_rows(spoke)
    fallback = centre_a + radius_a * _unit_rows(toward - centre_a)
    return np.where(valid[..., None], resolved, fallback)


def _carbonyl_rows(ca: np.ndarray, c: np.ndarray, n_next: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_place_carbonyl_oxygen`."""
    to_ca = _unit_rows(ca - c)
    to_n = _unit_rows(n_next - c)
    normal = np.cross(to_ca, to_n)
    normal = np.where(np.linalg.norm(normal, axis=-1, keepdims=True) < 1e-6, _Z_AXIS, normal)
    return c + C_O_BOND_LENGTH * _rotate_rows(to_ca, normal, -np.radians(CA_C_O_ANGLE))


def _terminal_nitrogen_rows(ca: np.ndarray, c: np.ndarray, ca_next: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_terminal_nitrogen`."""
    to_c = _unit_rows(c - ca)
    normal = np.cross(to_c, ca_next - ca)
    normal = np.where(np.linalg.norm(normal, axis=-1, keepdims=True) < 1e-6, _Z_AXIS, normal)
    return ca + N_CA_BOND_LENGTH * _rotate_rows(to_c, normal, np.radians(N_CA_C_ANGLE))


def _terminal_carbon_rows(ca: np.ndarray, n: np.ndarray, ca_prev: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_terminal_carbon`."""
    to_n = _unit_rows(n - ca)
    normal = np.cross(to_n, ca_prev - ca)
    normal = np.where(np.linalg.norm(normal, axis=-1, keepdims=True) < 1e-6, _Z_AXIS, normal)
    return ca + CA_C_BOND_LENGTH * _rotate_rows(to_n, normal, -np.radians(N_CA_C_ANGLE))


def _terminal_oxygen_rows(ca: np.ndarray, c: np.ndarray, ca_prev: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_terminal_oxygen`."""
    to_ca = _unit_rows(ca - c)
    normal = np.cross(to_ca, ca_prev - c)
    normal = np.where(np.linalg.norm(normal, axis=-1, keepdims=True) < 1e-6, _Z_AXIS, normal)
    return c + C_O_BOND_LENGTH * _rotate_rows(to_ca, normal, -np.radians(CA_C_O_ANGLE))


def _backbone_atoms(
    ca: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    """Place N, C and O for every peptide unit of one or many CA traces, with no per-residue loop.

    ``ca`` is ``(..., n_residues, 3)`` -- a single trace as ``(n, 3)`` or a batch of conformers as
    ``(B, n, 3)``. Every step runs over the residue axis (and any leading batch axes) at once: a
    peptide unit's placement depends only on the handful of alpha carbons around it, all known up
    front, so nothing here is sequential and there is nothing to loop over in Python.

    Returns the three atom arrays (same leading shape as ``ca``), the per-unit ``degenerate`` mask,
    and the per-unit predictor names. A unit's predictor is fixed by its position in the chain, so
    ``source`` is one tuple shared by every conformer in a batch.
    """
    n = ca.shape[-2]
    units = n - 1
    along = ca[..., 1:, :] - ca[..., :-1, :]
    if np.any(np.linalg.norm(along, axis=-1) < 1e-9):
        raise GeometryError(
            "Two consecutive alpha carbons coincide, giving a zero-length virtual bond."
        )

    # Reference axis fixing each frame's peptide plane. Past the first unit it is the incoming CA-CA
    # bond; the first unit has no predecessor, so it looks forward (or, for a two-residue chain with
    # nothing to look at, to a fixed axis).
    reference = np.empty_like(along)
    if units >= 2:
        reference[..., 1:, :] = ca[..., 1:units, :] - ca[..., : units - 1, :]
    # ca[1] - ca[2], NOT ca[2] - ca[1]. The forward marginal (_C/_N_FORWARD_MARGINAL) was MEASURED
    # in the frame this sign builds -- verified by re-deriving it in scripts/derive_peptide_table.py
    # -- and the frame's z is cross(reference, along), so the opposite sign flips z (and y) and
    # mis-places the first peptide unit of every region: measured C error 1.06 A against the 0.31 A
    # the marginal delivers in its own frame. Do not "simplify" this to ca[2] - ca[1].
    reference[..., 0, :] = ca[..., 1, :] - ca[..., 2, :] if n >= 3 else _Z_AXIS

    # Predictor per unit, fixed by position. An interior unit with a fifth alpha carbon (one beyond
    # the fourth) uses the 2D table, keyed on both adjacent pseudo-dihedrals; the last interior unit
    # has only four and falls back to the 1D table, keyed on tau_i alone; the chain's last unit uses
    # the trailing marginal and its first the leading one.
    index = np.arange(units)
    is_table_2d = (index >= 1) & (index <= n - 4)
    is_table_1d = (index == n - 3) & (index >= 1)
    is_marginal = (index == n - 2) & (n >= 3)

    local_c = np.broadcast_to(np.array(_C_FORWARD_MARGINAL), along.shape).copy()
    local_n = np.broadcast_to(np.array(_N_FORWARD_MARGINAL), along.shape).copy()
    if n >= 3:
        local_c[..., n - 2, :] = _C_MARGINAL
        local_n[..., n - 2, :] = _N_MARGINAL

    # Last interior unit (no fifth CA): the 1D table, keyed on tau_i over CA[i-1..i+2].
    units_1d = np.nonzero(is_table_1d)[0]
    if units_1d.size:
        tau = _dihedral_rows(
            ca[..., units_1d - 1, :],
            ca[..., units_1d, :],
            ca[..., units_1d + 1, :],
            ca[..., units_1d + 2, :],
        )
        bins = np.minimum(((tau + 180.0) // _BIN_WIDTH_DEGREES).astype(np.int64), len(_C_TABLE) - 1)
        local_c[..., units_1d, :] = _C_TABLE[bins]
        local_n[..., units_1d, :] = _N_TABLE[bins]

    # Interior units with a fifth CA: the 2D table, keyed on (tau_i, tau_{i+1}) -- tau_i over
    # CA[i-1..i+2] as in the 1D path, tau_{i+1} over CA[i..i+3].
    units_2d = np.nonzero(is_table_2d)[0]
    if units_2d.size:
        clip = _C_TABLE_2D.shape[0] - 1
        tau_i = _dihedral_rows(
            ca[..., units_2d - 1, :],
            ca[..., units_2d, :],
            ca[..., units_2d + 1, :],
            ca[..., units_2d + 2, :],
        )
        tau_j = _dihedral_rows(
            ca[..., units_2d, :],
            ca[..., units_2d + 1, :],
            ca[..., units_2d + 2, :],
            ca[..., units_2d + 3, :],
        )
        bins_i = np.minimum(((tau_i + 180.0) // _BIN_WIDTH_DEGREES).astype(np.int64), clip)
        bins_j = np.minimum(((tau_j + 180.0) // _BIN_WIDTH_DEGREES).astype(np.int64), clip)
        local_c[..., units_2d, :] = _C_TABLE_2D[bins_i, bins_j]
        local_n[..., units_2d, :] = _N_TABLE_2D[bins_i, bins_j]

    x, y, z, degenerate = _frames_rows(along, reference)
    apex = ca[..., :units, :]
    placed_c = apex + local_c[..., 0:1] * x + local_c[..., 1:2] * y + local_c[..., 2:3] * z
    placed_n = apex + local_n[..., 0:1] * x + local_n[..., 1:2] * y + local_n[..., 2:3] * z

    n_xyz = np.full(ca.shape, np.nan)
    c_xyz = np.full(ca.shape, np.nan)
    o_xyz = np.full(ca.shape, np.nan)

    # N belongs to residue i+1; anchor it to that residue's own alpha carbon so N-CA is exact.
    n_owner = ca[..., 1 : units + 1, :]
    n_xyz[..., 1 : units + 1, :] = n_owner + N_CA_BOND_LENGTH * _unit_rows(placed_n - n_owner)
    # C satisfies both its bonds; the table's prediction chooses where on the intersection circle.
    c_xyz[..., :units, :] = _two_spheres_rows(
        apex, CA_C_BOND_LENGTH, n_xyz[..., 1 : units + 1, :], C_N_PEPTIDE_BOND_LENGTH, placed_c
    )
    o_xyz[..., :units, :] = _carbonyl_rows(
        apex, c_xyz[..., :units, :], n_xyz[..., 1 : units + 1, :]
    )

    # The two atoms no peptide unit covers: N of the first residue, C and O of the last.
    n_xyz[..., 0, :] = _terminal_nitrogen_rows(ca[..., 0, :], c_xyz[..., 0, :], ca[..., 1, :])
    c_xyz[..., -1, :] = _terminal_carbon_rows(ca[..., -1, :], n_xyz[..., -1, :], ca[..., -2, :])
    o_xyz[..., -1, :] = _terminal_oxygen_rows(ca[..., -1, :], c_xyz[..., -1, :], ca[..., -2, :])

    sources = tuple(
        "table5"
        if is_table_2d[i]
        else "table"
        if is_table_1d[i]
        else "marginal"
        if is_marginal[i]
        else "forward"
        for i in range(units)
    )
    return n_xyz, c_xyz, o_xyz, degenerate, sources


def backbone_from_ca(ca_coords: np.ndarray) -> CABackboneResult:
    """Place N, C and O for every residue of a CA-only trace.

    Parameters
    ----------
    ca_coords
        ``(n_residues, 3)`` alpha-carbon coordinates of one continuous chain segment, in order.
        At least two residues.

    Returns
    -------
    CABackboneResult
        Backbone coordinates, all finite, with ideal bond lengths by construction.

    Raises
    ------
    ~dodo.exceptions.GeometryError
        If the input is not ``(n, 3)`` with ``n >= 2``, or contains a non-finite coordinate.
    """
    ca = np.asarray(ca_coords, dtype=np.float64)
    if ca.ndim != 2 or ca.shape[1] != 3:
        raise GeometryError(f"Expected CA coordinates shaped (n_residues, 3), got {ca.shape}.")
    if ca.shape[0] < 2:
        raise GeometryError(
            f"A backbone needs at least two alpha carbons to define one peptide unit, got "
            f"{ca.shape[0]}."
        )
    if not np.all(np.isfinite(ca)):
        raise GeometryError("CA coordinates contain a non-finite value.")

    n_res = ca.shape[0]
    n_xyz, c_xyz, o_xyz, degenerate, sources = _backbone_atoms(ca)
    collinear = int(np.count_nonzero(degenerate))

    notes = [
        f"Placed {n_res} residue(s) of backbone from alpha carbons: "
        f"{sources.count('table5')} peptide unit(s) from the five-CA (2D) table, "
        f"{sources.count('table')} from the four-CA (1D) table, "
        f"{sources.count('marginal')} from the trailing marginal, "
        f"{sources.count('forward')} from the leading marginal."
    ]
    if collinear:
        notes.append(f"{collinear} peptide unit(s) had collinear alpha carbons.")
    return CABackboneResult(
        n_xyz=n_xyz, c_xyz=c_xyz, o_xyz=o_xyz, source=tuple(sources), notes=tuple(notes)
    )


def _terminal_nitrogen(ca: np.ndarray, c: np.ndarray, ca_next: np.ndarray) -> np.ndarray:
    """N of the first residue, which no peptide unit covers.

    Fixed by the N-CA-C angle off the already-placed C. The remaining rotation about the CA-C axis
    is this residue's phi, which nothing in a CA trace constrains at a chain start, so it is
    resolved deterministically into the CA(0)/C(0)/CA(1) plane.
    """
    to_c = _unit(c - ca)
    reference = ca_next - ca
    normal = np.cross(to_c, reference)
    if float(np.linalg.norm(normal)) < 1e-6:
        normal = np.array([0.0, 0.0, 1.0])
    return ca + N_CA_BOND_LENGTH * _rotate(to_c, normal, np.radians(N_CA_C_ANGLE))


def _terminal_carbon(ca: np.ndarray, n: np.ndarray, ca_prev: np.ndarray) -> np.ndarray:
    """C of the last residue, the mirror case of :func:`_terminal_nitrogen`."""
    to_n = _unit(n - ca)
    reference = ca_prev - ca
    normal = np.cross(to_n, reference)
    if float(np.linalg.norm(normal)) < 1e-6:
        normal = np.array([0.0, 0.0, 1.0])
    return ca + CA_C_BOND_LENGTH * _rotate(to_n, normal, -np.radians(N_CA_C_ANGLE))


def _terminal_oxygen(ca: np.ndarray, c: np.ndarray, ca_prev: np.ndarray) -> np.ndarray:
    """O of the last residue, which has no following N to orient the carbonyl against.

    The obvious shortcut is to invent that nitrogen by continuing the chain straight out along
    CA->C, and this function exists because that is wrong in a way that is easy to miss: a notional
    N on the CA-C line is *collinear* with it, so the cross product defining the carbonyl plane
    degenerates to zero and the normal falls back to an arbitrary axis. On a 22-residue test
    sequence that put the final O 0.981 A from its own CA -- two unbonded atoms on top of each
    other, in a chain with no seams to blame.

    So the plane is taken from the previous alpha carbon instead, which is always available (this is
    only reached with at least two residues) and never collinear with CA-C in practice. Which side
    of that plane O falls on is this residue's psi, and nothing in a CA trace constrains psi at a
    free C-terminus, so the choice is arbitrary but consistent.
    """
    to_ca = _unit(ca - c)
    normal = np.cross(to_ca, ca_prev - c)
    if float(np.linalg.norm(normal)) < 1e-6:
        normal = np.array([0.0, 0.0, 1.0])
    placed: np.ndarray = c + C_O_BOND_LENGTH * _rotate(to_ca, normal, -np.radians(CA_C_O_ANGLE))
    return placed


def _existing_atom(structure: Structure, residue: int, name: str) -> np.ndarray | None:
    """Coordinates of a named atom already present in a residue, or None if it has none."""
    atoms = structure.atom_slice_for_residues(residue, residue + 1)
    for index in range(atoms.start, atoms.stop):
        if str(structure.atom_name[index]) == name:
            found: np.ndarray = structure.xyz[index]
            return found
    return None


#: Separation a seam placement must keep from the neighbouring residue's atoms, in Angstroms.
#:
#: MEASURED as the smallest value that stops the recurring seam defect. Real non-bonded backbone
#: contacts across a peptide junction run from about 2.7 A (a carbonyl O to the next residue's CB)
#: upward, so a 2.6 A floor rejects a collision without rejecting legitimate geometry. It is
#: checked rather than assumed -- see :func:`_place_on_cone` for the three successive
#: by-construction fixes this replaced.
SEAM_CLASH_DISTANCE: Final[float] = 2.6


def _place_on_cone(
    apex: np.ndarray,
    reference: np.ndarray,
    angle_degrees: float,
    length: float,
    *,
    prefer: np.ndarray | None = None,
    depart: np.ndarray | None = None,
    avoid: np.ndarray | None = None,
    clash: float = SEAM_CLASH_DISTANCE,
    samples: int = 48,
    hard_avoid: np.ndarray | None = None,
    hard_clash: float = 0.0,
    oxygen_toward: np.ndarray | None = None,
) -> np.ndarray:
    """Place an atom at ``length`` from ``apex`` and ``angle_degrees`` off ``apex``->``reference``.

    That leaves exactly one free parameter -- the azimuth about the reference axis -- and this picks
    it. Collisions are ruled out first, then ``prefer`` is approached and ``depart`` is escaped.

    ``hard_avoid``/``hard_clash`` is a second obstacle set with its own distance -- in practice a
    wide shell of every structure atom near the seam, held to the impossible-separation floor,
    where ``avoid`` stays the deliberately narrow neighbouring-residue set at the seam clash
    distance. They are separate because widening ``avoid`` itself would add soft gradients from
    atoms that were never a problem and shift seams that are currently fine. ``oxygen_toward``
    additionally scores, for each candidate carbon, the carbonyl oxygen that candidate would induce
    (:func:`_place_carbonyl_oxygen` toward that nitrogen) against ``hard_avoid`` -- because O is
    fully determined by the chosen C, an O collision can only be avoided by choosing a different C,
    which is exactly what this term makes the sweep do.

    This exists because the strained-seam fallback kept producing the same class of defect in a
    new place each time it was fixed by construction alone. Aiming the carbon straight at the
    neighbour's nitrogen made the three collinear and sent the carbonyl to an arbitrary axis
    (measured: O 0.6 A from its own alpha carbon on p300). Holding the N-CA-C angle instead fixed
    that, and let the oxygen land 0.975 A from the neighbour's nitrogen. Placing the oxygen trans
    to that nitrogen fixed *that*, and left it 0.96 A from the neighbour's CB. Every fix was
    locally right and blind to one more neighbour; the way to stop is to check instead of
    construct, which is what this does.

    The angle to ``reference`` is never traded away: it is the residue's own internal geometry, and
    damaging that is worse than a long bond to the next residue.
    """
    axis = _unit(reference - apex)
    fallback = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(fallback, axis))) > 0.9:
        fallback = np.array([0.0, 1.0, 0.0])
    first = _unit(fallback - float(np.dot(fallback, axis)) * axis)
    second = np.cross(axis, first)

    angle = np.radians(angle_degrees)
    azimuths = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    radial = np.cos(azimuths)[:, None] * first + np.sin(azimuths)[:, None] * second
    candidates = apex + length * (np.cos(angle) * axis + np.sin(angle) * radial)

    score = np.zeros(samples)
    if avoid is not None and avoid.shape[0]:
        gaps = np.linalg.norm(candidates[:, None, :] - avoid[None, :, :], axis=2)
        # A hard penalty, not a weighted term: a candidate that collides is not a candidate.
        score += 1e6 * np.square(np.clip(clash - gaps, 0.0, None)).sum(axis=1)
    if hard_avoid is not None and hard_avoid.shape[0]:
        gaps = np.linalg.norm(candidates[:, None, :] - hard_avoid[None, :, :], axis=2)
        score += 1e6 * np.square(np.clip(hard_clash - gaps, 0.0, None)).sum(axis=1)
        if oxygen_toward is not None:
            induced = np.array(
                [_place_carbonyl_oxygen(apex, c, oxygen_toward) for c in candidates]
            )
            gaps = np.linalg.norm(induced[:, None, :] - hard_avoid[None, :, :], axis=2)
            score += 1e6 * np.square(np.clip(hard_clash - gaps, 0.0, None)).sum(axis=1)
    if prefer is not None:
        score += np.linalg.norm(candidates - prefer, axis=1)
    if depart is not None:
        score -= np.linalg.norm(candidates - depart, axis=1)
    chosen: np.ndarray = candidates[int(np.argmin(score))]
    return chosen


def _seam_obstacles(structure: Structure, residues: Iterable[int]) -> np.ndarray:
    """Atoms of ``residues`` that a seam placement must not collide with.

    Deliberately narrow: the atoms that actually get hit are the neighbouring residue's own, and its
    side chain in particular -- an anchor's CB is the thing a carbonyl oxygen runs into.
    """
    collected: list[np.ndarray] = []
    for residue in residues:
        if not 0 <= residue < structure.n_residues:
            continue
        atoms = structure.atom_slice_for_residues(residue, residue + 1)
        for index in range(atoms.start, atoms.stop):
            collected.append(structure.xyz[index])
    if not collected:
        return np.empty((0, 3))
    return np.asarray(collected, dtype=np.float64)


def _would_overlap(point: np.ndarray, obstacles: np.ndarray, floor: float) -> bool:
    """Report whether ``point`` lands within ``floor`` of any obstacle -- no bond is that short.

    The exact two-sphere seam placement makes the peptide bond exact but is blind to every OTHER
    atom, so it can drop the placed atom onto one: measured on PTBP2, a rebuilt N landed 0.797 A
    from the folded residue's carbonyl oxygen. The atom the seam actually bonds to sits at the
    bond length, comfortably above ``floor`` (:data:`IMPOSSIBLE_SEPARATION`, 1.00 A), so this
    never rejects a good placement -- only one that would write something physically impossible,
    which is then handed to the obstacle-avoiding cone instead.
    """
    if obstacles.shape[0] == 0:
        return False
    return bool(np.min(np.linalg.norm(obstacles - point, axis=1)) < floor)


#: Radius, in A, of the shell of structure atoms a seam placement is vetoed against.
#:
#: DERIVED to be a superset. A seam carbon lies exactly ``CA_C_BOND_LENGTH`` (1.525 A) from its
#: anchor CA, the carbonyl oxygen it induces at most ~2.44 A from it, and a nitrogen exactly
#: ``N_CA_BOND_LENGTH`` (1.458 A) -- so every atom a seam can place lies within 2.44 A of the
#: anchor, and anything that could sit within the 1.00 A impossible floor of one lies within
#: 3.44 A. 4.0 covers that with margin. The narrow :func:`_seam_obstacles` set is NOT enough here:
#: measured on Q15642, the exact seam carbon's induced oxygen landed on a folded atom 56 residues
#: away in sequence -- nearby in space, invisible to a neighbouring-residues-only set.
_SEAM_VETO_RADIUS: Final[float] = 4.0


def _bond_angle(a: np.ndarray, apex: np.ndarray, b: np.ndarray) -> float:
    """Angle a-apex-b in degrees."""
    v1, v2 = a - apex, b - apex
    cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def _polish_coupled_clashes(
    placed: dict[int, dict[str, np.ndarray]],
    structure: Structure,
    generated: list[tuple[int, int]],
    *,
    rounds: int = 8,
    grid: int = 73,
    angle_window: tuple[float, float] = (80.0, 160.0),
) -> int:
    """Clear steric clashes coupled between two movable peptide units. Modifies ``placed``.

    The per-unit refinement (:func:`~dodo.construct.backbone_refine.refine_backbone`) moves one
    azimuth at a time, so a clash coupled between two movable units -- where clearing it by rotating
    either one alone makes another contact worse -- sits at a local minimum that more sweeps cannot
    escape (measured: 200 sweeps leaves the same clashes as 30, a wide window barely helps). For
    each residual clash this jointly searches the azimuths of the one or two *interior* units (both
    flanking alpha carbons rebuilt, so only rebuilt atoms move -- never a folded-domain or seam
    atom) that place the clashing atoms, keeping any combination that lowers the total van der Waals
    overlap without pushing a moved residue's N-CA-C angle out of ``angle_window``.

    Two units is enough: the residual clashes were measured to couple at most two movable units
    (p300's two clusters are {925,933} and {2190,2195}; the rest are a unit against a *fixed* CA),
    so a three-unit joint search has nothing to grip -- there is no three-unit cluster. What the
    borderline ones needed was azimuth *resolution*: at the old 15-degree grid (``grid=25``) the
    clearing orientation fell between samples, and 5 degrees (``grid=73``) clears one more p300
    contact (4 -> 3 introduced) with the angle distribution byte-identical and ~+1 s. The
    ``angle_window`` deliberately stays at ``(80, 160)`` rather than widening: (70, 170) clears one
    further contact but only by placing a 165.7-degree N-CA-C angle, trading a 0.05 A borderline
    clash for a badly distorted backbone angle. The two contacts that remain sit against fixed CAs
    that no azimuth can move.

    The overlap and its element limits are the validator's own
    (:func:`~dodo.validate.clashes.contact_limit` constants), and only pairs more than one residue
    apart are scored -- consecutive-residue contacts are backbone bonds or their 1-3/1-4 neighbours,
    and seam contacts belong to a different problem. Returns the number of joint moves applied.
    """
    from scipy.spatial import cKDTree

    from ..validate.clashes import (
        ALLOWED_OVERLAP,
        DEFAULT_VDW_RADIUS,
        HBOND_ALLOWED_OVERLAP,
        POLAR_ELEMENTS,
        VDW_RADII,
    )
    from .backbone_refine import (
        _angle,
        _angles_batch,
        _azimuth_frame,
        _place_oxygen,
        _place_oxygen_batch,
        _place_unit,
        _place_units_batch,
    )

    rebuilt: set[int] = set()
    for start, stop in generated:
        rebuilt.update(range(start, stop))
    ca = structure.ca_xyz
    n_res = structure.n_residues
    units = [
        u
        for u in range(n_res - 1)
        if u in rebuilt and (u + 1) in rebuilt and u in placed and (u + 1) in placed
    ]
    if not units:
        return 0
    unit_set = set(units)

    azimuth: dict[int, float] = {}
    for u in units:
        axis, first, second = _azimuth_frame(ca[u], ca[u + 1])
        offset = placed[u]["C"] - ca[u]
        radial = offset - float(np.dot(offset, axis)) * axis
        azimuth[u] = float(
            np.degrees(np.arctan2(float(np.dot(radial, second)), float(np.dot(radial, first))))
        )

    def limit(ea: str, eb: str) -> float:
        ra = VDW_RADII.get(ea, DEFAULT_VDW_RADIUS)
        rb = VDW_RADII.get(eb, DEFAULT_VDW_RADIUS)
        both_polar = ea in POLAR_ELEMENTS and eb in POLAR_ELEMENTS
        return ra + rb - (HBOND_ALLOWED_OVERLAP if both_polar else ALLOWED_OVERLAP)

    query_r = 2.0 * max(VDW_RADII.get(e, DEFAULT_VDW_RADIUS) for e in ("C", "N", "O", "S"))

    # The whole atom cloud: every fixed atom (folded residues all-atom + rebuilt alpha carbons --
    # exactly what structure.xyz holds while rebuilt regions are CA-only) plus the placed N/C/O.
    fixed_xyz = np.asarray(structure.xyz, dtype=np.float64)
    fixed_res = np.asarray(structure.residue_index)
    fixed_elem = np.array([str(e).upper() for e in structure.element], dtype=object)
    mov_keys = [
        (r, nm)
        for r in sorted(rebuilt)
        if r in placed
        for nm in ("N", "C", "O")
        if nm in placed[r]  # a seam atom may have been omitted as unsatisfiable
    ]
    row = {key: len(fixed_xyz) + i for i, key in enumerate(mov_keys)}
    all_xyz = np.vstack([fixed_xyz, np.array([placed[r][nm] for r, nm in mov_keys])])
    all_res = np.concatenate([fixed_res, np.array([r for r, _ in mov_keys], dtype=fixed_res.dtype)])
    all_elem = np.concatenate([fixed_elem, np.array([nm for _, nm in mov_keys], dtype=object)])

    def unit_of(index: int) -> int | None:
        if index < len(fixed_xyz):
            return None  # a fixed atom -- folded, or a rebuilt CA -- never moves
        r, nm = mov_keys[index - len(fixed_xyz)]
        u = r - 1 if nm == "N" else r
        return u if u in unit_set else None

    def place_into(u: int, azimuth_deg: float, target: np.ndarray) -> None:
        c, n_next = _place_unit(ca[u], ca[u + 1], azimuth_deg)
        target[row[(u, "C")]] = c
        target[row[(u + 1, "N")]] = n_next
        target[row[(u, "O")]] = _place_oxygen(ca[u], c, n_next)

    lo, hi = angle_window

    def angle_ok(u: int, positions: np.ndarray) -> bool:
        # Moving unit u sets residue u's C and residue u+1's N; check the N-CA-C at both, using the
        # N and C the move does not touch (placed by the neighbouring units, held fixed here). A
        # residue whose N or C was omitted at a seam has no angle to check, so it is skipped.
        for res in (u, u + 1):
            if (res, "N") in row and (res, "C") in row:
                angle = _angle(positions[row[(res, "N")]], ca[res], positions[row[(res, "C")]])
                if not lo <= angle <= hi:
                    return False
        return True

    def try_polish(unit_group: tuple[int, ...]) -> bool:
        # The three atoms each unit moves, each with the anchor and reach that bound where it can
        # go: C and O hang off CA(u), N off CA(u+1), and each travels on a circle no larger than its
        # bond(s), so a shell of ``query_r + reach`` around the anchor covers EVERY candidate
        # position the search will try. Ordering (C, N, O per unit) matches the old ``moved`` list.
        moved_specs = [
            spec
            for u in unit_group
            for spec in (
                (row[(u, "C")], ca[u], CA_C_BOND_LENGTH),
                (row[(u + 1, "N")], ca[u + 1], N_CA_BOND_LENGTH),
                (row[(u, "O")], ca[u], CA_C_BOND_LENGTH + C_O_BOND_LENGTH),
            )
        ]
        moved = [r for r, _, _ in moved_specs]
        static = np.ones(len(all_xyz), dtype=bool)
        static[moved] = False
        txyz, tres, telem = all_xyz[static], all_res[static], all_elem[static]
        static_tree = cKDTree(txyz)

        # Precompute each moved atom's static neighbours ONCE for the whole group search, not with a
        # ``query_ball_point`` per atom per candidate. The static cloud is fixed here and the shell
        # is a superset of any candidate's true neighbours; a neighbour past its own vdW limit adds
        # zero headroom, so scoring the superset gives the identical overlap sum. That per-candidate
        # tree query was measured at ~78% of --backbone wall time.
        nbr_xyz: dict[int, np.ndarray] = {}
        nbr_lim: dict[int, np.ndarray] = {}
        for r_idx, anchor, reach in moved_specs:
            residue = int(all_res[r_idx])
            element = str(all_elem[r_idx])
            found = np.asarray(
                static_tree.query_ball_point(anchor, query_r + reach), dtype=np.int64
            )
            if found.size:
                found = found[np.abs(tres[found] - residue) > 1]
            nbr_xyz[r_idx] = txyz[found]
            nbr_lim[r_idx] = np.array([limit(element, str(telem[j])) for j in found])

        # Overlap AMONG the group's own moved atoms -- excluded from the static cloud, but exactly
        # the clash a joint search is trying to clear. Their element pairs and residues are fixed,
        # so the vdW limit and within-one-residue exclusion are settled once here, not per score.
        inter_pairs: list[tuple[int, int, float]] = []
        for a_pos in range(len(moved)):
            for b_pos in range(a_pos + 1, len(moved)):
                ia, ib = moved[a_pos], moved[b_pos]
                if abs(int(all_res[ia]) - int(all_res[ib])) <= 1:
                    continue
                inter_pairs.append((ia, ib, limit(str(all_elem[ia]), str(all_elem[ib]))))

        def score(positions: np.ndarray) -> float:
            total = 0.0
            for r_idx in moved:
                pts = nbr_xyz[r_idx]
                if pts.shape[0]:
                    distance = np.sqrt(((pts - positions[r_idx]) ** 2).sum(axis=1))
                    headroom = nbr_lim[r_idx] - distance
                    total += float(headroom[headroom > 0.0].sum())
            for ia, ib, lim in inter_pairs:
                headroom = lim - float(np.linalg.norm(positions[ia] - positions[ib]))
                if headroom > 0.0:
                    total += headroom
            return total

        base = score(all_xyz)
        if base <= 1e-9:
            return False
        trial = all_xyz.copy()
        grid_offsets = np.linspace(-180.0, 180.0, grid)
        best = base
        best_choice: tuple[float, ...] | None = None

        u_rows_of = {u: (row[(u, "C")], row[(u + 1, "N")], row[(u, "O")]) for u in unit_group}

        def evaluate_last(u: int) -> tuple[float, float] | None:
            """Best ``(overlap, azimuth)`` for unit ``u`` over the whole grid at once.

            Every other unit of the group is held at its current ``trial`` placement. The whole
            grid is placed, angle-filtered and scored with array operations instead of a
            Python candidate loop -- which is what makes a fine grid affordable. It returns the same
            choice the scalar loop would: the total is ``const + var`` where ``const`` is the score
            of every atom this unit does not move (identical across candidates) and ``var`` is this
            unit's own contribution, so ``argmin`` over the valid candidates is the global minimum.
            """
            grid_az = azimuth[u] + grid_offsets
            c_b, n_b = _place_units_batch(ca[u], ca[u + 1], grid_az)
            o_b = _place_oxygen_batch(ca[u], c_b, n_b)

            valid = np.ones(grid_az.shape[0], dtype=bool)
            if (u, "N") in row and (u, "C") in row:
                ang = _angles_batch(trial[row[(u, "N")]], ca[u], c_b)
                valid &= (lo <= ang) & (ang <= hi)
            if (u + 1, "N") in row and (u + 1, "C") in row:
                ang = _angles_batch(n_b, ca[u + 1], trial[row[(u + 1, "C")]])
                valid &= (lo <= ang) & (ang <= hi)
            if not valid.any():
                return None

            moving = {u_rows_of[u][0]: c_b, u_rows_of[u][1]: n_b, u_rows_of[u][2]: o_b}
            var = np.zeros(grid_az.shape[0])
            for r_idx, pos in moving.items():
                pts = nbr_xyz[r_idx]
                if pts.shape[0]:
                    distance = np.sqrt(((pos[:, None, :] - pts[None, :, :]) ** 2).sum(-1))
                    headroom = nbr_lim[r_idx][None, :] - distance
                    var += np.where(headroom > 0.0, headroom, 0.0).sum(1)
            const = 0.0
            for r_idx in moved:
                if r_idx in moving:
                    continue
                pts = nbr_xyz[r_idx]
                if pts.shape[0]:
                    distance = np.sqrt(((pts - trial[r_idx]) ** 2).sum(1))
                    headroom = nbr_lim[r_idx] - distance
                    const += float(headroom[headroom > 0.0].sum())
            for ia, ib, lim in inter_pairs:
                a_moving, b_moving = ia in moving, ib in moving
                if a_moving and b_moving:
                    continue  # two atoms of unit u -- excluded from inter_pairs already
                if a_moving or b_moving:
                    pos = moving[ia] if a_moving else moving[ib]
                    other = trial[ib] if a_moving else trial[ia]
                    headroom = lim - np.sqrt(((pos - other) ** 2).sum(-1))
                    var += np.where(headroom > 0.0, headroom, 0.0)
                else:
                    headroom = lim - float(np.linalg.norm(trial[ia] - trial[ib]))
                    if headroom > 0.0:
                        const += headroom

            total = np.where(valid, const + var, np.inf)
            best_i = int(np.argmin(total))
            if not np.isfinite(total[best_i]):
                return None
            return float(total[best_i]), float(grid_az[best_i])

        def descend(depth: int, choice: tuple[float, ...]) -> None:
            nonlocal best, best_choice
            u = unit_group[depth]
            if depth + 1 == len(unit_group):
                result = evaluate_last(u)
                if result is not None and result[0] < best - 1e-9:
                    best = result[0]
                    best_choice = (*choice, result[1])
                return
            for delta in grid_offsets:
                z = azimuth[u] + float(delta)
                place_into(u, z, trial)
                if angle_ok(u, trial):
                    descend(depth + 1, (*choice, z))

        descend(0, ())
        if best_choice is None:
            return False
        for u, z in zip(unit_group, best_choice, strict=True):
            place_into(u, z, all_xyz)
            azimuth[u] = ((z + 180.0) % 360.0) - 180.0
        return True

    moves = 0
    for _ in range(rounds):
        tree = cKDTree(all_xyz)
        groups: set[tuple[int, ...]] = set()
        for i, j in tree.query_pairs(query_r, output_type="ndarray"):
            if abs(int(all_res[i]) - int(all_res[j])) <= 1:
                continue
            distance = float(np.linalg.norm(all_xyz[i] - all_xyz[j]))
            if distance >= limit(str(all_elem[i]), str(all_elem[j])):
                continue
            involved = tuple(sorted({u for u in (unit_of(i), unit_of(j)) if u is not None}))
            if involved:
                groups.add(involved)
        if not groups:
            break
        applied = False
        for group in groups:
            if try_polish(group):
                applied = True
                moves += 1
        if not applied:
            break

    for u in units:
        placed[u]["C"] = all_xyz[row[(u, "C")]].copy()
        placed[u + 1]["N"] = all_xyz[row[(u + 1, "N")]].copy()
        placed[u]["O"] = all_xyz[row[(u, "O")]].copy()
    return moves


def add_backbone_to_rebuilt(
    structure: Structure,
    *,
    refine: bool = True,
    on_region_done: Callable[[int], None] | None = None,
) -> BackboneResult:
    """Return a copy with N, C and O placed on the regions DODO rebuilt.

    Only the rebuilt regions gain atoms. Folded domains are returned exactly as they arrived,
    with every atom they already had -- DODO moves them rigidly and never regenerates them, and
    that is not negotiable here either.

    Parameters
    ----------
    structure
        A rebuilt structure. Regions DODO generated are identified through
        :meth:`~dodo.structure.Domain.generated_spans`, so a region whose build failed is left
        alone rather than given a backbone on top of untouched input coordinates.
    refine
        Run :func:`~dodo.construct.backbone_refine.refine_backbone` afterwards. On by default:
        measured over 3,643 held-out residues it improves every criterion at once -- N 0.188 to
        0.164 A, C 0.283 to 0.210, O 0.816 to 0.614, and the N-CA-C spread 11.70 to 3.47 degrees
        against a real 2.94.

    on_region_done
        Called with each region's residue count as its backbone is finished. This pass is slow
        enough on a large structure to look like a hang without it -- on p300 it is most of the
        wall time -- so a caller needs to report progress through it, not merely before it.

    Returns
    -------
    BackboneResult
        The new structure (the input is not modified) plus every seam left with a strained
        peptide bond, so the caller can report what could not be closed rather than leaving the
        validator to rediscover it.
    """
    from scipy.spatial import cKDTree

    from ..structure import Structure
    from ..validate.impossible import IMPOSSIBLE_SEPARATION

    # KD-tree over every structure atom (folded all-atom + rebuilt alpha carbons), for the seam
    # veto shells. Built lazily on the first seam: a from-sequence build has no seams and skips it.
    seam_tree: cKDTree | None = None

    generated: list[tuple[int, int]] = []
    for domain in structure.domains:
        for span in domain.generated_spans():
            generated.append((span.start, span.stop))
    if not generated:
        return BackboneResult(structure=structure)

    # Everything DODO did not generate: folded domains, and any region whose build failed and so
    # still holds its input coordinates. Also every alpha carbon, including the generated ones,
    # since a new backbone atom must not be placed on top of a neighbouring region's CA.
    inherited = structure.xyz[
        ~np.isin(
            structure.residue_index,
            np.concatenate([np.arange(a, b) for a, b in generated]),
        )
    ]
    obstacle_blocks = [inherited, structure.ca_xyz]

    placed: dict[int, dict[str, np.ndarray]] = {}
    strained: list[SeamStrain] = []
    for start, stop in generated:
        ca = structure.ca_xyz[start:stop]
        if ca.shape[0] < 2 or not np.all(np.isfinite(ca)):
            continue
        built = backbone_from_ca(ca)
        n_xyz, c_xyz, o_xyz = built.n_xyz, built.c_xyz, built.o_xyz
        if refine and ca.shape[0] >= 3:
            from .backbone_refine import refine_backbone

            # Regions are refined one at a time, so each one is given the backbone already placed
            # for the earlier ones. Without this they are mutually blind: refinement avoids the
            # folded domains and, within a region, itself, but two separate regions passing close
            # together would each place atoms into the other's space. Measured on p300, which has
            # six generated regions, that accounted for 4 of 19 remaining steric clashes.
            refined = refine_backbone(ca, n_xyz, c_xyz, obstacles=np.vstack(obstacle_blocks))
            n_xyz, c_xyz, o_xyz = refined.n_xyz, refined.c_xyz, refined.o_xyz
        obstacle_blocks.append(np.vstack([n_xyz, c_xyz, o_xyz]))
        if on_region_done is not None:
            on_region_done(stop - start)
        for offset, residue in enumerate(range(start, stop)):
            placed[residue] = {"N": n_xyz[offset], "C": c_xyz[offset], "O": o_xyz[offset]}

        # SEAMS. A rebuilt region is placed in isolation, so where it meets a residue DODO did not
        # rebuild, the peptide bond across the join is between an atom this function placed and one
        # that was already there -- and nothing has reconciled them. Measured on dnmt3a before this
        # was handled: C(A:PRO282)-N(A:GLU283) came out at 1.222 A against an ideal 1.329, and the
        # loop boundary at A:LYS382/A:LEU383 at 2.325 A.
        #
        # Both ends are fixed by construction against the REAL neighbour rather than a predicted
        # one: the outgoing C onto the next residue's existing N, and the incoming N onto the
        # previous residue's existing C.
        if stop < structure.n_residues:
            neighbour_n = _existing_atom(structure, stop, "N")
            if neighbour_n is not None:
                anchor_ca = structure.ca_xyz[stop - 1]
                obstacles = _seam_obstacles(structure, (stop, stop + 1))
                if seam_tree is None:
                    seam_tree = cKDTree(np.asarray(structure.xyz, dtype=np.float64))
                near = seam_tree.query_ball_point(anchor_ca, _SEAM_VETO_RADIUS)
                veto = np.asarray(structure.xyz, dtype=np.float64)[
                    np.asarray(near, dtype=np.int64)
                ]
                seam = _on_two_spheres(
                    anchor_ca,
                    CA_C_BOND_LENGTH,
                    neighbour_n,
                    C_N_PEPTIDE_BOND_LENGTH,
                    placed[stop - 1]["C"],
                )
                # The exact point holds both bond lengths and nothing else, so it is CHECKED before
                # it is accepted -- against everything nearby (measured on PTBP2, an exact seam atom
                # landed 0.797 A from a folded oxygen), against the carbonyl oxygen it induces (O is
                # fully determined by C, and on Q15642 the exact C's oxygen landed on a folded atom
                # 56 residues away in sequence), and against the residue's own N-CA-C angle (on
                # Q8N8A8 the exact C collapsed it to 78.9 degrees, putting the residue's own N and C
                # 1.90 A apart -- the geometry the bond validator flags as atoms on top of each
                # other).
                if seam is not None:
                    prev_n = placed[stop - 1].get("N")
                    angle_bad = prev_n is not None and not (
                        N_CA_C_WINDOW_MIN
                        <= _bond_angle(prev_n, anchor_ca, seam)
                        <= N_CA_C_WINDOW_MAX
                    )
                    induced_o = _place_carbonyl_oxygen(anchor_ca, seam, neighbour_n)
                    if (
                        angle_bad
                        or _would_overlap(seam, veto, IMPOSSIBLE_SEPARATION)
                        or _would_overlap(induced_o, veto, IMPOSSIBLE_SEPARATION)
                    ):
                        seam = None
                if seam is None:
                    # The exact placement is unavailable: the spheres do not intersect (on dnmt3a's
                    # loop CA(392) sits 3.741 A from N(393) where a peptide unit reaches ~2.43, so
                    # no C satisfies both bonds), or the exact point failed a check above. Hold this
                    # residue's own N-CA-C angle and spend the remaining freedom reaching toward the
                    # neighbour while clearing its atoms AND the wide shell, scoring the induced
                    # oxygen too: the seam bond stays long and is reported, but nothing impossible
                    # is written.
                    seam = _place_on_cone(
                        anchor_ca,
                        placed[stop - 1]["N"],
                        N_CA_C_ANGLE,
                        CA_C_BOND_LENGTH,
                        prefer=neighbour_n,
                        avoid=obstacles,
                        hard_avoid=veto,
                        hard_clash=IMPOSSIBLE_SEPARATION,
                        oxygen_toward=neighbour_n,
                    )
                    strained.append(
                        SeamStrain(
                            residue=stop - 1,
                            side="C",
                            bond_length=float(np.linalg.norm(seam - neighbour_n)),
                        )
                    )
                placed[stop - 1]["C"] = seam
                # The carbonyl needs the CA-C-N plane, which the two-spheres C (or the cone
                # fallback) makes well-defined -- it sits off the CA->N axis -- so O lands *trans*
                # to the neighbour's nitrogen rather than on an arbitrary axis.
                placed[stop - 1]["O"] = _place_carbonyl_oxygen(anchor_ca, seam, neighbour_n)
        if start > 0:
            neighbour_c = _existing_atom(structure, start - 1, "C")
            if neighbour_c is not None:
                first_ca = structure.ca_xyz[start]
                obstacles = _seam_obstacles(structure, (start - 1, start - 2))
                if seam_tree is None:
                    seam_tree = cKDTree(np.asarray(structure.xyz, dtype=np.float64))
                near = seam_tree.query_ball_point(first_ca, _SEAM_VETO_RADIUS)
                veto = np.asarray(structure.xyz, dtype=np.float64)[
                    np.asarray(near, dtype=np.int64)
                ]
                seam = _on_two_spheres(
                    first_ca,
                    N_CA_BOND_LENGTH,
                    neighbour_c,
                    C_N_PEPTIDE_BOND_LENGTH,
                    placed[start]["N"],
                )
                # Mirror of the C-side acceptance: the exact N is checked against the wide shell
                # (measured on PTBP2: an exact N 0.797 A from the folded carbonyl O) and against
                # the residue's own N-CA-C angle. No induced oxygen here -- this residue's O hangs
                # off its C, which the seam does not move.
                if seam is not None:
                    own_c = placed[start].get("C")
                    angle_bad = own_c is not None and not (
                        N_CA_C_WINDOW_MIN
                        <= _bond_angle(seam, first_ca, own_c)
                        <= N_CA_C_WINDOW_MAX
                    )
                    if angle_bad or _would_overlap(seam, veto, IMPOSSIBLE_SEPARATION):
                        seam = None
                if seam is None:
                    # Hold this residue's own N-CA-C angle and spend the remaining freedom reaching
                    # toward the neighbour's carbon while clearing its atoms and the wide shell.
                    seam = _place_on_cone(
                        first_ca,
                        placed[start]["C"],
                        N_CA_C_ANGLE,
                        N_CA_BOND_LENGTH,
                        prefer=neighbour_c,
                        avoid=obstacles,
                        hard_avoid=veto,
                        hard_clash=IMPOSSIBLE_SEPARATION,
                    )
                    strained.append(
                        SeamStrain(
                            residue=start,
                            side="N",
                            bond_length=float(np.linalg.norm(seam - neighbour_c)),
                        )
                    )
                placed[start]["N"] = seam

    # With every region placed and the seams reconciled, clear the steric clashes the per-region
    # refinement leaves behind -- the ones coupled between two movable units, which single-azimuth
    # descent cannot escape. Done here, over the whole assembled structure, so a clash between two
    # different rebuilt regions is visible where a per-region pass could not see it.
    if refine:
        _polish_coupled_clashes(placed, structure, generated)

    # Rebuild the flat atom records, inserting N before and C, O after each rebuilt CA so the
    # output keeps the conventional N-CA-C-O order within a residue.
    xyz: list[np.ndarray] = []
    names: list[str] = []
    elements: list[str] = []
    residue_names: list[str] = []
    residue_numbers: list[int] = []
    chain_ids: list[str] = []
    insertion_codes: list[str] = []
    b_factors: list[float] = []
    occupancies: list[float] = []

    for residue in range(structure.n_residues):
        atoms = structure.atom_slice_for_residues(residue, residue + 1)
        chain = structure.chains[int(structure.chain_index[residue])]
        existing = {str(name) for name in structure.atom_name[atoms]}
        extra = placed.get(residue, {})

        # Collect this residue's atoms in conventional N-CA-C-O order, then append once. Built as
        # a list rather than through a closure so nothing captures the loop variable.
        # b_factor and occupancy are per-RESIDUE on a Structure, not per-atom, so a placed atom
        # inherits its residue's values rather than carrying its own.
        b_value = float(structure.b_factor[residue])
        occupancy = float(structure.occupancy[residue])

        records: list[tuple[str, str, np.ndarray]] = []
        if not extra:
            # Nothing is being added to this residue, so it is passed through byte-for-byte -- same
            # atoms, same coordinates, SAME ORDER. The order matters: reordering into the
            # conventional N-CA-C-O below would rewrite 864 of dnmt3a's 912 residues, including
            # every folded-domain residue, and "folded domains are returned exactly as they
            # arrived" has to mean exactly.
            for index in range(atoms.start, atoms.stop):
                records.append(
                    (
                        str(structure.atom_name[index]),
                        str(structure.element[index]),
                        structure.xyz[index],
                    )
                )
            for name, element, point in records:
                xyz.append(np.asarray(point, dtype=np.float64))
                names.append(name)
                elements.append(element)
                residue_names.append(str(structure.residue_name[residue]))
                residue_numbers.append(int(structure.residue_number[residue]))
                chain_ids.append(chain.chain_id)
                insertion_codes.append(str(structure.insertion_code[residue]))
                b_factors.append(b_value)
                occupancies.append(occupancy)
            continue

        # A residue that gains atoms is emitted in conventional N-CA-C-O order, then everything else
        # in the order it arrived. A placed atom is only ever offered for a residue that lacks it,
        # so an existing atom of the same name wins -- which is what keeps this additive.
        emitted: set[str] = set()
        for name in ("N", "CA", "C", "O"):
            if name in extra:
                records.append((name, name[0], extra[name]))
                emitted.add(name)
            elif name in existing:
                for index in range(atoms.start, atoms.stop):
                    if str(structure.atom_name[index]) == name:
                        records.append((name, str(structure.element[index]), structure.xyz[index]))
                        emitted.add(name)
                        break
        # Side chains, and anything else the residue carries, in the order it arrived.
        for index in range(atoms.start, atoms.stop):
            name = str(structure.atom_name[index])
            if name in emitted:
                continue
            records.append((name, str(structure.element[index]), structure.xyz[index]))

        for name, element, point in records:
            xyz.append(np.asarray(point, dtype=np.float64))
            names.append(name)
            elements.append(element)
            residue_names.append(str(structure.residue_name[residue]))
            residue_numbers.append(int(structure.residue_number[residue]))
            chain_ids.append(chain.chain_id)
            insertion_codes.append(str(structure.insertion_code[residue]))
            b_factors.append(b_value)
            occupancies.append(occupancy)

    rebuilt = Structure.from_atom_records(
        xyz=np.asarray(xyz, dtype=np.float64),
        atom_name=names,
        element=elements,
        residue_name=residue_names,
        residue_number=residue_numbers,
        chain_id=chain_ids,
        insertion_code=insertion_codes,
        b_factor=b_factors,
        occupancy=occupancies,
        source=f"{structure.source} + placed backbone",
    )

    # Carry the region assignment across. from_atom_records builds a bare structure, so without
    # this the result has no domains at all -- every downstream consumer that asks what was
    # rebuilt, including the validator's provenance and the region annotation in the writer,
    # would silently see nothing. Residues map one-to-one, so the spans transfer unchanged.
    if rebuilt.n_residues != structure.n_residues:
        raise GeometryError(
            f"Placing a backbone changed the residue count from {structure.n_residues} to "
            f"{rebuilt.n_residues}, which it must never do."
        )
    for source_chain, target_chain in zip(structure.chains, rebuilt.chains, strict=True):
        target_chain.domains = [
            replace(domain, structure=rebuilt, rebuilt_loops=set(domain.rebuilt_loops))
            for domain in source_chain.domains
        ]
    return BackboneResult(structure=rebuilt, strained_seams=tuple(strained))
