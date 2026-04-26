from __future__ import annotations

import math

TWOPI = 2.0 * math.pi


def hz_to_rad_s(value_hz: float) -> float:
    return TWOPI * float(value_hz)


def rad_s_to_hz(value_rad_s: float) -> float:
    return float(value_rad_s) / TWOPI
