""" Registration methods """

from .demons import demons_registration
from .rigid import rigid_registration
from .bspline import bspline_registration
from .soft_demons import demons_with_soft_mask

__all__ = [
    "demons_registration",
    "rigid_registration",
    "bspline_registration",
    "demons_with_soft_mask",
]
