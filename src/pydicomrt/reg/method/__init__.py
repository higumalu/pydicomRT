""" Registration methods """

from .demons import demons_registration
from .rigid import rigid_registration

__all__ = [
    "demons_registration",
    "rigid_registration",
]
