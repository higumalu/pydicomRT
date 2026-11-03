"""
ctypes binding of the compiled core of ``rigid_registration(algorithm="stochastic")``.

The library is ``csrc/rigid_core.cpp``, built by CMake into this directory when pydicomrt is
installed. Importing this module raises ImportError when it is missing (an install from source
without a C++ compiler) or stale; ``_stochastic_rigid._core`` then falls back to numpy. Its functions
mirror ``_stochastic_rigid.evaluate`` / ``level`` / ``refine`` signature for signature.

A volume is passed as ``(values, origin, to_phys)``: float32 values in sitk.GetArrayFromImage order
(z, y, x), the physical origin, and the 3x3 matrix direction @ diag(spacing).
"""

import ctypes
import sys
from pathlib import Path

import numpy as np

# Bump together with rigid_core_api_version() in csrc/rigid_core.cpp on any change to the C ABI.
API_VERSION = 1


def _find_library():
    # Every directory of the package, not just this file's: an editable install leaves this file in
    # the source tree and the built library in site-packages, both on the package's __path__.
    directories = [Path(d) for d in sys.modules[__package__].__path__]
    for directory in directories:
        for name in ("librigid_core.so", "librigid_core.dylib", "rigid_core.dll", "librigid_core.dll"):
            if (directory / name).exists():
                return directory / name
    raise ImportError(f"pydicomrt: compiled rigid core not found in {', '.join(map(str, directories))}")


_lib = ctypes.CDLL(str(_find_library()))


class _Vol(ctypes.Structure):
    _fields_ = [("v", ctypes.POINTER(ctypes.c_float)), ("size", ctypes.c_int64 * 3),
                ("origin", ctypes.c_double * 3), ("to_phys", ctypes.c_double * 9),
                ("to_index", ctypes.c_double * 9), ("vmin", ctypes.c_double), ("vmax", ctypes.c_double)]


class _Level(ctypes.Structure):
    _fields_ = [("f", ctypes.POINTER(_Vol)), ("m", ctypes.POINTER(_Vol)), ("centre", ctypes.c_double * 3),
                ("bins", ctypes.c_int), ("seed", ctypes.c_uint64), ("threads", ctypes.c_int)]


_dp = ctypes.POINTER(ctypes.c_double)
_lp = ctypes.POINTER(_Level)
_lib.rigid_core_api_version.restype = ctypes.c_int
_lib.rigid_core_evaluate.restype = ctypes.c_int64
_lib.rigid_core_evaluate.argtypes = [_lp, _dp, _dp, _dp, ctypes.c_int64, _dp, _dp]
_lib.rigid_core_level.restype = ctypes.c_int
_lib.rigid_core_level.argtypes = [_lp, _dp, _dp, ctypes.c_int, ctypes.c_int64, ctypes.c_double, _dp]
_lib.rigid_core_refine.restype = ctypes.c_int
_lib.rigid_core_refine.argtypes = [_lp, _dp, _dp, ctypes.c_int64, ctypes.c_double, ctypes.c_double,
                                   ctypes.c_double, ctypes.c_int]

if _lib.rigid_core_api_version() != API_VERSION:
    raise ImportError("pydicomrt: compiled rigid core does not match its binding; rebuild pydicomrt")


def _f64(a, shape=None):
    a = np.ascontiguousarray(a, dtype=np.float64)
    return a if shape is None else a.reshape(shape)


def _ptr(a):
    return a.ctypes.data_as(_dp)


class _Pair:
    """ctypes view of a fixed/moving pair; holds the arrays alive for as long as it lives."""

    def __init__(self, fixed, moving, centre, bins, seed, threads):
        self.vols, self.keep = [], []
        for values, origin, to_phys in (fixed, moving):
            values = np.ascontiguousarray(values, dtype=np.float32)
            to_phys = _f64(to_phys, (3, 3))
            v = _Vol()
            v.v = values.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
            v.size[:] = list(values.shape[::-1])
            v.origin[:] = _f64(origin).tolist()
            v.to_phys[:] = to_phys.ravel().tolist()
            v.to_index[:] = np.linalg.inv(to_phys).ravel().tolist()
            v.vmin, v.vmax = float(values.min()), float(values.max())
            self.vols.append(v)
            self.keep.append(values)
        self.level = _Level()
        self.level.f, self.level.m = ctypes.pointer(self.vols[0]), ctypes.pointer(self.vols[1])
        self.level.centre[:] = _f64(centre).tolist()
        self.level.bins = int(bins)
        self.level.seed = int(seed) % 2 ** 64
        self.level.threads = max(1, int(threads))


def evaluate(fixed, moving, centre, bins, R, t, u, threads=1, gradient=True):
    """Mutual information (and its gradient) on fixed continuous indices u (n, 3): (mi, grad, n_inside)."""
    pair = _Pair(fixed, moving, centre, bins, 0, threads)
    R, t, u = _f64(R), _f64(t), _f64(u, (-1, 3))
    mi, grad = np.zeros(1), np.zeros(6)
    n = _lib.rigid_core_evaluate(ctypes.byref(pair.level), _ptr(R), _ptr(t), _ptr(u), len(u), _ptr(mi),
                                 _ptr(grad) if gradient else None)
    return float(mi[0]), (grad if gradient else None), int(n)


def level(fixed, moving, centre, bins, seed, threads, R, t, iterations, samples, step_mm):
    """One pyramid level of adaptive stochastic gradient ascent: (R, t, history (iterations, 12))."""
    pair = _Pair(fixed, moving, centre, bins, seed, threads)
    R, t = _f64(R).copy(), _f64(t).copy()
    history = np.zeros((int(iterations), 12))
    if _lib.rigid_core_level(ctypes.byref(pair.level), _ptr(R), _ptr(t), int(iterations), int(samples),
                             float(step_mm), _ptr(history)) != 0:
        raise RuntimeError("too few samples map inside the moving image")
    return R, t, history


def refine(fixed, moving, centre, bins, seed, threads, R, t, samples, step, relax, min_step, iterations):
    """Regular-step gradient ascent on one fixed sample set: (R, t, iterations taken)."""
    pair = _Pair(fixed, moving, centre, bins, seed, threads)
    R, t = _f64(R).copy(), _f64(t).copy()
    taken = _lib.rigid_core_refine(ctypes.byref(pair.level), _ptr(R), _ptr(t), int(samples), float(step),
                                   float(relax), float(min_step), int(iterations))
    if taken < 0:
        raise RuntimeError("too few samples map inside the moving image")
    return R, t, taken
