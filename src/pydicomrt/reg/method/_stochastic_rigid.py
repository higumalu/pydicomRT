"""
rigid_registration's default algorithm, algorithm="stochastic".

Mattes mutual information evaluated on a fresh set of random points every iteration, maximised by
adaptive stochastic gradient ascent over a four-level smoothing pyramid (the strategy elastix uses),
then finished by a deterministic regular-step ascent on one large fixed set of points, so the answer
is the optimum of a single MI estimate rather than wherever the stochastic walk stopped.

Every derivative is exact: the moving image is sampled with trilinear interpolation and the
derivative *of that interpolant*, not a precomputed gradient image. ITK's precomputed gradient
disagrees with the linear interpolant on thick slices, which is what walked the default algorithm
30 mm away from a correct start on a 4 mm MR.

Two implementations of the inner loops, chosen by ``backend``: numpy here, or the C++ core built with
the package (``_rigid_core``, a ctypes binding of ``csrc/rigid_core.cpp``), which is missing only
from an install from source without a C++ compiler. The random points come from a counter-based
generator -- point k of draw d is a pure function of (seed, d, k) -- which this module reproduces bit
for bit, so both backends see the same points. Single evaluations agree to rounding; over a whole
registration rounding occasionally tips one point into another histogram bin and the stochastic
steps carry it forward, so final answers differ by up to about 0.2 mm. Neither depends on the
thread count. A change here needs the same change
there, and ``test/reg/test_rigid_stochastic.py`` compares them.

The schedule below was chosen on public brain MR/CT (RIRE, SynthRAD2023) and checked on one pelvis
CT/CBCT pair; the research is in ``internal/rigid_opt/`` (not in version control).
"""

import numpy as np
import SimpleITK as sitk

LEVEL_SIGMAS = (4.0, 2.0, 1.0, 0.0)  # Gaussian sigma per level, in voxels of each image
LEVEL_SAMPLES = (2048, 4096, 16384, 65536)
LEVEL_ITERATIONS = 256
LEVEL_BINS = 50
AVERAGE_LAST = 0.5  # the last level's answer is the mean of its last half of iterates
FINISH_SAMPLES = 500_000
FINISH_BINS = 32
FINISH_STEP_MM, FINISH_RELAXATION, FINISH_MIN_STEP_MM, FINISH_ITERATIONS = 0.5, 0.7, 0.01, 200

ASGD_A, ASGD_ALPHA = 20.0, 0.602  # elastix's defaults
SIGMOID_MAX, SIGMOID_MIN = 1.0, -0.8
RHO_SAMPLES = 20000
# No step moves points by more than this many step_mm. The gain is set from the gradient at the
# start of a level, which is ~0 when that start is already the optimum (two identical images);
# the next, ordinary gradient times that gain would then throw the transform tens of mm away.
MAX_STEP_RATIO = 2.0
MIN_INSIDE = 64

_U64 = np.uint64


def _splitmix(x):
    x = x + _U64(0x9E3779B97F4A7C15)
    x = (x ^ (x >> _U64(30))) * _U64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> _U64(27))) * _U64(0x94D049BB133111EB)
    return x ^ (x >> _U64(31))


def _uniform(seed, draw, k):
    """Uniform [0, 1) numbers for counters k of draw `draw`, identical to the C++ core's."""
    with np.errstate(over="ignore"):
        inner = _splitmix(_U64(draw) * _U64(0xD1B54A32D192ED03) ^ _splitmix(k.astype(_U64)))
        h = _splitmix(_U64(seed) ^ inner)
    return (h >> _U64(11)).astype(np.float64) * (1.0 / 9007199254740992.0)


def level_seed(seed, level):
    return (int(seed) * 1000003 + level) % 2 ** 64


def volume(image):
    """(values, origin, to_phys) of a float image: the volume format both backends take."""
    values = sitk.GetArrayFromImage(image).astype(np.float32)
    to_phys = np.array(image.GetDirection()).reshape(3, 3) * np.array(image.GetSpacing())
    return values, np.array(image.GetOrigin()), to_phys


class _Volume:
    def __init__(self, vol):
        values, origin, to_phys = vol
        self.size = np.array(values.shape[::-1])  # x, y, z
        self.origin = np.asarray(origin, dtype=np.float64)
        self.to_phys = np.asarray(to_phys, dtype=np.float64).reshape(3, 3)
        self.to_index = np.linalg.inv(self.to_phys)
        self.flat = values.astype(np.float64).ravel()
        self.vmin, self.vmax = float(values.min()), float(values.max())

    def sample(self, u, grad=False):
        """Trilinear value (and exact index-space gradient) at continuous indices u (n, 3)."""
        i0 = np.clip(np.floor(u).astype(np.int64), 0, self.size - 2)
        f = u - i0
        fx, fy, fz = f[:, 0], f[:, 1], f[:, 2]
        nx, nxy = int(self.size[0]), int(self.size[0] * self.size[1])
        base = (i0[:, 2] * int(self.size[1]) + i0[:, 1]) * nx + i0[:, 0]
        v = self.flat
        c000, c100 = v[base], v[base + 1]
        c010, c110 = v[base + nx], v[base + nx + 1]
        c001, c101 = v[base + nxy], v[base + nxy + 1]
        c011, c111 = v[base + nxy + nx], v[base + nxy + nx + 1]
        c00, c10 = c000 + (c100 - c000) * fx, c010 + (c110 - c010) * fx
        c01, c11 = c001 + (c101 - c001) * fx, c011 + (c111 - c011) * fx
        c0, c1 = c00 + (c10 - c00) * fy, c01 + (c11 - c01) * fy
        value = c0 + (c1 - c0) * fz
        if not grad:
            return value
        gx = ((c100 - c000) * (1 - fy) + (c110 - c010) * fy) * (1 - fz) + \
             ((c101 - c001) * (1 - fy) + (c111 - c011) * fy) * fz
        gy = (c10 - c00) * (1 - fz) + (c11 - c01) * fz
        gz = c1 - c0
        return value, np.stack([gx, gy, gz], axis=1)


def _bspline3(u):
    a = np.abs(u)
    return np.where(a < 1, (4 - 6 * a ** 2 + 3 * a ** 3) / 6, np.where(a < 2, (2 - a) ** 3 / 6, 0.0))


def _bspline3_d(u):
    a = np.abs(u)
    outer = np.where(u > 0, -(2 - a) ** 2 / 2, (2 - a) ** 2 / 2)
    return np.where(a < 1, -2 * u + 1.5 * u * a, np.where(a < 2, outer, 0.0))


def _skew_exp(w):
    """Rotation matrix exp([w]x) (Rodrigues)."""
    theta = float(np.sqrt(w @ w))
    if theta < 1e-12:
        return np.eye(3)
    k = w / theta
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


def _apply_step(step, R, t):
    return _skew_exp(step[:3]) @ R, t + step[3:]


class _Problem:
    """Mutual information between a fixed and a moving volume, on sets of fixed points."""

    def __init__(self, fixed, moving, centre, bins, seed):
        self.f, self.m = _Volume(fixed), _Volume(moving)
        self.c = np.asarray(centre, dtype=np.float64)
        self.bins, self.seed = int(bins), int(seed) % 2 ** 64
        self.fscale = (bins - 1) / max(self.f.vmax - self.f.vmin, 1e-12)
        self.mscale = (bins - 5) / max(self.m.vmax - self.m.vmin, 1e-12)

    def draw(self, d, n):
        """Fixed points of draw d as continuous indices (n, 3)."""
        k = np.arange(3 * n, dtype=np.uint64).reshape(n, 3)
        return _uniform(self.seed, d, k) * (self.f.size - 1)

    def offsets(self, u):
        """Physical fixed points minus the centre."""
        return u @ self.f.to_phys.T + (self.f.origin - self.c)

    def fixed_part(self, u):
        """What depends only on the fixed points: offsets from the centre and fixed bins."""
        xc = self.offsets(u)
        fv = self.f.sample(u)
        fb = np.clip(np.floor((fv - self.f.vmin) * self.fscale + 0.5), 0, self.bins - 1).astype(np.int64)
        return xc, fb

    def rho(self, d):
        """RMS distance of draw d's first points from the centre: mm moved per radian of rotation."""
        xc = self.offsets(self.draw(d, RHO_SAMPLES))
        return float(np.sqrt((xc ** 2).sum(1).mean()))

    def evaluate(self, R, t, xc, fb, gradient=True):
        """(mi, gradient (dw, dt), number inside), or None when fewer than MIN_INSIDE map inside."""
        B = self.bins
        r = xc @ R.T
        u_m = (r + self.c + t - self.m.origin) @ self.m.to_index.T
        inside = np.all((u_m >= 0) & (u_m <= self.m.size - 1), axis=1)
        n = int(inside.sum())
        if n < MIN_INSIDE:
            return None
        u_m, r, fb = u_m[inside], r[inside], fb[inside]
        mv, g_idx = self.m.sample(u_m, grad=True)
        mb = np.clip((mv - self.m.vmin) * self.mscale + 2.0, 2.0, B - 3.0)
        tap = np.floor(mb).astype(np.int64)[:, None] - 1 + np.arange(4)[None, :]
        du = mb[:, None] - tap
        flat = fb[:, None] * B + tap
        joint = np.bincount(flat.ravel(), _bspline3(du).ravel(), B * B).reshape(B, B)
        joint = joint / joint.sum()
        pf, pm = joint.sum(1), joint.sum(0)
        nz = joint > 0
        safe = np.where(nz, joint, 1.0)
        mi = float((joint * np.where(nz, np.log(safe) - np.log(pf[:, None] * pm[None, :] + 1e-300), 0.0)).sum())
        if not gradient:
            return mi, None, n
        L = np.where(nz, np.log(safe) - np.log(pm[None, :] + 1e-300), 0.0).ravel()
        s = (_bspline3_d(du) * L[flat]).sum(1)
        g_phys = g_idx @ self.m.to_index  # d moving / d y
        d_t = (s[:, None] * g_phys).sum(0) * self.mscale
        d_w = (s[:, None] * np.cross(r, g_phys)).sum(0) * self.mscale
        return mi, np.concatenate([d_w, d_t]) / n, n


def evaluate(fixed, moving, centre, bins, R, t, u, threads=1, gradient=True):
    """Mutual information (and its gradient) on fixed continuous indices u (n, 3): (mi, grad, n_inside)."""
    p = _Problem(fixed, moving, centre, bins, 0)
    out = p.evaluate(np.asarray(R, dtype=np.float64), np.asarray(t, dtype=np.float64),
                     *p.fixed_part(np.asarray(u, dtype=np.float64).reshape(-1, 3)), gradient=gradient)
    return (out[0], out[1], out[2]) if out is not None else (0.0, None, 0)


def level(fixed, moving, centre, bins, seed, threads, R, t, iterations, samples, step_mm):
    """One pyramid level of adaptive stochastic gradient ascent: (R, t, history (iterations, 12))."""
    p = _Problem(fixed, moving, centre, bins, seed)
    R, t = np.array(R, dtype=np.float64), np.array(t, dtype=np.float64)
    rho = p.rho(0)
    scale = np.array([1 / rho] * 3 + [1.0] * 3)
    draw = [1]

    def gradient():
        for _ in range(8):
            out = p.evaluate(R, t, *p.fixed_part(p.draw(draw[0], samples)))
            draw[0] += 1
            if out is not None:
                return out[1]
        raise RuntimeError("too few samples map inside the moving image")

    # Gain: the first step moves points by about step_mm.
    g_norm = 0.0
    for _ in range(5):
        g_norm += float(np.linalg.norm(gradient() * scale)) / 5
    a = step_mm * (ASGD_A + 1) ** ASGD_ALPHA / max(g_norm, 1e-12)
    time, previous = 0.0, None
    history = np.zeros((iterations, 12))
    for k in range(iterations):
        g = gradient() * scale
        if previous is not None:
            cos = float(g @ previous) / (np.sqrt(g @ g) * np.sqrt(previous @ previous) + 1e-30)
            time = max(0.0, time + SIGMOID_MIN + (SIGMOID_MAX - SIGMOID_MIN) / (1 + np.exp(4 * cos)))
        previous = g
        gain = a / (ASGD_A + time + 1) ** ASGD_ALPHA
        moved = gain * np.sqrt(g @ g)  # mm, rotation counted at the RMS radius
        if moved > MAX_STEP_RATIO * step_mm:
            gain *= MAX_STEP_RATIO * step_mm / moved
        R, t = _apply_step(gain * g * scale, R, t)
        history[k, :9], history[k, 9:] = R.ravel(), t
    return R, t, history


def refine(fixed, moving, centre, bins, seed, threads, R, t, samples, step, relax, min_step, iterations):
    """Regular-step gradient ascent on one fixed point set (draw 0): (R, t, iterations taken)."""
    p = _Problem(fixed, moving, centre, bins, seed)
    R, t = np.array(R, dtype=np.float64), np.array(t, dtype=np.float64)
    rho = p.rho(0)
    scale = np.array([1 / rho] * 3 + [1.0] * 3)
    xc, fb = p.fixed_part(p.draw(0, samples))  # the point set never changes: compute its fixed side once
    previous = np.zeros(6)
    taken = 0
    for taken in range(iterations):
        out = p.evaluate(R, t, xc, fb)
        if out is None:
            raise RuntimeError("too few samples map inside the moving image")
        g = out[1] * scale
        if taken > 0 and g @ previous < 0:
            step *= relax
        if step < min_step:
            break
        previous = g
        R, t = _apply_step(step * g / (np.sqrt(g @ g) + 1e-30) * scale, R, t)
    else:
        taken = iterations
    return R, t, taken


def _core(backend):
    """The module that runs the inner loops for `backend`."""
    import importlib
    import sys
    this = sys.modules[__name__]
    if backend == "numpy":
        return this
    try:
        core = importlib.import_module("pydicomrt.reg.method._rigid_core")
    except ImportError as error:
        if backend == "cpp":
            raise ImportError(
                f"backend='cpp' needs the compiled core, which this installation lacks ({error}). "
                "Reinstall pydicomrt from a wheel, or from source with a C++17 compiler; "
                "backend='auto' runs on numpy instead"
            ) from None
        return this
    return core


def _smoothed(image, sigma_voxels):
    if sigma_voxels <= 0:
        return image
    return sitk.SmoothingRecursiveGaussian(image, [sigma_voxels * s for s in image.GetSpacing()])


def register(fixed_image, moving_image, *, seed=0, backend="auto", threads=1):
    """
    Rigid fixed -> moving resampling transform as (R, t, centre): y = R (x - centre) + centre + t.

    Starts from the geometric centres aligned, as CenteredTransformInitializer(GEOMETRY) does.
    """
    core = _core(backend)
    fixed = sitk.Cast(fixed_image, sitk.sitkFloat32)
    moving = sitk.Cast(moving_image, sitk.sitkFloat32)
    for name, image in (("fixed_image", fixed), ("moving_image", moving)):
        if image.GetDimension() != 3 or min(image.GetSize()) < 2:
            raise ValueError(f"{name} must be 3D with at least 2 voxels along every axis, "
                             f"got size {image.GetSize()}")

    def geometric_centre(image):
        return np.array(image.TransformContinuousIndexToPhysicalPoint([(n - 1) / 2 for n in image.GetSize()]))

    centre = geometric_centre(fixed)
    R, t = np.eye(3), geometric_centre(moving) - centre
    step_mm = float(np.mean(fixed.GetSpacing()))  # elastix's maximum step length rule

    history = None
    for k, (sigma, samples) in enumerate(zip(LEVEL_SIGMAS, LEVEL_SAMPLES)):
        R, t, history = core.level(volume(_smoothed(fixed, sigma)), volume(_smoothed(moving, sigma)), centre,
                                   LEVEL_BINS, level_seed(seed, k), threads, R, t, LEVEL_ITERATIONS,
                                   samples, step_mm)
    tail = history[int(len(history) * (1 - AVERAGE_LAST)):]
    u, _, vt = np.linalg.svd(tail[:, :9].reshape(-1, 3, 3).sum(0))
    R, t = u @ vt, tail[:, 9:].mean(0)

    R, t, _ = core.refine(volume(fixed), volume(moving), centre, FINISH_BINS,
                          level_seed(seed, len(LEVEL_SIGMAS)), threads, R, t, FINISH_SAMPLES,
                          FINISH_STEP_MM, FINISH_RELAXATION, FINISH_MIN_STEP_MM, FINISH_ITERATIONS)
    return R, t, centre
