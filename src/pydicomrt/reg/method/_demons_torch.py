"""
PyTorch implementation of the multiscale demons used by ``demons_registration``.

Selected with ``demons_registration(..., backend="torch")``. PyTorch is an optional
dependency (``pip install "pydicomrt[gpu]"``); nothing outside this module imports it, and
this module is only imported when that backend is requested.

The implementation follows the ITK 5.4 sources rather than the algorithm's description, so
that it reproduces the SimpleITK backend field for field (agreement is a few thousandths of a
millimetre on the test pairs):

* ``ESMDemonsRegistrationFunction`` -- symmetric forces, the update bounded to half a voxel
  through the mean squared spacing, the warped-moving gradient taken one-sided at image edges
  and next to samples that left the moving image, zero update for such samples;
* ``PDEDeformableRegistrationFilter`` -- update field smoothed at 1 voxel, displacement field
  at the regularisation width, both with ``itk::GaussianOperator`` (MaximumError 0.1,
  MaximumKernelWidth 30);
* ``FiniteDifferenceImageFilter::Halt`` -- a level stops early once the RMS of the raw update
  over the voxels still inside the moving image falls below MaximumRMSError (0.02). SimpleITK
  leaves this on by default, so the requested iteration count is a ceiling, not a count;
* ``sitk.DiscreteGaussian`` (MaximumError 0.01) for the pyramid, the impulse response of
  ``sitk.SmoothingRecursiveGaussian`` for the per-level field smoothing, and ITK's resampling
  convention -- the buffer reaches half a voxel past the outer voxel centres, sampled as the
  edge value, and anything beyond it gets the default value.

These details matter most at coarse pyramid levels, where most voxels sit near an edge: an
implementation that skips them drifts millimetres away from SimpleITK as iterations grow.

Geometry is handled in physical space with the full direction matrix, so oblique and
anisotropic images, and moving images on a different grid, behave as in the SimpleITK path.
Arrays are ``(z, y, x)``; displacement components are physical ``(x, y, z)`` in mm.

``metric="mind"`` (``_mind_demons_level``) is the one part that is not a port: demons on
MIND-SSC descriptors for images of different modalities, with no SimpleITK counterpart. It
reuses the geometry, sampling, kernels and update rule above but carries one field through
all levels instead of composing per-level refinements; see that function.
"""

import logging
from functools import lru_cache

import numpy as np
import SimpleITK as sitk
import torch
import torch.nn.functional as F
from scipy.special import ive

logger = logging.getLogger(__name__)

#: ITK defaults that FastSymmetricForcesDemonsRegistrationFilter runs with in this library.
MAXIMUM_UPDATE_STEP_LENGTH = 0.5
INTENSITY_DIFFERENCE_THRESHOLD = 0.001
DENOMINATOR_THRESHOLD = 1e-9
MAXIMUM_RMS_ERROR = 0.02

#: Peak device memory per voxel of the finest pyramid level (measured: ~100 B allocated,
#: ~175 B reserved by the caching allocator). Used for the out-of-memory message only.
BYTES_PER_FINEST_VOXEL = 180
UPDATE_FIELD_SIGMA_VOXELS = 1.0
DEMONS_MAXIMUM_ERROR = 0.1
DEMONS_MAXIMUM_KERNEL_WIDTH = 30
PYRAMID_MAXIMUM_ERROR = 0.01


def resolve_device(device):
    """Return a ``torch.device``; ``None`` picks CUDA when available, else the CPU."""
    if device is not None:
        return torch.device(device)
    if torch.cuda.is_available():
        return torch.device("cuda")
    logger.warning("backend='torch' found no CUDA device; running on the CPU")
    return torch.device("cpu")


class _Geometry:
    """Grid geometry in SimpleITK order: size/spacing/origin as (x, y, z), direction 3x3."""

    def __init__(self, size, spacing, origin, direction):
        self.size = tuple(int(s) for s in size)
        self.spacing = np.asarray(spacing, dtype=np.float64)
        self.origin = np.asarray(origin, dtype=np.float64)
        self.direction = np.asarray(direction, dtype=np.float64).reshape(3, 3)

    @classmethod
    def of(cls, image):
        return cls(image.GetSize(), image.GetSpacing(), image.GetOrigin(), image.GetDirection())

    @property
    def shape(self):
        return self.size[::-1]

    def level(self, resolution, isotropic):
        """The grid ``smooth_and_resample`` produces for one pyramid level."""
        if isotropic:
            scale = float(resolution) / self.spacing
            new_size = [int(sz / float(sf) + 0.5) for sz, sf in zip(self.size, scale)]
        else:
            new_size = [int(sz / float(resolution) + 0.5) for sz in self.size]
        for i, n in enumerate(new_size):
            if n < 2:
                raise ValueError(
                    f"Computed new size for dimension {i} is {n}; use a smaller shrink factor "
                    f"or a larger isotropic voxel size."
                )
        new_spacing = [(so - 1) * sp / (sn - 1) for so, sp, sn in zip(self.size, self.spacing, new_size)]
        return _Geometry(new_size, new_spacing, self.origin, self.direction)

    def matches(self, other):
        return (self.size == other.size and np.allclose(self.spacing, other.spacing, rtol=0, atol=1e-9)
                and np.allclose(self.origin, other.origin, rtol=0, atol=1e-9)
                and np.allclose(self.direction, other.direction, rtol=0, atol=1e-12))

    def to_image(self, array_zyx3):
        """A VectorFloat64 field image; converting float32 in SimpleITK is the fast path."""
        image = sitk.GetImageFromArray(np.ascontiguousarray(array_zyx3, dtype=np.float32), isVector=True)
        image.SetSpacing(self.spacing.tolist())
        image.SetOrigin(self.origin.tolist())
        image.SetDirection(self.direction.ravel().tolist())
        return sitk.Cast(image, sitk.sitkVectorFloat64)


def _source_index(source, target, device, displacement=None):
    """
    Continuous index into `source` of every `target` voxel centre, moved by `displacement`.

    p = o_t + D_t (s_t * i);  q = p + u(i);  j = (D_s^T (q - o_s)) / s_s.
    Returns (Z, Y, X, 3) with the last axis in (x, y, z) order. `displacement` is (3, Z, Y, X)
    in physical mm on the target grid.
    """
    to_source = np.diag(1.0 / source.spacing) @ source.direction.T
    matrix = to_source @ target.direction @ np.diag(target.spacing)
    offset = to_source @ (target.origin - source.origin)
    axes = [torch.arange(n, device=device, dtype=torch.float32) for n in target.size]
    ix = axes[0].view(1, 1, -1)
    iy = axes[1].view(1, -1, 1)
    iz = axes[2].view(-1, 1, 1)
    components = []
    for k in range(3):
        j = float(matrix[k, 0]) * ix + float(matrix[k, 1]) * iy + float(matrix[k, 2]) * iz + float(offset[k])
        j = j.expand(target.shape)
        if displacement is not None:
            j = j + (float(to_source[k, 0]) * displacement[0]
                     + float(to_source[k, 1]) * displacement[1]
                     + float(to_source[k, 2]) * displacement[2])
        components.append(j)
    return torch.stack(components, dim=-1)


def _sample(volume, index, source_size, pad, return_inside=False):
    """
    Linear interpolation of a (C, Z, Y, X) volume at continuous indices, as sitk.Resample
    does: inside the buffer, neighbours clamp to the edge; outside it, `pad`.
    """
    inside = None
    for k, n in enumerate(source_size):
        axis = (index[..., k] >= -0.5) & (index[..., k] < n - 0.5)
        inside = axis if inside is None else inside & axis
    for k, n in enumerate(source_size):
        index[..., k].mul_(2.0 / max(n - 1, 1)).sub_(1.0)  # in place: normalised coordinates
    out = F.grid_sample(volume.unsqueeze(0), index.unsqueeze(0), mode="bilinear",
                        padding_mode="border", align_corners=True).squeeze(0)
    out = torch.where(inside, out, torch.full_like(out, float(pad)))
    return (out, inside) if return_inside else out


@lru_cache(maxsize=64)
def _gaussian_operator(sigma_voxels, maximum_error, maximum_width):
    """itk::GaussianOperator coefficients: the discrete (modified Bessel) Gaussian."""
    variance = float(sigma_voxels) ** 2
    coeff = [float(ive(0, variance)), float(ive(1, variance))]
    total = coeff[0] + 2.0 * coeff[1]
    i = 2
    while total < 1.0 - maximum_error:
        coeff.append(float(ive(i, variance)))
        total += 2.0 * coeff[i]
        if coeff[i] <= 0.0 or len(coeff) > maximum_width:
            break
        i += 1
    half = np.asarray(coeff) / total
    return np.concatenate([half[:0:-1], half])


@lru_cache(maxsize=64)
def _recursive_gaussian(sigma_mm, spacing):
    """
    Impulse response of sitk.SmoothingRecursiveGaussian along one axis.

    The Deriche IIR filter is not a sampled Gaussian -- least of all at the sub-voxel widths
    the coarse pyramid levels use -- so its response is measured rather than approximated:
    smooth a plane impulse with SimpleITK itself and read the line back. ITK extends the
    signal with its edge value, which is exactly replicate padding of this kernel.
    """
    radius = int(np.ceil(6.0 * sigma_mm / spacing)) + 8
    length = 2 * radius + 1
    impulse = np.zeros((4, 4, length))
    impulse[..., radius] = 1.0
    image = sitk.GetImageFromArray(impulse)
    image.SetSpacing((float(spacing), 1.0, 1.0))
    smoothed = sitk.SmoothingRecursiveGaussian(image, [float(sigma_mm), 1.0, 1.0])
    return sitk.GetArrayFromImage(smoothed)[1, 1, :].copy()


def _operator_kernels(sigma_voxels_xyz, maximum_error, maximum_width):
    return [_gaussian_operator(float(s), maximum_error, maximum_width) if s > 0 else None
            for s in sigma_voxels_xyz]


def _recursive_kernels(sigma_mm_xyz, spacing_xyz):
    return [_recursive_gaussian(float(s), float(sp)) if s > 0 else None
            for s, sp in zip(sigma_mm_xyz, spacing_xyz)]


_KERNEL_CACHE = {}


def _kernel_weight(kernel, dim, channels, device, dtype):
    """The kernel as a depthwise conv3d weight on `device`, built once per configuration."""
    key = (kernel.tobytes(), len(kernel), dim, channels, str(device), dtype)
    weight = _KERNEL_CACHE.get(key)
    if weight is None:
        shape = [1, 1, 1, 1, 1]
        shape[dim] = len(kernel)
        weight = torch.tensor(kernel, dtype=dtype).view(shape).repeat(channels, 1, 1, 1, 1).to(device)
        if len(_KERNEL_CACHE) > 256:
            _KERNEL_CACHE.clear()
        _KERNEL_CACHE[key] = weight
    return weight


def _smooth(volume, kernels_xyz):
    """
    Separable filtering of a (C, Z, Y, X) volume along its index axes with one 1-D kernel per
    (x, y, z) axis (None skips the axis), edges replicated as in ITK's zero-flux Neumann
    boundary condition.
    """
    channels = volume.shape[0]
    out = volume.unsqueeze(0)
    for dim, kernel in zip((4, 3, 2), kernels_xyz):  # x, y, z -> tensor dims 4, 3, 2
        if kernel is None:
            continue
        radius = (len(kernel) - 1) // 2
        if radius == 0:
            continue
        weight = _kernel_weight(kernel, dim, channels, volume.device, volume.dtype)
        pad = [0, 0, 0, 0, 0, 0]
        pad[2 * (4 - dim)] = radius
        pad[2 * (4 - dim) + 1] = radius
        out = F.conv3d(F.pad(out, pad, mode="replicate"), weight, groups=channels)
    return out.squeeze(0)


def _central_gradient(volume, spacing):
    """CentralDifferenceImageFunction along index axes, zero on the border; (3, Z, Y, X)."""
    v = volume[0]
    grads = []
    for axis, sp in zip((2, 1, 0), spacing):
        g = torch.zeros_like(v)
        n = v.shape[axis]
        if n >= 3:
            centre, ahead, behind = [slice(None)] * 3, [slice(None)] * 3, [slice(None)] * 3
            centre[axis], ahead[axis], behind[axis] = slice(1, n - 1), slice(2, n), slice(0, n - 2)
            g[tuple(centre)] = (v[tuple(ahead)] - v[tuple(behind)]) / (2.0 * float(sp))
        grads.append(g)
    return torch.stack(grads)


def _warped_moving_gradient(warped, valid, spacing):
    """
    ESMDemonsRegistrationFunction's gradient of the warped moving image: central where both
    neighbours are valid, one-sided where only one is, zero where neither is.
    """
    v = warped[0]
    grads = []
    for axis, sp in zip((2, 1, 0), spacing):
        n = v.shape[axis]
        ahead = torch.zeros_like(v)
        behind = torch.zeros_like(v)
        has_ahead = torch.zeros_like(valid)
        has_behind = torch.zeros_like(valid)
        if n >= 2:
            lo, hi = [slice(None)] * 3, [slice(None)] * 3
            lo[axis], hi[axis] = slice(0, n - 1), slice(1, n)
            lo, hi = tuple(lo), tuple(hi)
            ahead[lo] = v[hi]
            has_ahead[lo] = valid[hi]
            behind[hi] = v[lo]
            has_behind[hi] = valid[lo]
        sp = float(sp)
        g = torch.where(has_ahead & has_behind, (ahead - behind) / (2.0 * sp),
            torch.where(has_ahead, (ahead - v) / sp,
            torch.where(has_behind, (v - behind) / sp, torch.zeros_like(v))))
        grads.append(g)
        del ahead, behind, has_ahead, has_behind
    return torch.stack(grads)


def _to_physical(vectors, direction):
    """Rotate (3, Z, Y, X) index-axis vectors into physical axes."""
    if np.allclose(direction, np.eye(3)):
        return vectors
    d = torch.tensor(direction, device=vectors.device, dtype=vectors.dtype)
    return torch.einsum("ij,j...->i...", d, vectors)


def _demons_level(fixed, moving, geometry, iterations, field_sigma_voxels, device):
    """FastSymmetricForcesDemonsRegistrationFilter on one level, starting from a zero field."""
    spacing = geometry.spacing
    normalizer = float(np.mean(spacing ** 2)) * MAXIMUM_UPDATE_STEP_LENGTH ** 2
    field_kernels = _operator_kernels(field_sigma_voxels, DEMONS_MAXIMUM_ERROR, DEMONS_MAXIMUM_KERNEL_WIDTH)
    update_kernels = _operator_kernels([UPDATE_FIELD_SIGMA_VOXELS] * 3, DEMONS_MAXIMUM_ERROR,
                                       DEMONS_MAXIMUM_KERNEL_WIDTH)
    fixed_gradient = _central_gradient(fixed, spacing)
    field = torch.zeros((3,) + geometry.shape, device=device, dtype=torch.float32)

    for iteration in range(iterations):
        if iteration > 0:
            # ITK smooths the field after every update; this is that smoothing, deferred.
            field = _smooth(field, field_kernels)
        index = _source_index(geometry, geometry, device, displacement=field)
        warped, inside = _sample(moving, index, geometry.size, 0.0, return_inside=True)
        del index
        gradient = fixed_gradient + _warped_moving_gradient(warped, inside, spacing)
        gradient = _to_physical(gradient, geometry.direction)
        speed = fixed[0] - warped[0]
        del warped
        denominator = (gradient * gradient).sum(0) + speed * speed / normalizer
        valid = inside & (speed.abs() >= INTENSITY_DIFFERENCE_THRESHOLD) & (denominator >= DENOMINATOR_THRESHOLD)
        scale = torch.where(valid, 2.0 * speed / denominator.clamp_min(DENOMINATOR_THRESHOLD),
                            torch.zeros_like(speed))
        del speed, denominator, valid
        update = gradient * scale
        del gradient, scale
        # ITK's RMSChange: raw update, averaged over the voxels whose warped sample stayed
        # inside the moving image (the others return before being counted).
        counted = int(inside.sum())
        rms = None
        if counted:
            squared = torch.where(inside, (update * update).sum(0), torch.zeros((), device=device))
            rms = float(torch.sqrt(squared.sum(dtype=torch.float64) / counted))
            del squared
        del inside
        field = field + _smooth(update, update_kernels)
        del update
        if rms is not None and rms < MAXIMUM_RMS_ERROR:
            break
    return _smooth(field, field_kernels)


# ---------------------------------------------------------------------------------------
# metric="mind": demons on MIND-SSC descriptors, for images of different modalities
# ---------------------------------------------------------------------------------------

#: MIND-SSC patch radius and neighbour distance, in voxels of the pyramid level.
MIND_RADIUS = 1
MIND_DILATION = 2
#: Descriptor channels warped per pass; bounds the per-iteration temporaries.
MIND_CHUNK = 4
#: Peak device memory per voxel of the finest level for metric="mind" (measured ~285 B
#: allocated, ~355 B reserved). For the out-of-memory message only.
MIND_BYTES_PER_FINEST_VOXEL = 360

# The twelve MIND-SSC comparisons: pairs of the six face neighbours at distance sqrt(2) from
# each other -- every edge of the octahedron they span. Offsets are (x, y, z) in units of the
# dilation.
_SIX = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
_MIND_PAIRS = [(a, b) for i, a in enumerate(_SIX) for b in _SIX[i + 1:]
               if sum((p - q) ** 2 for p, q in zip(a, b)) == 2]


def _shifted(volume, offset_xyz):
    """volume (Z, Y, X) sampled at x + offset, edges replicated."""
    out = volume
    for dim, step in zip((2, 1, 0), offset_xyz):
        if step == 0:
            continue
        n = out.shape[dim]
        index = torch.clamp(torch.arange(n, device=out.device) + step, 0, n - 1)
        out = out.index_select(dim, index)
    return out


def _mind_ssc(volume):
    """
    MIND-SSC descriptor (Heinrich et al., MICCAI 2013) of a (1, Z, Y, X) volume: (12, Z, Y, X).

    Each channel is how different two neighbouring patches are, normalised by the local mean
    of those differences and mapped through exp(-d). That pattern depends on where edges are,
    not on which side is brighter, so it is the same for MR and CT wherever both show the
    same boundary -- and unchanged by any gain or offset of either image.
    """
    v = volume[0]
    box = 2 * MIND_RADIUS + 1
    channels = []
    for a, b in _MIND_PAIRS:
        d = (_shifted(v, [MIND_DILATION * c for c in a]) - _shifted(v, [MIND_DILATION * c for c in b])) ** 2
        d = F.avg_pool3d(F.pad(d[None, None], (MIND_RADIUS,) * 6, mode="replicate"), box, stride=1)
        channels.append(d[0, 0])
        del d
    ssd = torch.stack(channels)
    del channels
    ssd -= ssd.min(0, keepdim=True)[0]
    variance = ssd.mean(0, keepdim=True)
    mean_variance = float(variance.mean())
    variance.clamp_(mean_variance * 1e-3, mean_variance * 1e3)
    return torch.exp(-ssd / variance)


def _index_gradient(volume, axis, spacing):
    """Central difference of a (C, Z, Y, X) volume along index axis 0 = x, 1 = y, 2 = z, in
    1/mm, one-sided at the edges."""
    dim = 3 - axis
    n = volume.shape[dim]
    if n < 2:
        return torch.zeros_like(volume)
    ahead = torch.cat([volume.narrow(dim, 1, n - 1), volume.narrow(dim, n - 1, 1)], dim)
    behind = torch.cat([volume.narrow(dim, 0, 1), volume.narrow(dim, 0, n - 1)], dim)
    steps = torch.full((n,), 2.0, device=volume.device, dtype=volume.dtype)
    steps[0] = steps[-1] = 1.0
    shape = [1, 1, 1, 1]
    shape[dim] = n
    return (ahead - behind) / (steps.view(shape) * float(spacing))


def _mind_demons_level(fixed, moving, moving_inside, geometry, iterations, field_sigma_voxels, device,
                       field):
    """
    The demons level of `_demons_level` with the intensity difference replaced by the
    difference of MIND-SSC descriptors -- multichannel symmetric-forces demons. For one
    channel the update is exactly ITK's ESM update, so MaximumUpdateStepLength, the field and
    update smoothing and the RMS early stop all keep their meaning.

    Unlike `_demons_level` it continues the *total* field (`field`, on this level's grid)
    instead of estimating a refinement from zero: every iteration's elastic smoothing then
    acts on the whole deformation, so a finer level can take back what a coarser one got
    wrong. Composed refinements leave each level's error frozen in, and with a descriptor
    metric that adds up -- 3 mm of invented motion in the brain on undeformed MR/CT pairs,
    against 0.7 mm this way.

    `moving` is the moving level resampled onto the fixed level grid *without* the field, so
    both descriptors see the same voxel neighbourhoods; `moving_inside` marks the voxels that
    came from inside the moving image. The moving descriptors are computed once and warped,
    rather than recomputed from a re-warped image every iteration.
    """
    spacing = geometry.spacing
    normalizer = float(np.mean(spacing ** 2)) * MAXIMUM_UPDATE_STEP_LENGTH ** 2
    field_kernels = _operator_kernels(field_sigma_voxels, DEMONS_MAXIMUM_ERROR, DEMONS_MAXIMUM_KERNEL_WIDTH)
    update_kernels = _operator_kernels([UPDATE_FIELD_SIGMA_VOXELS] * 3, DEMONS_MAXIMUM_ERROR,
                                       DEMONS_MAXIMUM_KERNEL_WIDTH)
    fixed_mind = _mind_ssc(fixed)
    moving_mind = torch.cat([_mind_ssc(moving), moving_inside.to(moving.dtype).reshape((1,) + geometry.shape)])
    channels = fixed_mind.shape[0]

    for iteration in range(iterations):
        if iteration > 0:
            field = _smooth(field, field_kernels)
        index = _source_index(geometry, geometry, device, displacement=field)
        inside = None
        for k, n in enumerate(geometry.size):
            axis = (index[..., k] >= -0.5) & (index[..., k] < n - 0.5)
            inside = axis if inside is None else inside & axis
            index[..., k].mul_(2.0 / max(n - 1, 1)).sub_(1.0)
        grid = index.unsqueeze(0)
        del index

        numerator = torch.zeros((3,) + geometry.shape, device=device, dtype=torch.float32)
        squared_gradient = torch.zeros(geometry.shape, device=device, dtype=torch.float32)
        squared_speed = torch.zeros(geometry.shape, device=device, dtype=torch.float32)
        for start in range(0, channels, MIND_CHUNK):
            stop = min(start + MIND_CHUNK, channels)
            warped = F.grid_sample(moving_mind[start:stop].unsqueeze(0), grid, mode="bilinear",
                                   padding_mode="border", align_corners=True).squeeze(0)
            speed = fixed_mind[start:stop] - warped
            squared_speed += (speed * speed).sum(0)
            for axis in range(3):
                gradient = (_index_gradient(fixed_mind[start:stop], axis, spacing[axis])
                            + _index_gradient(warped, axis, spacing[axis]))
                numerator[axis] += (speed * gradient).sum(0)
                squared_gradient += (gradient * gradient).sum(0)
                del gradient
            del warped, speed
        sampled_inside = F.grid_sample(moving_mind[channels:].unsqueeze(0), grid, mode="bilinear",
                                       padding_mode="border", align_corners=True)[0, 0] > 0.5
        del grid
        inside = inside & sampled_inside
        del sampled_inside

        # |D g| = |g| for a rotation D, so only the numerator needs turning into physical axes.
        numerator = _to_physical(numerator, geometry.direction)
        denominator = squared_gradient + squared_speed / normalizer
        del squared_gradient, squared_speed
        valid = inside & (denominator >= DENOMINATOR_THRESHOLD)
        scale = torch.where(valid, 2.0 / denominator.clamp_min(DENOMINATOR_THRESHOLD),
                            torch.zeros_like(denominator))
        update = numerator * scale
        del numerator, denominator, valid, scale

        counted = int(inside.sum())
        rms = None
        if counted:
            squared = torch.where(inside, (update * update).sum(0), torch.zeros((), device=device))
            rms = float(torch.sqrt(squared.sum(dtype=torch.float64) / counted))
            del squared
        del inside
        field = field + _smooth(update, update_kernels)
        del update
        if rms is not None and rms < MAXIMUM_RMS_ERROR:
            break
    return _smooth(field, field_kernels)


def _upload(image, device):
    array = sitk.GetArrayFromImage(sitk.Cast(image, sitk.sitkFloat32))
    return torch.from_numpy(array).to(device).unsqueeze(0)


def _pyramid(image, levels, sigmas_mm, device):
    """smooth_and_resample for every level; only one full-size copy is on the device at a time."""
    geometry = _Geometry.of(image)
    full = _upload(image, device)
    out = []
    for level, sigma in zip(levels, sigmas_mm):
        if sigma:
            sigma_voxels = [float(sigma) / sp for sp in geometry.spacing]
            width = int(max(8.0 * float(sigma) / sp for sp in geometry.spacing))
            smoothed = _smooth(full, _operator_kernels(sigma_voxels, PYRAMID_MAXIMUM_ERROR, width))
        else:
            smoothed = full
        out.append(_sample(smoothed, _source_index(geometry, level, device), geometry.size, 0.0))
        del smoothed
    del full
    return out


def multiscale_demons_torch(
    fixed_image,
    moving_image,
    resolution_staging,
    smoothing_sigmas,
    iteration_staging,
    isotropic_resample,
    regularization_kernel_mm,
    initial_displacement_field=None,
    default_value=None,
    device=None,
    metric="intensity",
):
    """
    Multiscale demons on PyTorch, mirroring ``multiscale_demons`` with the filter
    ``demons_registration`` configures.

    Parameters
    ----------
    fixed_image, moving_image : sitk.Image
        3D float images; the moving image may be on any grid.
    resolution_staging, smoothing_sigmas, iteration_staging, isotropic_resample
        As in ``multiscale_demons``.
    regularization_kernel_mm : array-like of float
        Displacement smoothing widths in mm, one per axis.
    initial_displacement_field : sitk.Image, optional
        Starting field, resampled onto the fixed grid.
    default_value : float, optional
        Padding for moving samples outside the moving image; inferred from raw intensities
        when omitted.
    device : str or torch.device, optional
        ``None`` uses CUDA when available, otherwise the CPU.
    metric : {"intensity", "mind"}
        What the demons force compares: intensities, as ITK does, or MIND-SSC descriptors,
        which makes it work across modalities.

    Returns
    -------
    sitk.Image
        Fixed -> moving displacement field (VectorFloat64) on the fixed grid.

    Raises
    ------
    torch.OutOfMemoryError
        When the device runs out of memory, with an estimate and the way out.
    """
    device = resolve_device(device)
    try:
        return _multiscale_demons_torch(
            fixed_image, moving_image, resolution_staging, smoothing_sigmas, iteration_staging,
            isotropic_resample, regularization_kernel_mm, initial_displacement_field,
            default_value, device, metric,
        )
    except torch.OutOfMemoryError as error:
        if device.type == "cuda":
            torch.cuda.empty_cache()
        finest = _Geometry.of(fixed_image).level(resolution_staging[-1], isotropic_resample)
        voxels = int(np.prod(finest.size))
        per_voxel = MIND_BYTES_PER_FINEST_VOXEL if metric == "mind" else BYTES_PER_FINEST_VOXEL
        raise torch.OutOfMemoryError(
            f"backend='torch' ran out of memory on {device}. The finest pyramid level has "
            f"{voxels / 1e6:.0f} M voxels and needs roughly {voxels * per_voxel / 1e9:.1f} GB. "
            f"Make that level coarser -- e.g. isotropic_resample=True with "
            f"resolution_staging=(8, 4, 1.5); the returned field is still on the fixed grid -- "
            f"or use backend='simpleitk'."
        ) from error


def _multiscale_demons_torch(fixed_image, moving_image, resolution_staging, smoothing_sigmas,
                             iteration_staging, isotropic_resample, regularization_kernel_mm,
                             initial_displacement_field, default_value, device, metric="intensity"):
    if default_value is None:
        default_value = -1000.0 if np.min(sitk.GetArrayViewFromImage(moving_image)) <= -1000 else 0.0
    kernel_mm = np.asarray(regularization_kernel_mm, dtype=np.float64)

    fixed_geometry, moving_geometry = _Geometry.of(fixed_image), _Geometry.of(moving_image)
    fixed_levels = [fixed_geometry.level(r, isotropic_resample) for r in resolution_staging]
    moving_levels = [moving_geometry.level(r, isotropic_resample) for r in resolution_staging]
    fixed_images = _pyramid(fixed_image, fixed_levels, smoothing_sigmas, device)
    moving_images = _pyramid(moving_image, moving_levels, smoothing_sigmas, device)

    if initial_displacement_field is not None:
        initial = sitk.Resample(sitk.Cast(initial_displacement_field, sitk.sitkVectorFloat64), fixed_image)
        total = torch.from_numpy(sitk.GetArrayFromImage(initial).astype(np.float32)).to(device)
        total = total.permute(3, 0, 1, 2).contiguous()
        total_geometry = fixed_geometry
    else:
        total, total_geometry = None, None

    for f_geometry, f_image, m_geometry, m_image, iterations in zip(
        fixed_levels, fixed_images, moving_levels, moving_images, iteration_staging
    ):
        field_sigma_voxels = (kernel_mm / f_geometry.spacing).tolist()
        if total is None:
            total = torch.zeros((3,) + f_geometry.shape, device=device, dtype=torch.float32)
        else:
            total = _sample(total, _source_index(total_geometry, f_geometry, device),
                            total_geometry.size, 0.0)
        total_geometry = f_geometry

        if metric == "mind":
            index = _source_index(m_geometry, f_geometry, device)
            moving_on_fixed, inside = _sample(m_image, index, m_geometry.size, default_value,
                                              return_inside=True)
            del index
            total = _mind_demons_level(f_image, moving_on_fixed, inside, f_geometry, int(iterations),
                                       field_sigma_voxels, device, total)
            del moving_on_fixed, inside
            continue

        index = _source_index(m_geometry, f_geometry, device, displacement=total)
        warped = _sample(m_image, index, m_geometry.size, default_value)
        del index
        refinement = _demons_level(f_image, warped, f_geometry, int(iterations), field_sigma_voxels, device)
        del warped

        # Moving was warped through T; the refinement D gives T(D(x)):
        # u_total(x) = v(x) + u(x + v(x)).
        index = _source_index(f_geometry, f_geometry, device, displacement=refinement)
        total = refinement + _sample(total, index, f_geometry.size, 0.0)
        del index, refinement
        # SmoothingRecursiveGaussian at the regularisation width, in the SimpleITK path.
        total = _smooth(total, _recursive_kernels(kernel_mm, f_geometry.spacing))

    field = total_geometry.to_image(total.permute(1, 2, 3, 0).cpu().numpy())
    del total, fixed_images, moving_images
    if device.type == "cuda":
        torch.cuda.empty_cache()  # return the cache to other processes sharing the card
    if total_geometry.matches(fixed_geometry):
        return field  # finest level is the fixed grid; resampling onto it changes nothing
    return sitk.Resample(field, fixed_image)
