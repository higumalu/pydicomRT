"""
B-spline deformable registration.

A B-spline transform has one coefficient per control point per axis, and the optimizer sees
only the voxels the metric samples. When the sampled voxels are fewer than the coefficients
the problem is underdetermined: the optimizer drives the metric to near zero by contorting
the grid, the result looks converged, and the deformation is nonsense everywhere the metric
did not look. Nothing raises. :func:`describe_sampling` reports that ratio per level and
:func:`bspline_registration` warns when a level is underdetermined.

Compared with ``demons_registration``, which regularises a dense field and is far more
forgiving, B-spline needs its resolution schedule, sampling and grid spacing chosen
together. Prefer demons unless you specifically need a parametric transform.

Measured on a synthetic phantom against a known B-spline deformation, recovering a 2.11 mm
mean displacement (lower is better), four threads:

===========================================================  =========  ========
configuration                                                error      time
===========================================================  =========  ========
previous defaults: staging (8, 4, 2), REGULAR sampling 10%    1.59 mm       1 s
only reaching full resolution: staging (4, 2, 1)              1.48 mm      17 s
only dense sampling                                           1.51 mm      16 s
both -- the current defaults                                  0.82 mm      17 s
``demons_registration`` defaults, for reference               0.80 mm       0.4 s
===========================================================  =========  ========

Neither change helps alone: the schedule and the sampling have to be fixed together.

**Sampling rate is not a speed dial.** It looks like one and it is not:

======================  =========  ========
sampling                error      time
======================  =========  ========
``None`` (dense path)    0.89 mm      37 s
``REGULAR`` at 100%      1.50 mm      47 s
``REGULAR`` at 50%       1.62 mm      19 s
``REGULAR`` at 10%       2.04 mm       4 s
``RANDOM`` at 10%        2.02 mm       5 s
``RANDOM`` at 1%         1.64 mm       1 s
======================  =========  ========

``REGULAR`` at 100% samples the same voxels as ``None`` and still scores 1.50 mm, so the
gap is the derivative code path rather than the sample count -- ITK evaluates a
local-support transform differently once a sampled point set is in play. Every sampled
configuration lands between 1.5 and 2.0 mm whatever the rate, and the ordering within that
band is noise. Sampling is a cliff, not a curve. Leave it at ``None``.

To trade accuracy for time, use the dials that do vary smoothly. Measured across three
different deformations, mean error by iteration cap: 5 -> 1.09 mm, 10 -> 1.01 mm,
20 -> 0.82 mm, 30 -> 0.88 mm, 50 -> 1.07 mm. Past about twenty the fit starts chasing
noise, so a higher cap costs time *and* accuracy. Coarsening ``initial_grid_spacing`` or
dropping a level trades similarly and predictably; threading is close to linear, taking the
default configuration from 111 s on one thread to 19 s on sixteen.
"""

import logging
from typing import List, Optional, Sequence, Tuple

import numpy as np
import SimpleITK as sitk

from pydicomrt.reg.pipeline.preprocessing import infer_default_pixel_value

from .common import (
    apply_transform,
    registration_command_iteration,
    smooth_and_resample,
    stage_iteration,
)

logger = logging.getLogger(__name__)

__all__ = [
    "bspline_registration",
    "control_point_spacing_distance_to_number",
    "describe_sampling",
    "OPTIMIZERS",
    "METRICS",
]

#: Optimizers this function knows how to configure.
#:
#: ``"LBFGSB"`` is bounded L-BFGS. ITK sizes its scales once, so it cannot be combined with
#: ``grid_scale_factors`` -- the coefficient count changes between levels and ITK raises
#: "Size of scales (N) must equal number of local parameters (M)". Pass
#: ``grid_scale_factors=None`` to use it.
OPTIMIZERS = ("LBFGS", "LBFGSB", "CGLS", "GRADIENT_DESCENT", "GRADIENT_DESCENT_LINE_SEARCH")

#: Similarity metrics this function knows how to configure.
#:
#: ITK's Demons metric is deliberately absent: it requires a displacement field transform
#: ("The moving transform must be a displacement field transform") and so can never pair
#: with a B-spline. It was listed as an option here and always threw. Use
#: ``demons_registration`` for that metric.
METRICS = ("mean_squares", "correlation", "mutual_information")

#: B-spline order used by ITK; a mesh of n cells carries n + 3 control points per axis.
_SPLINE_ORDER = 3


def control_point_spacing_distance_to_number(
    image: sitk.Image,
    grid_spacing,
    ) -> np.ndarray:
    """
    Convert a control point spacing in millimetres into a transform domain mesh size.

    Parameters
    ----------
    image : sitk.Image
        Image the grid will span.
    grid_spacing : float or sequence[float]
        Desired spacing between control points, in mm.

    Returns
    -------
    np.ndarray
        Mesh size per axis, at least 1 -- a mesh of 0 cells is not a transform.
    """
    image_spacing = np.array(image.GetSpacing())
    image_size = np.array(image.GetSize())
    number_points = image_size * image_spacing / np.array(grid_spacing)
    return np.maximum((number_points + 0.5).astype(int), 1)


def describe_sampling(
    fixed_image: sitk.Image,
    mesh_size,
    resolution_staging: Sequence[int],
    grid_scale_factors: Optional[Sequence[int]],
    sampling_rate: Optional[float],
    ) -> List[dict]:
    """
    Report metric samples against B-spline coefficients, per resolution level.

    Fewer samples than coefficients means the fit is underdetermined and the resulting
    deformation is unconstrained wherever the metric did not sample.

    Returns
    -------
    list[dict]
        One entry per level with ``level``, ``shrink_factor``, ``samples``, ``parameters``
        and ``underdetermined``.
    """
    size = np.array(fixed_image.GetSize())
    mesh_size = np.asarray(mesh_size)
    scales = list(grid_scale_factors) if grid_scale_factors else [1] * len(resolution_staging)
    rate = 1.0 if sampling_rate is None else float(sampling_rate)

    levels = []
    for index, (shrink, scale) in enumerate(zip(resolution_staging, scales)):
        voxels = float(np.prod(np.maximum(size // max(int(shrink), 1), 1)))
        samples = voxels * rate
        parameters = float(np.prod(mesh_size * scale + _SPLINE_ORDER) * fixed_image.GetDimension())
        levels.append({
            "level": index,
            "shrink_factor": int(shrink),
            "samples": samples,
            "parameters": parameters,
            "underdetermined": samples < parameters,
        })
    return levels


def _configure_optimizer(registration, optimizer, number_of_iterations, verbose):
    name = optimizer.upper()
    if name not in OPTIMIZERS:
        raise ValueError(f"optimizer must be one of {OPTIMIZERS}, got {optimizer!r}")

    if name == "LBFGSB":
        registration.SetOptimizerAsLBFGSB(
            gradientConvergenceTolerance=1e-5,
            numberOfIterations=number_of_iterations,
            maximumNumberOfCorrections=5,
            maximumNumberOfFunctionEvaluations=1024,
            costFunctionConvergenceFactor=1e7,
            trace=verbose,
        )
    elif name == "LBFGS":
        registration.SetOptimizerAsLBFGS2(numberOfIterations=number_of_iterations)
    elif name == "CGLS":
        registration.SetOptimizerAsConjugateGradientLineSearch(
            learningRate=0.05, numberOfIterations=number_of_iterations
        )
        registration.SetOptimizerScalesFromPhysicalShift()
    elif name == "GRADIENT_DESCENT":
        registration.SetOptimizerAsGradientDescent(
            learningRate=5.0,
            numberOfIterations=number_of_iterations,
            convergenceMinimumValue=1e-6,
            convergenceWindowSize=10,
        )
        registration.SetOptimizerScalesFromPhysicalShift()
    else:  # GRADIENT_DESCENT_LINE_SEARCH
        registration.SetOptimizerAsGradientDescentLineSearch(
            learningRate=1.0, numberOfIterations=number_of_iterations
        )
        registration.SetOptimizerScalesFromPhysicalShift()


def _configure_metric(registration, metric, histogram_bins):
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}, got {metric!r}")

    if metric == "correlation":
        registration.SetMetricAsCorrelation()
    elif metric == "mean_squares":
        registration.SetMetricAsMeanSquares()
    else:  # mutual_information
        registration.SetMetricAsMattesMutualInformation(numberOfHistogramBins=histogram_bins)


def bspline_registration(
    fixed_image: sitk.Image,
    moving_image: sitk.Image,
    fixed_mask: Optional[sitk.Image] = None,
    moving_mask: Optional[sitk.Image] = None,
    resolution_staging: Sequence[int] = (4, 2, 1),
    smoothing_sigmas: Sequence[float] = (4, 2, 0),
    sampling_rate: Optional[float] = None,
    optimizer: str = "LBFGS",
    metric: str = "mean_squares",
    initial_grid_spacing: float = 64,
    grid_scale_factors: Optional[Sequence[int]] = (1, 2, 4),
    interpolator: int = sitk.sitkBSpline,
    default_value: Optional[float] = None,
    number_of_iterations: int = 20,
    isotropic_resample: bool = False,
    initial_isotropic_size: float = 1,
    histogram_bins: int = 30,
    verbose: bool = False,
    ncores: int = 1,
    ) -> Tuple[sitk.Image, sitk.Transform, sitk.Image]:
    """
    Deformable registration with a B-spline transform.

    Returns the same triple as :func:`demons_registration` so the two are interchangeable.
    The returned transform is a *resampling* transform: it maps a point on the fixed grid
    back into the moving image.

    Parameters
    ----------
    fixed_image, moving_image : sitk.Image
        Images to register. Bring them into rough alignment first -- this stage refines, it
        does not search.
    fixed_mask, moving_mask : sitk.Image, optional
        Binary masks restricting where the metric is evaluated, in the fixed and moving
        image respectively.
    resolution_staging : sequence[int]
        Downsampling factor per level, coarsest first. The last entry should be 1 so the
        finest level sees the image at full resolution; the previous default stopped at 2,
        which meant it never did.
    smoothing_sigmas : sequence[float]
        Gaussian sigma per level, in **millimetres** -- not multiples of the shrink factor.
    sampling_rate : float, optional
        Fraction of voxels the metric samples. ``None`` (the default) selects ITK's dense
        path and is the only setting that registers well here -- any sampled configuration
        loses most of the accuracy regardless of rate. Reduce ``number_of_iterations`` or
        coarsen ``initial_grid_spacing`` for speed instead; see the module docstring.
    optimizer : str
        One of :data:`OPTIMIZERS`. ``"LBFGSB"`` requires ``grid_scale_factors=None``.
    metric : str
        One of :data:`METRICS`.
    initial_grid_spacing : float
        Control point spacing at the coarsest level, in mm.
    grid_scale_factors : sequence[int], optional
        Mesh multiplier per level, so the grid refines as the resolution does. ``None``
        keeps one grid throughout, which is what ``"LBFGSB"`` needs.
    interpolator : int
        Interpolator for the final resampling.
    default_value : float, optional
        Padding value. Inferred from the moving image when omitted (CT-like images use
        -1000).
    number_of_iterations : int
        Maximum optimizer iterations per level. Past roughly twenty the fit begins chasing
        noise, so raising this costs time and accuracy together -- see the module docstring.
    isotropic_resample : bool
        Resample both images to isotropic voxels before registering.
    initial_isotropic_size : float
        Voxel size in mm when ``isotropic_resample`` is set.
    histogram_bins : int
        Bins for ``metric="mutual_information"``.
    verbose : bool
        Log per-iteration metric values.
    ncores : int
        Registration threads. Speedup is close to linear, but ITK reduces the metric in
        thread completion order, so only ``ncores=1`` is bit-reproducible.

    Returns
    -------
    registered_image : sitk.Image
        Moving image resampled onto the fixed grid.
    output_transform : sitk.Transform
        Displacement field transform, fixed -> moving.
    deformation_field : sitk.Image
        The displacement field itself, in mm.

    Raises
    ------
    ValueError
        If the optimizer or metric name is unknown, the resolution schedules disagree, or
        ``"LBFGSB"`` is combined with ``grid_scale_factors``.
    """
    if len(resolution_staging) != len(smoothing_sigmas):
        raise ValueError(
            "resolution_staging and smoothing_sigmas must describe the same number of levels "
            f"(got {len(resolution_staging)} and {len(smoothing_sigmas)})"
        )
    if grid_scale_factors is not None and len(grid_scale_factors) != len(resolution_staging):
        raise ValueError(
            "grid_scale_factors must have one entry per resolution level "
            f"(got {len(grid_scale_factors)} and {len(resolution_staging)})"
        )
    if ncores <= 0:
        raise ValueError("`ncores` must be a positive integer.")
    if optimizer.upper() == "LBFGSB" and grid_scale_factors is not None:
        raise ValueError(
            "optimizer='LBFGSB' cannot be combined with grid_scale_factors: ITK sizes the "
            "optimizer scales once, but a scaled grid changes the coefficient count between "
            "levels. Pass grid_scale_factors=None, or use optimizer='LBFGS'."
        )

    default_value = infer_default_pixel_value(moving_image, default_value)

    fixed_image = sitk.Cast(fixed_image, sitk.sitkFloat32)
    moving_image_type = moving_image.GetPixelID()
    moving_image = sitk.Cast(moving_image, sitk.sitkFloat32)

    # (Optional) isotropic resample. Changes the sampling geometry, so it is off by default.
    fixed_image_original = fixed_image
    if isotropic_resample:
        fixed_image_original.MakeUnique()
        fixed_image = smooth_and_resample(
            fixed_image, isotropic_voxel_size_mm=initial_isotropic_size
        )
        moving_image = smooth_and_resample(
            moving_image, isotropic_voxel_size_mm=initial_isotropic_size
        )

    mesh_size = control_point_spacing_distance_to_number(fixed_image, initial_grid_spacing)

    levels = describe_sampling(
        fixed_image, mesh_size, resolution_staging, grid_scale_factors, sampling_rate
    )
    for level in levels:
        if level["underdetermined"]:
            logger.warning(
                "B-spline level %d (shrink %d) samples about %.0f voxels for %.0f "
                "coefficients: the fit is underdetermined and the deformation will be "
                "unconstrained away from the sampled voxels. Raise sampling_rate, coarsen "
                "initial_grid_spacing, or drop the coarsest level.",
                level["level"], level["shrink_factor"], level["samples"], level["parameters"],
            )

    registration = sitk.ImageRegistrationMethod()
    registration.SetNumberOfThreads(int(ncores))
    registration.SetShrinkFactorsPerLevel(list(resolution_staging))
    registration.SetSmoothingSigmasPerLevel(list(smoothing_sigmas))
    registration.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    _configure_optimizer(registration, optimizer, number_of_iterations, verbose)
    _configure_metric(registration, metric, histogram_bins)
    registration.SetInterpolator(sitk.sitkLinear)

    if sampling_rate is None:
        registration.SetMetricSamplingStrategy(sitk.ImageRegistrationMethod.NONE)
    else:
        registration.SetMetricSamplingStrategy(sitk.ImageRegistrationMethod.REGULAR)
        if isinstance(sampling_rate, (int, float)):
            registration.SetMetricSamplingPercentage(float(sampling_rate))
        else:
            registration.SetMetricSamplingPercentagePerLevel(list(sampling_rate))

    if moving_mask is not None:
        registration.SetMetricMovingMask(moving_mask)
    if fixed_mask is not None:
        registration.SetMetricFixedMask(fixed_mask)

    if verbose:
        logger.debug("Initial B-spline mesh size: %s", mesh_size)

    initial_transform = sitk.BSplineTransformInitializer(
        fixed_image, transformDomainMeshSize=[int(i) for i in mesh_size]
    )
    if grid_scale_factors is None:
        registration.SetInitialTransform(initial_transform, inPlace=True)
    else:
        registration.SetInitialTransformAsBSpline(
            initial_transform, inPlace=True, scaleFactors=list(grid_scale_factors)
        )

    if verbose:
        registration.AddCommand(
            sitk.sitkIterationEvent,
            lambda: registration_command_iteration(registration),
        )
        registration.AddCommand(
            sitk.sitkMultiResolutionIterationEvent,
            lambda: stage_iteration(registration),
        )

    output_transform = registration.Execute(fixed=fixed_image, moving=moving_image)

    registered_image = apply_transform(
        input_image=moving_image,
        reference_image=fixed_image_original,
        transform=output_transform,
        default_value=default_value,
        interpolator=interpolator,
    )
    registered_image = sitk.Cast(registered_image, moving_image_type)

    deformation_field = sitk.TransformToDisplacementField(
        output_transform,
        sitk.sitkVectorFloat64,
        fixed_image_original.GetSize(),
        fixed_image_original.GetOrigin(),
        fixed_image_original.GetSpacing(),
        fixed_image_original.GetDirection(),
    )
    # DisplacementFieldTransform takes ownership of the image it is given, leaving the
    # original empty, so it is handed a copy and the caller keeps a usable field.
    output_transform = sitk.DisplacementFieldTransform(sitk.Image(deformation_field))

    return registered_image, output_transform, deformation_field
