"""
Rigid registration accuracy and optimizer sanity.

The failure mode these guard against is silent: the optimizer overshoots, the convergence
window declares success, and a plausible-looking transform comes back tens of millimetres
off. Nothing raises, so only a ground-truth comparison catches it.
"""

import numpy as np
import pytest
import SimpleITK as sitk

from synthetic import make_textured_phantom, phantom_centre
from pydicomrt.reg.method.rigid import rigid_registration, to_centre_free_affine

CENTRE = phantom_centre()
# Sample the mapping at the corners, not just the centre: a rotation error is invisible at
# the centre of rotation and largest at the periphery.
PROBE_POINTS = [
    (0.0, 0.0, 0.0),
    (94.0, 0.0, 0.0),
    (0.0, 94.0, 0.0),
    (0.0, 0.0, 93.0),
    (94.0, 94.0, 93.0),
    tuple(CENTRE),
]


def known_rigid_transform(angle_deg=0.0, shift=(0.0, 0.0, 0.0)) -> sitk.Euler3DTransform:
    transform = sitk.Euler3DTransform()
    transform.SetCenter([float(c) for c in CENTRE])
    transform.SetRotation(0.0, 0.0, np.deg2rad(angle_deg))
    transform.SetTranslation([float(v) for v in shift])
    return transform


def warp(image: sitk.Image, transform: sitk.Transform) -> sitk.Image:
    """Produce a moving image for which ``transform`` is the fixed -> moving mapping."""
    return sitk.Resample(image, image, transform.GetInverse(), sitk.sitkLinear, -1000.0)


def mapping_error_mm(got: sitk.Transform, expected: sitk.Transform) -> float:
    return float(np.mean([
        np.linalg.norm(np.array(got.TransformPoint(p)) - np.array(expected.TransformPoint(p)))
        for p in PROBE_POINTS
    ]))


class TestCentreFreeAffine:
    """The DICOM REG matrix has no field for a rotation centre, so it must be folded in."""

    def test_conversion_preserves_the_mapping_exactly(self):
        centred = known_rigid_transform(angle_deg=20.0, shift=(5.0, -3.0, 2.0))

        affine = to_centre_free_affine(centred)

        assert affine.GetCenter() == (0.0, 0.0, 0.0)
        for point in PROBE_POINTS:
            assert affine.TransformPoint(point) == pytest.approx(
                centred.TransformPoint(point), abs=1e-9
            )

    def test_dropping_the_centre_would_be_wrong_under_rotation(self):
        """Confirms the test above is actually load-bearing."""
        centred = known_rigid_transform(angle_deg=20.0, shift=(5.0, -3.0, 2.0))

        naive = sitk.AffineTransform(3)
        naive.SetMatrix(centred.GetMatrix())
        naive.SetTranslation(centred.GetTranslation())

        assert mapping_error_mm(naive, centred) > 1.0

    def test_identity_rotation_is_unaffected(self):
        centred = known_rigid_transform(angle_deg=0.0, shift=(5.0, -3.0, 2.0))

        affine = to_centre_free_affine(centred)

        assert affine.GetTranslation() == pytest.approx((5.0, -3.0, 2.0))


class TestRigidAccuracy:
    @pytest.mark.parametrize(
        "angle_deg, shift",
        [
            (0.0, (-8.0, 0.0, 0.0)),
            (0.0, (-20.0, 10.0, 5.0)),
            (5.0, (-4.0, 3.0, 0.0)),
            (-8.0, (15.0, -9.0, 7.0)),
        ],
    )
    def test_recovers_a_known_rigid_transform(self, textured_phantom, angle_deg, shift):
        truth = known_rigid_transform(angle_deg, shift)
        moving = warp(textured_phantom, truth)

        recovered = rigid_registration(textured_phantom, moving)

        # Worst measured case is 0.78 mm; 2 mm leaves room for interpolation differences
        # without letting a real regression through.
        assert mapping_error_mm(recovered, truth) < 2.0

    def test_identical_images_give_an_identity_transform(self, textured_phantom):
        recovered = rigid_registration(textured_phantom, textured_phantom)

        assert mapping_error_mm(recovered, sitk.Transform(3, sitk.sitkIdentity)) < 2.0

    def test_does_not_walk_away_from_an_exact_start(self, textured_phantom):
        """
        When the two images differ only by a whole-grid shift, the centred initializer hands
        the optimizer the exact answer. Plain gradient descent took it from there to 71 mm
        out and still reported convergence -- the failure that motivated the regular-step
        default, and one an inexact starting point hides.
        """
        shift = (8.0, -6.0, 4.0)
        shifted = sitk.Image(textured_phantom)
        shifted.SetOrigin(tuple(o + s for o, s in zip(textured_phantom.GetOrigin(), shift)))
        truth = sitk.TranslationTransform(3, shift)

        recovered = rigid_registration(textured_phantom, shifted)

        assert mapping_error_mm(recovered, truth) < 2.0

    def test_works_when_both_images_share_a_grid(self, textured_phantom):
        """
        The centred initializer only aligns *grid* centres, so it contributes nothing once
        both images sit on the same grid -- everything then rests on the optimizer.
        """
        truth = known_rigid_transform(shift=(-12.0, 6.0, -4.0))
        moving = warp(textured_phantom, truth)
        assert moving.GetOrigin() == textured_phantom.GetOrigin()

        recovered = rigid_registration(textured_phantom, moving)

        assert mapping_error_mm(recovered, truth) < 5.0


class TestOptimizerDoesNotRegress:
    def test_result_is_never_worse_than_its_initialisation(self, textured_phantom):
        """
        The regression that motivated the bounded multi-resolution defaults: starting from a
        perfect alignment, the optimizer walked away and still reported convergence.
        """
        truth = known_rigid_transform(shift=(-8.0, 0.0, 0.0))
        moving = warp(textured_phantom, truth)

        metric = sitk.ImageRegistrationMethod()
        metric.SetMetricAsMattesMutualInformation(numberOfHistogramBins=100)
        metric.SetInterpolator(sitk.sitkLinear)

        def score(transform):
            resampled = sitk.Resample(moving, textured_phantom, transform, sitk.sitkLinear, -1000.0)
            return metric.MetricEvaluate(
                sitk.Cast(textured_phantom, sitk.sitkFloat32),
                sitk.Cast(resampled, sitk.sitkFloat32),
            )

        initial_score = score(sitk.CenteredTransformInitializer(
            sitk.Cast(textured_phantom, sitk.sitkFloat32),
            sitk.Cast(moving, sitk.sitkFloat32),
            sitk.VersorRigid3DTransform(),
            sitk.CenteredTransformInitializerFilter.GEOMETRY,
        ))
        final_score = score(rigid_registration(textured_phantom, moving))

        # Mattes MI is returned negated, so lower is better.
        assert final_score <= initial_score + 1e-3


class TestOptimizerChoice:
    def test_gradient_descent_is_still_selectable(self, textured_phantom):
        truth = known_rigid_transform(shift=(-8.0, 0.0, 0.0))
        moving = warp(textured_phantom, truth)

        recovered = rigid_registration(
            textured_phantom, moving, optimizer="gradient_descent"
        )

        # The fast path is offered for speed, not accuracy -- a loose bound on purpose.
        assert mapping_error_mm(recovered, truth) < 10.0

    def test_unknown_optimizer_is_rejected(self, textured_phantom):
        with pytest.raises(ValueError, match="regular_step"):
            rigid_registration(textured_phantom, textured_phantom, optimizer="adam")


class TestParameterValidation:
    def test_mismatched_resolution_schedules_are_rejected(self, textured_phantom):
        with pytest.raises(ValueError, match="resolution levels"):
            rigid_registration(
                textured_phantom,
                textured_phantom,
                shrink_factors=(4, 2, 1),
                smoothing_sigmas=(2.0, 0.0),
            )

    def test_single_resolution_is_still_available(self, textured_phantom):
        truth = known_rigid_transform(shift=(-4.0, 0.0, 0.0))
        moving = warp(textured_phantom, truth)

        recovered = rigid_registration(
            textured_phantom, moving, shrink_factors=None, smoothing_sigmas=None
        )

        assert isinstance(recovered, sitk.AffineTransform)
