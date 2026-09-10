"""
B-spline deformable registration.

The failure this module is prone to is silent: with fewer metric samples than B-spline
coefficients the optimizer drives the metric to near zero by contorting the grid, reports
convergence, and returns a deformation that is unconstrained wherever the metric did not
look. The sampling guard and the argument validation are what make that visible.

A small phantom and few iterations keep these fast; accuracy itself is checked once, loosely,
against a known deformation.
"""

import logging

import numpy as np
import pytest
import SimpleITK as sitk

from synthetic import make_textured_phantom
from pydicomrt.reg.method import bspline_registration, demons_registration
from pydicomrt.reg.method.bspline import (
    METRICS,
    OPTIMIZERS,
    control_point_spacing_distance_to_number,
    describe_sampling,
)

def quick(**overrides):
    """
    Few iterations on a coarse schedule: these check wiring and guards, not convergence.
    The default dense sampling makes even a five-iteration run costly at full resolution.
    """
    # (4, 2) rather than anything coarser: ITK's Gaussian smoothing needs at least four
    # voxels along each axis, and shrinking this phantom by 8 leaves fewer.
    settings = dict(
        number_of_iterations=5,
        resolution_staging=(4, 2),
        smoothing_sigmas=(4, 2),
        grid_scale_factors=(1, 2),
    )
    settings.update(overrides)
    return settings


def make_deformed_pair(size=(24, 24, 16), spacing=(4.0, 4.0, 6.0), amplitude=6.0):
    """
    A phantom and a copy of it warped by a known B-spline transform.

    ``fixed`` is built by resampling ``moving`` through ``truth``, so ``truth`` is exactly
    the fixed -> moving mapping a registration should recover. Warping the fixed image
    rather than the moving one avoids needing an inverse, which a B-spline does not have.
    """
    moving = make_textured_phantom(size=size, spacing=spacing)

    rng = np.random.default_rng(0)
    truth = sitk.BSplineTransformInitializer(moving, [4, 4, 4])
    truth.SetParameters(
        (np.array(truth.GetParameters()) + rng.uniform(-amplitude, amplitude, len(truth.GetParameters()))).tolist()
    )
    fixed = sitk.Resample(moving, moving, truth, sitk.sitkLinear, -1000.0)
    return fixed, moving, truth


@pytest.fixture(scope="module")
def deformed_pair():
    """Small pair, for the wiring and guard checks."""
    return make_deformed_pair()


def mapping_error_mm(transform, truth, reference):
    size = reference.GetSize()
    points = [
        reference.TransformIndexToPhysicalPoint((x, y, z))
        for x in (size[0] // 5, size[0] // 2, 4 * size[0] // 5)
        for y in (size[1] // 5, size[1] // 2, 4 * size[1] // 5)
        for z in (size[2] // 5, size[2] // 2, 4 * size[2] // 5)
    ]
    return float(np.mean([
        np.linalg.norm(np.array(transform.TransformPoint(p)) - np.array(truth.TransformPoint(p)))
        for p in points
    ]))


class TestReturnContract:
    def test_returns_the_same_triple_as_demons(self, deformed_pair):
        """
        The two deformable methods have to be interchangeable. bspline_registration used to
        return two values against demons' three, so it could not be dropped in.
        """
        fixed, moving, _ = deformed_pair

        bspline_result = bspline_registration(fixed, moving, **quick())
        demons_result = demons_registration(
            fixed, moving, resolution_staging=(4, 2, 1), iteration_staging=(1, 1, 1)
        )

        assert len(bspline_result) == len(demons_result) == 3
        for bspline_item, demons_item in zip(bspline_result, demons_result):
            assert type(bspline_item) is type(demons_item)

    def test_registered_image_lands_on_the_fixed_grid(self, deformed_pair):
        fixed, moving, _ = deformed_pair

        registered, _, _ = bspline_registration(fixed, moving, **quick())

        assert registered.GetSize() == fixed.GetSize()
        assert registered.GetOrigin() == pytest.approx(fixed.GetOrigin())

    def test_deformation_field_covers_the_fixed_grid(self, deformed_pair):
        fixed, moving, _ = deformed_pair

        _, _, field = bspline_registration(fixed, moving, **quick())

        assert field.GetSize() == fixed.GetSize()
        assert field.GetNumberOfComponentsPerPixel() == 3


class TestSamplingGuard:
    """
    Samples versus coefficients is what governs whether the fit means anything, and it is
    not visible from the arguments -- it comes out of the shrink factor, the grid scale and
    the sampling rate together.
    """

    def test_reports_a_ratio_per_level(self, deformed_pair):
        fixed, _, _ = deformed_pair
        mesh = control_point_spacing_distance_to_number(fixed, 64)

        levels = describe_sampling(fixed, mesh, (4, 2, 1), (1, 2, 4), None)

        assert [level["level"] for level in levels] == [0, 1, 2]
        assert all(level["samples"] > 0 for level in levels)
        assert all(level["parameters"] > 0 for level in levels)

    def test_flags_the_configuration_that_used_to_be_the_default(self, deformed_pair):
        """staging (8, 4, 2) at 10% sampling left ~14 samples for 375 coefficients."""
        fixed, _, _ = deformed_pair
        mesh = control_point_spacing_distance_to_number(fixed, 64)

        levels = describe_sampling(fixed, mesh, (8, 4, 2), (1, 2, 4), 0.1)

        assert all(level["underdetermined"] for level in levels)
        assert levels[0]["samples"] < 20
        assert levels[0]["parameters"] > 300

    def test_dense_sampling_at_full_resolution_is_not_underdetermined(self, deformed_pair):
        fixed, _, _ = deformed_pair
        mesh = control_point_spacing_distance_to_number(fixed, 64)

        levels = describe_sampling(fixed, mesh, (4, 2, 1), (1, 2, 4), None)

        assert levels[-1]["underdetermined"] is False

    def test_warns_when_a_level_is_underdetermined(self, deformed_pair, caplog):
        fixed, moving, _ = deformed_pair

        with caplog.at_level(logging.WARNING, logger="pydicomrt.reg.method.bspline"):
            bspline_registration(fixed, moving, **quick(
                resolution_staging=(4, 2, 1), smoothing_sigmas=(4, 2, 1),
                grid_scale_factors=(1, 2, 4), sampling_rate=0.1,
            ))

        assert "underdetermined" in caplog.text

    def test_stays_quiet_when_the_fit_is_determined(self, deformed_pair, caplog):
        """
        Well-determined depends on image size as much as on settings: this phantom is 96 mm
        across, so the default 64 mm grid leaves a coarse level underdetermined however the
        sampling is set. One full-resolution level with a coarse grid is comfortably
        determined, and the guard should then say nothing.
        """
        fixed, moving, _ = deformed_pair

        with caplog.at_level(logging.WARNING, logger="pydicomrt.reg.method.bspline"):
            bspline_registration(
                fixed, moving,
                number_of_iterations=5,
                resolution_staging=(1,), smoothing_sigmas=(0,), grid_scale_factors=None,
                initial_grid_spacing=96,
            )

        assert "underdetermined" not in caplog.text


class TestArgumentValidation:
    """These all used to fall through silently or fail deep inside ITK."""

    def test_unknown_optimizer_is_rejected(self, deformed_pair):
        fixed, moving, _ = deformed_pair

        with pytest.raises(ValueError, match="optimizer must be one of"):
            bspline_registration(fixed, moving, **quick(optimizer="lbgfs"))

    def test_unknown_metric_is_rejected(self, deformed_pair):
        fixed, moving, _ = deformed_pair

        with pytest.raises(ValueError, match="metric must be one of"):
            bspline_registration(fixed, moving, **quick(metric="mutual_info"))

    def test_lbfgsb_with_grid_scale_factors_is_rejected_with_an_explanation(self, deformed_pair):
        """
        ITK sizes the optimizer scales once, so a scaled grid makes LBFGSB fail with
        "Size of scales (375) must equal number of local parameters (1029)". The option was
        documented and always threw.
        """
        fixed, moving, _ = deformed_pair

        with pytest.raises(ValueError, match="LBFGSB.*grid_scale_factors"):
            bspline_registration(fixed, moving, **quick(optimizer="LBFGSB"))

    def test_lbfgsb_works_without_grid_scale_factors(self, deformed_pair):
        fixed, moving, _ = deformed_pair

        registered, transform, _ = bspline_registration(
            fixed, moving, **quick(optimizer="LBFGSB", grid_scale_factors=None)
        )

        assert registered.GetSize() == fixed.GetSize()

    def test_mismatched_resolution_schedules_are_rejected(self, deformed_pair):
        fixed, moving, _ = deformed_pair

        with pytest.raises(ValueError, match="same number of levels"):
            bspline_registration(
                fixed, moving,
                **quick(resolution_staging=(4, 2, 1), smoothing_sigmas=(4, 2),
                        grid_scale_factors=(1, 2, 4))
            )

    def test_mismatched_grid_scale_factors_are_rejected(self, deformed_pair):
        fixed, moving, _ = deformed_pair

        with pytest.raises(ValueError, match="one entry per resolution level"):
            bspline_registration(
                fixed, moving,
                **quick(resolution_staging=(4, 2, 1), smoothing_sigmas=(4, 2, 0),
                        grid_scale_factors=(1, 2))
            )

    @pytest.mark.parametrize("optimizer", OPTIMIZERS)
    def test_every_documented_optimizer_runs(self, deformed_pair, optimizer):
        """A documented option that always throws is worse than an undocumented one."""
        fixed, moving, _ = deformed_pair
        # quick() runs two levels, and grid_scale_factors has to match that length.
        scale_factors = None if optimizer == "LBFGSB" else (1, 2)

        registered, _, _ = bspline_registration(
            fixed, moving, **quick(optimizer=optimizer, grid_scale_factors=scale_factors)
        )

        assert registered.GetSize() == fixed.GetSize()

    @pytest.mark.parametrize("metric", METRICS)
    def test_every_documented_metric_runs(self, deformed_pair, metric):
        fixed, moving, _ = deformed_pair

        registered, _, _ = bspline_registration(fixed, moving, **quick(metric=metric))

        assert registered.GetSize() == fixed.GetSize()


class TestMeshSize:
    def test_spacing_converts_to_a_mesh(self, deformed_pair):
        fixed, _, _ = deformed_pair  # 96 mm on every axis

        assert list(control_point_spacing_distance_to_number(fixed, 96)) == [1, 1, 1]
        assert list(control_point_spacing_distance_to_number(fixed, 32)) == [3, 3, 3]
        assert list(control_point_spacing_distance_to_number(fixed, 16)) == [6, 6, 6]

    def test_never_returns_a_degenerate_mesh(self, deformed_pair):
        """A grid coarser than the image gave a mesh of 0, which is not a transform."""
        fixed, _, _ = deformed_pair

        assert list(control_point_spacing_distance_to_number(fixed, 10_000)) == [1, 1, 1]


@pytest.mark.slow
class TestAccuracy:
    """Full-size and full defaults, which the dense metric makes slow."""

    def test_recovers_most_of_a_known_deformation(self):
        """
        Runs the shipped defaults, so it measures what a caller actually gets. Loose on
        purpose: the point is the direction of travel, since the previous defaults left
        1.59 mm of a 2.11 mm deformation -- barely better than doing nothing.
        """
        fixed, moving, truth = make_deformed_pair(size=(48, 48, 32), spacing=(2.0, 2.0, 3.0))

        _, transform, _ = bspline_registration(fixed, moving, ncores=4)

        identity = sitk.Transform(3, sitk.sitkIdentity)
        before = mapping_error_mm(identity, truth, moving)
        after = mapping_error_mm(transform, truth, moving)

        assert after < before / 2


class TestSpeedDials:
    """
    Which knobs actually trade accuracy for time, measured rather than assumed.

    The sampling rate looks like the obvious dial and is not one: every sampled
    configuration loses most of the accuracy whatever the rate, because ITK evaluates a
    local-support transform through a different derivative path once a point set is in
    play. These pin the shape of that finding so a future "let's sample less, it'll be
    faster" change has to confront it.
    """

    def test_dense_sampling_is_the_default(self):
        import inspect

        default = inspect.signature(bspline_registration).parameters["sampling_rate"].default

        assert default is None

    def test_iteration_cap_defaults_below_the_overfitting_point(self):
        """
        Measured across three deformations: 20 iterations averaged 0.82 mm, 50 averaged
        1.07 mm. A higher cap costs time and accuracy together, so the default sits at the
        turn rather than as high as it will go.
        """
        import inspect

        default = inspect.signature(bspline_registration).parameters["number_of_iterations"].default

        assert default <= 30

    def test_thread_count_is_selectable_and_validated(self, deformed_pair):
        """
        Threading is close to linear here, and ITK reduces the metric in thread completion
        order, so the count has to be the caller's to choose.
        """
        fixed, moving, _ = deformed_pair

        registered, _, _ = bspline_registration(fixed, moving, **quick(ncores=2))
        assert registered.GetSize() == fixed.GetSize()

        with pytest.raises(ValueError, match="positive integer"):
            bspline_registration(fixed, moving, **quick(ncores=0))

    def test_single_thread_is_reproducible(self, deformed_pair):
        fixed, moving, _ = deformed_pair

        first = bspline_registration(fixed, moving, **quick(ncores=1))[1]
        second = bspline_registration(fixed, moving, **quick(ncores=1))[1]

        probe = fixed.TransformIndexToPhysicalPoint((12, 12, 8))
        assert first.TransformPoint(probe) == pytest.approx(second.TransformPoint(probe))
