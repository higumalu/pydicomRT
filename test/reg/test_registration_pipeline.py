"""
Registration pipeline contracts.

The pipeline hands back both an image and the transforms that produced it. If those two
disagree, downstream code that exports the transform to DICOM REG and code that reads the
returned image quietly diverge -- so the self-consistency test below is the important one.
"""

import numpy as np
import pytest
import SimpleITK as sitk

from synthetic import make_textured_phantom, phantom_centre
from pydicomrt.reg.method import rigid_registration
from pydicomrt.reg.pipeline import compose_transforms, registration_pipeline

CENTRE = phantom_centre()


def known_shift(shift) -> sitk.Euler3DTransform:
    transform = sitk.Euler3DTransform()
    transform.SetCenter([float(c) for c in CENTRE])
    transform.SetTranslation([float(v) for v in shift])
    return transform


def warp(image: sitk.Image, transform: sitk.Transform) -> sitk.Image:
    return sitk.Resample(image, image, transform.GetInverse(), sitk.sitkLinear, -1000.0)


class TestComposeTransforms:
    """
    SimpleITK evaluates a CompositeTransform back-to-front, so a stage-ordered list composes
    without reversing. Getting this backwards produces a transform that is wrong only when
    the stages do not commute -- i.e. as soon as there is any rotation.
    """

    def test_composition_applies_the_last_stage_first(self):
        first = sitk.TranslationTransform(3, (10.0, 0.0, 0.0))
        second = sitk.AffineTransform(3)
        second.SetMatrix((2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

        composed = compose_transforms([first, second])

        point = (1.0, 0.0, 0.0)
        assert composed.TransformPoint(point) == pytest.approx(
            first.TransformPoint(second.TransformPoint(point))
        )

    def test_composing_preserves_content_that_per_stage_resampling_discards(
        self, textured_phantom
    ):
        """
        Resampling once per stage crops to the fixed grid after *every* stage, so anything a
        stage pushes out of view is replaced by padding and can never come back -- even when
        a later stage would have brought it in. Composing first avoids that entirely.

        The two shifts here cancel, so the composed result is the original image and the
        difference is purely the data the staged path threw away.
        """
        out_of_view = known_shift((30.0, 0.0, 0.0))
        back_again = known_shift((-30.0, 0.0, 0.0))

        stage_by_stage = textured_phantom
        for transform in [out_of_view, back_again]:
            stage_by_stage = sitk.Resample(
                stage_by_stage, textured_phantom, transform, sitk.sitkLinear, -1000.0
            )
        in_one_go = sitk.Resample(
            textured_phantom,
            textured_phantom,
            compose_transforms([out_of_view, back_again]),
            sitk.sitkLinear,
            -1000.0,
        )

        original = sitk.GetArrayFromImage(textured_phantom)
        staged = sitk.GetArrayFromImage(stage_by_stage)
        composed = sitk.GetArrayFromImage(in_one_go)

        # The transforms cancel, so composing must reproduce the original image.
        assert np.abs(composed - original).max() < 1e-3

        # The staged path lost a band of anatomy to padding; the composed path did not.
        padding = -1000.0
        assert (staged <= padding).sum() > (composed <= padding).sum()
        assert np.abs(staged - original).max() > 100.0

    def test_skips_missing_stages(self):
        only = sitk.TranslationTransform(3, (1.0, 2.0, 3.0))

        assert compose_transforms([None, None]) is None
        assert compose_transforms([only, None]) is only
        assert compose_transforms([None, only]) is only

    def test_wraps_multiple_stages(self):
        composed = compose_transforms([
            sitk.TranslationTransform(3, (1.0, 0.0, 0.0)),
            sitk.TranslationTransform(3, (0.0, 2.0, 0.0)),
        ])

        assert isinstance(composed, sitk.CompositeTransform)
        assert composed.TransformPoint((0.0, 0.0, 0.0)) == pytest.approx((1.0, 2.0, 0.0))


class TestPipelineOutputContract:
    @pytest.fixture(scope="class")
    def rigid_only_result(self):
        fixed = make_textured_phantom()
        moving = warp(fixed, known_shift((-8.0, 4.0, 0.0)))
        result = registration_pipeline(
            fixed, moving, perform_rigid=True, perform_deformable=False
        )
        return (fixed, moving) + result

    def test_output_is_a_single_resample_of_the_untouched_moving_image(self, rigid_only_result):
        """
        The image handed back must be exactly what a caller gets by applying the returned
        transforms themselves. Pre-resampling the moving image onto the fixed grid before
        transforming it, or resampling once per stage, both break this.
        """
        fixed, moving, registered, rigid, deformable, _ = rigid_only_result

        expected = sitk.Resample(
            moving,
            fixed,
            compose_transforms([rigid, deformable]),
            sitk.sitkLinear,
            -1000.0,
            moving.GetPixelID(),
        )

        assert np.array_equal(
            sitk.GetArrayFromImage(registered), sitk.GetArrayFromImage(expected)
        )

    def test_output_sits_on_the_fixed_grid(self, rigid_only_result):
        fixed, _, registered, _, _, _ = rigid_only_result

        assert registered.GetSize() == fixed.GetSize()
        assert registered.GetOrigin() == pytest.approx(fixed.GetOrigin())
        assert registered.GetSpacing() == pytest.approx(fixed.GetSpacing())
        assert registered.GetDirection() == pytest.approx(fixed.GetDirection())

    def test_output_keeps_the_moving_pixel_type(self, rigid_only_result):
        _, moving, registered, _, _, _ = rigid_only_result

        assert registered.GetPixelID() == moving.GetPixelID()

    def test_rigid_stage_recovers_the_known_shift(self, rigid_only_result):
        _, _, _, rigid, _, _ = rigid_only_result
        truth = known_shift((-8.0, 4.0, 0.0))

        probes = [(0.0, 0.0, 0.0), (94.0, 94.0, 93.0), tuple(CENTRE)]
        error = np.mean([
            np.linalg.norm(np.array(rigid.TransformPoint(p)) - np.array(truth.TransformPoint(p)))
            for p in probes
        ])
        assert error < 5.0

    def test_deformable_stage_is_reported_as_absent(self, rigid_only_result):
        _, _, _, _, deformable, deformation_field = rigid_only_result

        assert deformable is None
        assert deformation_field is None


class TestPreprocessingDoesNotLeak:
    def test_window_clipping_only_affects_the_registration_not_the_output(self):
        """
        Clipping is a registration aid. If it reached the returned image the caller would
        silently lose the original intensities -- HU values they may go on to export.
        """
        fixed = make_textured_phantom()
        moving = warp(fixed, known_shift((-6.0, 0.0, 0.0)))

        registered, _, _, _ = registration_pipeline(
            fixed,
            moving,
            perform_rigid=True,
            perform_deformable=False,
            preprocess_config={"window_clip": [-10, 500]},
        )

        # Air and the +900 HU insert both survive; clipping would have crushed them to the
        # [-10, 500] window.
        array = sitk.GetArrayViewFromImage(registered)
        assert array.min() < -500
        assert array.max() > 600


class TestLargeInitialOffset:
    """
    The pipeline used to resample both images onto a shared grid before the rigid stage --
    a union extent, or the fixed grid after a pre-alignment step. Either one reduces
    ``CenteredTransformInitializer`` to a no-op, because aligning the centres of two
    identical grids yields zero, and the optimizer then has to find the whole offset alone.

    An offset larger than the images themselves is what makes that visible: it is exactly
    the situation the initializer exists for.
    """

    OFFSET = (0.0, 150.0, 0.0)

    @staticmethod
    def shifted(image, offset):
        moved = sitk.Image(image)
        moved.SetOrigin(tuple(o + s for o, s in zip(image.GetOrigin(), offset)))
        return moved

    def test_recovers_an_offset_larger_than_the_image(self, textured_phantom):
        moving = self.shifted(textured_phantom, self.OFFSET)

        _, rigid, _, _ = registration_pipeline(
            textured_phantom, moving, perform_rigid=True, perform_deformable=False
        )

        recovered = np.array(rigid.TransformPoint(CENTRE)) - np.array(CENTRE)
        assert recovered == pytest.approx(self.OFFSET, abs=2.0)

    def test_matches_calling_rigid_registration_directly(self, textured_phantom):
        """
        The pipeline is a convenience wrapper. Whatever preprocessing it arranges, it must
        not do worse on the registration itself than the function it wraps.
        """
        moving = self.shifted(textured_phantom, self.OFFSET)

        _, pipeline_rigid, _, _ = registration_pipeline(
            textured_phantom, moving, perform_rigid=True, perform_deformable=False
        )
        direct = rigid_registration(textured_phantom, moving)

        probes = [(0.0, 0.0, 0.0), (94.0, 94.0, 93.0), tuple(CENTRE)]
        difference = np.mean([
            np.linalg.norm(
                np.array(pipeline_rigid.TransformPoint(p))
                - np.array(direct.TransformPoint(p))
            )
            for p in probes
        ])
        assert difference < 1.0

    def test_output_still_lands_on_the_fixed_grid(self, textured_phantom):
        moving = self.shifted(textured_phantom, self.OFFSET)

        registered, _, _, _ = registration_pipeline(
            textured_phantom, moving, perform_rigid=True, perform_deformable=False
        )

        assert registered.GetSize() == textured_phantom.GetSize()
        assert registered.GetOrigin() == pytest.approx(textured_phantom.GetOrigin())


class TestPipelineValidation:
    def test_requires_at_least_one_stage(self, textured_phantom):
        with pytest.raises(ValueError, match="At least one"):
            registration_pipeline(
                textured_phantom,
                textured_phantom,
                perform_rigid=False,
                perform_deformable=False,
            )

    def test_accepts_the_nested_preprocess_config(self):
        fixed = make_textured_phantom()
        moving = warp(fixed, known_shift((-4.0, 0.0, 0.0)))

        registered, rigid, _, _ = registration_pipeline(
            fixed,
            moving,
            perform_rigid=True,
            perform_deformable=False,
            preprocess_config={"rigid": {"window_clip": [-200, 800]}},
        )

        assert rigid is not None
        assert registered.GetSize() == fixed.GetSize()
