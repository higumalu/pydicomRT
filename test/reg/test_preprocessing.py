"""
Registration preprocessing helpers.

The padding value is the subtle one: pad a CT with 0 instead of -1000 and you have invented
a shell of soft tissue around the patient, which the similarity metric will happily try to
match. Nothing errors, the registration just pulls toward the wrong answer.
"""

import numpy as np
import pytest
import SimpleITK as sitk

from synthetic import make_textured_phantom
from pydicomrt.reg.pipeline.preprocessing import (
    CT_AIR_HU,
    align_image_extents,
    crop_image_to_extent,
    get_image_center,
    get_image_physical_extent,
    get_images_distance,
    get_intersection_extent,
    infer_default_pixel_value,
    preprocess_image,
    window_clip,
)


def uniform_image(value, size=(10, 10, 10), spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    image = sitk.GetImageFromArray(np.full(size[::-1], value, dtype=np.float32))
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    return image


class TestInferDefaultPixelValue:
    def test_explicit_value_always_wins(self):
        ct_like = uniform_image(-1000.0)

        assert infer_default_pixel_value(ct_like, default_value=0.0) == 0.0
        assert infer_default_pixel_value(ct_like, default_value=-500.0) == -500.0

    def test_ct_like_image_pads_with_air(self, textured_phantom):
        assert infer_default_pixel_value(textured_phantom) == CT_AIR_HU

    def test_non_ct_image_pads_with_zero(self):
        mr_like = uniform_image(250.0)

        assert infer_default_pixel_value(mr_like) == 0.0

    def test_clipped_ct_no_longer_reads_as_ct(self, textured_phantom):
        """
        Documents the limit of the inference: window clipping raises the minimum above the
        air value, so a preprocessed CT is indistinguishable from an MR here. This is why
        the pipeline resolves the padding value on the raw images first.
        """
        clipped = window_clip(textured_phantom, [-10, 500])

        assert infer_default_pixel_value(clipped) == 0.0

    def test_vector_images_have_no_air_value(self):
        field = sitk.Image(4, 4, 4, sitk.sitkVectorFloat64)

        assert infer_default_pixel_value(field) == 0.0


class TestWindowClip:
    def test_clamps_to_the_requested_range(self, textured_phantom):
        clipped = window_clip(textured_phantom, [-10, 500])
        array = sitk.GetArrayViewFromImage(clipped)

        assert array.min() >= -10
        assert array.max() <= 500

    def test_range_order_does_not_matter(self, textured_phantom):
        ascending = sitk.GetArrayFromImage(window_clip(textured_phantom, [-10, 500]))
        descending = sitk.GetArrayFromImage(window_clip(textured_phantom, [500, -10]))

        assert np.array_equal(ascending, descending)

    def test_geometry_is_untouched(self, textured_phantom):
        clipped = window_clip(textured_phantom, [-10, 500])

        assert clipped.GetSize() == textured_phantom.GetSize()
        assert clipped.GetOrigin() == textured_phantom.GetOrigin()
        assert clipped.GetSpacing() == textured_phantom.GetSpacing()

    @pytest.mark.parametrize("bad_range", [[-10], [-10, 0, 10], "not a range", 5])
    def test_malformed_range_is_rejected(self, textured_phantom, bad_range):
        with pytest.raises(ValueError, match="two elements"):
            window_clip(textured_phantom, bad_range)

    def test_preprocess_image_without_config_is_a_no_op(self, textured_phantom):
        assert preprocess_image(textured_phantom, None) is textured_phantom
        assert preprocess_image(textured_phantom, {}) is textured_phantom


class TestPhysicalExtent:
    def test_extent_covers_the_voxel_centres(self):
        image = uniform_image(0.0, size=(10, 10, 10), spacing=(2.0, 3.0, 4.0), origin=(1.0, 2.0, 3.0))

        minimum, maximum = get_image_physical_extent(image)

        assert minimum == pytest.approx([1.0, 2.0, 3.0])
        assert maximum == pytest.approx([1.0 + 9 * 2.0, 2.0 + 9 * 3.0, 3.0 + 9 * 4.0])

    def test_centre_is_midway_between_the_corners(self):
        image = uniform_image(0.0, size=(10, 10, 10), spacing=(2.0, 2.0, 2.0), origin=(0.0, 0.0, 0.0))

        assert get_image_center(image) == pytest.approx([9.0, 9.0, 9.0])

    def test_distance_between_image_centres(self):
        first = uniform_image(0.0, origin=(0.0, 0.0, 0.0))
        second = uniform_image(0.0, origin=(3.0, 4.0, 0.0))

        assert get_images_distance(first, second) == pytest.approx(5.0)


class TestIntersection:
    def test_overlapping_images_report_their_shared_box(self):
        first = uniform_image(0.0, origin=(0.0, 0.0, 0.0))
        second = uniform_image(0.0, origin=(5.0, 0.0, 0.0))

        minimum, maximum, overlaps = get_intersection_extent(first, second)

        assert overlaps
        assert minimum == pytest.approx([5.0, 0.0, 0.0])
        assert maximum == pytest.approx([9.0, 9.0, 9.0])

    def test_disjoint_images_report_no_overlap(self):
        first = uniform_image(0.0, origin=(0.0, 0.0, 0.0))
        second = uniform_image(0.0, origin=(500.0, 0.0, 0.0))

        minimum, maximum, overlaps = get_intersection_extent(first, second)

        assert not overlaps
        assert minimum is None and maximum is None

    def test_align_extents_refuses_disjoint_images(self):
        first = uniform_image(0.0, origin=(0.0, 0.0, 0.0))
        second = uniform_image(0.0, origin=(500.0, 0.0, 0.0))

        with pytest.raises(ValueError, match="no intersection"):
            align_image_extents(first, second, use_initial_rigid=False)

    def test_align_extents_puts_both_images_on_one_grid(self):
        fixed = make_textured_phantom(origin=(0.0, 0.0, 0.0))
        moving = make_textured_phantom(origin=(10.0, 6.0, 0.0))

        fixed_aligned, moving_aligned, _ = align_image_extents(
            fixed, moving, use_initial_rigid=False
        )

        assert fixed_aligned.GetSize() == moving_aligned.GetSize()
        assert fixed_aligned.GetOrigin() == pytest.approx(moving_aligned.GetOrigin())
        assert fixed_aligned.GetSpacing() == pytest.approx(moving_aligned.GetSpacing())


class TestCropPadding:
    def test_crop_beyond_the_source_pads_with_the_inferred_air_value(self):
        ct_like = make_textured_phantom()
        minimum, maximum = get_image_physical_extent(ct_like)

        # Reach 20 mm past the far edge so the crop must invent voxels
        cropped = crop_image_to_extent(ct_like, minimum, maximum + np.array([20.0, 0.0, 0.0]))

        assert sitk.GetArrayViewFromImage(cropped)[:, :, -1].max() == pytest.approx(CT_AIR_HU)

    def test_explicit_padding_value_is_honoured(self):
        ct_like = make_textured_phantom()
        minimum, maximum = get_image_physical_extent(ct_like)

        cropped = crop_image_to_extent(
            ct_like, minimum, maximum + np.array([20.0, 0.0, 0.0]), default_value=0.0
        )

        assert sitk.GetArrayViewFromImage(cropped)[:, :, -1].max() == pytest.approx(0.0)
