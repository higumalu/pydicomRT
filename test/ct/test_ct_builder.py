"""
CT Image builder.

Geometry is the whole job here: a slice written with the row and column spacings swapped, or
with the plane normal on the wrong axis, still loads and still displays. It is only wrong in
patient coordinates, which is where the dose gets delivered.

Every geometry assertion uses anisotropic spacing or a rotated direction, because square
voxels and an axis-aligned volume pass whether or not the code is correct.
"""

import numpy as np
import pytest
import SimpleITK as sitk

from synthetic import make_ct_series
from pydicomrt.ct import CTBuilder, check_ct_iod
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder

# Anisotropic on every axis so nothing can be confused with anything else.
SPACING = (1.5, 3.0, 5.0)
SIZE = (10, 8, 6)  # sitk (x, y, z)
ORIGIN = (10.0, -20.0, 30.0)


def make_volume(spacing=SPACING, origin=ORIGIN, direction=None):
    rng = np.random.default_rng(0)
    array = rng.integers(-1000, 2000, SIZE[::-1]).astype(np.float32)  # numpy (z, y, x)
    image = sitk.GetImageFromArray(array)
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    if direction is not None:
        image.SetDirection(direction)
    return image


@pytest.fixture
def reference_series():
    return make_ct_series(shape=(3, 8, 8))


@pytest.fixture
def builder(reference_series):
    return CTBuilder(reference_series)


class TestAxialGeometry:
    def test_round_trips_through_the_loader_exactly(self, builder):
        """The strongest single check: build a series, read it back, compare geometry."""
        source = make_volume()

        recovered = SimpleITKImageBuilder().from_image_series(
            builder.set_volume(source).set_plane("AXIAL").build()
        )

        assert recovered.GetSize() == source.GetSize()
        assert recovered.GetSpacing() == pytest.approx(source.GetSpacing())
        assert recovered.GetOrigin() == pytest.approx(source.GetOrigin())
        assert recovered.GetDirection() == pytest.approx(source.GetDirection())

    def test_hounsfield_values_survive(self, builder):
        source = make_volume()

        recovered = SimpleITKImageBuilder().from_image_series(
            builder.set_volume(source).set_plane("AXIAL").build()
        )

        assert np.array_equal(
            sitk.GetArrayFromImage(recovered), sitk.GetArrayFromImage(source)
        )

    def test_pixel_spacing_is_row_then_column(self, builder):
        """
        DICOM orders PixelSpacing [between rows, between columns], the reverse of the index
        order. Writing (x, y) straight through corrupted the geometry on a round trip.
        """
        slices = builder.set_volume(make_volume()).set_plane("AXIAL").build()

        # rows advance along y (3.0), columns along x (1.5)
        assert list(slices[0].PixelSpacing) == pytest.approx([3.0, 1.5])
        assert slices[0].Rows == SIZE[1]
        assert slices[0].Columns == SIZE[0]

    def test_slices_step_along_the_normal(self, builder):
        slices = builder.set_volume(make_volume()).set_plane("AXIAL").build()

        first = np.array(slices[0].ImagePositionPatient, dtype=float)
        second = np.array(slices[1].ImagePositionPatient, dtype=float)

        assert first == pytest.approx(ORIGIN)
        assert (second - first) == pytest.approx([0.0, 0.0, SPACING[2]])
        assert slices[0].SliceThickness == pytest.approx(SPACING[2])


class TestPlaneViews:
    @pytest.mark.parametrize(
        "plan_view, normal_axis, slice_count",
        [("AXIAL", 2, SIZE[2]), ("CORONAL", 1, SIZE[1]), ("SAGITTAL", 0, SIZE[0])],
    )
    def test_each_plane_stacks_along_the_right_axis(
        self, builder, plan_view, normal_axis, slice_count
    ):
        """
        Named for the plane the slices lie in, so the stacking axis is the one *not* named:
        axial stacks along z, coronal along y, sagittal along x. CORONAL and SAGITTAL used
        to be the other way round.
        """
        slices = builder.set_volume(make_volume()).set_plane(plan_view).build()

        assert len(slices) == slice_count

        step = (np.array(slices[1].ImagePositionPatient, dtype=float)
                - np.array(slices[0].ImagePositionPatient, dtype=float))
        assert int(np.argmax(np.abs(step))) == normal_axis
        assert np.abs(step).max() == pytest.approx(SPACING[normal_axis])

    @pytest.mark.parametrize(
        "plan_view, expected_rows, expected_columns",
        [("AXIAL", SIZE[1], SIZE[0]), ("CORONAL", SIZE[2], SIZE[0]), ("SAGITTAL", SIZE[2], SIZE[1])],
    )
    def test_in_plane_dimensions(self, builder, plan_view, expected_rows, expected_columns):
        slices = builder.set_volume(make_volume()).set_plane(plan_view).build()

        assert slices[0].Rows == expected_rows
        assert slices[0].Columns == expected_columns

    def test_pixel_data_matches_the_declared_dimensions(self, builder):
        """A Rows/Columns mismatch renders as diagonal tearing rather than an error."""
        for plan_view in ("AXIAL", "CORONAL", "SAGITTAL"):
            slice_ds = builder.set_volume(make_volume()).set_plane(plan_view).build()[0]

            assert slice_ds.pixel_array.shape == (slice_ds.Rows, slice_ds.Columns)

    def test_unknown_plane_is_rejected(self, builder):
        with pytest.raises(ValueError, match="Invalid plan view"):
            builder.set_plane("OBLIQUE")

    def test_building_without_a_volume_is_rejected(self, reference_series):
        with pytest.raises(ValueError, match="no volume set"):
            CTBuilder(reference_series).build()

    def test_plane_name_is_case_insensitive(self, builder):
        assert len(builder.set_volume(make_volume()).set_plane("axial").build()) == SIZE[2]


class TestRotatedVolumes:
    ROTATED = (0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)

    def test_orientation_uses_direction_columns(self, builder):
        """
        ImageOrientationPatient is the first two *columns* of the direction matrix, not the
        first six entries of its row-major flattening. The two agree only for a symmetric
        matrix, which is why an axis-aligned volume cannot catch this.
        """
        source = make_volume(direction=self.ROTATED)

        slice_ds = builder.set_volume(source).set_plane("AXIAL").build()[0]

        matrix = np.array(self.ROTATED).reshape(3, 3)
        expected = matrix[:, 0].tolist() + matrix[:, 1].tolist()
        assert list(slice_ds.ImageOrientationPatient) == pytest.approx(expected)
        assert list(slice_ds.ImageOrientationPatient) != pytest.approx(list(self.ROTATED[:6]))

    def test_rotated_volume_round_trips(self, builder):
        source = make_volume(direction=self.ROTATED)

        recovered = SimpleITKImageBuilder().from_image_series(
            builder.set_volume(source).set_plane("AXIAL").build()
        )

        assert recovered.GetDirection() == pytest.approx(source.GetDirection())
        assert recovered.GetOrigin() == pytest.approx(source.GetOrigin())
        assert np.array_equal(
            sitk.GetArrayFromImage(recovered), sitk.GetArrayFromImage(source)
        )


class TestPixelEncoding:
    def test_stored_values_are_declared_unsigned(self, builder):
        """
        The data is written as uint16. Declaring PixelRepresentation signed made every
        stored value above 32767 read back negative.
        """
        slice_ds = builder.set_volume(make_volume()).set_plane("AXIAL").build()[0]

        assert slice_ds.PixelRepresentation == 0
        assert slice_ds.BitsAllocated == 16

    def test_high_density_values_do_not_wrap(self, builder):
        """A stored value beyond the signed range: metal implants reach these numbers."""
        array = np.full(SIZE[::-1], 40000.0, dtype=np.float32)
        image = sitk.GetImageFromArray(array)
        image.SetSpacing(SPACING)

        recovered = SimpleITKImageBuilder().from_image_series(
            builder.set_volume(image).set_plane("AXIAL").build()
        )

        assert sitk.GetArrayViewFromImage(recovered).min() == pytest.approx(40000.0)

    def test_values_above_the_stored_range_are_clipped_not_wrapped(self, builder):
        """The upper clip used to be applied before the intercept offset, so it overflowed."""
        array = np.full(SIZE[::-1], 1e6, dtype=np.float32)
        image = sitk.GetImageFromArray(array)
        image.SetSpacing(SPACING)

        recovered = SimpleITKImageBuilder().from_image_series(
            builder.set_volume(image).set_plane("AXIAL").build()
        )
        values = sitk.GetArrayViewFromImage(recovered)

        assert values.max() == pytest.approx(2 ** 16 - 1 - 1024)
        assert values.min() > 0  # not wrapped to a negative or tiny value

    def test_padding_value_is_an_integer(self, builder):
        """It used to be written as raw bytes; the VR here is US."""
        slice_ds = builder.set_volume(make_volume()).set_plane("AXIAL").build()[0]

        assert isinstance(slice_ds.PixelPaddingValue, int)
        # air, in the stored representation
        assert slice_ds.PixelPaddingValue + slice_ds.RescaleIntercept == -1024


class TestDatasetIdentity:
    def test_passes_the_iod_check(self, builder):
        for slice_ds in builder.set_volume(make_volume()).set_plane("AXIAL").build():
            result = check_ct_iod(slice_ds)
            assert result["result"] is True, result["content"]

    def test_slices_share_a_series_and_differ_by_instance(self, builder):
        slices = builder.set_volume(make_volume()).set_plane("AXIAL").build()

        assert len({s.SeriesInstanceUID for s in slices}) == 1
        assert len({s.SOPInstanceUID for s in slices}) == len(slices)
        assert [s.InstanceNumber for s in slices] == list(range(1, len(slices) + 1))

    def test_sop_instance_uid_matches_file_meta(self, builder):
        for slice_ds in builder.set_volume(make_volume()).set_plane("AXIAL").build():
            assert slice_ds.SOPInstanceUID == slice_ds.file_meta.MediaStorageSOPInstanceUID

    def test_uid_prefix_is_honoured(self, builder):
        prefix = "1.2.826.0.1.3680043.2.1125."
        builder.set_uid_prefix(prefix)

        slices = builder.set_volume(make_volume()).set_plane("AXIAL").build()

        assert slices[0].SeriesInstanceUID.startswith(prefix)
        assert slices[0].SOPInstanceUID.startswith(prefix)
        assert slices[0].file_meta.ImplementationClassUID.startswith(prefix)

    def test_marked_as_derived(self, builder):
        slice_ds = builder.set_volume(make_volume()).set_plane("AXIAL").build()[0]

        assert list(slice_ds.ImageType)[:2] == ["DERIVED", "SECONDARY"]

    def test_inherits_patient_and_frame_of_reference(self, builder, reference_series):
        slice_ds = builder.set_volume(make_volume()).set_plane("AXIAL").build()[0]

        assert slice_ds.PatientID == reference_series[0].PatientID
        assert slice_ds.StudyInstanceUID == reference_series[0].StudyInstanceUID
        assert slice_ds.FrameOfReferenceUID == reference_series[0].FrameOfReferenceUID

    def test_empty_reference_series_is_rejected(self):
        with pytest.raises(ValueError, match="image_series"):
            CTBuilder([])


class TestWriteAndRead:
    def test_series_reads_back_through_gdcm(self, builder, tmp_path):
        """
        SimpleITK's own reader is an independent implementation of the same standard, so it
        catches a convention we have misread rather than merely regressed against.
        """
        source = make_volume()
        for index, slice_ds in enumerate(builder.set_volume(source).set_plane("AXIAL").build()):
            slice_ds.save_as(tmp_path / f"slice_{index:03d}.dcm", enforce_file_format=True)

        reader = sitk.ImageSeriesReader()
        reader.SetFileNames(reader.GetGDCMSeriesFileNames(str(tmp_path)))
        recovered = reader.Execute()

        assert recovered.GetSize() == source.GetSize()
        assert recovered.GetSpacing() == pytest.approx(source.GetSpacing())
        assert recovered.GetOrigin() == pytest.approx(source.GetOrigin())
        assert np.array_equal(
            sitk.GetArrayFromImage(sitk.Cast(recovered, sitk.sitkFloat32)),
            sitk.GetArrayFromImage(source),
        )


class TestBuildFromNumpy:
    def test_matches_the_sitk_path(self, builder):
        source = make_volume()

        from_array = (
            builder.set_volume_from_array(
                volume=sitk.GetArrayFromImage(source),
                origin=source.GetOrigin(),
                spacing=source.GetSpacing(),
                direction=source.GetDirection(),
            )
            .build()
        )

        recovered = SimpleITKImageBuilder().from_image_series(from_array)
        assert recovered.GetSpacing() == pytest.approx(source.GetSpacing())
        assert recovered.GetOrigin() == pytest.approx(source.GetOrigin())
