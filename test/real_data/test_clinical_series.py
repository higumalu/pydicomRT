"""
Series loading and REG parsing against real clinical DICOM.

These are the checks synthetic fixtures cannot make: the fixtures are built by the same
understanding of the standard that the library uses, so they agree with it by construction.
Real files written by a scanner and a TPS do not.
"""

import numpy as np
import pytest
import SimpleITK as sitk

from pydicomrt.reg import get_spatial_registrations
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder


class TestAgainstSimpleITKReader:
    """
    Cross-check the DICOM -> SimpleITK conversion against SimpleITK's own GDCM reader.

    This is the strongest available check on the geometry conventions: an independent
    implementation of the same standard, reading the same files.
    """

    @staticmethod
    def _gdcm_reference(directory):
        reader = sitk.ImageSeriesReader()
        reader.SetFileNames(reader.GetGDCMSeriesFileNames(str(directory)))
        return reader.Execute()

    @pytest.mark.parametrize(
        "series_name, fixture_name",
        [("CTSIM", "planning_ct_series"), ("CBCT", "cbct_series")],
    )
    def test_geometry_matches_gdcm(self, series_name, fixture_name, series_dir, request):
        directory = series_dir(series_name)
        ours = SimpleITKImageBuilder().from_image_series(request.getfixturevalue(fixture_name))
        reference = self._gdcm_reference(directory)

        assert ours.GetSize() == reference.GetSize()
        assert ours.GetOrigin() == pytest.approx(reference.GetOrigin())
        assert ours.GetSpacing() == pytest.approx(reference.GetSpacing())
        assert ours.GetDirection() == pytest.approx(reference.GetDirection())

    @pytest.mark.parametrize(
        "series_name, fixture_name",
        [("CTSIM", "planning_ct_series"), ("CBCT", "cbct_series")],
    )
    def test_voxel_values_match_gdcm(self, series_name, fixture_name, series_dir, request):
        directory = series_dir(series_name)
        ours = SimpleITKImageBuilder().from_image_series(request.getfixturevalue(fixture_name))
        reference = self._gdcm_reference(directory)

        # GetArrayFromImage, not GetArrayView: a view over a temporary Cast() result
        # outlives the image it points into.
        ours_array = sitk.GetArrayFromImage(ours).astype(np.float64)
        reference_array = sitk.GetArrayFromImage(sitk.Cast(reference, sitk.sitkFloat64))

        # Covers slice ordering and the rescale slope/intercept together
        assert np.array_equal(ours_array, reference_array)


class TestSeriesGeometry:
    def test_planning_ct_and_cbct_are_distinct_frames_of_reference(
        self, planning_ct_series, cbct_series
    ):
        assert (
            planning_ct_series[0].FrameOfReferenceUID
            != cbct_series[0].FrameOfReferenceUID
        )

    @pytest.mark.parametrize("fixture_name", ["planning_ct_series", "cbct_series"])
    def test_slices_are_ordered_and_evenly_spaced(self, fixture_name, request):
        ds_list = request.getfixturevalue(fixture_name)

        z_positions = np.array([float(ds.ImagePositionPatient[2]) for ds in ds_list])
        steps = np.diff(z_positions)

        assert np.all(steps > 0), "slices are not sorted along +z"
        assert steps.std() < 1e-6, "slice spacing is not uniform"

    @pytest.mark.parametrize("fixture_name", ["planning_ct_image", "cbct_image"])
    def test_reaches_air_so_padding_inference_works(self, fixture_name, request):
        """The -1000 padding inference keys off the intensity minimum reaching air."""
        from pydicomrt.reg.pipeline import infer_default_pixel_value

        image = request.getfixturevalue(fixture_name)

        assert infer_default_pixel_value(image) == -1000.0


class TestClinicalRegParsing:
    def test_parser_finds_the_registration(
        self, clinical_reg_ds, planning_ct_series, cbct_series
    ):
        reg_dict = get_spatial_registrations(clinical_reg_ds)

        fixed_uid = planning_ct_series[0].FrameOfReferenceUID
        moving_uid = cbct_series[0].FrameOfReferenceUID

        assert fixed_uid in reg_dict, "planning CT frame is not the registered RCS"
        assert moving_uid in reg_dict[fixed_uid]

    def test_parsed_matrix_is_a_valid_rigid_transform(
        self, clinical_reg_ds, planning_ct_series, cbct_series
    ):
        reg_dict = get_spatial_registrations(clinical_reg_ds)
        fixed_uid = planning_ct_series[0].FrameOfReferenceUID
        moving_uid = cbct_series[0].FrameOfReferenceUID

        matrix = np.array(reg_dict[fixed_uid][moving_uid], dtype=float).reshape(4, 4)

        rotation = matrix[:3, :3]
        assert rotation @ rotation.T == pytest.approx(np.identity(3), abs=1e-6)
        assert np.linalg.det(rotation) == pytest.approx(1.0, abs=1e-6)
        assert matrix[3] == pytest.approx([0.0, 0.0, 0.0, 1.0])

    def test_registration_is_declared_rigid(self, clinical_reg_ds):
        for item in clinical_reg_ds.RegistrationSequence:
            for matrix_registration in item.MatrixRegistrationSequence:
                for matrix in matrix_registration.MatrixSequence:
                    assert matrix.FrameOfReferenceTransformationMatrixType == "RIGID"


class TestClinicalTransformDirection:
    """
    Pin down which way the stored matrix points.

    A REG matrix maps the item's frame into the registered RCS (moving -> fixed), while a
    SimpleITK resampling transform maps fixed -> moving. Both directions "work" numerically,
    so the only way to tell them apart is to check which one actually aligns anatomy.
    """

    @staticmethod
    def _body_overlap(fixed_image, moving_image, transform, shrink=4):
        # Register on a coarse grid: this is a direction check, not an accuracy check.
        coarse = sitk.Shrink(fixed_image, [shrink] * 3)
        moved = sitk.Resample(moving_image, coarse, transform, sitk.sitkLinear, -1000.0)

        fixed_body = sitk.GetArrayFromImage(coarse) > -300
        moved_body = sitk.GetArrayFromImage(moved) > -300
        total = fixed_body.sum() + moved_body.sum()
        return 2.0 * np.logical_and(fixed_body, moved_body).sum() / total

    def test_clinical_transform_aligns_anatomy(
        self, planning_ct_image, cbct_image, clinical_transform
    ):
        aligned = self._body_overlap(planning_ct_image, cbct_image, clinical_transform)
        unregistered = self._body_overlap(planning_ct_image, cbct_image, sitk.Transform())

        # The CBCT covers only part of the CT, which caps the achievable overlap well
        # below 1.0 -- the point is that registering more than triples it.
        assert aligned > 0.7
        assert aligned > 3 * unregistered

    def test_the_inverse_direction_is_clearly_wrong(
        self, planning_ct_image, cbct_image, clinical_transform
    ):
        """Guards the guard: confirm the overlap measure can tell the directions apart."""
        backwards = clinical_transform.GetInverse()

        assert self._body_overlap(planning_ct_image, cbct_image, backwards) < 0.2
