"""
CTBuilder against a real clinical CT series.

The synthetic CT tests use volumes this library also built, so they agree with its own
conventions by construction. This starts from DICOM a scanner wrote: load it, rebuild it,
and check that the geometry and Hounsfield values come back.
"""

import numpy as np
import pytest
import SimpleITK as sitk

from pydicomrt.ct import CTBuilder, check_ct_iod
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder


@pytest.fixture(scope="module")
def rebuilt_series(planning_ct_series, planning_ct_image):
    """The real CT, taken through SimpleITK and written back out as a new series."""
    return CTBuilder(planning_ct_series).set_volume(planning_ct_image).set_plane("AXIAL").build()


class TestRebuildingARealSeries:
    def test_slice_count_is_preserved(self, rebuilt_series, planning_ct_series):
        assert len(rebuilt_series) == len(planning_ct_series)

    def test_geometry_survives_the_rebuild(self, rebuilt_series, planning_ct_image):
        recovered = SimpleITKImageBuilder().from_image_series(rebuilt_series)

        assert recovered.GetSize() == planning_ct_image.GetSize()
        assert recovered.GetSpacing() == pytest.approx(planning_ct_image.GetSpacing())
        assert recovered.GetOrigin() == pytest.approx(planning_ct_image.GetOrigin())
        assert recovered.GetDirection() == pytest.approx(planning_ct_image.GetDirection())

    def test_hounsfield_values_survive_the_rebuild(self, rebuilt_series, planning_ct_image):
        """
        Clinical CT reaches about -1024 to 3071 HU, comfortably inside the stored range, so
        this should be lossless rather than merely close.
        """
        recovered = SimpleITKImageBuilder().from_image_series(rebuilt_series)

        assert np.array_equal(
            sitk.GetArrayFromImage(recovered),
            sitk.GetArrayFromImage(planning_ct_image),
        )

    def test_slice_positions_match_the_original(self, rebuilt_series, planning_ct_series):
        original = np.array(
            [ds.ImagePositionPatient for ds in planning_ct_series], dtype=float
        )
        rebuilt = np.array(
            [ds.ImagePositionPatient for ds in rebuilt_series], dtype=float
        )

        assert rebuilt == pytest.approx(original, abs=1e-6)

    def test_pixel_spacing_matches_the_original(self, rebuilt_series, planning_ct_series):
        assert list(rebuilt_series[0].PixelSpacing) == pytest.approx(
            [float(v) for v in planning_ct_series[0].PixelSpacing]
        )

    def test_orientation_matches_the_original(self, rebuilt_series, planning_ct_series):
        assert list(rebuilt_series[0].ImageOrientationPatient) == pytest.approx(
            [float(v) for v in planning_ct_series[0].ImageOrientationPatient]
        )

    def test_every_slice_passes_the_iod_check(self, rebuilt_series):
        for slice_ds in rebuilt_series:
            result = check_ct_iod(slice_ds)
            assert result["result"] is True, result["content"]

    def test_patient_context_is_carried_over(self, rebuilt_series, planning_ct_series):
        source = planning_ct_series[0]
        rebuilt = rebuilt_series[0]

        assert rebuilt.PatientID == source.PatientID
        assert rebuilt.StudyInstanceUID == source.StudyInstanceUID
        assert rebuilt.FrameOfReferenceUID == source.FrameOfReferenceUID

    def test_rebuild_is_a_new_series(self, rebuilt_series, planning_ct_series):
        """Derived images must not masquerade as instances of the acquired series."""
        source_sop_uids = {ds.SOPInstanceUID for ds in planning_ct_series}

        assert rebuilt_series[0].SeriesInstanceUID != planning_ct_series[0].SeriesInstanceUID
        assert not source_sop_uids & {ds.SOPInstanceUID for ds in rebuilt_series}
        assert list(rebuilt_series[0].ImageType)[:2] == ["DERIVED", "SECONDARY"]


class TestReadableByGdcm:
    def test_written_series_reads_back_through_gdcm(
        self, rebuilt_series, planning_ct_image, tmp_path_factory
    ):
        """
        SimpleITK's reader is an independent implementation, so this catches a convention
        misread rather than a regression against our own reader.
        """
        directory = tmp_path_factory.mktemp("rebuilt_ct")
        for index, slice_ds in enumerate(rebuilt_series):
            slice_ds.save_as(directory / f"slice_{index:04d}.dcm", enforce_file_format=True)

        reader = sitk.ImageSeriesReader()
        reader.SetFileNames(reader.GetGDCMSeriesFileNames(str(directory)))
        recovered = reader.Execute()

        assert recovered.GetSize() == planning_ct_image.GetSize()
        assert recovered.GetSpacing() == pytest.approx(planning_ct_image.GetSpacing())
        assert recovered.GetOrigin() == pytest.approx(planning_ct_image.GetOrigin())
        assert np.array_equal(
            sitk.GetArrayFromImage(sitk.Cast(recovered, sitk.sitkFloat64)),
            sitk.GetArrayFromImage(sitk.Cast(planning_ct_image, sitk.sitkFloat64)),
        )
