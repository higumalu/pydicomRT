"""
Deformable REG builder tests.

Grew out of a print-only demo script; the assertions now pin down the IOD fields that a
deformable REG consumer actually reads back.
"""

import numpy as np
import pytest
import SimpleITK as sitk
from pydicom.dataset import Dataset
from pydicom.uid import generate_uid

from pydicomrt.reg import check_deformable_reg_iod
from pydicomrt.reg.builder import DeformableSpatialRegistrationBuilder

DEFORMABLE_SPATIAL_REGISTRATION_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.66.3"
CT_IMAGE_STORAGE_SOP_CLASS_UID = "1.2.840.10008.5.1.4.1.1.2"


# -----------------------------------------------------------
# Build fake data
# -----------------------------------------------------------
def make_fake_fixed_ref_ds() -> Dataset:
    ds = Dataset()
    ds.PatientName = "Test^Patient"
    ds.PatientID = "TP001"
    ds.PatientBirthDate = "19700101"
    ds.PatientSex = "O"
    ds.StudyInstanceUID = generate_uid()
    ds.StudyID = "1"
    ds.StudyDate = "20250101"
    ds.StudyTime = "120000"
    ds.AccessionNumber = ""
    ds.ReferringPhysicianName = ""
    ds.FrameOfReferenceUID = generate_uid()
    ds.PositionReferenceIndicator = "SCI"
    ds.SeriesInstanceUID = generate_uid()
    ds.SOPClassUID = CT_IMAGE_STORAGE_SOP_CLASS_UID
    ds.SOPInstanceUID = generate_uid()
    return ds


def make_fake_moving_instance(frame_of_ref_uid: str, series_uid: str = None) -> Dataset:
    ds = Dataset()
    ds.SOPClassUID = CT_IMAGE_STORAGE_SOP_CLASS_UID
    ds.SOPInstanceUID = generate_uid()
    ds.FrameOfReferenceUID = frame_of_ref_uid
    ds.SeriesInstanceUID = series_uid or generate_uid()
    return ds


def make_displacement_field_transform(
    size=(8, 9, 10),
    spacing=(1.5, 1.5, 2.0),
    origin=(0.0, 0.0, 0.0),
    direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    shift=(1.0, 0.0, 0.0)
) -> sitk.DisplacementFieldTransform:
    img_x = sitk.Image(size, sitk.sitkFloat64)
    img_x += float(shift[0])
    img_y = sitk.Image(size, sitk.sitkFloat64)
    img_y += float(shift[1])
    img_z = sitk.Image(size, sitk.sitkFloat64)
    img_z += float(shift[2])

    vec = sitk.Compose(img_x, img_y, img_z)
    vec = sitk.Cast(vec, sitk.sitkVectorFloat64)

    vec.SetSpacing(spacing)
    vec.SetOrigin(origin)
    vec.SetDirection(direction)

    return sitk.DisplacementFieldTransform(vec)


def identity_4x4_flat_f32():
    return [float(x) for x in np.eye(4, dtype=np.float32).reshape(-1)]


@pytest.fixture
def deformable_reg_ds():
    fixed_ref = make_fake_fixed_ref_ds()
    moving_series_uid = generate_uid()
    moving_list = [
        make_fake_moving_instance(fixed_ref.FrameOfReferenceUID, moving_series_uid),
        make_fake_moving_instance(fixed_ref.FrameOfReferenceUID, moving_series_uid),
    ]
    displacement_transform = make_displacement_field_transform(shift=(1.0, 0.0, 0.0))

    builder = DeformableSpatialRegistrationBuilder([fixed_ref])
    builder.set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    builder.add_registration(
        moving_series=moving_list,
        vectorial_field_transform=displacement_transform,
        pre_transform=identity_4x4_flat_f32(),
        post_transform=identity_4x4_flat_f32(),
    )
    return builder.build(), fixed_ref, moving_list


def test_builds_deformable_registration_iod(deformable_reg_ds):
    reg_ds, fixed_ref, _ = deformable_reg_ds

    assert reg_ds.SOPClassUID == DEFORMABLE_SPATIAL_REGISTRATION_SOP_CLASS_UID
    assert reg_ds.SOPClassUID == reg_ds.file_meta.MediaStorageSOPClassUID
    assert reg_ds.SOPInstanceUID == reg_ds.file_meta.MediaStorageSOPInstanceUID
    assert reg_ds.Modality == "REG"
    # Patient and study context must be inherited so the REG files alongside its images
    assert reg_ds.PatientID == fixed_ref.PatientID
    assert reg_ds.StudyInstanceUID == fixed_ref.StudyInstanceUID
    assert reg_ds.FrameOfReferenceUID == fixed_ref.FrameOfReferenceUID


def test_generated_uids_honour_the_prefix(deformable_reg_ds):
    reg_ds, _, _ = deformable_reg_ds
    prefix = "1.2.826.0.1.3680043.2.1125."

    assert reg_ds.SOPInstanceUID.startswith(prefix)
    assert reg_ds.SeriesInstanceUID.startswith(prefix)


def test_references_every_moving_instance(deformable_reg_ds):
    reg_ds, _, moving_list = deformable_reg_ds
    item = reg_ds.DeformableRegistrationSequence[0]

    assert item.SourceFrameOfReferenceUID == moving_list[0].FrameOfReferenceUID
    referenced_uids = [ref.ReferencedSOPInstanceUID for ref in item.ReferencedImageSequence]
    assert referenced_uids == [ds.SOPInstanceUID for ds in moving_list]


def test_grid_geometry_matches_the_displacement_field(deformable_reg_ds):
    reg_ds, _, _ = deformable_reg_ds
    grid = reg_ds.DeformableRegistrationSequence[0].DeformableRegistrationGridSequence[0]

    assert list(grid.GridDimensions) == [8, 9, 10]
    assert list(grid.GridResolution) == pytest.approx([1.5, 1.5, 2.0])
    assert list(grid.ImagePositionPatient) == pytest.approx([0.0, 0.0, 0.0])
    # ImageOrientationPatient is VM 6: row direction followed by column direction
    assert len(grid.ImageOrientationPatient) == 6
    assert list(grid.ImageOrientationPatient) == pytest.approx([1.0, 0.0, 0.0, 0.0, 1.0, 0.0])


def test_vector_grid_data_is_float32_xyz_triplets(deformable_reg_ds):
    reg_ds, _, _ = deformable_reg_ds
    grid = reg_ds.DeformableRegistrationSequence[0].DeformableRegistrationGridSequence[0]

    voxel_count = int(np.prod(grid.GridDimensions))
    assert len(grid.VectorGridData) == voxel_count * 3 * 4  # 3 components, float32

    vectors = np.frombuffer(grid.VectorGridData, dtype="<f4").reshape(-1, 3)
    # The field was built as a uniform +1 mm shift along x
    assert np.allclose(vectors[:, 0], 1.0)
    assert np.allclose(vectors[:, 1:], 0.0)


def test_orientation_uses_direction_columns_not_flattened_prefix():
    """
    A rotated grid separates the two: for a non-symmetric direction matrix the first six
    entries of the row-major flattening are not the row and column direction vectors.
    """
    fixed_ref = make_fake_fixed_ref_ds()
    rotated_direction = (0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    displacement_transform = make_displacement_field_transform(direction=rotated_direction)

    builder = DeformableSpatialRegistrationBuilder([fixed_ref])
    builder.add_registration(
        moving_series=[make_fake_moving_instance(fixed_ref.FrameOfReferenceUID)],
        vectorial_field_transform=displacement_transform,
        pre_transform=identity_4x4_flat_f32(),
        post_transform=identity_4x4_flat_f32(),
    )
    reg_ds = builder.build()
    grid = reg_ds.DeformableRegistrationSequence[0].DeformableRegistrationGridSequence[0]

    matrix = np.array(rotated_direction).reshape(3, 3)
    expected = matrix[:, 0].tolist() + matrix[:, 1].tolist()
    assert list(grid.ImageOrientationPatient) == pytest.approx(expected)
    assert list(grid.ImageOrientationPatient) != pytest.approx(list(rotated_direction[:6]))


def test_pre_and_post_deformation_matrices_are_recorded(deformable_reg_ds):
    reg_ds, _, _ = deformable_reg_ds
    item = reg_ds.DeformableRegistrationSequence[0]

    pre = item.PreDeformationMatrixRegistrationSequence[0]
    post = item.PostDeformationMatrixRegistrationSequence[0]
    assert pre.FrameOfReferenceTransformationMatrixType == "RIGID"
    assert post.FrameOfReferenceTransformationMatrixType == "RIGID"
    assert list(pre.FrameOfReferenceTransformationMatrix) == pytest.approx(identity_4x4_flat_f32())
    assert list(post.FrameOfReferenceTransformationMatrix) == pytest.approx(identity_4x4_flat_f32())


class TestIodConformance:
    def test_passes_the_iod_check(self, deformable_reg_ds):
        reg_ds, _, _ = deformable_reg_ds

        result = check_deformable_reg_iod(reg_ds)

        assert result["result"] is True, result["content"]

    def test_content_identification_is_present(self, deformable_reg_ds):
        """Instance Number, Content Date and Content Time are Type 1 and were all missing."""
        reg_ds, _, _ = deformable_reg_ds

        assert reg_ds.InstanceNumber == 1
        assert reg_ds.ContentDate
        assert reg_ds.ContentTime
        assert len(reg_ds.ContentLabel) <= 16  # VR CS

    def test_registration_type_code_sits_at_the_item_level(self, deformable_reg_ds):
        """
        The deformable IOD has no Matrix Registration Sequence, so unlike the rigid one the
        code sequence belongs directly in the Deformable Registration Sequence item.
        """
        reg_ds, _, _ = deformable_reg_ds
        item = reg_ds.DeformableRegistrationSequence[0]

        assert item.RegistrationTypeCodeSequence[0].CodeValue == "125024"

    def test_lists_the_referenced_series(self, deformable_reg_ds):
        reg_ds, fixed_ref, moving_list = deformable_reg_ds

        series_uids = {item.SeriesInstanceUID for item in reg_ds.ReferencedSeriesSequence}
        assert series_uids == {fixed_ref.SeriesInstanceUID, moving_list[0].SeriesInstanceUID}

    def test_survives_write_and_read(self, deformable_reg_ds, tmp_path):
        from pydicom import dcmread

        reg_ds, _, _ = deformable_reg_ds
        path = tmp_path / "deformable.dcm"
        reg_ds.save_as(path, enforce_file_format=True)

        reloaded = dcmread(path)

        assert check_deformable_reg_iod(reloaded)["result"] is True
        grid = (reloaded.DeformableRegistrationSequence[0]
                .DeformableRegistrationGridSequence[0])
        assert list(grid.GridDimensions) == [8, 9, 10]


class TestBuilderValidation:
    def test_building_without_a_registration_is_rejected(self):
        builder = DeformableSpatialRegistrationBuilder([make_fake_fixed_ref_ds()])

        with pytest.raises(ValueError, match="no registration added"):
            builder.build()

    def test_empty_moving_series_is_rejected(self):
        builder = DeformableSpatialRegistrationBuilder([make_fake_fixed_ref_ds()])

        with pytest.raises(ValueError, match="moving_series"):
            builder.add_registration(
                moving_series=[],
                vectorial_field_transform=make_displacement_field_transform(),
                pre_transform=identity_4x4_flat_f32(),
                post_transform=identity_4x4_flat_f32(),
            )

    def test_malformed_pre_transform_is_rejected(self):
        fixed_ref = make_fake_fixed_ref_ds()
        builder = DeformableSpatialRegistrationBuilder([fixed_ref])

        with pytest.raises(ValueError, match="16 elements"):
            builder.add_registration(
                moving_series=[make_fake_moving_instance(fixed_ref.FrameOfReferenceUID)],
                vectorial_field_transform=make_displacement_field_transform(),
                pre_transform=[1, 2, 3],
                post_transform=identity_4x4_flat_f32(),
            )

    def test_affine_pre_transform_must_be_declared(self):
        """A scaled matrix silently declared RIGID is exactly the kind of thing that loads."""
        fixed_ref = make_fake_fixed_ref_ds()
        scaled = np.diag([1.5, 1.0, 1.0, 1.0]).ravel().tolist()

        builder = DeformableSpatialRegistrationBuilder([fixed_ref])
        with pytest.raises(ValueError, match="orthonormal"):
            builder.add_registration(
                moving_series=[make_fake_moving_instance(fixed_ref.FrameOfReferenceUID)],
                vectorial_field_transform=make_displacement_field_transform(),
                pre_transform=scaled,
                post_transform=identity_4x4_flat_f32(),
            )

        builder.add_registration(
            moving_series=[make_fake_moving_instance(fixed_ref.FrameOfReferenceUID)],
            vectorial_field_transform=make_displacement_field_transform(),
            pre_transform=scaled,
            post_transform=identity_4x4_flat_f32(),
            pre_transform_type="RIGID_SCALE",
        )
        item = builder.build().DeformableRegistrationSequence[0]
        assert (item.PreDeformationMatrixRegistrationSequence[0]
                .FrameOfReferenceTransformationMatrixType) == "RIGID_SCALE"
