"""Physical mappings at the registration/DICOM boundary, using synthetic data only."""

from itertools import product

import numpy as np
import pytest
import SimpleITK as sitk
from pydicom import dcmread
from pydicom.dataset import Dataset

from synthetic import make_ct_series
from pydicomrt.reg import (
    SpatialRegistrationBuilder, DeformableSpatialRegistrationBuilder,
    affine_to_homogeneous_matrix, get_spatial_registrations,
    get_deformable_registrations, check_spatial_reg_iod, check_deformable_reg_iod,
)
from pydicomrt.reg.method.demons import multiscale_demons, demons_registration
from pydicomrt.reg.pipeline.preprocessing import (
    crop_image_to_extent, get_image_physical_extent, align_image_extents,
    create_reference_image_from_extent,
)


def field_image(size=(16, 16, 16)):
    """Nonuniform field on a rotated, anisotropic grid."""
    z, y, x = np.indices(size[::-1])
    values = np.stack((0.03 * y, 0.02 * x, 0.01 * z), axis=-1)
    image = sitk.GetImageFromArray(values, isVector=True)
    image.SetSpacing((1.0, 1.5, 2.0))
    image.SetOrigin((-10.0, 3.0, -5.0))
    angle = 0.3
    image.SetDirection((np.cos(angle), -np.sin(angle), 0, np.sin(angle), np.cos(angle), 0, 0, 0, 1))
    return image


def matrix_point(matrix, point):
    return (np.asarray(matrix).reshape(4, 4) @ np.r_[point, 1])[:3]


def test_centred_rigid_roundtrip_preserves_point_mapping(tmp_path):
    fixed, moving = make_ct_series(), make_ct_series()
    transform = sitk.Euler3DTransform()
    transform.SetCenter((12.0, -7.0, 4.0))
    transform.SetRotation(0.1, -0.2, 0.3)
    transform.SetTranslation((3.0, 5.0, -2.0))
    matrix = affine_to_homogeneous_matrix(transform.GetInverse())
    ds = SpatialRegistrationBuilder(fixed).add_registration(moving, matrix).build()
    path = tmp_path / 'rigid.dcm'
    ds.save_as(path, enforce_file_format=True)
    stored = get_spatial_registrations(dcmread(path))[fixed[0].FrameOfReferenceUID][moving[0].FrameOfReferenceUID]
    for point in [(12., -7., 4.), (0., 0., 0.), (-9., 13., 7.)]:
        assert matrix_point(stored, transform.TransformPoint(point)) == pytest.approx(point, abs=1e-6)


def test_deformable_roundtrip_matches_rigid_after_residual(tmp_path):
    fixed, moving = make_ct_series(), make_ct_series()
    field = field_image()
    residual = sitk.DisplacementFieldTransform(sitk.Image(field))
    rigid = sitk.Euler3DTransform()
    rigid.SetCenter((5., -2., 3.))
    rigid.SetRotation(0.1, -0.1, 0.2)
    rigid.SetTranslation((10., -4., 2.))
    ds = (DeformableSpatialRegistrationBuilder(fixed)
          .add_registration(moving, residual, post_transform=affine_to_homogeneous_matrix(rigid))
          .build())
    path = tmp_path / 'deformable.dcm'
    ds.save_as(path, enforce_file_format=True)
    reloaded = dcmread(path)
    assert check_deformable_reg_iod(reloaded)['result']
    entry = get_deformable_registrations(reloaded)[0]
    assert np.array(entry['PreDeformationMatrixRegistration']).reshape(4, 4) == pytest.approx(np.eye(4))
    grid = entry['DeformableRegistrationGrid']
    # Independently decode DICOM geometry and evaluate its equation at grid nodes.
    iop = np.array(grid['ImageOrientationPatient']).reshape(2, 3)
    axes = np.column_stack((iop[0], iop[1], np.cross(*iop)))
    for index in [(3, 4, 5), (8, 7, 6), (12, 10, 9)]:
        point = np.asarray(grid['ImagePositionPatient']) + axes @ (np.asarray(grid['GridResolution']) * index)
        x, y, z = index
        displaced = matrix_point(entry['PreDeformationMatrixRegistration'], point) + grid['VectorGridData'][z, y, x]
        actual = matrix_point(entry['PostDeformationMatrixRegistration'], displaced)
        expected = rigid.TransformPoint(residual.TransformPoint(tuple(point)))
        assert actual == pytest.approx(expected, abs=1e-6)


def test_spatial_parser_composes_all_matrices_in_dicom_order():
    fixed, moving = make_ct_series(), make_ct_series()
    first = np.eye(4)
    first[:3, 3] = (10, 0, 0)
    second = np.eye(4)
    second[:2, :2] = [[0, -1], [1, 0]]
    ds = SpatialRegistrationBuilder(fixed).add_registration(moving, first).build()
    item = Dataset()
    item.FrameOfReferenceTransformationMatrix = second.ravel().tolist()
    item.FrameOfReferenceTransformationMatrixType = 'RIGID'
    ds.RegistrationSequence[0].MatrixRegistrationSequence[0].MatrixSequence.append(item)
    parsed = get_spatial_registrations(ds, calc_inverse_matrix=True)
    forward = parsed[fixed[0].FrameOfReferenceUID][moving[0].FrameOfReferenceUID]
    reverse = parsed[moving[0].FrameOfReferenceUID][fixed[0].FrameOfReferenceUID]
    assert matrix_point(forward, (1, 2, 3)) == pytest.approx((-2, 11, 3))
    assert matrix_point(reverse, (-2, 11, 3)) == pytest.approx((1, 2, 3))


def deformable_dataset():
    return (DeformableSpatialRegistrationBuilder(make_ct_series())
            .add_registration(make_ct_series(), sitk.DisplacementFieldTransform(field_image()))
            .build())


def test_grid_only_deformable_is_valid_and_parses_identity_matrices():
    ds = deformable_dataset()
    assert check_deformable_reg_iod(ds)['result']
    parsed = get_deformable_registrations(ds)
    assert len(parsed) == 1
    for name in ['PreDeformationMatrixRegistration', 'PostDeformationMatrixRegistration']:
        assert np.asarray(parsed[0][name]).reshape(4, 4) == pytest.approx(np.eye(4))


@pytest.mark.parametrize('kind', ['rigid', 'deformable'])
def test_empty_required_registration_sequence_is_rejected(kind):
    if kind == 'rigid':
        ds = SpatialRegistrationBuilder(make_ct_series()).add_registration(make_ct_series(), np.eye(4)).build()
        ds.RegistrationSequence = []
        result = check_spatial_reg_iod(ds)
    else:
        ds = deformable_dataset()
        ds.DeformableRegistrationSequence = []
        result = check_deformable_reg_iod(ds)
    assert not result['result']


def test_bad_grid_byte_count_is_rejected():
    ds = deformable_dataset()
    ds.DeformableRegistrationSequence[0].DeformableRegistrationGridSequence[0].VectorGridData = b'1234'
    assert not check_deformable_reg_iod(ds)['result']
    assert get_deformable_registrations(ds) == []


def test_builder_rejects_missing_fixed_frame():
    series = make_ct_series()
    for ds in series:
        del ds.FrameOfReferenceUID
    with pytest.raises(ValueError, match='FrameOfReferenceUID'):
        SpatialRegistrationBuilder(series)


class ConstantRefinement:
    """Known refinement isolates composition from optimizer behaviour."""
    def SetNumberOfIterations(self, count):
        pass

    def GetStandardDeviations(self):
        return (0.1, 0.1, 0.1)

    def Execute(self, fixed, moving):
        values = np.zeros(tuple(reversed(fixed.GetSize())) + (3,))
        values[..., 1] = 1
        result = sitk.GetImageFromArray(values, isVector=True)
        result.CopyInformation(fixed)
        return result


def test_demons_refines_existing_mapping_in_correct_order():
    image = sitk.Image([16] * 3, sitk.sitkFloat32)
    image.SetSpacing((1., 2., 3.))
    values = np.zeros((16, 16, 16, 3))
    values[..., 0] = np.indices((16, 16, 16))[1] * 0.2  # u_x = 0.1 * physical y
    initial = sitk.GetImageFromArray(values, isVector=True)
    initial.CopyInformation(image)
    result = multiscale_demons(
        ConstantRefinement(), image, image, initial_displacement_field=initial,
        resolution_staging=(1,), smoothing_sigmas=(0,), iteration_staging=(1,),
    )
    transform = sitk.DisplacementFieldTransform(result)
    # D shifts y by 1, then T shifts x by 0.1 * the updated y.
    assert transform.TransformPoint((8., 16., 24.)) == pytest.approx((9.7, 17., 24.), abs=1e-5)


def test_demons_handles_different_input_grids():
    fixed = sitk.Image([16] * 3, sitk.sitkFloat32) + 1
    moving = sitk.Image([18, 17, 16], sitk.sitkFloat32) + 1
    moving.SetOrigin((-1., -1., 0.))
    moving.SetSpacing((1.1, 1., 1.))
    registered, transform, field = demons_registration(
        fixed, moving, resolution_staging=(1,), iteration_staging=(1,), smoothing_sigmas=0,
    )
    assert registered.GetSize() == fixed.GetSize()
    assert field.GetOrigin() == fixed.GetOrigin()
    assert transform.TransformPoint((8., 8., 8.)) == pytest.approx((8., 8., 8.), abs=1e-6)


def test_full_extent_crop_preserves_all_voxel_centres():
    values = np.arange(16 * 12 * 10, dtype=np.float32).reshape(10, 12, 16)
    image = sitk.GetImageFromArray(values)
    image.SetSpacing((0.7, 1.5, 2.5))
    image.SetOrigin((-10., 15., -4.))
    cropped = crop_image_to_extent(image, *get_image_physical_extent(image))
    assert cropped.GetSize() == image.GetSize()
    assert sitk.GetArrayFromImage(cropped) == pytest.approx(values)


def test_oblique_crop_keeps_patient_coordinates():
    image = sitk.Image([16, 12, 10], sitk.sitkFloat32) + 17
    image.SetSpacing((1., 2., 3.))
    image.SetDirection((0., -1., 0., 1., 0., 0., 0., 0., 1.))
    cropped = crop_image_to_extent(image, *get_image_physical_extent(image))
    for index in [(2, 3, 4), (10, 8, 6)]:
        point = image.TransformIndexToPhysicalPoint(index)
        assert cropped[cropped.TransformPhysicalPointToIndex(point)] == pytest.approx(17)


def test_rotated_reference_covers_requested_box():
    angle = 0.4
    direction = (np.cos(angle), -np.sin(angle), 0, np.sin(angle), np.cos(angle), 0, 0, 0, 1)
    image = create_reference_image_from_extent(np.array([-2., 1., 3.]), np.array([8., 12., 20.]),
                                               (1., 2., 3.), direction)
    for corner in product((-2., 8.), (1., 12.), (3., 20.)):
        index = np.array(image.TransformPhysicalPointToContinuousIndex(corner))
        assert np.all(index >= -1e-8)
        assert np.all(index <= np.array(image.GetSize()) - 1 + 1e-8)


def test_align_extents_matches_direction_and_preserves_moving_type():
    fixed = sitk.Image([12] * 3, sitk.sitkInt16)
    moving = sitk.Image([12] * 3, sitk.sitkFloat32) + 1.75
    moving.SetOrigin((11., 0., 0.))
    moving.SetDirection((0., -1., 0., 1., 0., 0., 0., 0., 1.))
    a, b, _ = align_image_extents(fixed, moving, use_initial_rigid=False)
    assert a.GetDirection() == b.GetDirection()
    assert b.GetPixelID() == sitk.sitkFloat32
    assert b[5, 5, 5] == pytest.approx(1.75)


def test_deformable_builder_rejects_left_handed_grid():
    image = field_image()
    image.SetDirection((1., 0., 0., 0., 1., 0., 0., 0., -1.))
    with pytest.raises(ValueError, match='right-handed'):
        DeformableSpatialRegistrationBuilder(make_ct_series()).add_registration(
            make_ct_series(), sitk.DisplacementFieldTransform(image))


def test_builder_checks_dataset_before_returning():
    builder = SpatialRegistrationBuilder(make_ct_series()).add_registration(make_ct_series(), np.eye(4))
    builder.set_instance_number(None)
    with pytest.raises(ValueError, match='InstanceNumber'):
        builder.build()


def test_real_demons_result_survives_dicom_export(tmp_path):
    from synthetic import make_textured_phantom
    from pydicomrt.reg.method import rigid_registration

    fixed = make_textured_phantom(size=(24, 24, 24), spacing=(1., 1.5, 2.))
    moving = sitk.Image(fixed)
    moving.SetOrigin(tuple(a + b for a, b in zip(fixed.GetOrigin(), (3., -2., 1.))))
    rigid = rigid_registration(fixed, moving, iterations=2, shrink_factors=(1,), smoothing_sigmas=(0,))
    moving_rigid = sitk.Resample(moving, fixed, rigid, sitk.sitkLinear, -1000.)
    _, residual, field = demons_registration(fixed, moving_rigid, resolution_staging=(1,),
                                             iteration_staging=(2,), smoothing_sigmas=0)
    ds = (DeformableSpatialRegistrationBuilder(make_ct_series())
          .add_registration(make_ct_series(), residual, post_transform=affine_to_homogeneous_matrix(rigid))
          .build())
    path = tmp_path / 'computed_deformable.dcm'
    ds.save_as(path, enforce_file_format=True)
    parsed = get_deformable_registrations(dcmread(path))[0]
    for index in [(5, 6, 7), (12, 10, 11)]:
        point = fixed.TransformIndexToPhysicalPoint(index)
        x, y, z = index
        vector = parsed['DeformableRegistrationGrid']['VectorGridData'][z, y, x]
        actual = matrix_point(parsed['PostDeformationMatrixRegistration'], np.asarray(point) + vector)
        expected = rigid.TransformPoint(residual.TransformPoint(point))
        assert actual == pytest.approx(expected, abs=1e-5)


def test_demons_regularization_is_in_mm_at_every_level():
    class RecordingFilter(ConstantRefinement):
        def __init__(self):
            self.widths = []

        def SetStandardDeviations(self, sigma):
            self.sigma = np.asarray(sigma)

        def GetStandardDeviations(self):
            return self.sigma.tolist()

        def Execute(self, fixed, moving):
            self.widths.append(self.sigma * fixed.GetSpacing())
            return super().Execute(fixed, moving)

    image = sitk.Image([16] * 3, sitk.sitkFloat32)
    image.SetSpacing((0.7, 1.5, 2.5))
    algorithm = RecordingFilter()
    multiscale_demons(algorithm, image, image, resolution_staging=(2, 1),
                     smoothing_sigmas=(0, 0), iteration_staging=(1, 1),
                     regularization_kernel_mm=(1.5, 1.5, 1.5))
    assert len(algorithm.widths) == 2
    for width in algorithm.widths:
        assert width == pytest.approx((1.5, 1.5, 1.5))


def test_demons_accepts_float32_initial_field_without_consuming_it():
    fixed = sitk.Image([16] * 3, sitk.sitkFloat32) + 1
    initial = sitk.Image([16] * 3, sitk.sitkVectorFloat32)
    _, transform, _ = demons_registration(fixed, fixed, initial_displacement_field=initial,
                                          resolution_staging=(1,), iteration_staging=(1,),
                                          smoothing_sigmas=0)
    assert initial.GetSize() == (16, 16, 16)
    assert transform.TransformPoint((8., 8., 8.)) == pytest.approx((8., 8., 8.))


def test_rigid_scale_does_not_accept_shear():
    matrix = np.eye(4)
    matrix[0, 1] = 0.3
    with pytest.raises(ValueError, match='orthogonal'):
        SpatialRegistrationBuilder(make_ct_series()).add_registration(
            make_ct_series(), matrix, matrix_type='RIGID_SCALE')


def test_registration_references_cannot_mix_frames():
    moving = make_ct_series()
    moving[-1].FrameOfReferenceUID = make_ct_series()[0].FrameOfReferenceUID
    with pytest.raises(ValueError, match='FrameOfReferenceUID'):
        SpatialRegistrationBuilder(make_ct_series()).add_registration(moving, np.eye(4))
