"""Spatial and deformable registration module"""

from .builder import (
    SpatialRegistrationBuilder,
    DeformableSpatialRegistrationBuilder,
)

from .parser import (
    get_spatial_registrations,
    get_deformable_registrations,
)

from .check import (
    check_spatial_reg_iod,
    check_deformable_reg_iod,
)

from .type_transform import (
    affine_to_homogeneous_matrix,
    sitk_displacement_field_to_deformable_registration_grid,
)

__all__ = [
    'SpatialRegistrationBuilder',
    'DeformableSpatialRegistrationBuilder',

    'get_spatial_registrations',
    'get_deformable_registrations',

    'check_spatial_reg_iod',
    'check_deformable_reg_iod',

    'affine_to_homogeneous_matrix',
    'sitk_displacement_field_to_deformable_registration_grid',
]
