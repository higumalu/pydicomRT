"""
Registration pipeline module.

Provides a packaged registration pipeline integrating rigid and deformable registration workflows.
"""

from .registration_pipeline import registration_pipeline, compose_transforms
from .preprocessing import (
    preprocess_image,
    infer_default_pixel_value,
    window_clip,
    align_image_extents,
    get_image_physical_extent,
    get_intersection_extent,
    crop_image_to_extent,
    get_initial_rigid_transform,
    get_image_center,
    get_images_distance,
    create_reference_image_from_extent,
)

__all__ = [
    "registration_pipeline",
    "compose_transforms",
    "preprocess_image",
    "infer_default_pixel_value",
    "window_clip",
    "align_image_extents",
    "get_image_physical_extent",
    "get_intersection_extent",
    "crop_image_to_extent",
    "get_initial_rigid_transform",
    "get_image_center",
    "get_images_distance",
    "create_reference_image_from_extent",
]
