

class SpatialRegistrationBuilder:
    def __init__(self, fixed_ds_list, moving_ds_list, transform_matrix):
        self.fixed_ds_list = fixed_ds_list
        self.moving_ds_list = moving_ds_list
        self.transform_matrix = transform_matrix


class DeformableSpatialRegistrationBuilder:
    def __init__(
        self, fixed_ds_list, moving_ds_list,
        pre_transform,
        vectorial_field,
        post_transform
        ):
        self.fixed_ds_list = fixed_ds_list
        self.moving_ds_list = moving_ds_list
        self.pre_transform = pre_transform
        self.vectorial_field = vectorial_field
        self.post_transform = post_transform