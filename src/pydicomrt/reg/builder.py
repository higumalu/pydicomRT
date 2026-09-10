"""
DICOM Spatial Registration and Deformable Spatial Registration builders.

Both objects are described in PS3.3:

- Spatial Registration IOD          A.39.1, Spatial Registration Module C.20.2
- Deformable Spatial Registration   A.39.2, Deformable Spatial Registration Module C.20.3

The matrices written here run **moving -> fixed**: a Frame of Reference Transformation
Matrix maps points in the item's Frame of Reference into the registered RCS. SimpleITK
resampling transforms run the other way, so invert before handing one to these builders.
"""

import datetime
import numbers
from importlib.metadata import PackageNotFoundError, version
from typing import TypeVar

import numpy as np

from pydicomrt.reg.type_transform import sitk_displacement_field_to_deformable_registration_grid

from abc import ABC, abstractmethod
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid, ImplicitVRLittleEndian

#: Fallback UID root (pydicom's), used when no prefix is set via set_uid_prefix().
DEFAULT_UID_PREFIX = "1.2.826.0.1.3680043.8.498."

try:
    _PYDICOMRT_VERSION = version("pydicomrt")
except PackageNotFoundError:  # pragma: no cover - running from a source tree
    _PYDICOMRT_VERSION = "unknown"

#: Values PS3.3 C.20.2 allows for Frame of Reference Transformation Matrix Type (0070,030C).
TRANSFORMATION_MATRIX_TYPES = ("RIGID", "RIGID_SCALE", "AFFINE")

# Registration Type Code Sequence entries, PS3.16 context group CID 7100.
IDENTITY_REGISTRATION_CODE = {
    "CodeValue": "125021",
    "CodingSchemeDesignator": "DCM",
    "CodeMeaning": "Frame of Reference Identity",
}
CONTENT_BASED_REGISTRATION_CODE = {
    "CodeValue": "125024",
    "CodingSchemeDesignator": "DCM",
    "CodeMeaning": "Image Content-based Alignment",
}


#: Self type for the chained setters, so a subclass keeps its own type through them.
#: typing.Self would say this more directly but arrives in 3.11, past this package's floor.
BuilderT = TypeVar("BuilderT", bound="BaseRegistrationBuilder")


def _code_sequence(code: dict) -> Sequence:
    item = Dataset()
    item.CodeValue = code["CodeValue"]
    item.CodingSchemeDesignator = code["CodingSchemeDesignator"]
    item.CodeMeaning = code["CodeMeaning"]
    return Sequence([item])


def validate_transformation_matrix(matrix, matrix_type: str = "RIGID") -> list:
    """
    Check a Frame of Reference Transformation Matrix before it is written.

    A malformed matrix does not fail at write time -- it produces a file that loads and
    misplaces the image, so the check happens here rather than downstream.

    Parameters
    ----------
    matrix : sequence or np.ndarray
        16 numbers in row-major order, or a 4x4 array.
    matrix_type : str
        One of ``TRANSFORMATION_MATRIX_TYPES``. ``"RIGID"`` additionally requires the upper
        3x3 to be a rotation (orthonormal, determinant +1).

    Returns
    -------
    list
        The matrix as 16 floats, row-major, ready for the DICOM element.

    Raises
    ------
    ValueError
        If the type is unknown, the shape is wrong, the values are not numeric, the bottom
        row is not ``[0, 0, 0, 1]``, or a RIGID matrix is not a rotation.
    """
    if matrix_type not in TRANSFORMATION_MATRIX_TYPES:
        raise ValueError(
            f"matrix_type must be one of {TRANSFORMATION_MATRIX_TYPES}, got {matrix_type!r}"
        )

    array = np.asarray(matrix, dtype=object).ravel()
    if array.size != 16:
        raise ValueError(
            f"transformation matrix must have 16 elements (4x4 row-major), got {array.size}"
        )
    if not all(isinstance(value, numbers.Real) for value in array.tolist()):
        raise ValueError("transformation matrix must contain only real numbers")

    values = np.asarray(array.tolist(), dtype=float).reshape(4, 4)
    if not np.isfinite(values).all():
        raise ValueError("transformation matrix contains NaN or infinity")
    if not np.allclose(values[3], [0.0, 0.0, 0.0, 1.0], atol=1e-6):
        raise ValueError(
            f"the bottom row of the matrix must be [0, 0, 0, 1], got {values[3].tolist()}. "
            "The matrix is row-major: translation belongs in the last column."
        )

    if matrix_type == "RIGID":
        rotation = values[:3, :3]
        if not np.allclose(rotation @ rotation.T, np.identity(3), atol=1e-4):
            raise ValueError(
                "matrix_type is RIGID but the upper 3x3 is not orthonormal; use "
                "'RIGID_SCALE' or 'AFFINE' if scaling or shear is intended"
            )
        if not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-4):
            raise ValueError(
                "matrix_type is RIGID but the upper 3x3 has determinant "
                f"{np.linalg.det(rotation):.6f}; a rotation has determinant +1"
            )

    if matrix_type == "RIGID_SCALE":
        gram = values[:3, :3].T @ values[:3, :3]
        if np.any(np.diag(gram) <= 0) or not np.allclose(gram, np.diag(np.diag(gram)), atol=1e-4):
            raise ValueError("RIGID_SCALE requires orthogonal, non-zero axes; use AFFINE for shear")
    return values.ravel().tolist()


def _matrix_registration_item(matrix, matrix_type: str, registration_code: dict) -> Dataset:
    """
    Build one Matrix Registration Sequence (0070,0309) item.

    Registration Type Code Sequence (0070,030D) sits *inside* this item, alongside Matrix
    Sequence (0070,030A) -- not one level up in the Registration Sequence item.
    """
    matrix_item = Dataset()
    matrix_item.FrameOfReferenceTransformationMatrixType = matrix_type
    matrix_item.FrameOfReferenceTransformationMatrix = validate_transformation_matrix(
        matrix, matrix_type
    )

    registration_item = Dataset()
    registration_item.MatrixSequence = Sequence([matrix_item])
    registration_item.RegistrationTypeCodeSequence = _code_sequence(registration_code)
    return registration_item


def _referenced_image_sequence(ds_list) -> Sequence:
    sequence = Sequence()
    for ds in ds_list:
        referenced_image = Dataset()
        referenced_image.ReferencedSOPClassUID = ds.SOPClassUID
        referenced_image.ReferencedSOPInstanceUID = ds.SOPInstanceUID
        sequence.append(referenced_image)
    return sequence


# --------------------- Basic Registration Builder ---------------------
class BaseRegistrationBuilder(ABC):
    """
    Shared scaffolding for the two registration IODs.

    Patient, study and frame-of-reference context is copied from the *fixed* series, so the
    object lands in the same study as the images it describes. The fixed series' Frame of
    Reference becomes the registered RCS, which is what every matrix is expressed relative to.
    """

    def __init__(self, fixed_series):
        if not fixed_series:
            raise ValueError("fixed_series must be non-empty")
        self.fixed_series = list(fixed_series)
        self.uid_prefix = None
        self.ref_ds = self.fixed_series[0]
        self.instance_number = 1
        self.registration_dataset_list = []
        # Every series this object references, keyed by Series Instance UID, for the
        # Common Instance Reference module.
        self._referenced_series = {}
        self._track_referenced_series(self.fixed_series)

    def set_uid_prefix(self: BuilderT, uid_prefix: str) -> BuilderT:
        """
        Set the root under which generated UIDs are minted.

        Parameters
        ----------
        uid_prefix : str
            Organisation UID root, dot-terminated, e.g. ``"1.2.826.0.1.3680043.2.1125."``.
            Unlike the RTSTRUCT and dose builders, the registration builders do **not**
            read ``DICOM_UID_PREFIX``; unset, they fall back to pydicom's default root.

        Returns
        -------
        SpatialRegistrationBuilder or DeformableSpatialRegistrationBuilder
            self, so calls chain.
        """
        self.uid_prefix = uid_prefix
        return self

    def set_instance_number(self: BuilderT, instance_number: int) -> BuilderT:
        """
        Set InstanceNumber (0020,0013).

        Parameters
        ----------
        instance_number : int
            Position within the series. Default 1. Give each object a distinct number when
            writing several registrations into one series -- InstanceNumber is Type 1 here,
            and some archives reject duplicates within a series.

        Returns
        -------
        SpatialRegistrationBuilder or DeformableSpatialRegistrationBuilder
            self, so calls chain.
        """
        self.instance_number = instance_number
        return self

    @property
    def _effective_uid_prefix(self) -> str:
        return self.uid_prefix or DEFAULT_UID_PREFIX

    @property
    @abstractmethod
    def _sop_class_uid(self) -> str:
        raise NotImplementedError

    @property
    @abstractmethod
    def _base_file_name(self) -> str:
        raise NotImplementedError

    def _track_referenced_series(self, ds_list):
        frames = {str(getattr(ds, "FrameOfReferenceUID", "")) for ds in ds_list}
        if len(frames) != 1 or "" in frames:
            raise ValueError("reference images must share one non-empty FrameOfReferenceUID")
        for ds in ds_list:
            series_uid = getattr(ds, "SeriesInstanceUID", None)
            if series_uid is None:
                continue
            instances = self._referenced_series.setdefault(series_uid, {})
            instances[ds.SOPInstanceUID] = ds.SOPClassUID

    def _generate_file_meta(self):
        file_meta = FileMetaDataset()
        file_meta.FileMetaInformationVersion = b"\x00\x01"
        file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = self._sop_class_uid
        file_meta.MediaStorageSOPInstanceUID = generate_uid(prefix=self.uid_prefix)
        # A bare "1" is not a usable Implementation Class UID; fall back to a real root.
        file_meta.ImplementationClassUID = self._effective_uid_prefix + "1"
        file_meta.ImplementationVersionName = "pydicomRT"
        return file_meta

    def _generate_base_dataset(self):
        file_meta = self._generate_file_meta()
        ds = FileDataset(self._base_file_name, {}, file_meta=file_meta, preamble=b"\0" * 128)
        return ds

    def _add_required_elements(self, ds: Dataset):
        dt = datetime.datetime.now()
        ds.SpecificCharacterSet = "ISO_IR 192"
        ds.InstanceCreationDate = dt.strftime("%Y%m%d")
        ds.InstanceCreationTime = dt.strftime("%H%M%S")
        ds.Modality = "REG"
        ds.Manufacturer = ""
        ds.ManufacturerModelName = "pydicomRT"
        ds.InstitutionName = ""
        # Which build produced this object; the only equipment attribute we can honestly
        # fill in, since the rest describe an acquisition device this is not.
        ds.SoftwareVersions = f"pydicomRT {_PYDICOMRT_VERSION}"
        ds.SOPClassUID = ds.file_meta.MediaStorageSOPClassUID
        ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID
        return ds

    def _add_patient_information_from_ref_ds(self, ds: Dataset, ref_ds: Dataset):
        ds.PatientID = getattr(ref_ds, "PatientID", "")
        ds.PatientName = getattr(ref_ds, "PatientName", "")
        ds.PatientBirthDate = getattr(ref_ds, "PatientBirthDate", "")
        ds.PatientBirthTime = getattr(ref_ds, "PatientBirthTime", "")
        ds.PatientSex = getattr(ref_ds, "PatientSex", "")
        ds.PatientAge = getattr(ref_ds, "PatientAge", "")
        ds.PatientWeight = getattr(ref_ds, "PatientWeight", "")
        ds.PatientPosition = getattr(ref_ds, "PatientPosition", "")
        return ds

    def _add_study_information_from_ref_ds(self, ds: Dataset, ref_ds: Dataset):
        ds.StudyInstanceUID = getattr(ref_ds, "StudyInstanceUID", "")
        ds.StudyID = getattr(ref_ds, "StudyID", "")
        ds.StudyDescription = getattr(ref_ds, "StudyDescription", "")
        ds.StudyDate = getattr(ref_ds, "StudyDate", "")
        ds.StudyTime = getattr(ref_ds, "StudyTime", "")
        ds.AccessionNumber = getattr(ref_ds, "AccessionNumber", "")
        ds.ReferringPhysicianName = getattr(ref_ds, "ReferringPhysicianName", "")
        ds.PhysiciansOfRecord = getattr(ref_ds, "PhysiciansOfRecord", "")
        return ds

    def _add_series_information(
        self,
        ds: Dataset,
        ref_ds: Dataset,
        series_number: str = "50",
        series_desc_prefix: str = "Spatial Registration",
        content_label: str = "SPATIAL_REG",
        ):
        dt = datetime.datetime.now()
        ds.SeriesInstanceUID = generate_uid(prefix=self.uid_prefix)
        ds.SeriesNumber = series_number
        ds.SeriesDescription = f"{series_desc_prefix} {dt.strftime('%Y%m%d%H%M')}"
        ds.SeriesDate = dt.strftime("%Y%m%d")
        ds.SeriesTime = dt.strftime("%H%M%S")

        # Content identification, all Type 1 in the registration modules.
        ds.ContentDate = dt.strftime("%Y%m%d")
        ds.ContentTime = dt.strftime("%H%M%S")
        ds.InstanceNumber = self.instance_number
        # ContentLabel is VR CS, capped at 16 characters -- "SPATIAL REGISTRATION"
        # overflowed it and made the element non-conformant.
        ds.ContentLabel = content_label
        ds.ContentDescription = series_desc_prefix
        ds.ContentCreatorName = "pydicomRT"

        # The registered RCS: Spatial REG maps into this frame; Deformable REG
        # samples from this frame into the Source Frame of Reference.
        ds.FrameOfReferenceUID = getattr(ref_ds, "FrameOfReferenceUID", "")
        ds.PositionReferenceIndicator = getattr(ref_ds, "PositionReferenceIndicator", "")
        return ds

    def _add_referenced_series_sequence(self, ds: Dataset):
        """
        Common Instance Reference module: which series this object was built from.

        Optional in the IOD, but it is how a viewer finds the images a registration belongs
        to without opening every series in the study.
        """
        sequence = Sequence()
        for series_uid, instances in self._referenced_series.items():
            series_item = Dataset()
            series_item.SeriesInstanceUID = series_uid
            series_item.ReferencedInstanceSequence = Sequence()
            for sop_instance_uid, sop_class_uid in instances.items():
                instance_item = Dataset()
                instance_item.ReferencedSOPClassUID = sop_class_uid
                instance_item.ReferencedSOPInstanceUID = sop_instance_uid
                series_item.ReferencedInstanceSequence.append(instance_item)
            sequence.append(series_item)
        ds.ReferencedSeriesSequence = sequence
        return ds


# --------------------- Spatial Registration Builder ---------------------

class SpatialRegistrationBuilder(BaseRegistrationBuilder):
    """
    Build a Spatial Registration (rigid/affine) object.

    Serialises a 4x4 matrix as DICOM Spatial Registration (SOP class
    ``1.2.840.10008.5.1.4.1.1.66.1``). Several moving series can be registered to one
    fixed series in a single object.

    Parameters
    ----------
    fixed_series : list of Dataset
        The series everything is registered *to*. Patient and study context is copied
        from it, and its Frame of Reference becomes the registered RCS -- the frame every
        stored matrix is expressed relative to.

    Methods
    -------
    add_registration(moving_series, matrix, matrix_type)
        Register one moving series. Call repeatedly for several.
    set_instance_number(instance_number)
        Set InstanceNumber.
    set_uid_prefix(uid_prefix)
        Set the UID root.
    build(include_identity=True)
        Assemble the ``FileDataset``.

    Raises
    ------
    ValueError
        If ``fixed_series`` is empty.

    See Also
    --------
    get_spatial_registrations : Read the result back.
    DeformableSpatialRegistrationBuilder : When the transform is a displacement field.
    check_spatial_reg_iod : Validate before writing.

    Notes
    -----
    **The matrix runs moving to fixed**, which is the opposite of what the registration
    functions in :mod:`pydicomrt.reg.method` return. Those return *resampling* transforms
    (fixed to moving, the direction ``sitk.Resample`` wants), so invert before passing one
    here. Getting this backwards produces a valid file that misplaces the image by twice
    the offset.

    Examples
    --------
    >>> transform = rigid_registration(fixed_image, moving_image)     # doctest: +SKIP
    >>> matrix = affine_to_homogeneous_matrix(transform.GetInverse()) # doctest: +SKIP
    >>> reg_ds = (                                                    # doctest: +SKIP
    ...     SpatialRegistrationBuilder(fixed_series)
    ...     .add_registration(moving_series, matrix.astype("float32").ravel().tolist())
    ...     .build()
    ... )
    """

    def __init__(self, fixed_series):
        super().__init__(fixed_series)

    @property
    def _sop_class_uid(self) -> str:
        # Spatial Registration (rigid/affine) SOP Class UID
        return "1.2.840.10008.5.1.4.1.1.66.1"

    @property
    def _base_file_name(self) -> str:
        return "spatial_registration"

    def add_registration(
        self,
        moving_series,
        rigid_transform_matrix,
        matrix_type: str = "RIGID",
    ) -> "SpatialRegistrationBuilder":
        """
        Register a moving series to the fixed series' Frame of Reference.

        Parameters
        ----------
        moving_series : list[Dataset]
            The series being registered.
        rigid_transform_matrix : sequence
            16 numbers, row-major 4x4, mapping **moving -> fixed**. A SimpleITK resampling
            transform runs fixed -> moving, so invert it first.
        matrix_type : str
            ``"RIGID"``, ``"RIGID_SCALE"`` or ``"AFFINE"``. RIGID is validated as a rotation.

        Returns
        -------
        SpatialRegistrationBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If the moving series is empty or the matrix does not validate.
        """
        if not moving_series:
            raise ValueError("moving_series must be non-empty")

        moving_series = list(moving_series)
        self._track_referenced_series(moving_series)

        registration_item = Dataset()
        registration_item.FrameOfReferenceUID = getattr(
            moving_series[0], "FrameOfReferenceUID", ""
        )
        registration_item.ReferencedImageSequence = _referenced_image_sequence(moving_series)
        registration_item.MatrixRegistrationSequence = Sequence([
            _matrix_registration_item(
                rigid_transform_matrix, matrix_type, CONTENT_BASED_REGISTRATION_CODE
            )
        ])

        self.registration_dataset_list.append(registration_item)
        return self

    def _identity_registration_item(self):
        """
        The item registering the fixed Frame of Reference to itself.

        Real systems emit this (the reference file from a Varian TPS does), and it is what
        marks which Frame of Reference the other matrices are expressed relative to. Without
        it a reader has only the top-level Frame of Reference UID to go on.
        """
        registration_item = Dataset()
        registration_item.FrameOfReferenceUID = getattr(
            self.ref_ds, "FrameOfReferenceUID", ""
        )
        registration_item.ReferencedImageSequence = _referenced_image_sequence(
            self.fixed_series
        )
        registration_item.MatrixRegistrationSequence = Sequence([
            _matrix_registration_item(
                np.identity(4).ravel().tolist(), "RIGID", IDENTITY_REGISTRATION_CODE
            )
        ])
        return registration_item

    def build(self, include_identity: bool = True) -> FileDataset:
        """
        Assemble the dataset.

        Parameters
        ----------
        include_identity : bool
            Append the identity item for the fixed Frame of Reference. Leave it on unless a
            downstream system objects to it.

        Raises
        ------
        ValueError
            If no registration has been added or the assembled dataset fails validation.
        """
        if not self.registration_dataset_list:
            raise ValueError(
                "no registration added; call add_registration() before build()"
            )

        ds = self._generate_base_dataset()
        self._add_required_elements(ds)
        self._add_patient_information_from_ref_ds(ds, self.ref_ds)
        self._add_study_information_from_ref_ds(ds, self.ref_ds)
        self._add_series_information(ds, self.ref_ds)
        self._add_referenced_series_sequence(ds)

        ds.RegistrationSequence = Sequence(list(self.registration_dataset_list))
        if include_identity:
            ds.RegistrationSequence.append(self._identity_registration_item())
        from .check import check_spatial_reg_iod
        result = check_spatial_reg_iod(ds)
        if not result["result"]:
            raise ValueError("Invalid Spatial Registration: " + "; ".join(result["content"]))
        return ds


# --------------------- Deformable Spatial Registration Builder ---------------------

class DeformableSpatialRegistrationBuilder(BaseRegistrationBuilder):
    """
    Build a Deformable Spatial Registration object.

    Serialises a displacement field, optionally bracketed by two affine matrices, as
    DICOM Deformable Spatial Registration (SOP class ``1.2.840.10008.5.1.4.1.1.66.3``).

    The three parts compose in a fixed order: **pre-matrix, then displacement grid, then
    post-matrix**. That is what makes the usual rigid-then-deformable pipeline
    expressible: for a resampling pipeline ``rigid(deform(point))``, use identity
    for ``pre_transform``, the residual field for the grid, and the forward rigid
    resampling matrix for ``post_transform``.

    Parameters
    ----------
    fixed_series : list of Dataset
        The series everything is registered to. Supplies patient and study context, and
        its Frame of Reference becomes the registered RCS.

    Methods
    -------
    add_registration(moving_series, vectorial_field_transform, pre_transform, post_transform)
        Register one moving series.
    set_instance_number(instance_number)
        Set InstanceNumber.
    set_uid_prefix(uid_prefix)
        Set the UID root.
    build()
        Assemble the ``FileDataset``.

    Raises
    ------
    ValueError
        If ``fixed_series`` is empty.

    See Also
    --------
    get_deformable_registrations : Read the result back.
    SpatialRegistrationBuilder : When the transform is a single matrix.
    demons_registration, bspline_registration : Produce the displacement transform.

    Notes
    -----
    Unlike Spatial Registration, this object maps **fixed to moving** (Registered RCS
    to Source RCS), as specified by DICOM PS3.3 C.20.3.1.1. Do not invert the
    resampling transform. The displacement grid is located in the fixed frame.
    At a fixed grid point ``x``, the DICOM equation is
    ``source = post(pre(x) + displacement(x))``: the vector is indexed on the original
    registered grid, even when a pre-matrix is present.

    A displacement grid is bulky -- three ``float32`` per voxel, so roughly 12 bytes where
    the CT stores 2. Full-resolution fields on a whole CT run to hundreds of megabytes.

    Examples
    --------
    >>> rigid = rigid_registration(fixed_image, moving_image)          # doctest: +SKIP
    >>> moving_rigid = sitk.Resample(                                  # doctest: +SKIP
    ...     moving_image, fixed_image, rigid, sitk.sitkLinear, -1000.0,
    ...     moving_image.GetPixelID())
    >>> _, deform, _ = demons_registration(fixed_image, moving_rigid)  # doctest: +SKIP
    >>> reg_ds = (                                                     # doctest: +SKIP
    ...     DeformableSpatialRegistrationBuilder(fixed_series)
    ...     .add_registration(
    ...         moving_series=moving_series,
    ...         vectorial_field_transform=deform,
    ...         pre_transform=np.eye(4).ravel().tolist(),
    ...         post_transform=affine_to_homogeneous_matrix(rigid).ravel().tolist(),
    ...     )
    ...     .build()
    ... )
    """

    def __init__(self, fixed_series):
        super().__init__(fixed_series)

    @property
    def _sop_class_uid(self) -> str:
        # Deformable Spatial Registration SOP Class UID
        return "1.2.840.10008.5.1.4.1.1.66.3"

    @property
    def _base_file_name(self) -> str:
        return "deformable_spatial_registration"

    def add_registration(
        self,
        moving_series,
        vectorial_field_transform,
        pre_transform=None,
        post_transform=None,
        pre_transform_type: str = "RIGID",
        post_transform_type: str = "RIGID",
    ) -> "DeformableSpatialRegistrationBuilder":
        """
        Register a moving series to the fixed series' Frame of Reference.

        Parameters
        ----------
        moving_series : list[Dataset]
            The series being registered.
        vectorial_field_transform : sitk.DisplacementFieldTransform
            The displacement field, converted to a Deformable Registration Grid.
        pre_transform, post_transform : sequence or None, optional
            16 numbers each, row-major 4x4, applied before and after the grid. Pass an
            identity matrix or None for either if unused (None omits the sequence).
        pre_transform_type, post_transform_type : str
            ``"RIGID"``, ``"RIGID_SCALE"`` or ``"AFFINE"``.

        Returns
        -------
        DeformableSpatialRegistrationBuilder
            self, so calls chain.
        """
        if not moving_series:
            raise ValueError("moving_series must be non-empty")

        moving_series = list(moving_series)
        self._track_referenced_series(moving_series)

        registration_item = Dataset()
        registration_item.SourceFrameOfReferenceUID = getattr(
            moving_series[0], "FrameOfReferenceUID", ""
        )
        registration_item.ReferencedImageSequence = _referenced_image_sequence(moving_series)

        deformable_registration_grid = sitk_displacement_field_to_deformable_registration_grid(
            vectorial_field_transform
        )
        registration_item.DeformableRegistrationGridSequence = Sequence([
            deformable_registration_grid
        ])

        for key, matrix, matrix_type in (
            ("PreDeformationMatrixRegistrationSequence", pre_transform, pre_transform_type),
            ("PostDeformationMatrixRegistrationSequence", post_transform, post_transform_type),
        ):
            if matrix is not None:
                item = Dataset()
                item.FrameOfReferenceTransformationMatrixType = matrix_type
                item.FrameOfReferenceTransformationMatrix = validate_transformation_matrix(matrix, matrix_type)
                setattr(registration_item, key, Sequence([item]))

        # For the deformable IOD this sits at the Deformable Registration Sequence item
        # level -- unlike the rigid IOD, where it lives inside Matrix Registration Sequence.
        registration_item.RegistrationTypeCodeSequence = _code_sequence(
            CONTENT_BASED_REGISTRATION_CODE
        )

        self.registration_dataset_list.append(registration_item)
        return self

    def build(self) -> FileDataset:
        """
        Assemble the Deformable Spatial Registration object.

        Returns
        -------
        FileDataset
            Ready to ``save_as(path, enforce_file_format=True)``.

        Raises
        ------
        ValueError
            If no registration has been added or the assembled dataset fails validation.

        Notes
        -----
        Unlike :meth:`SpatialRegistrationBuilder.build` there is no ``include_identity``
        parameter. The deformable IOD carries the fixed frame in
        the top-level FrameOfReferenceUID rather than needing an identity item to mark it.
        """
        if not self.registration_dataset_list:
            raise ValueError("no registration added; call add_registration() before build()")

        ds = self._generate_base_dataset()
        self._add_required_elements(ds)
        self._add_patient_information_from_ref_ds(ds, self.ref_ds)
        self._add_study_information_from_ref_ds(ds, self.ref_ds)
        self._add_series_information(
            ds,
            self.ref_ds,
            series_desc_prefix="Deformable Spatial Registration",
            content_label="DEFORMABLE_REG",
        )
        self._add_referenced_series_sequence(ds)

        ds.DeformableRegistrationSequence = Sequence(list(self.registration_dataset_list))
        from .check import check_deformable_reg_iod
        result = check_deformable_reg_iod(ds)
        if not result["result"]:
            raise ValueError("Invalid Deformable Registration: " + "; ".join(result["content"]))
        return ds
