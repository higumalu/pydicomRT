import os
import numpy as np
import SimpleITK as sitk

from datetime import datetime
from pydicom.tag import Tag
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid
from pydicom.uid import ImplicitVRLittleEndian, RTDoseStorage, RTPlanStorage

from pydicomrt.utils.sitk_transform import (
    sitk_direction_to_image_orientation_patient,
    sitk_spacing_to_pixel_spacing,
)


DICOM_UID_PREFIX = os.getenv('DICOM_UID_PREFIX', "1.2.826.0.1.3680043.8.498.")  # PYDICOM_ROOT_UID

#: Stored dose values are unsigned 32-bit, so this is the largest one.
_STORED_MAX = 2 ** 32 - 1

#: Default Dose Grid Scaling. With unsigned 32-bit storage this covers up to ~429 Gy at
#: 0.1 microGy resolution; add_dose_grid_to_ds() scales up automatically beyond that.
DEFAULT_DOSE_GRID_SCALING = 1e-7

#: Dose Summation Types that make Referenced RT Plan Sequence Type 1C (PS3.3 C.8.8.3.2).
SUMMATION_TYPES_REQUIRING_PLAN = frozenset({
    "PLAN", "MULTI_PLAN", "PLAN_OVERVIEW", "FRACTION", "BEAM", "BRACHY",
    "FRACTION_SESSION", "BEAM_SESSION", "BRACHY_SESSION", "CONTROL_POINT",
})


def generate_base_dataset(uid_prefix: str = None) -> FileDataset:
    """
    Create an empty RT Dose dataset with its required elements.

    File meta, SOP class, the mandatory Type 1 and Type 2 elements, and empty sequences.
    Carries no patient, study or dose data. Used by :class:`RTDoseBuilder`.

    Parameters
    ----------
    uid_prefix : str, optional
        UID root. Defaults to the ``DICOM_UID_PREFIX`` environment variable.

    Returns
    -------
    FileDataset
        An RT Dose skeleton.

    See Also
    --------
    RTDoseBuilder : The supported entry point.
    cp_information_from_ds : Add patient and study context to the result.
    """
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    file_name = "bear_dose"
    file_meta = get_file_meta(uid_prefix)
    ds = FileDataset(file_name, {}, file_meta=file_meta, preamble=b"\0" * 128)
    add_required_elements_to_ds(ds)
    add_sequence_lists_to_ds(ds)
    return ds

def get_file_meta(uid_prefix: str = None) -> FileMetaDataset:
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    # Group length is recomputed on write, so it is not set here.
    file_meta = FileMetaDataset()
    file_meta.FileMetaInformationVersion = b"\x00\x01"
    file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
    file_meta.MediaStorageSOPClassUID = RTDoseStorage
    file_meta.MediaStorageSOPInstanceUID = (
        generate_uid(prefix=uid_prefix)
    )
    file_meta.ImplementationClassUID = uid_prefix + "1"
    file_meta.ImplementationVersionName = "pydicomRT"
    return file_meta

def add_required_elements_to_ds(ds: FileDataset):
    dt = datetime.now()
    # Append data elements required by the DICOM standarad
    ds.SpecificCharacterSet = "ISO_IR 192"
    ds.InstanceCreationDate = dt.strftime("%Y%m%d")
    ds.InstanceCreationTime = dt.strftime("%H%M%S")
    ds.ContentDate = dt.strftime("%Y%m%d")
    ds.ContentTime = dt.strftime("%H%M%S")
    ds.Modality = "RTDOSE"
    ds.Manufacturer = "pydicomRT"
    ds.ManufacturerModelName = "modelv1"
    ds.InstitutionName = "pydicomRT"
    # General Image module, Type 2 -- was missing entirely.
    ds.InstanceNumber = 1
    # Set values already defined in the file meta
    ds.SOPClassUID = ds.file_meta.MediaStorageSOPClassUID
    ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID

def add_sequence_lists_to_ds(ds: FileDataset):
    ds.ReferencedRTPlanSequence = Sequence()

# ---------------------------------------------------------------------------------------- #

def add_patient_information(ds: FileDataset, reference_ds: Dataset):
    ds.PatientName = getattr(reference_ds, "PatientName", "Unknown")
    ds.PatientID = getattr(reference_ds, "PatientID", "Unknown")
    ds.PatientBirthDate = getattr(reference_ds, "PatientBirthDate", "")
    ds.PatientSex = getattr(reference_ds, "PatientSex", "")
    ds.PatientAge = getattr(reference_ds, "PatientAge", "")
    ds.PatientSize = getattr(reference_ds, "PatientSize", "")
    ds.PatientWeight = getattr(reference_ds, "PatientWeight", "")
    return ds

def add_study_information(ds: FileDataset, reference_ds: Dataset, uid_prefix: str = None):
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    dt = datetime.now()
    ds.StudyInstanceUID = getattr(reference_ds, "StudyInstanceUID", generate_uid(prefix=uid_prefix))
    ds.StudyDate = getattr(reference_ds, "StudyDate", dt.strftime("%Y%m%d"))
    ds.StudyTime = getattr(reference_ds, "StudyTime", dt.strftime("%H%M%S"))
    ds.StudyID = getattr(reference_ds, "StudyID", "")
    ds.AccessionNumber = getattr(reference_ds, "AccessionNumber", "")
    ds.ReferringPhysicianName = getattr(reference_ds, "ReferringPhysicianName", "")
    return ds

def add_series_information(ds: FileDataset, reference_ds: Dataset, series_number: int = 1, uid_prefix: str = None):
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    dt = datetime.now()
    ds.Modality = "RTDOSE"
    ds.SeriesInstanceUID = generate_uid(prefix=uid_prefix)
    ds.SeriesDate = dt.strftime("%Y%m%d")
    ds.SeriesTime = dt.strftime("%H%M%S")
    ds.SeriesDescription = getattr(reference_ds, "SeriesDescription", "") + "_dose" + dt.strftime("%Y%m%d%H%M%S")
    # A new series gets its own number; copying the reference image's used to make the dose
    # claim membership of the series it was merely derived from.
    ds.SeriesNumber = series_number
    ds.OperatorsName = getattr(reference_ds, "OperatorsName", "")
    return ds

def add_frame_of_reference_information(ds: FileDataset, reference_ds: Dataset, uid_prefix: str = None):
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    ds.FrameOfReferenceUID = getattr(reference_ds, "FrameOfReferenceUID", generate_uid(prefix=uid_prefix))
    ds.PositionReferenceIndicator = getattr(reference_ds, "PositionReferenceIndicator", "")
    return ds

def add_rf_rt_plan_seq_from_dose_ds(ds: FileDataset, reference_ds: Dataset):
    sop_class_uid = getattr(getattr(reference_ds, "file_meta", None), "MediaStorageSOPClassUID", None)
    if sop_class_uid != RTDoseStorage:
        raise ValueError("reference_ds is not an RT Dose (SOPClassUID does not match RT Dose storage).")
    rf_rt_plan_seq = getattr(reference_ds, "ReferencedRTPlanSequence", None)
    if rf_rt_plan_seq is None or len(rf_rt_plan_seq) == 0:
        raise ValueError("reference_ds does not have a ReferencedRTPlanSequence.")
    ds.ReferencedRTPlanSequence = rf_rt_plan_seq
    return ds

def add_rf_rt_plan_seq_from_plan_ds(ds: FileDataset, reference_ds: Dataset):
    # 1.2.840.10008.5.1.4.1.1.481.3 is RT Structure Set, not RT Plan; this used to check for
    # that one, so it rejected every real plan and accepted a structure set instead.
    sop_class_uid = getattr(getattr(reference_ds, "file_meta", None), "MediaStorageSOPClassUID", None)
    if sop_class_uid != RTPlanStorage:
        raise ValueError(
            "reference_ds is not an RT Plan (SOPClassUID does not match RT Plan storage)."
        )


    rf_rt_plan_seq_block = Dataset()
    rf_rt_plan_seq_block.ReferencedSOPClassUID = reference_ds.SOPClassUID
    rf_rt_plan_seq_block.ReferencedSOPInstanceUID = reference_ds.SOPInstanceUID
    # TODO ReferencedFractionGroupSequence, ReferencedPlanOverviewIndex
    '''
    (300C,0020) ReferencedFractionGroupSequence
    Required if Dose Summation Type (3004,000A) is 
    FRACTION, BEAM, BRACHY, FRACTION_SESSION, BEAM_SESSION, BRACHY_SESSION or CONTROL_POINT.

    (300C,0118) ReferencedPlanOverviewIndex
    The value of Plan Overview Index (300C,0117) from 
    the Plan Overview Sequence (300C,0116) to which this RT Plan corresponds.
    Shall be unique, i.e., not be duplicated within another Item of this Referenced RT Plan Sequence (300C,0002).
    Required if Plan Overview Sequence (300C,0116) is present.
    '''
    ds.ReferencedRTPlanSequence.append(rf_rt_plan_seq_block)
    return ds


# ---------------------------------------------------------------------------------------- #

def cp_information_from_ds(
    ds: FileDataset,
    reference_ds: Dataset,
    series_number: int = 1,
    uid_prefix: str = None,
    ):
    """
    Copy patient, study and frame-of-reference context into a dose dataset.

    What makes the dose file into the same study as the images it was computed on, rather
    than appearing as an unrelated study in the archive.

    Parameters
    ----------
    ds : FileDataset
        The dose dataset to populate, from :func:`generate_base_dataset`. Modified in
        place.
    reference_ds : Dataset
        A planning CT slice or the RT Plan. Patient identity, StudyInstanceUID and
        FrameOfReferenceUID are read from it.
    series_number : int, optional
        SeriesNumber for the dose. Default 1.
    uid_prefix : str, optional
        UID root for the newly minted Series and SOP Instance UIDs. Defaults to the
        ``DICOM_UID_PREFIX`` environment variable.

    Returns
    -------
    FileDataset
        ``ds``, modified in place.

    See Also
    --------
    RTDoseBuilder : The supported entry point.

    Notes
    -----
    Study identity is inherited; series identity is new. A dose is a new series within an
    existing study, so reusing the reference's SeriesInstanceUID would file it as more
    slices of the CT.
    """
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    ds = add_patient_information(ds, reference_ds)
    ds = add_study_information(ds, reference_ds, uid_prefix)
    ds = add_series_information(ds, reference_ds, series_number, uid_prefix)
    ds = add_frame_of_reference_information(ds, reference_ds, uid_prefix)
    return ds

# ---------------------------------------------------------------------------------------- #

def add_dose_grid_to_ds(
    ds: FileDataset,
    dose_sitk_image: sitk.Image,
    dose_grid_scaling: float = None,
    dose_type: str = "PHYSICAL",
    dose_summation_type: str = "PLAN",
    dose_units: str = "GY",
    ):
    """
    Attach a dose grid to an RT Dose dataset.

    Parameters
    ----------
    ds : FileDataset
        Dataset from :func:`generate_base_dataset`, with patient context already copied.
    dose_sitk_image : sitk.Image
        The dose grid. Values are in ``dose_units``.
    dose_grid_scaling : float, optional
        Stored value -> dose factor. Defaults to 1e-7, raised automatically if the maximum
        dose would not otherwise fit in unsigned 32-bit.
    dose_type : str
        ``"PHYSICAL"``, ``"EFFECTIVE"`` or ``"ERROR"``. The previous default was EFFECTIVE,
        which means biologically weighted -- rarely what a plain dose grid holds.
    dose_summation_type : str
        ``"PLAN"``, ``"BEAM"``, ``"FRACTION"`` and so on. Several values make Referenced RT
        Plan Sequence required; see :func:`check_dose_plan_reference`.
    dose_units : str
        ``"GY"`` or ``"RELATIVE"``.

    Raises
    ------
    ValueError
        If the grid holds negative dose, or is not 3D.
    """
    array = sitk.GetArrayFromImage(dose_sitk_image).astype(np.float64)
    if dose_sitk_image.GetDimension() != 3:
        raise ValueError(
            f"dose grid must be 3D, got {dose_sitk_image.GetDimension()}D"
        )
    if array.size and array.min() < 0:
        raise ValueError(
            f"dose grid holds negative values (min {array.min():g}); stored dose is "
            "unsigned and negatives would wrap around to very large doses"
        )

    max_dose = float(array.max()) if array.size else 0.0
    if dose_grid_scaling is None:
        dose_grid_scaling = DEFAULT_DOSE_GRID_SCALING
        # Beyond ~429 Gy the default no longer fits; scale up rather than wrap silently.
        if max_dose / dose_grid_scaling > _STORED_MAX:
            dose_grid_scaling = max_dose / _STORED_MAX
    elif dose_grid_scaling <= 0:
        raise ValueError(f"dose_grid_scaling must be positive, got {dose_grid_scaling}")

    stored = np.rint(array / dose_grid_scaling)
    if stored.size and stored.max() > _STORED_MAX:
        raise ValueError(
            f"dose up to {max_dose:g} does not fit at dose_grid_scaling={dose_grid_scaling:g}; "
            f"use at least {max_dose / _STORED_MAX:g}"
        )

    size = list(dose_sitk_image.GetSize())
    spacing = list(dose_sitk_image.GetSpacing())
    direction = list(dose_sitk_image.GetDirection())
    origin = list(dose_sitk_image.GetOrigin())

    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 32
    ds.BitsStored = 32
    ds.HighBit = 31
    # The stored values are unsigned; declaring them signed made every dose above
    # ~215 Gy read back negative.
    ds.PixelRepresentation = 0
    ds.DoseUnits = dose_units
    ds.DoseType = dose_type
    ds.DoseSummationType = dose_summation_type
    ds.NumberOfFrames = size[2]
    ds.FrameIncrementPointer = Tag(0x3004, 0x000C)

    ds.ImagePositionPatient = origin
    # ImageOrientationPatient is VM 6 (row direction + column direction) and PixelSpacing is
    # ordered (row, column); neither matches SimpleITK's layout directly.
    ds.ImageOrientationPatient = sitk_direction_to_image_orientation_patient(direction)
    ds.PixelSpacing = sitk_spacing_to_pixel_spacing(spacing)
    ds.SliceThickness = spacing[2]
    ds.Columns = size[0]
    ds.Rows = size[1]
    ds.GridFrameOffsetVector = [ov * spacing[2] for ov in range(size[2])]
    ds.DoseGridScaling = dose_grid_scaling
    ds.PixelData = stored.astype(np.uint32).tobytes()
    return ds


def check_dose_plan_reference(ds: Dataset) -> list:
    """
    Check the Type 1C condition on Referenced RT Plan Sequence.

    Most Dose Summation Types -- PLAN among them -- require a plan reference. The builder
    creates the sequence empty, so a dose written without calling
    :func:`add_rf_rt_plan_seq_from_plan_ds` claims a summation it cannot substantiate.

    Returns
    -------
    list
        Problem descriptions; empty when the dataset is consistent.
    """
    summation_type = getattr(ds, "DoseSummationType", None)
    if summation_type is None:
        return ["Missing DoseSummationType"]

    if summation_type not in SUMMATION_TYPES_REQUIRING_PLAN:
        return []

    plan_sequence = getattr(ds, "ReferencedRTPlanSequence", None)
    if not plan_sequence:
        return [
            f"DoseSummationType is {summation_type}, which makes ReferencedRTPlanSequence "
            "Type 1C, but it is absent or empty"
        ]
    return []


class RTDoseBuilder:
    """
    Build an RT Dose object from a dose grid.

    Patient, study and frame-of-reference context is taken from a reference dataset -- a
    slice of the planning CT, or the plan -- so the dose files into the same study.

    Parameters
    ----------
    reference_ds : Dataset
        A single dataset, not a series: one planning CT slice, or the RT Plan. Supplies
        patient identity, study and frame of reference.

    Methods
    -------
    set_dose_grid(dose_image, ...)
        The dose values and their geometry. Required.
    add_referenced_plan(plan_ds)
        Reference the RT Plan. Required for most Dose Summation Types.
    set_series_number(series_number)
        Set SeriesNumber.
    set_uid_prefix(uid_prefix)
        Set the UID root.
    build()
        Assemble the ``FileDataset``.

    Raises
    ------
    ValueError
        If ``reference_ds`` is None.

    See Also
    --------
    get_dose_image : Read the result back.
    check_rtdose_iod : Validate, including the conditional plan reference.

    Notes
    -----
    Dose is stored as scaled unsigned 32-bit integers. Negative values are rejected --
    unsigned storage would turn them into very large positive ones -- so clip or offset
    beforehand if your grid dips below zero numerically.

    Examples
    --------
    >>> dose_ds = (                                                 # doctest: +SKIP
    ...     RTDoseBuilder(planning_ct_series[0])
    ...     .set_dose_grid(dose_image)
    ...     .add_referenced_plan(plan_ds)
    ...     .build()
    ... )
    >>> dose_ds.save_as("rtdose.dcm", enforce_file_format=True)     # doctest: +SKIP
    """

    def __init__(self, reference_ds: Dataset):
        if reference_ds is None:
            raise ValueError("reference_ds must be a Dataset")
        self.reference_ds = reference_ds
        self.uid_prefix = DICOM_UID_PREFIX
        self.series_number = 1
        self._dose_image = None
        self._grid_kwargs = {}
        self._plan_datasets = []

    def set_uid_prefix(self, uid_prefix: str) -> "RTDoseBuilder":
        """
        Set the root under which generated UIDs are minted.

        Parameters
        ----------
        uid_prefix : str
            Organisation UID root, dot-terminated. Defaults to the ``DICOM_UID_PREFIX``
            environment variable.

        Returns
        -------
        RTDoseBuilder
            self, so calls chain.
        """
        self.uid_prefix = uid_prefix
        return self

    def set_series_number(self, series_number: int) -> "RTDoseBuilder":
        """
        Set SeriesNumber (0020,0011).

        Parameters
        ----------
        series_number : int
            Default 1. Give each dose its own number when writing several into one study
            -- beam doses alongside a plan dose, say -- so a viewer can tell them apart.

        Returns
        -------
        RTDoseBuilder
            self, so calls chain.
        """
        self.series_number = series_number
        return self

    def set_dose_grid(
        self,
        dose_image: sitk.Image,
        dose_grid_scaling: float = None,
        dose_type: str = "PHYSICAL",
        dose_summation_type: str = "PLAN",
        dose_units: str = "GY",
        ) -> "RTDoseBuilder":
        """
        Set the dose grid and how it is described.

        Parameters
        ----------
        dose_image : sitk.Image
            Dose values in ``dose_units``. Geometry -- origin, spacing, direction --
            becomes the grid's; it need not match the CT grid, and typically is coarser.
        dose_grid_scaling : float, optional
            Factor relating stored integers to real dose. Default 1e-7, covering about
            429 Gy. Raised automatically if the grid exceeds what the current factor can
            represent, so the default rarely needs changing; lower it for finer
            quantisation over a smaller range.
        dose_type : str, optional
            DoseType: ``"PHYSICAL"`` (default), ``"EFFECTIVE"`` or ``"ERROR"``.
        dose_summation_type : str, optional
            DoseSummationType: ``"PLAN"`` (default), ``"FRACTION"``, ``"BEAM"``,
            ``"BRACHY"`` and so on. Everything except ``"FRACTION"`` makes Referenced RT
            Plan Sequence mandatory -- add one with :meth:`add_referenced_plan`.
        dose_units : str, optional
            DoseUnits: ``"GY"`` (default) or ``"RELATIVE"``.

        Returns
        -------
        RTDoseBuilder
            self, so calls chain.

        See Also
        --------
        add_dose_grid_to_ds : The underlying function.
        add_referenced_plan : Required for most summation types.

        Notes
        -----
        Calling this twice replaces the grid rather than accumulating; the last call wins.
        """
        self._dose_image = dose_image
        self._grid_kwargs = dict(
            dose_grid_scaling=dose_grid_scaling,
            dose_type=dose_type,
            dose_summation_type=dose_summation_type,
            dose_units=dose_units,
        )
        return self

    def add_referenced_plan(self, plan_ds: Dataset) -> "RTDoseBuilder":
        """
        Reference the RT Plan this dose belongs to.

        Parameters
        ----------
        plan_ds : Dataset
            An RT Plan (SOP class ``1.2.840.10008.5.1.4.1.1.481.5``). Only its SOP Class
            and Instance UIDs are stored.

        Returns
        -------
        RTDoseBuilder
            self, so calls chain.

        Notes
        -----
        Referenced RT Plan Sequence is Type 1C, required for every Dose Summation Type
        except ``"FRACTION"`` -- including the ``"PLAN"`` default. Omitting it produces a
        file :func:`check_rtdose_iod` reports and most TPSs refuse.

        Can be called more than once; each plan is appended.

        Raises
        ------
        ValueError
            If ``plan_ds`` is not an RT Plan.
        """
        self._plan_datasets.append(plan_ds)
        return self

    def build(self) -> FileDataset:
        """
        Assemble the RT Dose object.

        Returns
        -------
        FileDataset
            Ready to ``save_as(path, enforce_file_format=True)``.

        Raises
        ------
        ValueError
            If no dose grid has been set, or a referenced plan is not an RT Plan.
        """
        if self._dose_image is None:
            raise ValueError("no dose grid set; call set_dose_grid() before build()")

        ds = generate_base_dataset(self.uid_prefix)
        cp_information_from_ds(
            ds, self.reference_ds,
            series_number=self.series_number, uid_prefix=self.uid_prefix,
        )
        add_dose_grid_to_ds(ds, self._dose_image, **self._grid_kwargs)
        for plan_ds in self._plan_datasets:
            add_rf_rt_plan_seq_from_plan_ds(ds, plan_ds)
        return ds
