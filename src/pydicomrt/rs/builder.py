"""
RT Structure Set construction.

:class:`RTStructBuilder` is the entry point. The module-level functions below are the
building blocks it is assembled from; use them directly only if you need to deviate from
what the builder does.
"""

import datetime
import os
from typing import Dict, List, Optional, Sequence as TypingSequence

import numpy as np
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ImplicitVRLittleEndian, generate_uid

from .add_new_roi import create_roi_into_rs_ds
from .make_contour_sequence import (
    add_contour_sequence_from_dcm_ctr_dict,
    add_contour_sequence_from_mask3d,
)


DICOM_UID_PREFIX = os.getenv('DICOM_UID_PREFIX', "1.2.826.0.1.3680043.8.498.")  # PYDICOM_ROOT_UID

def create_rtstruct_dataset(image_series, uid_prefix: str = None) -> FileDataset:
    """
    Create an empty RT Structure Set bound to an image series.

    Everything except the ROIs: file meta, patient and study context copied from the
    images, the referenced frame-of-reference sequence listing every slice, and the three
    empty ROI sequences ready to be filled.

    Called by :class:`RTStructBuilder`. Use it directly when you need to interleave your
    own dataset edits with ROI creation.

    Parameters
    ----------
    image_series : list of Dataset
        Slices sorted along the slice normal. Patient identity, study, frame of reference
        and the per-slice references all come from here.
    uid_prefix : str, optional
        UID root. Defaults to the ``DICOM_UID_PREFIX`` environment variable.

    Returns
    -------
    FileDataset
        A conformant but ROI-less RT Structure Set. It will not pass
        :func:`check_rtstruct_iod` until at least one ROI is added.

    See Also
    --------
    RTStructBuilder : The supported entry point.
    create_roi_into_rs_ds : Add an ROI to the result.
    """
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    ds = generate_base_dataset(uid_prefix)
    add_study_and_series_information(ds, image_series, uid_prefix)
    add_patient_information(ds, image_series)
    add_refd_frame_of_ref_sequence(ds, image_series, uid_prefix)
    add_rs_series_information(ds, image_series)
    return ds

def generate_base_dataset(uid_prefix: str = None) -> FileDataset:
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    file_name = "rs_test"
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
    file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3"
    file_meta.MediaStorageSOPInstanceUID = (
        generate_uid(prefix=uid_prefix)
    )
    file_meta.ImplementationClassUID = uid_prefix + "1"
    file_meta.ImplementationVersionName = "pydicomRT"
    return file_meta

def add_required_elements_to_ds(ds: FileDataset):
    dt = datetime.datetime.now()
    # Append data elements required by the DICOM standarad
    ds.SpecificCharacterSet = "ISO_IR 192"
    ds.InstanceCreationDate = dt.strftime("%Y%m%d")
    ds.InstanceCreationTime = dt.strftime("%H%M%S")
    ds.StructureSetLabel = "RS_" + dt.strftime("%Y%m%d")
    ds.StructureSetDate = dt.strftime("%Y%m%d")
    ds.StructureSetTime = dt.strftime("%H%M%S.%f")
    ds.Modality = "RTSTRUCT"
    ds.Manufacturer = "pydicomRT"
    ds.ManufacturerModelName = "modelv1"
    ds.InstitutionName = ""
    # Structure Set module, Type 2 -- was missing entirely.
    ds.InstanceNumber = 1
    # Set values already defined in the file meta
    ds.SOPClassUID = ds.file_meta.MediaStorageSOPClassUID
    ds.SOPInstanceUID = ds.file_meta.MediaStorageSOPInstanceUID

    # Approval module: Type 2C fields are only required once the status is APPROVED.
    ds.ApprovalStatus = "UNAPPROVED"

def add_sequence_lists_to_ds(ds: FileDataset):
    ds.StructureSetROISequence = Sequence()
    ds.ROIContourSequence = Sequence()
    ds.RTROIObservationsSequence = Sequence()

# ---------------------------------------------------------------------------------------- #

def add_study_and_series_information(ds: FileDataset, image_series, uid_prefix: str = None):
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    reference_ds = image_series[0]  # All elements in series should have the same data
    ds.StudyDate = reference_ds.StudyDate
    ds.SeriesDate = getattr(reference_ds, "SeriesDate", "")
    ds.StudyTime = reference_ds.StudyTime
    ds.SeriesTime = getattr(reference_ds, "SeriesTime", "")
    ds.StudyDescription = getattr(reference_ds, "StudyDescription", "")
    ds.SeriesDescription = getattr(reference_ds, "SeriesDescription", "")
    ds.StudyInstanceUID = reference_ds.StudyInstanceUID
    ds.SeriesInstanceUID = generate_uid(prefix=uid_prefix)  # TODO: find out if random generation is ok
    ds.StudyID = reference_ds.StudyID
    ds.AccessionNumber = getattr(reference_ds, "AccessionNumber", "")
    ds.SeriesNumber = "1"  # TODO: find out if we can just use 1 (Should be fine since its a new series)
    ds.OperatorsName = getattr(reference_ds, "OperatorsName", "")
    ds.ReferringPhysicianName = getattr(reference_ds, "ReferringPhysicianName", "")
    ds.PhysiciansOfRecord = getattr(reference_ds, "PhysiciansOfRecord", "")

def add_patient_information(ds: FileDataset, image_series):
    reference_ds = image_series[0]  # All elements in series should have the same data
    ds.PatientName = getattr(reference_ds, "PatientName", "")
    ds.PatientID = getattr(reference_ds, "PatientID", "")
    ds.PatientBirthDate = getattr(reference_ds, "PatientBirthDate", "")
    ds.PatientSex = getattr(reference_ds, "PatientSex", "")
    ds.PatientAge = getattr(reference_ds, "PatientAge", "")
    ds.PatientSize = getattr(reference_ds, "PatientSize", "")
    ds.PatientWeight = getattr(reference_ds, "PatientWeight", "")

def add_refd_frame_of_ref_sequence(ds: FileDataset, image_series, uid_prefix: str = None):
    uid_prefix = uid_prefix or DICOM_UID_PREFIX
    refd_frame_of_ref = Dataset()
    refd_frame_of_ref.FrameOfReferenceUID = getattr(image_series[0], 'FrameOfReferenceUID', generate_uid(prefix=uid_prefix))
    refd_frame_of_ref.RTReferencedStudySequence = create_frame_of_ref_study_sequence(image_series)

    ds.ReferencedFrameOfReferenceSequence = Sequence()
    ds.ReferencedFrameOfReferenceSequence.append(refd_frame_of_ref)
    ds.FrameOfReferenceUID = getattr(image_series[0], 'FrameOfReferenceUID', generate_uid(prefix=uid_prefix))

def create_frame_of_ref_study_sequence(image_series) -> Sequence:
    reference_ds = image_series[0]
    rt_refd_series = Dataset()
    rt_refd_series.SeriesInstanceUID = reference_ds.SeriesInstanceUID
    rt_refd_series.ContourImageSequence = create_contour_image_sequence(image_series)
    rt_refd_series_sequence = Sequence()
    rt_refd_series_sequence.append(rt_refd_series)
    rt_refd_study = Dataset()
    rt_refd_study.ReferencedSOPClassUID = "1.2.840.10008.3.1.2.3.1"
    rt_refd_study.ReferencedSOPInstanceUID = reference_ds.StudyInstanceUID
    rt_refd_study.RTReferencedSeriesSequence = rt_refd_series_sequence
    rt_refd_study_sequence = Sequence()
    rt_refd_study_sequence.append(rt_refd_study)
    return rt_refd_study_sequence

def create_contour_image_sequence(image_series) -> Sequence:
    contour_image_sequence = Sequence()
    for series in image_series:
        contour_image = Dataset()
        contour_image.ReferencedSOPClassUID = series.SOPClassUID
        contour_image.ReferencedSOPInstanceUID = series.SOPInstanceUID
        contour_image_sequence.append(contour_image)
    return contour_image_sequence

def add_rs_series_information(ds: FileDataset, image_series):
    dt = datetime.datetime.now()
    ds.StructureSetDescription = dt.strftime("%Y%m%d") + "MBB"


DEFAULT_CONTOUR_CONFIG = {
    "ex_noise_size": 10,
    "in_noise_size": 10,
    "lowpass_ratio": 10,
    "ctr_precision": 8,
}


class RTStructBuilder:
    """
    Build an RT Structure Set from masks or contours.

    Patient, study and frame-of-reference context is taken from the image series the
    structures were drawn on, so the output files into the same study as those images.

    Parameters
    ----------
    image_series : list of Dataset
        The slices the structures were drawn on, sorted along the slice normal. Every mask
        added must match this series in shape, and the contours are placed using its
        geometry.

    Methods
    -------
    add_roi(mask, name, ...)
        Add an ROI from a 3D binary mask.
    add_roi_from_contours(contours, name, ...)
        Add an ROI from points already in patient coordinates.
    set_structure_set_label(label)
        Override the generated StructureSetLabel.
    set_uid_prefix(uid_prefix)
        Override the UID root.
    build()
        Assemble the ``FileDataset``.

    Raises
    ------
    ValueError
        If ``image_series`` is empty.

    See Also
    --------
    rtstruct_to_masks : The reverse direction, RTSTRUCT to voxels.
    check_rtstruct_iod : Validate the result before writing it.

    Notes
    -----
    The UID root defaults to the ``DICOM_UID_PREFIX`` environment variable. Set it, or
    call :meth:`set_uid_prefix`, before writing files that leave your machine.

    Examples
    --------
    >>> rs_ds = (                                                   # doctest: +SKIP
    ...     RTStructBuilder(image_series)
    ...     .add_roi(mask=ctv_mask, name="CTV", color=[0, 255, 0], interpreted_type="CTV")
    ...     .add_roi(mask=cord_mask, name="SpinalCord", interpreted_type="ORGAN")
    ...     .build()
    ... )
    >>> rs_ds.save_as("rtstruct.dcm", enforce_file_format=True)     # doctest: +SKIP
    """

    def __init__(self, image_series: List[Dataset]):
        if not image_series:
            raise ValueError("image_series must be non-empty")
        self.image_series = list(image_series)
        self.uid_prefix = DICOM_UID_PREFIX
        self.structure_set_label = None
        self._roi_specs: List[dict] = []

    def set_uid_prefix(self, uid_prefix: str) -> "RTStructBuilder":
        """
        Set the root under which generated UIDs are minted.

        Parameters
        ----------
        uid_prefix : str
            Organisation UID root, dot-terminated, e.g. ``"1.2.826.0.1.3680043.2.1125."``.
            Defaults to the ``DICOM_UID_PREFIX`` environment variable.

        Returns
        -------
        RTStructBuilder
            self, so calls chain.

        Notes
        -----
        UIDs generated under a root you do not own can collide with another
        organisation's. That matters for anything leaving your machine and not at all for
        scratch files.
        """
        self.uid_prefix = uid_prefix
        return self

    def set_structure_set_label(self, label: str) -> "RTStructBuilder":
        """
        Override the generated StructureSetLabel.

        Parameters
        ----------
        label : str
            VR SH, so 16 characters. Longer labels are what TPSs display, and a
            non-conformant length is a common reason for an import to be refused.

        Returns
        -------
        RTStructBuilder
            self, so calls chain.
        """
        self.structure_set_label = label
        return self

    def _next_roi_number(self) -> int:
        return max((spec["number"] for spec in self._roi_specs), default=0) + 1

    def _add(self, source: str, payload, name, number, color, description,
             interpreted_type, contour_config) -> "RTStructBuilder":
        if not name:
            raise ValueError("an ROI needs a name")

        number = self._next_roi_number() if number is None else int(number)
        if any(spec["number"] == number for spec in self._roi_specs):
            raise ValueError(f"roi_number {number} was already added")

        self._roi_specs.append({
            "source": source,
            "payload": payload,
            "name": name,
            "number": number,
            "color": list(color),
            "description": description,
            "interpreted_type": interpreted_type,
            "contour_config": contour_config or DEFAULT_CONTOUR_CONFIG,
        })
        return self

    def add_roi(
        self,
        mask: np.ndarray,
        name: str,
        number: Optional[int] = None,
        color: TypingSequence[int] = (255, 0, 0),
        description: str = "",
        interpreted_type: str = "ORGAN",
        contour_config: Optional[Dict] = None,
        ) -> "RTStructBuilder":
        """
        Add an ROI from a 3D mask.

        Parameters
        ----------
        mask : np.ndarray
            ``(slice, row, column)`` -- the same axis order as the image series stacks to.
            Non-zero voxels are inside the ROI.
        name : str
            ROIName, as it appears in the TPS.
        number : int, optional
            ROINumber. Assigned in sequence when omitted.
        color : sequence[int]
            ROIDisplayColor, RGB.
        description : str
            ROIDescription.
        interpreted_type : str
            RTROIInterpretedType, e.g. ``"ORGAN"``, ``"PTV"``, ``"EXTERNAL"``.
        contour_config : dict, optional
            Contour extraction tuning; see :data:`DEFAULT_CONTOUR_CONFIG`.

        Returns
        -------
        RTStructBuilder
            self, so calls chain.
        """
        expected = (len(self.image_series), self.image_series[0].Rows,
                    self.image_series[0].Columns)
        if tuple(np.shape(mask)) != expected:
            raise ValueError(
                f"mask shape {tuple(np.shape(mask))} does not match the image series "
                f"{expected} (slice, row, column)"
            )
        return self._add("mask", mask, name, number, color, description,
                         interpreted_type, contour_config)

    def add_roi_from_contours(
        self,
        contours: Dict,
        name: str,
        number: Optional[int] = None,
        color: TypingSequence[int] = (255, 0, 0),
        description: str = "",
        interpreted_type: str = "ORGAN",
        ) -> "RTStructBuilder":
        """
        Add an ROI from contour points already in patient coordinates.

        Parameters
        ----------
        contours : dict
            Keyed by the SOP Instance UID of the image each contour lies on::

                {sop_instance_uid: {"sop_class_uid": str,
                                    "contours": [[x1, y1, z1, x2, y2, z2, ...], ...]}}

            Points are patient coordinates in mm, interleaved into one flat list per
            contour -- the DICOM ContourData layout, and what :func:`get_contours` returns
            under ``"dcm_contour"``. The UIDs must belong to the builder's image series.
        name : str
            ROIName, as it appears in the TPS.
        number : int, optional
            ROINumber. Assigned in sequence when omitted.
        color : sequence of int
            ROIDisplayColor, RGB 0-255. Default ``(255, 0, 0)``.
        description : str
            ROIDescription. Default empty.
        interpreted_type : str
            RTROIInterpretedType, e.g. ``"ORGAN"``, ``"PTV"``, ``"EXTERNAL"``. Default
            ``"ORGAN"``.

        Returns
        -------
        RTStructBuilder
            self, so calls chain.

        Raises
        ------
        ValueError
            If ``name`` is empty, or ``number`` was already added.

        See Also
        --------
        add_roi : When you have a mask rather than points.
        get_contours : Produces this structure from an existing RTSTRUCT.

        Notes
        -----
        The points are written as given -- no contour filtering, no resampling, no
        planarity check. This is the path to use when contours must round-trip unchanged.
        """
        return self._add("contours", contours, name, number, color, description,
                         interpreted_type, None)

    def build(self) -> FileDataset:
        """
        Assemble the RT Structure Set.

        Contours are extracted here, not when ``add_roi`` is called, so the cost of
        ``cv2.findContours`` over every mask lands on this call.

        Returns
        -------
        FileDataset
            Ready to ``save_as(path, enforce_file_format=True)``.

        Raises
        ------
        ValueError
            If no ROI has been added.

        Notes
        -----
        The builder is not consumed. Calling ``build()`` twice produces two datasets with
        *different* SOP Instance UIDs and the same content, which is usually not what you
        want -- keep the first result rather than rebuilding.
        """
        if not self._roi_specs:
            raise ValueError("no ROI added; call add_roi() before build()")

        ds = create_rtstruct_dataset(self.image_series, uid_prefix=self.uid_prefix)
        if self.structure_set_label is not None:
            ds.StructureSetLabel = self.structure_set_label

        for spec in self._roi_specs:
            create_roi_into_rs_ds(
                ds,
                spec["color"],
                spec["number"],
                spec["name"],
                spec["description"],
                spec["interpreted_type"],
            )
            if spec["source"] == "mask":
                add_contour_sequence_from_mask3d(
                    ds, self.image_series, spec["number"], spec["payload"],
                    spec["contour_config"],
                )
            else:
                add_contour_sequence_from_dcm_ctr_dict(
                    ds, self.image_series, spec["number"], spec["payload"],
                )
        return ds
