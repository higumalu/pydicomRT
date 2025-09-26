from pydicom.dataset import Dataset
from pydicom.uid import generate_uid
import SimpleITK as sitk
import numpy as np
from pydicomrt.reg.builder import DeformableSpatialRegistrationBuilder

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
    return ds

def make_fake_moving_instance(frame_of_ref_uid: str) -> Dataset:
    ds = Dataset()
    ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"  # CT Image Storage
    ds.SOPInstanceUID = generate_uid()
    ds.FrameOfReferenceUID = frame_of_ref_uid
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

    tfm = sitk.DisplacementFieldTransform(vec)

    tfm.GetOrigin = lambda: vec.GetOrigin()
    tfm.GetSpacing = lambda: vec.GetSpacing()
    tfm.GetSize = lambda: vec.GetSize()
    tfm.GetDirection = lambda: vec.GetDirection()

    return tfm

def identity_4x4_flat_f32():
    M = np.eye(4, dtype=np.float32)
    return [float(x) for x in M.reshape(-1)]

# -----------------------------------------------------------
# Demo
# -----------------------------------------------------------
def demo():
    fixed_ref = make_fake_fixed_ref_ds()

    moving_list = [
        make_fake_moving_instance(fixed_ref.FrameOfReferenceUID),
        make_fake_moving_instance(fixed_ref.FrameOfReferenceUID),
    ]

    disp_tfm = make_displacement_field_transform(
        size=(8, 9, 10),
        spacing=(1.5, 1.5, 2.0),
        origin=(0.0, 0.0, 0.0),
        direction=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        shift=(1.0, 0.0, 0.0)
    )
    preM = identity_4x4_flat_f32()
    postM = identity_4x4_flat_f32()

    # 4) Builder: Add a deformable registration block and build the REG dataset
    builder = DeformableSpatialRegistrationBuilder([fixed_ref])
    builder.set_uid_prefix("1.2.826.0.1.3680043.2.1125.")
    builder.add_deformable_registration(
        moving_ds_list=moving_list,
        vectorial_field_transform=disp_tfm,
        pre_transform=preM,
        post_transform=postM
    )
    reg_ds = builder.build()

    # 5) Verify the output
    print("=== REG Dataset Summary ===")
    print("SOPClassUID: ", reg_ds.SOPClassUID)
    print("SOPInstanceUID:", reg_ds.SOPInstanceUID)
    print("Modality   : ", reg_ds.Modality)
    print("SeriesDesc : ", getattr(reg_ds, "SeriesDescription", ""))
    print("FoR UID    : ", getattr(reg_ds, "FrameOfReferenceUID", ""))

    seq = reg_ds.DeformableRegistrationSequence
    print("\nDeformableRegistrationSequence items:", len(seq))
    item = seq[0]
    print("  SourceFoR UID:", item.SourceFrameOfReferenceUID)
    print("  #ReferencedImages:", len(item.ReferencedImageSequence))

    grid = item.DeformableRegistrationGridSequence[0]
    print("  GridDimensions : ", grid.GridDimensions)
    print("  GridResolution : ", grid.GridResolution)
    print("  IOP            : ", list(grid.ImageOrientationPatient))
    print("  IPP            : ", list(grid.ImagePositionPatient))
    print("  VectorGridData bytes:", len(grid.VectorGridData))
    code = item.RegistrationTypeCodeSequence[0]
    print("  Code (Value,Meaning):", code.CodeValue, code.CodeMeaning)

    # 6)（Optional）：
    # from pydicom.filewriter import write_file
    # write_file("REG_Deformable.dcm", reg_ds)
    # print("\nSaved to REG_Deformable.dcm")


if __name__ == "__main__":
    demo()
