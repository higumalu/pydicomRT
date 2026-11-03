"""
demons_registration(backend="torch"): optional, and the same registration as SimpleITK.

PyTorch is an optional extra, so the first group of tests runs without it: the default
backend must never import torch, and asking for the torch backend without it must say how to
install it. The parity tests need torch and run on the CPU everywhere, plus CUDA when present.

Parity is asserted on the mapping -- where the transform sends probe points -- and on the
field, against the SimpleITK backend on anisotropic, rotated and mismatched grids. The
synthetic pairs are small enough that demons is well conditioned throughout; on large
volumes the two backends agree to the same degree SimpleITK agrees with itself under input
noise of a thousandth of a HU, which near truncated image edges can be millimetres.

metric="mind" has no SimpleITK counterpart, so it is held to what it is for: recovering a
known deformation between images of different "modalities" -- the moving phantom passed
through a non-monotonic intensity map, as CT -> MR is -- inventing no motion between two
such images that are already aligned, and not caring about gain or offset.
"""

import subprocess
import sys

import numpy as np
import pytest
import SimpleITK as sitk

from synthetic import make_textured_phantom
from pydicomrt.reg.method import demons_registration


def rotation(about_z=0.3, about_y=0.4):
    cz, sz = np.cos(about_z), np.sin(about_z)
    cy, sy = np.cos(about_y), np.sin(about_y)
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    return tuple((rz @ ry).ravel())


def deformed_pair(spacing=(1.0, 1.5, 2.5), direction=None, amplitude=4.0):
    """A phantom and a copy warped by a known B-spline: fixed -> moving is `truth`."""
    moving = make_textured_phantom(size=(40, 36, 28), spacing=spacing)
    if direction is not None:
        moving.SetDirection(direction)
    rng = np.random.default_rng(0)
    truth = sitk.BSplineTransformInitializer(moving, [3, 3, 3])
    params = np.array(truth.GetParameters())
    truth.SetParameters((params + rng.uniform(-amplitude, amplitude, params.size)).tolist())
    fixed = sitk.Resample(moving, moving, truth, sitk.sitkLinear, -1000.0)
    return fixed, moving, truth


def probe_points(image):
    size = image.GetSize()
    return [
        image.TransformIndexToPhysicalPoint((int(a * size[0]), int(b * size[1]), int(c * size[2])))
        for a in (0.3, 0.5, 0.7) for b in (0.3, 0.5, 0.7) for c in (0.3, 0.5, 0.7)
    ]


def assert_same_registration(reference, candidate, fixed):
    for point in probe_points(fixed):
        assert candidate[1].TransformPoint(point) == pytest.approx(reference[1].TransformPoint(point), abs=1e-3)
    difference = np.linalg.norm(
        sitk.GetArrayFromImage(candidate[2]) - sitk.GetArrayFromImage(reference[2]), axis=-1
    )
    assert difference.max() < 0.05
    assert candidate[2].GetPixelID() == reference[2].GetPixelID() == sitk.sitkVectorFloat64
    for a, b in ((candidate[0], reference[0]), (candidate[2], reference[2])):
        assert a.GetSize() == b.GetSize()
        assert a.GetSpacing() == pytest.approx(b.GetSpacing())
        assert a.GetOrigin() == pytest.approx(b.GetOrigin())
        assert a.GetDirection() == pytest.approx(b.GetDirection())


# ---------------------------------------------------------------------------------------
# Without torch
# ---------------------------------------------------------------------------------------

def test_default_backend_never_imports_torch():
    """Importing the library and running a default registration must not pull in torch."""
    code = (
        "import sys, SimpleITK as sitk\n"
        "import pydicomrt.reg, pydicomrt.reg.method, pydicomrt.reg.pipeline\n"
        "import pydicomrt.rs, pydicomrt.dose, pydicomrt.ct, pydicomrt.utils\n"
        "from pydicomrt.reg.method import demons_registration\n"
        "image = sitk.Image([16] * 3, sitk.sitkFloat32) + 1\n"
        "demons_registration(image, image, resolution_staging=(1,), iteration_staging=(1,),"
        " smoothing_sigmas=0)\n"
        "print('torch' in sys.modules)\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stdout.strip() == "False"


def test_missing_torch_says_how_to_install_it(monkeypatch):
    monkeypatch.setitem(sys.modules, "torch", None)  # makes `import torch` raise ImportError
    monkeypatch.delitem(sys.modules, "pydicomrt.reg.method._demons_torch", raising=False)
    image = sitk.Image([16] * 3, sitk.sitkFloat32) + 1
    with pytest.raises(ImportError, match=r"pydicomrt\[gpu\]"):
        demons_registration(image, image, resolution_staging=(1,), iteration_staging=(1,), backend="torch")


def test_unknown_backend_is_rejected():
    image = sitk.Image([16] * 3, sitk.sitkFloat32)
    with pytest.raises(ValueError, match="backend"):
        demons_registration(image, image, backend="cuda")


def test_device_belongs_to_the_torch_backend():
    image = sitk.Image([16] * 3, sitk.sitkFloat32)
    with pytest.raises(ValueError, match="device"):
        demons_registration(image, image, device="cuda")


def test_torch_backend_interpolates_linearly_only():
    image = sitk.Image([16] * 3, sitk.sitkFloat32)
    with pytest.raises(ValueError, match="interp_order"):
        demons_registration(image, image, backend="torch", interp_order=sitk.sitkBSpline)


def test_unknown_metric_is_rejected():
    image = sitk.Image([16] * 3, sitk.sitkFloat32)
    with pytest.raises(ValueError, match="metric"):
        demons_registration(image, image, metric="mutual_information")


def test_mind_metric_needs_the_torch_backend():
    image = sitk.Image([16] * 3, sitk.sitkFloat32)
    with pytest.raises(ValueError, match="backend='torch'"):
        demons_registration(image, image, metric="mind")


def test_intensity_defaults_are_unchanged():
    """None-valued defaults resolve per metric; for intensity they are the old literals."""
    fixed, moving, _ = deformed_pair()
    implicit = demons_registration(fixed, moving)[2]
    explicit = demons_registration(fixed, moving, resolution_staging=(8, 4, 1), iteration_staging=(10, 10, 10),
                                   regularization_kernel_mm=1.5, smoothing_sigma_factor=1.0)[2]
    assert np.array_equal(sitk.GetArrayFromImage(implicit), sitk.GetArrayFromImage(explicit))


def test_torch_backend_is_3d_only():
    image = sitk.Image([16, 16], sitk.sitkFloat32)
    with pytest.raises(ValueError, match="3D"):
        demons_registration(image, image, backend="torch")


# ---------------------------------------------------------------------------------------
# With torch
# ---------------------------------------------------------------------------------------

def _devices():
    torch = pytest.importorskip("torch")
    devices = ["cpu"]
    if torch.cuda.is_available():
        devices.append("cuda")
    return devices


@pytest.fixture(params=["cpu", "cuda"])
def device(request):
    if request.param not in _devices():
        pytest.skip("no CUDA device")
    return request.param


def run_both(fixed, moving, device, **kwargs):
    reference = demons_registration(fixed, moving, default_value=-1000.0, **kwargs)
    candidate = demons_registration(fixed, moving, default_value=-1000.0, backend="torch",
                                    device=device, **kwargs)
    return reference, candidate


def test_matches_simpleitk_on_anisotropic_grid(device):
    fixed, moving, _ = deformed_pair()
    assert_same_registration(*run_both(fixed, moving, device), fixed)


def test_matches_simpleitk_on_rotated_grid(device):
    """A direction that mixes all three axes: row-major and column readings differ."""
    fixed, moving, _ = deformed_pair(direction=rotation())
    assert_same_registration(*run_both(fixed, moving, device), fixed)


def test_matches_simpleitk_with_moving_on_another_grid(device):
    fixed, moving, _ = deformed_pair(direction=rotation())
    regridded = sitk.Resample(
        moving, [44, 40, 30], sitk.Transform(), sitk.sitkLinear,
        [o - 3.0 for o in moving.GetOrigin()], [0.9, 1.3, 2.2], (1, 0, 0, 0, 1, 0, 0, 0, 1), -1000.0,
    )
    assert_same_registration(*run_both(fixed, regridded, device), fixed)


def test_matches_simpleitk_with_isotropic_staging(device):
    fixed, moving, _ = deformed_pair()
    assert_same_registration(
        *run_both(fixed, moving, device, isotropic_resample=True, resolution_staging=(6.0, 3.0, 2.0)),
        fixed,
    )


def test_matches_simpleitk_from_an_initial_field(device):
    fixed, moving, _ = deformed_pair(direction=rotation())
    initial = sitk.TransformToDisplacementField(
        sitk.TranslationTransform(3, (0.5, -0.3, 0.2)), sitk.sitkVectorFloat64,
        fixed.GetSize(), fixed.GetOrigin(), fixed.GetSpacing(), fixed.GetDirection(),
    )
    assert_same_registration(*run_both(fixed, moving, device, initial_displacement_field=initial), fixed)


def test_recovers_a_known_deformation(device):
    fixed, moving, truth = deformed_pair(direction=rotation())
    _, transform, _ = demons_registration(fixed, moving, default_value=-1000.0, backend="torch", device=device)
    points = probe_points(fixed)
    before = np.mean([np.linalg.norm(np.subtract(truth.TransformPoint(p), p)) for p in points])
    after = np.mean([np.linalg.norm(np.subtract(truth.TransformPoint(p), transform.TransformPoint(p)))
                     for p in points])
    assert after < before / 2


def test_repeated_runs_are_identical(device):
    fixed, moving, _ = deformed_pair()
    first = demons_registration(fixed, moving, backend="torch", device=device)[2]
    second = demons_registration(fixed, moving, backend="torch", device=device)[2]
    assert np.array_equal(sitk.GetArrayFromImage(first), sitk.GetArrayFromImage(second))


def test_pipeline_forwards_the_backend(device):
    from pydicomrt.reg.pipeline import registration_pipeline

    fixed, moving, _ = deformed_pair()
    _, _, reference, _ = registration_pipeline(fixed, moving, perform_rigid=False)
    _, _, candidate, _ = registration_pipeline(
        fixed, moving, perform_rigid=False, deformable_kwargs={"backend": "torch", "device": device}
    )
    for point in probe_points(fixed):
        assert candidate.TransformPoint(point) == pytest.approx(reference.TransformPoint(point), abs=1e-3)


def test_out_of_memory_says_how_to_fit(monkeypatch):
    torch = pytest.importorskip("torch")
    from pydicomrt.reg.method import _demons_torch

    def exhausted(*args, **kwargs):
        raise torch.OutOfMemoryError("CUDA out of memory")

    monkeypatch.setattr(_demons_torch, "_multiscale_demons_torch", exhausted)
    fixed, moving, _ = deformed_pair()
    with pytest.raises(torch.OutOfMemoryError, match="isotropic_resample=True"):
        demons_registration(fixed, moving, backend="torch", device="cpu")


# ---------------------------------------------------------------------------------------
# metric="mind": across modalities
# ---------------------------------------------------------------------------------------

def mr_like(image):
    """A non-monotonic intensity map, the way tissue goes from CT to T1: air stays dark, the
    -700 insert turns brightest, the 900 one darkest, body and the others in between."""
    array = sitk.GetArrayFromImage(image)
    mapped = np.interp(array, [-1000, -700, 0, 150, 300, 900], [0, 900, 500, 300, 650, 40])
    out = sitk.GetImageFromArray(mapped.astype(np.float32))
    out.CopyInformation(image)
    return out


def mean_error(transform, truth, fixed):
    points = probe_points(fixed)
    return np.mean([np.linalg.norm(np.subtract(truth.TransformPoint(p), transform.TransformPoint(p)))
                    for p in points])


def mean_motion(transform, fixed):
    return np.mean([np.linalg.norm(np.subtract(transform.TransformPoint(p), p)) for p in probe_points(fixed)])


@pytest.mark.parametrize("direction", [None, rotation()], ids=["axial", "rotated"])
def test_mind_recovers_a_deformation_across_modalities(device, direction):
    fixed, moving, truth = deformed_pair(direction=direction)
    moving = mr_like(moving)
    before = mean_error(sitk.Transform(3, sitk.sitkIdentity), truth, fixed)
    _, transform, _ = demons_registration(fixed, moving, backend="torch", device=device, metric="mind")
    assert mean_error(transform, truth, fixed) < before / 2


def test_mind_beats_intensity_across_modalities(device):
    """The reason metric="mind" exists: intensity demons follows the wrong intensities."""
    fixed, moving, truth = deformed_pair()
    moving = mr_like(moving)
    _, mind, _ = demons_registration(fixed, moving, backend="torch", device=device, metric="mind")
    _, intensity, _ = demons_registration(fixed, moving, backend="torch", device=device)
    assert mean_error(mind, truth, fixed) < mean_error(intensity, truth, fixed) / 2


def test_mind_invents_no_motion_between_aligned_modalities(device):
    fixed = make_textured_phantom(size=(40, 36, 28), spacing=(1.0, 1.5, 2.5))
    fixed.SetDirection(rotation())
    _, transform, _ = demons_registration(fixed, mr_like(fixed), backend="torch", device=device, metric="mind")
    assert mean_motion(transform, fixed) < 0.3


def test_mind_ignores_gain_and_offset(device):
    fixed, moving, _ = deformed_pair()
    moving = mr_like(moving)
    rescaled = sitk.Cast(moving * 3.0 + 500.0, sitk.sitkFloat32)
    kwargs = dict(backend="torch", device=device, metric="mind", default_value=0.0)
    a = sitk.GetArrayFromImage(demons_registration(fixed, moving, **kwargs)[2])
    b = sitk.GetArrayFromImage(demons_registration(fixed, rescaled, **dict(kwargs, default_value=500.0))[2])
    assert np.abs(a - b).max() < 0.01


def covering_grid(image, spacing, margin=2.0):
    """An axis-aligned grid with its own spacing that covers all of an oblique image."""
    size = image.GetSize()
    corners = np.array([image.TransformContinuousIndexToPhysicalPoint((i, j, k))
                        for i in (0, size[0] - 1) for j in (0, size[1] - 1) for k in (0, size[2] - 1)])
    low, high = corners.min(0) - margin, corners.max(0) + margin
    return [int(np.ceil((h - l) / s)) + 1 for l, h, s in zip(low, high, spacing)], low.tolist()


def test_mind_with_moving_on_another_grid(device):
    """MR and CT never share a grid -- here an axial MR against an oblique CT, with its own
    spacing -- so the moving image is sampled onto the fixed one before descriptors are taken."""
    fixed, moving, truth = deformed_pair(direction=rotation())
    spacing = [0.9, 1.3, 2.2]
    size, origin = covering_grid(moving, spacing)
    regridded = sitk.Resample(mr_like(moving), size, sitk.Transform(), sitk.sitkLinear,
                              origin, spacing, (1, 0, 0, 0, 1, 0, 0, 0, 1), 0.0)
    before = mean_error(sitk.Transform(3, sitk.sitkIdentity), truth, fixed)
    _, transform, field = demons_registration(fixed, regridded, backend="torch", device=device, metric="mind")
    assert mean_error(transform, truth, fixed) < before / 2
    assert field.GetSize() == fixed.GetSize()
    assert field.GetDirection() == pytest.approx(fixed.GetDirection())


def test_mind_repeated_runs_are_identical(device):
    fixed, moving, _ = deformed_pair()
    moving = mr_like(moving)
    first = demons_registration(fixed, moving, backend="torch", device=device, metric="mind")[2]
    second = demons_registration(fixed, moving, backend="torch", device=device, metric="mind")[2]
    assert np.array_equal(sitk.GetArrayFromImage(first), sitk.GetArrayFromImage(second))


def test_pipeline_forwards_the_metric(device):
    from pydicomrt.reg.pipeline import registration_pipeline

    fixed, moving, truth = deformed_pair()
    _, _, transform, _ = registration_pipeline(
        fixed, mr_like(moving), perform_rigid=False,
        deformable_kwargs={"backend": "torch", "device": device, "metric": "mind"},
    )
    before = mean_error(sitk.Transform(3, sitk.sitkIdentity), truth, fixed)
    assert mean_error(transform, truth, fixed) < before / 2
