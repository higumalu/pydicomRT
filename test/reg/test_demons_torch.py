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
