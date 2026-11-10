"""
rigid_registration's default algorithm ("stochastic"): accuracy, and the numpy / C++ backends
agreeing. The previous SimpleITK path, algorithm="itk", is in test_rigid_registration.py.

The compiled core (``_rigid_core``, built from ``csrc/rigid_core.cpp`` at install time) is missing
only from an install from source without a C++ compiler. Then the parity and thread tests skip and
everything else runs on numpy; ``backend="cpp"`` must then say why it cannot run.
"""

import sys

import numpy as np
import pytest
import SimpleITK as sitk

from synthetic import make_textured_phantom, phantom_centre
from pydicomrt.reg.method import _stochastic_rigid as stochastic
from pydicomrt.reg.method.rigid import rigid_registration
from reg.test_rigid_registration import known_rigid_transform, mapping_error_mm, warp


@pytest.fixture(scope="module")
def phantom():
    return make_textured_phantom()


@pytest.fixture(scope="module")
def thick_rotated_case(phantom):
    """
    A moving image on 6 mm slices, on a rotated grid, with a non-monotonic intensity map, and
    the fixed -> moving truth. Thick slices are where ITK's precomputed gradient went wrong;
    the rotated grid is what catches direction-matrix mistakes; the intensity map needs MI.
    """
    size, spacing = np.array([40, 40, 12]), np.array([2.5, 2.5, 6.0])
    rotation = sitk.Euler3DTransform()
    rotation.SetRotation(0.1, -0.05, 0.3)
    direction = np.array(rotation.GetMatrix()).reshape(3, 3)
    grid = sitk.Image(size.tolist(), sitk.sitkFloat32)
    grid.SetSpacing(spacing.tolist())
    grid.SetDirection(direction.ravel().tolist())
    grid.SetOrigin((np.array(phantom_centre()) - direction @ ((size - 1) * spacing / 2)).tolist())
    truth = known_rigid_transform(5.0, (-4.0, 3.0, 2.0))
    moving = sitk.Resample(phantom, grid, truth.GetInverse(), sitk.sitkLinear, -1000.0)
    return sitk.Abs(moving - 100.0), truth


CORE = "pydicomrt.reg.method._rigid_core"


def compiled_core():
    return pytest.importorskip(CORE)


class TestAccuracy:
    def test_recovers_a_known_rigid_transform(self, phantom):
        truth = known_rigid_transform(-8.0, (15.0, -9.0, 7.0))

        recovered = rigid_registration(phantom, warp(phantom, truth))

        # Measured 0.19 mm.
        assert mapping_error_mm(recovered, truth) < 1.0

    def test_thick_rotated_multimodal_on_numpy(self, phantom, thick_rotated_case):
        """Also the end-to-end run of the numpy backend, whether or not the core is installed."""
        moving, truth = thick_rotated_case

        recovered = rigid_registration(phantom, moving, backend="numpy")

        # Measured 0.22 mm; the "itk" default lands 1.2 mm out on the same case.
        assert mapping_error_mm(recovered, truth) < 1.0

    @pytest.mark.parametrize("shift", [(8.0, -6.0, 4.0), (0.0, 0.0, 0.0), (0.0, 150.0, 0.0)])
    def test_does_not_walk_away_from_an_exact_start(self, phantom, shift):
        """
        A whole-grid shift: the geometric-centre start is already the answer, and the two
        images agree exactly there, so the gradient that sets the step gain is ~0. Before
        steps were capped, the next ordinary gradient times that gain threw the transform
        53 mm away (150 mm case). The offset larger than the image is the one the
        initialiser exists for.
        """
        moving = sitk.Image(phantom)
        moving.SetOrigin(tuple(o + s for o, s in zip(phantom.GetOrigin(), shift)))

        recovered = rigid_registration(phantom, moving)

        # Measured below 0.02 mm on both backends.
        assert mapping_error_mm(recovered, sitk.TranslationTransform(3, shift)) < 0.5

    def test_returns_a_centre_free_rigid_affine(self, phantom):
        recovered = rigid_registration(phantom, phantom)

        rotation = np.array(recovered.GetMatrix()).reshape(3, 3)
        assert isinstance(recovered, sitk.AffineTransform)
        assert recovered.GetCenter() == (0.0, 0.0, 0.0)
        assert rotation @ rotation.T == pytest.approx(np.eye(3), abs=1e-12)
        assert np.linalg.det(rotation) == pytest.approx(1.0)


class TestBackends:
    def test_random_points_match_the_cores_generator(self):
        """Values printed by the C++ generator; numpy must reproduce them bit for bit."""
        got = stochastic._uniform(5, 3, np.arange(4, dtype=np.uint64))

        assert got.tolist() == [0.17944501125524803, 0.65760126145369791,
                                0.93854684261697296, 0.87877132867957486]

    @pytest.fixture
    def level_inputs(self, phantom):
        moving = warp(phantom, known_rigid_transform(5.0, (-4.0, 3.0, 0.0)))
        fixed = stochastic.volume(sitk.Cast(phantom, sitk.sitkFloat32))
        moving = stochastic.volume(sitk.Cast(moving, sitk.sitkFloat32))
        return fixed, moving, np.array(phantom_centre())

    def test_one_evaluation_agrees_to_rounding(self, level_inputs):
        core = compiled_core()
        fixed, moving, centre = level_inputs
        rotation = stochastic._skew_exp(np.array([0.03, -0.02, 0.05]))
        points = np.random.default_rng(1).random((20000, 3)) * (np.array(fixed[0].shape[::-1]) - 1)

        mi_np, grad_np, n_np = stochastic.evaluate(fixed, moving, centre, 32, rotation, (2.0, -1.0, 1.0), points)
        mi_c, grad_c, n_c = core.evaluate(fixed, moving, centre, 32, rotation, (2.0, -1.0, 1.0), points, threads=4)

        assert n_np == n_c
        assert mi_np == pytest.approx(mi_c, rel=1e-12)
        assert grad_np == pytest.approx(grad_c, rel=1e-9, abs=1e-12 * np.abs(grad_np).max())

    def test_a_pyramid_level_follows_the_same_path(self, level_inputs):
        """
        Same random points, same steps. Only the first iterations can be compared: stochastic
        ascent amplifies a rounding difference two- to threefold per step on this phantom
        (1e-15 at the first iterate, millimetres by the thirtieth), which is also why whole
        registrations differ by tenths of a millimetre between backends.
        """
        core = compiled_core()
        fixed, moving, centre = level_inputs
        args = (fixed, moving, centre, 50, 7, 4, np.eye(3), np.zeros(3), 4, 4096, 2.0)

        _, _, path_np = stochastic.level(*args)
        _, _, path_c = core.level(*args)

        assert np.abs(path_np - path_c).max() < 1e-10

    def test_answer_does_not_depend_on_the_thread_count(self, phantom):
        compiled_core()
        moving = warp(phantom, known_rigid_transform(5.0, (-4.0, 3.0, 0.0)))

        one = rigid_registration(phantom, moving, backend="cpp", threads=1)
        four = rigid_registration(phantom, moving, backend="cpp", threads=4)

        assert one.GetParameters() == four.GetParameters()

    def test_cpp_without_the_core_says_how_to_get_it(self, monkeypatch, phantom):
        monkeypatch.setitem(sys.modules, CORE, None)  # `import` now raises ImportError

        with pytest.raises(ImportError, match="C\\+\\+17 compiler"):
            rigid_registration(phantom, phantom, backend="cpp")

    def test_auto_falls_back_to_numpy(self, monkeypatch):
        monkeypatch.setitem(sys.modules, CORE, None)

        assert stochastic._core("auto") is stochastic


class TestValidation:
    def test_unknown_algorithm_is_rejected(self, phantom):
        with pytest.raises(ValueError, match="stochastic"):
            rigid_registration(phantom, phantom, algorithm="elastix")

    def test_unknown_backend_is_rejected(self, phantom):
        with pytest.raises(ValueError, match="numpy"):
            rigid_registration(phantom, phantom, backend="torch")

    def test_2d_images_are_rejected(self):
        image = sitk.Image([32, 32], sitk.sitkFloat32)
        with pytest.raises(ValueError, match="3D"):
            rigid_registration(image, image)
