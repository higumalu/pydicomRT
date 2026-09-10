"""
Fixtures for tests that run against real clinical DICOM.

The data itself is never committed -- ``testdata/`` is gitignored. Point
``PYDICOMRT_TESTDATA`` at a directory laid out as::

    <testdata>/CTSIM/*.dcm          planning CT series
    <testdata>/CBCT/*.dcm           on-treatment CBCT series
    <testdata>/registration.dcm     Spatial Registration produced by the TPS

Every test here skips cleanly when that directory is absent, so the suite still runs on a
fresh checkout and in CI.
"""

import os
from pathlib import Path

import numpy as np
import pytest
import SimpleITK as sitk

pydicom = pytest.importorskip("pydicom")

from pydicomrt.utils import load_sorted_image_series
from pydicomrt.utils.sitk_transform import SimpleITKImageBuilder


def testdata_root() -> Path:
    return Path(os.environ.get("PYDICOMRT_TESTDATA", "testdata"))


def _require(path: Path) -> Path:
    if not path.exists():
        pytest.skip(f"real test data not found at {path} (set PYDICOMRT_TESTDATA)")
    return path


@pytest.fixture(scope="session")
def series_dir():
    """Resolve a series directory by name, skipping when the data is not present."""
    def resolve(name: str) -> Path:
        return _require(testdata_root() / name)
    return resolve


@pytest.fixture(scope="session")
def planning_ct_series():
    return load_sorted_image_series(str(_require(testdata_root() / "CTSIM")))


@pytest.fixture(scope="session")
def cbct_series():
    return load_sorted_image_series(str(_require(testdata_root() / "CBCT")))


@pytest.fixture(scope="session")
def clinical_reg_ds():
    return pydicom.dcmread(str(_require(testdata_root() / "registration.dcm")))


@pytest.fixture(scope="session")
def planning_ct_image(planning_ct_series):
    return SimpleITKImageBuilder().from_image_series(planning_ct_series)


@pytest.fixture(scope="session")
def cbct_image(cbct_series):
    return SimpleITKImageBuilder().from_image_series(cbct_series)


@pytest.fixture(scope="session")
def clinical_transform(clinical_reg_ds, cbct_series):
    """
    The TPS registration, expressed the way ``rigid_registration`` returns transforms.

    A DICOM REG matrix maps points from the item's Frame of Reference into the registered
    RCS, i.e. moving -> fixed. A SimpleITK resampling transform runs the other way, mapping
    a point on the fixed grid back into the moving image, so the stored matrix is inverted
    here. Getting this backwards still produces a plausible-looking number, so the
    direction is asserted against image content in ``test_clinical_transform.py``.
    """
    moving_frame_of_reference = cbct_series[0].FrameOfReferenceUID

    for item in clinical_reg_ds.RegistrationSequence:
        if item.FrameOfReferenceUID != moving_frame_of_reference:
            continue
        matrix = np.array(
            item.MatrixRegistrationSequence[0].MatrixSequence[0]
            .FrameOfReferenceTransformationMatrix,
            dtype=float,
        ).reshape(4, 4)

        inverse = np.linalg.inv(matrix)
        transform = sitk.AffineTransform(3)
        transform.SetMatrix(inverse[:3, :3].ravel().tolist())
        transform.SetTranslation(inverse[:3, 3].tolist())
        return transform

    pytest.fail(
        f"no registration item for the CBCT frame of reference {moving_frame_of_reference}"
    )
