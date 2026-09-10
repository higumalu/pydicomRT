"""
Shared test fixtures.

Everything here is synthesised in-memory: no patient data on disk, no network, and a fixed
RNG seed so the geometry assertions are reproducible. The builders themselves live in
synthetic.py. Tests that need real clinical DICOM are under test/real_data/.
"""

import pytest
import SimpleITK as sitk

from synthetic import make_blob_image, make_ct_series, make_textured_phantom

# ITK reduces metric values across threads in whatever order they finish, so a registration
# run is only bit-reproducible when it is single-threaded. Left at the default, the same
# registration here varied by ~2.7 mm between runs, which is enough to make an accuracy
# assertion flap depending on machine load.
sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)


@pytest.fixture
def textured_phantom():
    return make_textured_phantom()


@pytest.fixture
def ct_series():
    """Isotropic-in-plane synthetic CT series (the common real-world case)."""
    return make_ct_series()


@pytest.fixture
def anisotropic_ct_series():
    """
    CT series with distinguishable row and column spacing.

    Square pixels hide axis-ordering mistakes, so anything asserting which spacing lands on
    which axis must use this fixture.
    """
    return make_ct_series(shape=(4, 8, 12), pixel_spacing=(3.0, 1.5), slice_spacing=5.0)


@pytest.fixture
def blob_image():
    return make_blob_image()
