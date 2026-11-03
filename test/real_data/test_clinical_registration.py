"""
Rigid registration measured against a TPS-produced registration (opt in with -m slow).

Runs full-size CT and CBCT volumes single-threaded; the default (stochastic) rigid takes about
15 s here. Each registration is shared through a module- or class-scoped fixture. The numbers
in the tolerances below were measured on this dataset; see ``conftest.py`` for the expected
data layout.

What this dataset can and cannot check:

- It CAN check the DICOM plumbing (loading, ordering, rescale, REG parsing) and rigid
  registration accuracy against a clinically produced answer.
- It CANNOT check the row/column spacing or direction-matrix conventions: both series are
  axial with square pixels, which is precisely the case where those bugs are invisible.
  Anisotropic and oblique coverage stays with the synthetic fixtures.
"""

import numpy as np
import pytest
import SimpleITK as sitk

from pydicomrt.reg.method.rigid import rigid_registration
from pydicomrt.reg.pipeline import registration_pipeline

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module", autouse=True)
def single_threaded():
    """
    The default algorithm's answer does not depend on the thread count, but the pipeline's
    preprocessing and any ``algorithm="itk"`` run do, so the module pins one thread.
    """
    previous = sitk.ProcessObject.GetGlobalDefaultNumberOfThreads()
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)
    yield
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(previous)


@pytest.fixture(scope="module")
def recovered_transform(planning_ct_image, cbct_image):
    return rigid_registration(planning_ct_image, cbct_image)


def translation_delta(recovered, clinical) -> np.ndarray:
    return np.array(recovered.GetTranslation()) - np.array(clinical.GetTranslation())


class TestAgainstClinicalRegistration:
    """
    Measured offset from the TPS registration (+1.6, -4.7, -1.1) mm, 5.1 mm in all (seeds
    0-2: 5.09-5.16 mm; numpy and C++ backends within 0.01 mm). elastix's default rigid
    finishes 16.7 mm out on this pair, 16 mm of it superior-inferior. The bounds lock in
    that result, not clinical adequacy -- one case is not a validation.
    """

    def test_recovers_the_left_right_offset(self, recovered_transform, clinical_transform):
        delta = translation_delta(recovered_transform, clinical_transform)

        assert abs(delta[0]) < 3.0

    def test_recovers_the_anterior_posterior_offset(
        self, recovered_transform, clinical_transform
    ):
        """
        The worst axis, at about 5 mm. Part of that gap is not an optimizer failure:
        sweeping y around the clinical answer shows mutual information over the full field
        of view peaking away from it, so the metric and the TPS genuinely disagree.
        """
        delta = translation_delta(recovered_transform, clinical_transform)

        assert abs(delta[1]) < 6.0

    def test_recovers_the_superior_inferior_offset(
        self, recovered_transform, clinical_transform
    ):
        """
        The axis where the CBCT, about a fifth of the planning CT's length, leaves the most
        room to slide -- and where elastix ends up 16 mm out.
        """
        delta = translation_delta(recovered_transform, clinical_transform)

        assert abs(delta[2]) < 3.0

    def test_residual_rotation_is_small(self, recovered_transform, clinical_transform):
        rotation = np.array(recovered_transform.GetMatrix()).reshape(3, 3)
        angle = np.degrees(np.arccos(np.clip((np.trace(rotation) - 1) / 2, -1, 1)))

        assert angle < 1.0

    def test_agrees_with_the_clinical_registration_overall(
        self, recovered_transform, clinical_transform
    ):
        residual = np.linalg.norm(translation_delta(recovered_transform, clinical_transform))
        clinical_offset = np.linalg.norm(clinical_transform.GetTranslation())

        assert residual < 6.5
        assert residual < clinical_offset / 10.0


class TestPipelineOnRealData:
    """
    The pipeline used to land ~194 mm from the clinical answer. It resampled both images onto
    a shared union-extent grid before registering, which left the initialiser nothing to
    work with -- aligning the centres of two identical grids yields zero -- so the ~180 mm
    couch offset between the CT and the CBCT was never estimated. Registering in the images'
    own grids fixed it.
    """

    @pytest.fixture(scope="class")
    def pipeline_result(self, planning_ct_image, cbct_image):
        return registration_pipeline(
            planning_ct_image,
            cbct_image,
            perform_rigid=True,
            perform_deformable=False,
        )

    def test_pipeline_matches_direct_registration(
        self, pipeline_result, clinical_transform, recovered_transform
    ):
        _, pipeline_rigid, _, _ = pipeline_result

        # With no preprocessing configured the pipeline does exactly what the bare function
        # does, so the two agree to the last digit rather than merely closely.
        assert np.array(pipeline_rigid.GetTranslation()) == pytest.approx(
            np.array(recovered_transform.GetTranslation()), abs=1e-9
        )

        pipeline_error = np.linalg.norm(translation_delta(pipeline_rigid, clinical_transform))
        assert pipeline_error < 6.5

    def test_output_lands_on_the_planning_ct_grid(self, pipeline_result, planning_ct_image):
        registered = pipeline_result[0]

        assert registered.GetSize() == planning_ct_image.GetSize()
        assert registered.GetOrigin() == pytest.approx(planning_ct_image.GetOrigin())
        assert registered.GetSpacing() == pytest.approx(planning_ct_image.GetSpacing())
