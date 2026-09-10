"""
Rigid registration measured against a TPS-produced registration (opt in with -m slow).

Runs full-size CT and CBCT volumes single-threaded, so one registration takes roughly nine
minutes and the module about twenty. Each registration is shared through a module- or
class-scoped fixture; re-running them per test cost an hour. The numbers in the tolerances
below were measured on this dataset; see ``conftest.py`` for the expected data layout.

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

# Library defaults. A four-level pyramid was measured too and came out slightly worse here
# (8.5 mm against 7.5 mm) as well as slower, so there is nothing to override.
DEFAULTS = {}


@pytest.fixture(scope="module", autouse=True)
def single_threaded():
    """
    ITK reduces the metric in thread completion order. On this pair the MI optimum along z
    is shallow enough that the thread count changes the answer by more than 15 mm, so these
    measurements are only meaningful single-threaded.
    """
    previous = sitk.ProcessObject.GetGlobalDefaultNumberOfThreads()
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)
    yield
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(previous)


@pytest.fixture(scope="module")
def recovered_transform(planning_ct_image, cbct_image):
    return rigid_registration(planning_ct_image, cbct_image, **DEFAULTS)


def translation_delta(recovered, clinical) -> np.ndarray:
    return np.array(recovered.GetTranslation()) - np.array(clinical.GetTranslation())


class TestAgainstClinicalRegistration:
    def test_recovers_the_left_right_offset(self, recovered_transform, clinical_transform):
        """x is the axis the metric constrains best here."""
        delta = translation_delta(recovered_transform, clinical_transform)

        assert abs(delta[0]) < 5.0

    def test_recovers_the_anterior_posterior_offset(
        self, recovered_transform, clinical_transform
    ):
        """
        y lands within about 10 mm. Part of that gap is not an optimizer failure: sweeping y
        around the clinical answer shows mutual information over the full field of view
        peaking roughly 10 mm away from it, so the metric and the TPS genuinely disagree.
        """
        delta = translation_delta(recovered_transform, clinical_transform)

        assert abs(delta[1]) < 12.0

    def test_residual_rotation_is_small(self, recovered_transform, clinical_transform):
        rotation = np.array(recovered_transform.GetMatrix()).reshape(3, 3)
        angle = np.degrees(np.arccos(np.clip((np.trace(rotation) - 1) / 2, -1, 1)))

        assert angle < 2.0

    def test_recovers_the_superior_inferior_offset(
        self, recovered_transform, clinical_transform
    ):
        """
        This was the worst axis by a wide margin under the old gradient-descent optimizer
        (~37 mm out) and looked like a field-of-view limitation: the CBCT covers about a
        fifth of the planning CT, and z is where the two differ most. It was not. Switching
        to a step that shrinks on gradient reversal brought it to under 1 mm, which means
        the optimizer had been walking away from the answer rather than failing to see it.
        """
        delta = translation_delta(recovered_transform, clinical_transform)

        assert abs(delta[2]) < 5.0

    def test_agrees_with_the_clinical_registration_overall(
        self, recovered_transform, clinical_transform
    ):
        """
        Measured at 7.5 mm, essentially all of it the y disagreement above. 12 mm is a
        regression bound rather than a statement of clinical adequacy -- one case is not a
        validation, and the residual should be checked against your own tolerance.
        """
        residual = np.linalg.norm(translation_delta(recovered_transform, clinical_transform))
        clinical_offset = np.linalg.norm(clinical_transform.GetTranslation())

        assert residual < 12.0
        assert residual < clinical_offset / 10.0


class TestConfigurationRegression:
    """
    Locks in the improvement over the original configuration, so a future change to the
    defaults cannot quietly undo it on real data.

    Measured on this pair: 41.6 mm for the original unbounded single-resolution gradient
    descent, 37.5 mm once it was made bounded and multi-resolution, and 7.5 mm after
    switching to a regular step. Multi-resolution alone barely helped -- the optimizer was
    leaving good positions, not failing to reach them.
    """

    def test_beats_the_unbounded_single_resolution_configuration(
        self, planning_ct_image, cbct_image, clinical_transform, recovered_transform
    ):
        def previous_configuration(fixed, moving):
            fixed = sitk.Cast(fixed, sitk.sitkFloat32)
            moving = sitk.Cast(moving, sitk.sitkFloat32)
            method = sitk.ImageRegistrationMethod()
            method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=100)
            method.SetInterpolator(sitk.sitkLinear)
            method.SetOptimizerAsGradientDescent(
                learningRate=2.0, numberOfIterations=100,
                convergenceMinimumValue=1e-6, convergenceWindowSize=10,
            )
            method.SetOptimizerScalesFromPhysicalShift()
            initial = sitk.CenteredTransformInitializer(
                fixed, moving, sitk.VersorRigid3DTransform(),
                sitk.CenteredTransformInitializerFilter.GEOMETRY,
            )
            initial.SetCenter((0.0, 0.0, 0.0))
            method.SetInitialTransform(initial, inPlace=False)
            result = method.Execute(fixed, moving)
            result = (result.GetNthTransform(0)
                      if isinstance(result, sitk.CompositeTransform) else result)
            affine = sitk.AffineTransform(3)
            affine.SetMatrix(result.GetMatrix())
            affine.SetTranslation(result.GetTranslation())
            return affine

        previous = previous_configuration(planning_ct_image, cbct_image)

        current_error = np.linalg.norm(
            translation_delta(recovered_transform, clinical_transform)
        )
        previous_error = np.linalg.norm(translation_delta(previous, clinical_transform))

        assert current_error < previous_error


class TestPipelineOnRealData:
    """
    The pipeline used to land ~194 mm from the clinical answer against ~7.5 mm for the bare
    function. It resampled both images onto a shared union-extent grid before registering,
    which reduced CenteredTransformInitializer to a no-op -- aligning the centres of two
    identical grids yields zero -- so the ~180 mm couch offset between the CT and the CBCT
    was never estimated. Registering in the images' own grids fixed it.
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

        # With no preprocessing configured the pipeline now does exactly what the bare
        # function does, so the two agree to the last digit rather than merely closely.
        assert np.array(pipeline_rigid.GetTranslation()) == pytest.approx(
            np.array(recovered_transform.GetTranslation()), abs=1e-9
        )

        pipeline_error = np.linalg.norm(translation_delta(pipeline_rigid, clinical_transform))
        assert pipeline_error < 12.0

    def test_output_lands_on_the_planning_ct_grid(self, pipeline_result, planning_ct_image):
        registered = pipeline_result[0]

        assert registered.GetSize() == planning_ct_image.GetSize()
        assert registered.GetOrigin() == pytest.approx(planning_ct_image.GetOrigin())
        assert registered.GetSpacing() == pytest.approx(planning_ct_image.GetSpacing())
