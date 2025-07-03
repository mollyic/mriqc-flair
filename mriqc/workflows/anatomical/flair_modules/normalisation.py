
# Class to run report functionality 
import niworkflows.interfaces.reportlets.base as nrb
from nipype.interfaces.ants import (RegistrationSynQuick)
from nipype.interfaces.ants.registration import (RegistrationSynQuickOutputSpec, RegistrationSynQuickInputSpec)
from nipype.interfaces.mixins import reporting
from mriqc import config


class _RegistrationSynQuickInputSpecRPT(nrb._SVGReportCapableInputSpec, RegistrationSynQuickInputSpec):
    pass

class _RegistrationSynQuickOutputSpecRPT(reporting.ReportCapableOutputSpec, RegistrationSynQuickOutputSpec):
    pass
class RegistrationSynQuickRPT(nrb.RegistrationRC, RegistrationSynQuick):
    """Report generating version of ANTs RegistrationSynQuick interface: 
        *derived from niworkflows.interfaces.reportlets.registration."""
    input_spec = _RegistrationSynQuickInputSpecRPT
    output_spec = _RegistrationSynQuickOutputSpecRPT
    def _post_run_hook(self, runtime):
        # Get arguments from ANTS
        self._fixed_image = self.inputs.fixed_image
        if isinstance(self._fixed_image, (list, tuple)):
            self._fixed_image = self.inputs.fixed_image[0]
        
        self._moving_image = self.aggregate_outputs(runtime=runtime).warped_image
        config.loggers.workflow.info(
            "Report - setting fixed (%s) and moving (%s) images",
            self._fixed_image,
            self._moving_image,
        )

        return super()._post_run_hook(runtime)