import shutil
import urllib.request
from typing import List

from nipype.interfaces import utility as niu
from nipype.interfaces.ants import RegistrationSynQuick
from nipype.interfaces.ants.registration import (
    RegistrationSynQuickInputSpec,
    RegistrationSynQuickOutputSpec,
)
from nipype.interfaces.base import traits
from nipype.interfaces.mixins import reporting
from nipype.pipeline import engine as pe
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from niworkflows.interfaces.norm import (
    NIWORKFLOWS_LOG,
    SpatialNormalization,
    _SpatialNormalizationInputSpec,
)
from niworkflows.interfaces.reportlets import base as nrb
from niworkflows.interfaces.reportlets.registration import (
    _SpatialNormalizationOutputSpecRPT,
)

from mriqc import config


class _RegistrationSynQuickInputSpecRPT(
    nrb._SVGReportCapableInputSpec, RegistrationSynQuickInputSpec
):
    pass


class _RegistrationSynQuickOutputSpecRPT(
    reporting.ReportCapableOutputSpec, RegistrationSynQuickOutputSpec
):
    pass


class RegistrationSynQuickRPT(nrb.RegistrationRC, RegistrationSynQuick):
    """
    Custom report-capable interface for ANTs RegistrationSynQuick.

    This class adds SVG report generation to the lightweight RegistrationSynQuick
    interface, using B-spline SyN registration [Tustison2013].

    .. [Tustison2013] Tustison NJ, Avants BB., *Explicit B-spline regularization in
        diffeomorphic image registration*, Front Neuroinform 7(39), 2013
        doi:`10.3389/fninf.2013.00039 <http://dx.doi.org/10.3389/fninf.2013.00039>`_.

    Created by Molly Ireland
    """

    input_spec = _RegistrationSynQuickInputSpecRPT
    output_spec = _RegistrationSynQuickOutputSpecRPT

    def _post_run_hook(self, runtime):
        from mriqc import config

        # Get arguments from ANTS
        self._fixed_image = self.inputs.fixed_image
        if isinstance(self._fixed_image, (list, tuple)):
            self._fixed_image = self.inputs.fixed_image[0]

        self._moving_image = self.aggregate_outputs(runtime=runtime).warped_image
        config.loggers.workflow.info(
            'Report - setting fixed (%s) and moving (%s) images',
            self._fixed_image,
            self._moving_image,
        )

        return super()._post_run_hook(runtime)


def quicksyn_normalisation(name='QuickSynSpatialNormalization'):
    """Create a simplified workflow to perform spatial normalization with Ants QuickSyn."""
    from templateflow.api import get as get_template

    from mriqc.workflows.anatomical.flair_modules.normalisation import (
        RegistrationSynQuickRPT as RobustMNINormalization,
    )

    # Define workflow interface
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'moving_image',
                'moving_mask',
                'modality',
                'tpl_resolution',
                'tpl_id',
                'tpl_tissues',
                'reference_image',
                'reference_mask',
            ]
        ),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['out_tpms', 'out_report', 'ind2std_xfm', 'hmask_mni2nat']),
        name='outputnode',
    )

    # Spatial normalization
    norm = pe.Node(
        RobustMNINormalization(
            dimension=3,
            num_threads=config.nipype.omp_nthreads,
            transform_type='b',
            generate_report=True,
        ),
        name='SpatialNormalization',
        # Request all MultiProc processes when ants_nthreads > n_procs
        num_threads=config.nipype.omp_nthreads,
        mem_gb=3,
    )

    # Project standard TPMs into T1w space
    tpms_std2t1w = pe.MapNode(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            interpolation='Linear',
            float=config.execution.ants_float,
        ),
        iterfield=['input_image'],
        name='tpms_std2t1w',
    )
    tpms_std2t1w.inputs.invert_transform_flags = [True]

    # Project standard headmask into native space
    hmask_mni2nat = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            interpolation='Linear',
            float=config.execution.ants_float,
        ),
        name='hmask_mni2nat',
    )

    hmask_mni2nat.inputs.input_image = get_template(
        template='MNI152Lin', desc='head', suffix='mask', resolution='2'
    )
    hmask_mni2nat.inputs.invert_transform_flags = [True]

    # fmt: off
    workflow.connect([
        (inputnode, norm, [('moving_image', 'moving_image'),
                           ('reference_image', 'fixed_image')]),
        (inputnode, tpms_std2t1w, [('moving_image', 'reference_image'),
                                   ('tpl_tissues', 'input_image')]),
        (norm, tpms_std2t1w, [
            ('out_matrix', 'transforms')
        ]),
        (norm, outputnode, [
            ('out_matrix', 'ind2std_xfm'),
            ('out_report', 'out_report'),
        ]),
        (tpms_std2t1w, outputnode, [('output_image', 'out_tpms')]),
        (inputnode, hmask_mni2nat, [('moving_image', 'reference_image')]),
        (norm, hmask_mni2nat, [('out_matrix', 'transforms')]),
        (hmask_mni2nat, outputnode, [('output_image', 'hmask_mni2nat')]),
    ])
    # fmt: on

    return workflow


class _WrapSpatialNormalizationInputSpec(_SpatialNormalizationInputSpec):
    reference = traits.Enum(
        'T1w',
        'T2w',
        'boldref',
        'PDw',
        'FLAIR',
        mandatory=True,
        usedefault=True,
        desc='set the reference modality for registration (extended)',
    )


class WrapSpatialNormalization(SpatialNormalization):
    input_spec = _WrapSpatialNormalizationInputSpec


class _WrapSpatialNormalizationInputSpecRPT(
    nrb._SVGReportCapableInputSpec, _WrapSpatialNormalizationInputSpec
):
    pass


class WrapSpatialNormalizationRPT(nrb.RegistrationRC, WrapSpatialNormalization):
    """
    Wrapper around NiWorkflows' SpatialNormalization and its reporting variant
    to extend support for FLAIR.
    This class adds SVG report generation to the SpatialNormalization interface
    """

    input_spec = _WrapSpatialNormalizationInputSpecRPT
    output_spec = _SpatialNormalizationOutputSpecRPT

    def _post_run_hook(self, runtime):
        # We need to dig into the internal ants.Registration interface
        self._fixed_image = self._get_ants_args()['fixed_image']
        if isinstance(self._fixed_image, (list, tuple)):
            self._fixed_image = self._fixed_image[0]  # get first item if list

        if self._get_ants_args().get('fixed_image_mask') is not None:
            self._fixed_image_mask = self._get_ants_args().get('fixed_image_mask')
        self._moving_image = self.aggregate_outputs(runtime=runtime).warped_image
        NIWORKFLOWS_LOG.info(
            'Report - setting fixed (%s) and moving (%s) images',
            self._fixed_image,
            self._moving_image,
        )

        return super()._post_run_hook(runtime)


def _download_file(url: str, output_path: str) -> None:
    urllib.request.urlretrieve(url, output_path)
    print(f'\t * Downloaded: {output_path}')
    return output_path


def _get_custom_templates(
    modality: str = 'FLAIR', template_id: str = 'GG853', template_res: int = 1
) -> List:
    """
    Wrapper function to handle FLAIR templates not available on TemplateFlow Database. Downloads Braindr files from S3

    Parameters
    ----------
    modality:  Modality of the template
    template_id: Template id (FLAIR: 'GG853').
    template_res: Resolution level of the template (default: 1).
    """

    import json
    from pathlib import Path

    from templateflow import conf

    license_text = """\
    Copyright (c) 2012 Anderson M. Winkler, Peter Kochunov, David C. Glahn.

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    """

    tpl_desc = {
        'Authors': ['Winkler AM', 'Kochunov P', 'Glahn DC'],
        'BIDSVersion': '1.1.0',
        'HowToAcknowledge': 'If you use these templates in your research, please cite it as: Winkler AM, Kochunov P, Glahn DC. FLAIR Templates. Available at http://brainder.org.',
        'Identifier': 'GG853',
        'License': 'See LICENSE file',
        'Name': 'GG853 - FLAIR templates for the Genetics of Brain Structure and Function Study (GOBS) 2012',
        'ReferencesAndLinks': [
            'https://brainder.org/download/flair/',
            'https://s3.us-east-2.amazonaws.com/brainder/publications/2009/HBM2009_flair_poster.pdf',
        ],
        'TemplateFlowVersion': '1.0.0',
        'res': {
            '01': {
                'origin': [-91.0, -126.0, -72.0],
                'shape': [182, 218, 182],
                'zooms': [1.0, 1.0, 1.0],
            }
        },
    }

    tf_tpl_dir: Path = conf.TF_HOME / ('tpl-' + template_id)
    tpl_license = tf_tpl_dir / 'LICENSE'
    tpl_desc_json = tf_tpl_dir / 'template_description.json'

    tf_tpl_dir.mkdir(parents=True, exist_ok=True)

    # Synthstripped brain mask of GG853
    brain_mask = Path.cwd() / 'mriqc' / 'templates' / 'tpl-GG853_res-01_desc-brain_mask.nii.gz'
    if not (tf_tpl_dir / brain_mask.name).exists():
        shutil.copy2(brain_mask, tf_tpl_dir / brain_mask.name)

    # Files to download from Braindr
    dwnld_head = 'https://s3.us-east-2.amazonaws.com/brainder/software/flair/templates/'
    dwnld_lst = [
        ('GG-853-FLAIR-1.0mm.nii.gz', f'tpl-{template_id}_res-0{template_res}_{modality}.nii.gz'),
        (
            'GG-853-GM-1.0mm.nii.gz',
            f'tpl-{template_id}_res-0{template_res}_label-GM_probseg.nii.gz',
        ),
        (
            'GG-853-WM-1.0mm.nii.gz',
            f'tpl-{template_id}_res-0{template_res}_label-WM_probseg.nii.gz',
        ),
        (
            'GG-853-CSF-1.0mm.nii.gz',
            f'tpl-{template_id}_res-0{template_res}_label-CSF_probseg.nii.gz',
        ),
    ]

    for url_name, tf_name in dwnld_lst:
        out_path = tf_tpl_dir / tf_name
        if not out_path.exists():
            _download_file(f'{dwnld_head}{url_name}', out_path)

    if not tpl_desc_json.exists():
        with tpl_desc_json.open('w') as f:
            json.dump(tpl_desc, f, indent=4)

    if not tpl_license.exists():
        tpl_license.write_text(license_text)


if __name__ == '__main__':
    """Run the FLAIR template check/download script"""

    print('Running FLAIR template check/download...')
    _get_custom_templates()
