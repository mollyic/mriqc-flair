# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2021 The NiPreps Developers <nipreps@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# We support and encourage derived works from this project, please read
# about our expectations at
#
#     https://www.nipreps.org/community/licensing/
#
""" Testing version for spatial Normalisation with reports with normal spatial normalization"""

# general purpose
from pkg_resources import resource_filename as pkgr_fn
from warnings import warn

# nipype
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.ants import AI, ImageMath, ThresholdImage, RegistrationSynQuick

# niworkflows
from niworkflows.interfaces.fixes import (
    FixHeaderRegistration as Registration,
    FixHeaderApplyTransforms as ApplyTransforms,
)
from niworkflows.interfaces.nibabel import ApplyMask, RegridToZooms
from mriqc import config
from niworkflows.interfaces.reportlets.base import RegistrationRC

from mriqc.workflows.anatomical.base import init_brain_tissue_segmentation


def spatial_normalization(name='SpatialNormalization'):
    """Create a simplified workflow to perform fast spatial normalization."""
    from mriqc.workflows.anatomical.flair_modules.normalisation import (
        WrapSpatialNormalizationRPT as RobustMNINormalization,
    )
    from pathlib import Path

    # Have the template id handy
    tpl_id = config.workflow.template_id

    # Define workflow interface
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['moving_image', 
                                      'moving_mask', 
                                      'modality', 
                                      'reference_mask', 
                                      'reference_image', 
                                      'tpl_resolution', 
                                      'tpl_id', 
                                      'tpl_tissues']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['out_tpms', 'out_report', 'ind2std_xfm', 'hmask_mni2nat']),
        name='outputnode',
    )

    # Spatial normalization
    norm = pe.Node(
        RobustMNINormalization(
            num_threads=config.nipype.omp_nthreads,
            float=config.execution.ants_float,
            generate_report=True,
        ),
        name='SpatialNormalization',
        # Request all MultiProc processes when ants_nthreads > n_procs
        num_threads=config.nipype.omp_nthreads,
        mem_gb=3,
    )
    if config.workflow.species.lower() == 'human':
        workflow.connect([
        (inputnode, norm, [('reference_mask', 'reference_mask'),
                           ('reference_image', 'reference_image'),
                           ('tpl_resolution', 'flavor'), 
                           ('tpl_id', 'template'), 
                           ])
                           ])
    else:
        norm.inputs.reference_image = str(get_template(tpl_id, suffix='T2w'))
        norm.inputs.reference_mask = str(get_template(tpl_id, desc='brain', suffix='mask')[0])
        norm.inputs.flavor = ['testing', 'fast'][config.execution.debug] 
        norm.inputs.template = tpl_id

    # Project standard TPMs into T1w space
    tpms_std2t1w = pe.MapNode(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            interpolation='Gaussian',
            float=config.execution.ants_float,
        ),
        iterfield=['input_image'],
        name='tpms_std2t1w',
    )

    #tpms_std2t1w.inputs.invert_transform_flags = [True]

    # Project standard headmask to native space
    hmask_std2t1w = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            interpolation='Gaussian',  # Choose the appropriate interpolation method
            float=config.execution.ants_float,
        ),
        name='syn_hmask_mni2nat',
    )
    
    hmask_std2t1w.inputs.input_image = Path(config.workflow.hmask_MNI)
    hmask_std2t1w.inputs.invert_transform_flags = [True]


    # fmt: off
    workflow.connect([
        (inputnode, norm, [('moving_image', 'moving_image'),
                           ('moving_mask', 'moving_mask'),
                           ('modality', 'reference')]),
        (inputnode, tpms_std2t1w, [('moving_image', 'reference_image'), 
                                   ('tpl_tissues', 'input_image')]),
        (norm, tpms_std2t1w, [
            ('inverse_composite_transform', 'transforms'),
        ]),
        (norm, outputnode, [
            ('composite_transform', 'ind2std_xfm'),
            ('out_report', 'out_report'),
        ]),
        (tpms_std2t1w, outputnode, [('output_image', 'out_tpms')]),
    ])
    # fmt: on

    return workflow