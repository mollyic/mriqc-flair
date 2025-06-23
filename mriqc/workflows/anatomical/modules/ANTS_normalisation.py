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
"""Nipype translation of ANTs' workflows."""

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
    from nipype.interfaces.ants import RegistrationSynQuick

    # Define workflow interface
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=[
                'in_files',
                'in_mask',
                'tissue_tpls',
                'tpl_target_path',
                'tpl_mask_path',
            ]
        ),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['out_tpms', 'ind2std_xfm','hmask_mni2nat']),
        name='outputnode',
    )

    # Spatial normalization
    norm = pe.Node(
        RegistrationSynQuick(
            dimension=3,
            transform_type='b',
        ),
        name='SpatialNormalization',
        # Request all MultiProc processes when ants_nthreads > n_procs
        num_threads=config.nipype.omp_nthreads,
        mem_gb=3,
    )

    def _get_templates(modality):
        if config.workflow.species.lower() == 'human':
            tpd_id = 'GG853' if modality == 'FLAIR' else config.workflow.template_id
    
    if config.workflow.species.lower() == 'human':
        #template function here?
        template = 'find me'


    # Project standard TPMs into Native space
    tpms_std2ntv = pe.MapNode(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            interpolation='Gaussian',
            float=config.execution.ants_float,
        ),
        iterfield=['input_image'],
        name='tpms_std2ntv',
    )
    tpms_std2ntv.inputs.invert_transform_flags = [True]

    tpms_std2ntv.inputs.input_image = [
        str(p)
        for p in get_template(
            config.workflow.template_id,
            suffix='probseg',
            resolution=(1 if config.workflow.species.lower() == 'human' else None),
            label=['CSF', 'GM', 'WM'],
        )
    ]
    # d. MNI hmask to native space
    syn_hmask_mni2nat = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            interpolation='Gaussian',  # Choose the appropriate interpolation method
            float=config.execution.ants_float,
        ),
        name='syn_hmask_mni2nat',
    )
    syn_hmask_mni2nat.inputs.input_image = config.workflow.hmask_MNI
    syn_hmask_mni2nat.inputs.invert_transform_flags = [True]


    # fmt: off
    workflow.connect([
        (inputnode, norm, [('tpl_target_path', 'fixed_image'),
                           ('in_files', 'moving_image')]),
            # convert tissues from MNI to NORm
        (norm, tpms_std2ntv, [('out_matrix', 'transforms')]),
        (inputnode, tpms_std2ntv,[('in_files', 'reference_image'), 
                                  ('tissue_tpls', 'input_image')]),
            # SEGMENTATION
            #(inputnode, syn_bts, [('in_files', 'inputnode.in_file')]),
            #(tpms_std2ntv, syn_bts, [('output_image', 'inputnode.std_tpms')]),
            #(inputnode, syn_bts, [('in_mask', 'inputnode.brainmask')]),
            # Headmask MNI-NAT
            (inputnode, syn_hmask_mni2nat, [('in_files', 'reference_image')]),
            (norm, syn_hmask_mni2nat, [('out_matrix', 'transforms')]),
            # OUTPUTS----------------------------------------------------------------------------
            (tpms_std2ntv,outputnode, [('output_image', 'out_tpms')],
            ),  # Tissue tpls from MNI-NAT
            (norm, outputnode, [('out_matrix', 'ind2std_xfm')]),  # registration transformation
            #(syn_bts, outputnode,
            #    [
            #        ('outputnode.out_segm', 'out_segm'),  # fsl fast segmentation
            #       ('outputnode.out_pvms', 'out_pvms'),
            #    ],
            #),  # prior probability maps
            (syn_hmask_mni2nat, outputnode, [('output_image', 'hmask_mni2nat')]),  # hmask MNI-NAT
        ]
    )


    return workflow


