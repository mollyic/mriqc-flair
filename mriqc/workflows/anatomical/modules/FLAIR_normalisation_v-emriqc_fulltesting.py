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


class MNItransform_AI(RegistrationRC, Registration):
    def __init__(self, *args, **kwargs):
        # Initialize parent classes
        RegistrationRC.__init__(self, *args, **kwargs)
        Registration.__init__(self, *args, **kwargs)


def ants_brain_extraction_wf(
    name='ants_brain_extraction_wf',
    use_float=True,
    normalization_quality='testing',
    omp_nthreads=None,
    mem_gb=3.0,
    use_laplacian=False,
):
    wf = pe.Workflow(name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
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
        niu.IdentityInterface(
            fields=[
                'out_tpms',
                'ind2std_xfm',
                'out_segm',
                'out_pvms',
                'hmask_mni2nat',
                'ai_out_tpms',
                'ai_ind2std_xfm',
                'ai_out_segm',
                'ai_out_pvms',
                'ai_hmask_mni2nat',
            ]
        ),
        name='outputnode',
    )

    res_tmpl = pe.Node(RegridToZooms(zooms=(4, 4, 4), smooth=True), name='res_tmpl')
    res_target = pe.Node(RegridToZooms(zooms=(4, 4, 4), smooth=True), name='res_target')

    """
    Normalisation QuickSyn (winner over  BTS with antsAI and previous method  )
    2.  QuickSyn
        - initial affine transform passed to Registration function
        - Registration with rigid, affine and syn transform
    """

    # a. RegistrationSynQuick
    syn_norm = pe.Node(
        RegistrationSynQuick(
            dimension=3,
            transform_type='b',
        ),
        name='syn_norm',
    )

    # b. tissue templates MNI to native space
    syn_tpms_mni2nat = pe.MapNode(
        ApplyTransforms(
            interpolation='Gaussian'),
        name='syn_tpms_mni2nat',
        mem_gb=1,
        iterfield=['input_image'],
    )
    syn_tpms_mni2nat.inputs.invert_transform_flags = [True]

    # c. tissue segmentation
    #syn_bts = init_brain_tissue_segmentation(name='syn_bts')

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


    wf.connect(
        [
            # resample scans at lower resolution
            # (inputnode, res_target, [('in_files', 'in_file')]),
            # (inputnode, res_tmpl, [('tpl_target_path', 'in_file')]),
            (inputnode, syn_norm, [('tpl_target_path', 'fixed_image')]),
            (inputnode, syn_norm, [('in_files', 'moving_image')]),
            # convert tissues from MNI to NORm
            (syn_norm, syn_tpms_mni2nat, [('out_matrix', 'transforms')]),
            (inputnode,syn_tpms_mni2nat, [
                ('in_files', 'reference_image'), 
                ('tissue_tpls', 'input_image')],
            ),
            # SEGMENTATION
            #(inputnode, syn_bts, [('in_files', 'inputnode.in_file')]),
            #(syn_tpms_mni2nat, syn_bts, [('output_image', 'inputnode.std_tpms')]),
            #(inputnode, syn_bts, [('in_mask', 'inputnode.brainmask')]),
            # Headmask MNI-NAT
            (inputnode, syn_hmask_mni2nat, [('in_files', 'reference_image')]),
            (syn_norm, syn_hmask_mni2nat, [('out_matrix', 'transforms')]),
            # OUTPUTS----------------------------------------------------------------------------
            (syn_tpms_mni2nat,outputnode, [('output_image', 'out_tpms')],
            ),  # Tissue tpls from MNI-NAT
            (syn_norm, outputnode, [('out_matrix', 'ind2std_xfm')]),  # registration transformation
            #(syn_bts, outputnode,
            #    [
            #        ('outputnode.out_segm', 'out_segm'),  # fsl fast segmentation
            #       ('outputnode.out_pvms', 'out_pvms'),
            #    ],
            #),  # prior probability maps
            (syn_hmask_mni2nat, outputnode, [('output_image', 'hmask_mni2nat')]),  # hmask MNI-NAT
        ]
    )


    return wf


