# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2023 The NiPreps Developers <nipreps@gmail.com>
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
# Modified by Molly Ireland on 2025-03-13
#
"""Shared workflows."""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe


def synthstrip_wf(name='synthstrip_wf', omp_nthreads=None):
    """Create a brain-extraction workflow using SynthStrip."""
    from nipype.interfaces.ants import N4BiasFieldCorrection
    from nipype.interfaces.fsl import BET
    from niworkflows.interfaces.nibabel import ApplyMask, IntensityClip

    from mriqc.interfaces.synthstrip import SynthStrip

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_files', 'bspline']), name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['out_corrected', 'out_brain', 'bias_image', 'out_mask', 'out_skin_mask']),
        name='outputnode',
    )

    # truncate target intensity for N4 correction
    pre_clip = pe.Node(IntensityClip(p_min=10, p_max=99.9), name='pre_clip')

    pre_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            num_threads=omp_nthreads,
            rescale_intensities=True,
            copy_header=True,
        ),
        name='pre_n4',
    )

    post_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            save_bias=True,
            num_threads=omp_nthreads,
            n_iterations=[50] * 4,
            copy_header=True,
        ),
        name='post_n4',
    )

    synthstrip = pe.Node(
        SynthStrip(num_threads=omp_nthreads),
        name='synthstrip',
        num_threads=omp_nthreads,
    )

    betted_skin = pe.Node(BET(frac =0.2), name='betted_skin', num_threads=omp_nthreads)
    betted_skin.inputs.args = '-A'
    betted_skin.inputs.surfaces = True

    final_masked = pe.Node(ApplyMask(), name='final_masked')

    workflow = pe.Workflow(name=name)
    # fmt: off
    workflow.connect([
        (inputnode, pre_clip, [('in_files', 'in_file')]),
        (inputnode, pre_n4, [('bspline', 'bspline_fitting_distance')]),
        (inputnode, post_n4, [('bspline', 'bspline_fitting_distance')]),
        (pre_clip, pre_n4, [('out_file', 'input_image')]),
        (pre_n4, synthstrip, [('output_image', 'in_file')]),
        (synthstrip, post_n4, [('out_mask', 'weight_image')]),
        (synthstrip, final_masked, [('out_mask', 'in_mask')]),
        (pre_clip, post_n4, [('out_file', 'input_image')]),
        (post_n4, final_masked, [('output_image', 'in_file')]),
        (final_masked, outputnode, [('out_file', 'out_brain')]),
        (post_n4, outputnode, [('bias_image', 'bias_image')]),
        (synthstrip, outputnode, [('out_mask', 'out_mask')]),
        (post_n4, outputnode, [('output_image', 'out_corrected')]),
        (pre_n4, betted_skin, [('output_image', 'in_file')]),
        (betted_skin, outputnode, [('outskin_mask_file', 'out_skin_mask')]),
    ])
    # fmt: on
    return workflow
