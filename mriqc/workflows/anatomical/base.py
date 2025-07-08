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
# Modified by Molly Ireland
#
"""
Anatomical workflow
===================

.. image :: _static/anatomical_workflow_source.svg

The anatomical workflow follows the following steps:

#. Conform (reorientations, revise data types) input data and read
   associated metadata.
#. Skull-stripping (AFNI).
#. Calculate head mask -- :py:func:`headmsk_wf`.
#. Spatial Normalization to MNI (ANTs)
#. Calculate air mask above the nasial-cerebelum plane -- :py:func:`airmsk_wf`.
#. Brain tissue segmentation (FAST).
#. Extraction of IQMs -- :py:func:`compute_iqms`.
#. Individual-reports generation --
   :py:func:`~mriqc.workflows.anatomical.output.init_anat_report_wf`.

This workflow is orchestrated by :py:func:`anat_qc_workflow`.

For the skull-stripping, we use ``afni_wf`` from ``niworkflows.anat.skullstrip``:

.. workflow::

    from niworkflows.anat.skullstrip import afni_wf
    from mriqc.testing import mock_config
    with mock_config():
        wf = afni_wf()

"""
from mriqc import config
from mriqc.interfaces import (
    ArtifactMask,
    ComputeQI2,
    ConformImage,
    IQMFileSink,
    RotationMask,
    HeadMask_review,
    StructuralQC,
)
from mriqc.interfaces.reports import AddProvenance
from mriqc.interfaces.datalad import DataladIdentityInterface
from mriqc.messages import BUILDING_WORKFLOW
from mriqc.workflows.utils import get_fwhmx
from mriqc.workflows.anatomical.output import init_anat_report_wf
from nipype.interfaces import utility as niu
from mriqc.workflows.anatomical.flair_modules.clean_segs_ants import clean_tissue_segs
from nipype.pipeline import engine as pe

from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from templateflow.api import get as get_template


def anat_qc_workflow(name="anatMRIQC"):
    """
    One-subject-one-session-one-run pipeline to extract the NR-IQMs from
    anatomical images

    .. workflow::

        import os.path as op
        from mriqc.workflows.anatomical.base import anat_qc_workflow
        from mriqc.testing import mock_config
        with mock_config():
            wf = anat_qc_workflow()

    """

    dataset = config.workflow.inputs.get("T1w", []) + config.workflow.inputs.get("T2w", []) + config.workflow.inputs.get("FLAIR", [])

    message = BUILDING_WORKFLOW.format(
        modality="anatomical",
        detail=(
            f"for {len(dataset)} NIfTI files."
            if len(dataset) > 2
            else f"({' and '.join(('<%s>' % v for v in dataset))})."
        ),
    )
    config.loggers.workflow.info(message)

    # Initialize workflow
    workflow = pe.Workflow(name=name)

    # Define workflow, inputs and outputs
    # 0. Get data
    inputnode = pe.Node(niu.IdentityInterface(fields=["in_file"]), name="inputnode")
    inputnode.iterables = [("in_file", dataset)]

    get_info = pe.Node(
            niu.Function(
                function = _get_info, 
                input_names=["in_file"],
                output_names=[
                    'modality', 
                    'bspline', 
                    'tpl_target_path', 
                    'tpl_mask_path',
                    'wm_tpl', 
                    'tissue_tpls', 
                    'likelihood_model'
                    ],
                ), 
            name = 'get_info')
    datalad_get = pe.Node(
        DataladIdentityInterface(fields=["in_file"], dataset_path=config.execution.bids_dir),
        name="datalad_get",
    )

    outputnode = pe.Node(niu.IdentityInterface(fields=["out_json"]), name="outputnode")

    # 1. Reorient anatomical image
    to_ras = pe.Node(ConformImage(check_dtype=False), name="conform")
    # 2. species specific skull-stripping
    if config.workflow.species.lower() == "human":
        skull_stripping = synthstrip_wf(omp_nthreads=config.nipype.omp_nthreads)
        ss_bias_field = "outputnode.bias_image"
    else:
        from nirodents.workflows.brainextraction import init_rodent_brain_extraction_wf

        skull_stripping = init_rodent_brain_extraction_wf(template_id=config.workflow.template_id)
        ss_bias_field = "final_n4.bias_image"
    # 3. Head mask
    hmsk = headmsk_wf(omp_nthreads=config.nipype.omp_nthreads)
    # 4. Spatial Normalization, using ANTs
    norm = spatial_normalization()
    # 5. Air mask (with and without artifacts)
    amw = airmsk_wf()
    # 6. Brain tissue segmentation
    bts = init_brain_tissue_segmentation()
    # 7. Compute IQMs
    iqmswf = compute_iqms()
    # Reports
    anat_report_wf = init_anat_report_wf()

    # Additional morphological changess
    clean_segs = clean_tissue_segs()
    # Connect all nodes
    # fmt: off
    workflow.connect([
        (inputnode, get_info, [("in_file", "in_file")]),
        (get_info, skull_stripping, [('bspline', "inputnode.bspline")]),
        (get_info, hmsk, [('modality', "inputnode.modality")]),
        (get_info, clean_segs, [('modality', "inputnode.modality")]),
        (get_info, iqmswf, [("modality", "inputnode.modality")]),
        (get_info, bts, [("likelihood_model", "inputnode.likelihood_model")]),
        (inputnode, datalad_get, [("in_file", "in_file")]),
        (inputnode, anat_report_wf, [
            ("in_file", "inputnode.name_source"),
        ]),
        (datalad_get, to_ras, [("in_file", "in_file")]),
        (datalad_get, iqmswf, [("in_file", "inputnode.in_file")]),
        (get_info, norm, [
            ('tpl_target_path', "inputnode.tpl_target_path"), 
            ('tpl_mask_path', "inputnode.tpl_mask_path"), 
            ('tissue_tpls', "inputnode.tissue_tpls")]),
        (to_ras, skull_stripping, [("out_file", "inputnode.in_files")]),
        (skull_stripping, hmsk, [
            ("outputnode.out_skin_mask", "inputnode.skinmask"),
            ("outputnode.out_corrected", "inputnode.in_file"),
            ("outputnode.out_mask", "inputnode.brainmask"),
        ]),
        (skull_stripping, bts, [("outputnode.out_mask", "inputnode.brainmask")]),
        (skull_stripping, norm, [
            ("outputnode.out_corrected", "inputnode.moving_image"),
            ("outputnode.out_mask", "inputnode.moving_mask")]),
        (norm, bts, [("outputnode.out_tpms", "inputnode.std_tpms")]),
        (bts, clean_segs, [
            ("outputnode.out_segm", "inputnode.segmentation"), 
            ("outputnode.out_pvms", "inputnode.pvms")]),
        (norm, amw, [
            ("outputnode.ind2std_xfm", "inputnode.ind2std_xfm")]),
        (norm, iqmswf, [
            ("outputnode.out_tpms", "inputnode.std_tpms")]),
        (norm, anat_report_wf, ([
            ("outputnode.out_report", "inputnode.mni_report")])),
        (norm, hmsk, [("outputnode.out_tpms", "inputnode.in_tpms"), 
                      ("outputnode.ind2std_xfm", "inputnode.ind2std_xfm"), 
                      ("outputnode.hmask_mni2nat", "inputnode.mask_tmpl")]),
        (to_ras, amw, [("out_file", "inputnode.in_file")]),
        (skull_stripping, amw, [("outputnode.out_mask", "inputnode.in_mask")]),
        (hmsk, amw, [("outputnode.out_file", "inputnode.head_mask")]),
        (to_ras, iqmswf, [("out_file", "inputnode.in_ras")]),
        (skull_stripping, iqmswf, [("outputnode.out_corrected", "inputnode.inu_corrected"),
                                   (ss_bias_field, "inputnode.in_inu"),
                                   ("outputnode.out_mask", "inputnode.brainmask")]),
        (amw, iqmswf, [("outputnode.air_mask", "inputnode.airmask"),
                       ("outputnode.hat_mask", "inputnode.hatmask"),
                       ("outputnode.art_mask", "inputnode.artmask"),
                       ("outputnode.rot_mask", "inputnode.rotmask")]),
        (hmsk, bts, [("outputnode.out_denoised", "inputnode.in_file")]),
        (clean_segs, iqmswf, [("outputnode.out_segm", "inputnode.segmentation"),
                       ("outputnode.out_pvms", "inputnode.pvms")]),
        (hmsk, iqmswf, [("outputnode.out_file", "inputnode.headmask")]),
        (to_ras, anat_report_wf, [("out_file", "inputnode.in_ras")]),
        (skull_stripping, anat_report_wf, [
            ("outputnode.out_corrected", "inputnode.inu_corrected"),
            ("outputnode.out_mask", "inputnode.brainmask")]),
        (hmsk, anat_report_wf, [("outputnode.out_file", "inputnode.headmask")]),
        (amw, anat_report_wf, [
            ("outputnode.air_mask", "inputnode.airmask"),
            ("outputnode.art_mask", "inputnode.artmask"),
            ("outputnode.rot_mask", "inputnode.rotmask"),
        ]),
        (bts, anat_report_wf, [("outputnode.out_segm", "inputnode.segmentation")]),
        (iqmswf, anat_report_wf, [("outputnode.noisefit", "inputnode.noisefit")]),
        (iqmswf, anat_report_wf, [("outputnode.out_file", "inputnode.in_iqms")]),
        (iqmswf, outputnode, [("outputnode.out_file", "out_json")]),
    ])
    # fmt: on

    # Upload metrics
    if not config.execution.no_sub:
        from mriqc.interfaces.webapi import UploadIQMs

        upldwf = pe.Node(
            UploadIQMs(
                endpoint=config.execution.webapi_url,
                auth_token=config.execution.webapi_token,
                strict=config.execution.upload_strict,
            ),
            name="UploadMetrics",
        )

        # fmt: off
        workflow.connect([
            (iqmswf, upldwf, [("outputnode.out_file", "in_iqms")]),
            (upldwf, anat_report_wf, [("api_id", "inputnode.api_id")]),
        ])
        # fmt: on

    return workflow


def spatial_normalization(name="SpatialNormalization"):
    """Create a simplied workflow to perform spatial normalization with Ants QuickSyn."""
    from mriqc.workflows.anatomical.flair_modules.normalisation import (
        RegistrationSynQuickRPT as RobustMNINormalization
    )

    from templateflow.api import get as get_template
    # Define workflow interface
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=[
            "moving_image", 
            'moving_mask', 
            'tissue_tpls',
            'tpl_target_path',
            'tpl_mask_path',]), 
            name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["out_tpms", "out_report", "ind2std_xfm", 'hmask_mni2nat']),
        name="outputnode",
    )

    # Spatial normalization
    norm = pe.Node(
        RobustMNINormalization(
            dimension=3,
            num_threads=config.nipype.omp_nthreads,
            transform_type = "b",
            generate_report=True,
        ),
        name="SpatialNormalization",
        # Request all MultiProc processes when ants_nthreads > n_procs
        num_threads=config.nipype.omp_nthreads,
        mem_gb=3,
    )

    # Project standard TPMs into T1w space
    tpms_std2t1w = pe.MapNode(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            interpolation="Linear",
            float=config.execution.ants_float,
        ),
        iterfield=["input_image"],
        name="tpms_std2t1w",
    )
    tpms_std2t1w.inputs.invert_transform_flags = [True]

    # Project standard headmask into native space
    hmask_mni2nat = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            interpolation="Linear",
            float=config.execution.ants_float,
        ),
        name="hmask_mni2nat"
        )

    hmask_mni2nat.inputs.input_image = get_template(
        template='MNI152Lin',
        desc='head',
        suffix='mask',
        resolution ='2'
        )
    hmask_mni2nat.inputs.invert_transform_flags = [True]

    # fmt: off
    workflow.connect([
        (inputnode, norm, [("moving_image", "moving_image"), 
                           ("tpl_target_path", "fixed_image")]),
        (inputnode, tpms_std2t1w, [("moving_image", "reference_image"),
                                   ("tissue_tpls", "input_image")]),
        (norm, tpms_std2t1w, [
            ("out_matrix", "transforms")
        ]),
        (norm, outputnode, [
            ("out_matrix", "ind2std_xfm"),
            ("out_report", "out_report"),
        ]),
        (tpms_std2t1w, outputnode, [("output_image", "out_tpms")]),
        (inputnode, hmask_mni2nat, [("moving_image", "reference_image")]),
        (norm, hmask_mni2nat, [("out_matrix", "transforms")]),
        (hmask_mni2nat, outputnode, [("output_image", "hmask_mni2nat")]), 
    ])
    # fmt: on

    return workflow


def init_brain_tissue_segmentation(name="brain_tissue_segmentation"):
    """
    Setup a workflow for brain tissue segmentation.

    .. workflow::

        from mriqc.workflows.anatomical.base import init_brain_tissue_segmentation
        from mriqc.testing import mock_config
        with mock_config():
            wf = init_brain_tissue_segmentation()

    """
    from nipype.interfaces.ants import Atropos

    def _format_tpm_names(in_files, fname_string=None):
        from pathlib import Path
        import nibabel as nb
        import glob

        out_path = Path.cwd().absolute()

        # copy files to cwd and rename iteratively
        for count, fname in enumerate(in_files):
            img = nb.load(fname)
            extension = "".join(Path(fname).suffixes)
            out_fname = f"priors_{1 + count:02}{extension}"
            nb.save(img, Path(out_path, out_fname))

        if fname_string is None:
            fname_string = f"priors_%02d{extension}"

        out_files = [str(prior) for prior in glob.glob(str(Path(out_path, f"priors*{extension}")))]

        # return path with c-style format string for Atropos
        file_format = str(Path(out_path, fname_string))
        return file_format, out_files

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["in_file", "brainmask", "std_tpms", 'likelihood_model']),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["out_segm", "out_pvms"]),
        name="outputnode",
    )

    format_tpm_names = pe.Node(
        niu.Function(
            input_names=["in_files"],
            output_names=["file_format"],
            function=_format_tpm_names,
            execution={"keep_inputs": True, "remove_unnecessary_outputs": False},
        ),
        name="format_tpm_names",
    )

    segment = pe.Node(
        Atropos(
            initialization="PriorProbabilityImages",
            number_of_tissue_classes=3,
            prior_weighting=0.1,
            mrf_radius=[1, 1, 1],
            mrf_smoothing_factor=0.01,
            save_posteriors=True,
            out_classified_image_name="segment.nii.gz",
            output_posteriors_name_template="segment_%02d.nii.gz",
            num_threads=config.nipype.omp_nthreads,
        ),
        name="segmentation",
        mem_gb=5,
        num_threads=config.nipype.omp_nthreads,
    )

    # fmt: off
    workflow.connect([
        (inputnode, segment, [("in_file", "intensity_images"),
                              ("brainmask", "mask_image"), 
                              ('likelihood_model', 'likelihood_model')]),
        (inputnode, format_tpm_names, [('std_tpms', 'in_files')]),
        (format_tpm_names, segment, [(('file_format', _pop), 'prior_image')]),
        (segment, outputnode, [("classified_image", "out_segm"),
                               ("posteriors", "out_pvms")]),
    ])
    # fmt: on
    return workflow


def compute_iqms(name="ComputeIQMs"):
    """
    Setup the workflow that actually computes the IQMs.

    .. workflow::

        from mriqc.workflows.anatomical.base import compute_iqms
        from mriqc.testing import mock_config
        with mock_config():
            wf = compute_iqms()

    """
    from niworkflows.interfaces.bids import ReadSidecarJSON

    from mriqc.interfaces.anatomical import Harmonize
    from mriqc.workflows.utils import _tofloat

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "in_file",
                "modality",
                "in_ras",
                "brainmask",
                "airmask",
                "artmask",
                "headmask",
                "rotmask",
                "hatmask",
                "segmentation",
                "inu_corrected",
                "in_inu",
                "pvms",
                "metadata",
                "std_tpms",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["out_file", "noisefit"]),
        name="outputnode",
    )

    # Extract metadata
    meta = pe.Node(ReadSidecarJSON(index_db=config.execution.bids_database_dir), name="metadata")

    # Add provenance
    addprov = pe.Node(AddProvenance(), name="provenance", run_without_submitting=True)

    # AFNI check smoothing
    fwhm_interface = get_fwhmx()

    fwhm = pe.Node(fwhm_interface, name="smoothness")

    # Harmonize
    homog = pe.Node(Harmonize(), name="harmonize")
    if config.workflow.species.lower() != "human":
        homog.inputs.erodemsk = False
        homog.inputs.thresh = 0.8

    # Mortamet's QI2
    getqi2 = pe.Node(ComputeQI2(), name="ComputeQI2")

    # Compute python-coded measures
    measures = pe.Node(StructuralQC(human=config.workflow.species.lower() == "human"), "measures")

    datasink = pe.Node(
        IQMFileSink(
            out_dir=config.execution.output_dir,
            dataset=config.execution.dsname,
        ),
        name="datasink",
        run_without_submitting=True,
    )

    def _getwm(inlist):
        return inlist[-1]

    # fmt: off
    workflow.connect([
        (inputnode, meta, [("in_file", "in_file")]),
        (inputnode, datasink, [("in_file", "in_file"),
                              ('modality', "modality")]),
        (inputnode, addprov, [('modality', "modality")]),
        (meta, datasink, [("subject", "subject_id"),
                          ("session", "session_id"),
                          ("task", "task_id"),
                          ("acquisition", "acq_id"),
                          ("reconstruction", "rec_id"),
                          ("run", "run_id"),
                          ("out_dict", "metadata")]),
        (inputnode, addprov, [("in_file", "in_file"),
                              ("airmask", "air_msk"),
                              ("rotmask", "rot_msk")]),
        (inputnode, getqi2, [("in_ras", "in_file"),
                             ("hatmask", "air_msk")]),
        (inputnode, homog, [("inu_corrected", "in_file"),
                            (("pvms", _getwm), "wm_mask"), 
                            ('modality', "modality")]),
        (inputnode, measures, [("in_inu", "in_bias"),
                               ("in_ras", "in_file"),
                               ("airmask", "air_msk"),
                               ("headmask", "head_msk"),
                               ("artmask", "artifact_msk"),
                               ("rotmask", "rot_msk"),
                               ("segmentation", "in_segm"),
                               ("pvms", "in_pvms"),
                               ("std_tpms", "mni_tpms")]),
        (inputnode, fwhm, [("in_ras", "in_file"),
                           ("brainmask", "mask")]),
        (homog, measures, [("out_file", "in_noinu")]),
        (fwhm, measures, [(("fwhm", _tofloat), "in_fwhm")]),
        (measures, datasink, [("out_qc", "root")]),
        (addprov, datasink, [("out_prov", "provenance")]),
        (getqi2, datasink, [("qi2", "qi_2")]),
        (getqi2, outputnode, [("out_file", "noisefit")]),
        (datasink, outputnode, [("out_file", "out_file")]),
    ])
    # fmt: on

    return workflow


def headmsk_wf(name="HeadMaskWorkflow", omp_nthreads=1):
    """
    Computes a head mask as in [Mortamet2009]_.

    .. workflow::

        from mriqc.testing import mock_config
        from mriqc.workflows.anatomical.base import headmsk_wf
        with mock_config():
            wf = headmsk_wf()

    """

    from niworkflows.interfaces.nibabel import ApplyMask

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["in_file", "brainmask", "skinmask", 'modality', 
                                      "in_tpms", "mask_tmpl", "ind2std_xfm"]), name="inputnode"
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["out_file", "out_denoised"]), name="outputnode"
    )

    def _select_wm(inlist):
        return [f for f in inlist if "WM" in f][0]

    enhance = pe.Node(
        niu.Function(
            input_names=["in_file", "wm_tpm", "modality"],
            output_names=["out_file"],
            function=_enhance,
        ),
        name="Enhance",
        num_threads=omp_nthreads,
    )

    gradient = pe.Node(
        niu.Function(
            input_names=["in_file", "brainmask", "sigma"],
            output_names=["out_file"],
            function=image_gradient,
        ),
        name="Grad",
        num_threads=omp_nthreads,
    )
    thresh = pe.Node(
        niu.Function(
            input_names=["in_file", "brainmask", "aniso", "thresh"],
            output_names=["out_file"],
            function=gradient_threshold,
        ),
        name="GradientThreshold",
        num_threads=omp_nthreads,
    )
    if config.workflow.species != "human":
        gradient.inputs.sigma = 3.0
        thresh.inputs.aniso = True
        thresh.inputs.thresh = 4.0

    apply_mask = pe.Node(ApplyMask(), name="apply_mask")
    review = pe.Node(
        HeadMask_review(),
        name = "ReviewMask",            
        execution={"keep_inputs": True, "remove_unnecessary_outputs": False}
        )

    # fmt: off
    workflow.connect([
        (inputnode, enhance, [("in_file", "in_file"),
                              ('modality', "modality"),
                              (("in_tpms", _select_wm), "wm_tpm")]),
        (inputnode, thresh, [("skinmask", "brainmask")
                             ('modality', "modality"),]),
        (inputnode, gradient, [("skinmask", "brainmask")]),
        (inputnode, apply_mask, [("brainmask", "in_mask")]),
        (enhance, gradient, [("out_file", "in_file")]),
        (gradient, thresh, [("out_file", "in_file")]),
        (enhance, apply_mask, [("out_file", "in_file")]),
        (thresh, review, [("out_file", "hmask")]),
        (inputnode, review, [("mask_tmpl", "hmask_mni2nat"),
                             ("ind2std_xfm", "ind2std_xfm")]),   
        (review, outputnode, [("out_file", "out_file"),]),
        (apply_mask, outputnode,[("out_file", "out_denoised")]),
    ])
    # fmt: on

    return workflow


def airmsk_wf(name="AirMaskWorkflow"):
    """
    Calculate air, artifacts and "hat" masks to evaluate noise in the background.

    This workflow mostly addresses the implementation of Step 1 in [Mortamet2009]_.
    This work proposes to look at the signal distribution in the background, where
    no signals are expected, to evaluate the spread of the noise.
    It is in the background where [Mortamet2009]_ proposed to also look at the presence
    of ghosts and artifacts, where they are very easy to isolate.

    However, [Mortamet2009]_ proposes not to look at the background around the face
    because of the likely signal leakage through the phase-encoding axis sourcing from
    eyeballs (and their motion).
    To avoid that, [Mortamet2009]_ proposed atlas-based identification of two landmarks
    (nasion and cerebellar projection on to the occipital bone).
    MRIQC, for simplicity, used a such a mask created in MNI152NLin2009cAsym space and
    projected it on to the individual.
    Such a solution is inadequate because it doesn't drop full in-plane slices as there
    will be a large rotation of the individual's tilt of the head with respect to the
    template.
    The new implementation (23.1.x series) follows [Mortamet2009]_ more closely,
    projecting the two landmarks from the template space and leveraging
    *NiTransforms* to do that.

    .. workflow::

        from mriqc.testing import mock_config
        from mriqc.workflows.anatomical.base import airmsk_wf
        with mock_config():
            wf = airmsk_wf()

    """
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "in_file",
                "in_mask",
                "head_mask",
                "ind2std_xfm",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["hat_mask", "air_mask", "art_mask", "rot_mask"]),
        name="outputnode",
    )

    rotmsk = pe.Node(RotationMask(), name="RotationMask")
    qi1 = pe.Node(ArtifactMask(), name="ArtifactMask")

    # fmt: off
    workflow.connect([
        (inputnode, rotmsk, [("in_file", "in_file")]),
        (inputnode, qi1, [("in_file", "in_file"),
                          ("head_mask", "head_mask"),
                          ("ind2std_xfm", "ind2std_xfm")]),
        (qi1, outputnode, [("out_hat_msk", "hat_mask"),
                           ("out_air_msk", "air_mask"),
                           ("out_art_msk", "art_mask")]),
        (rotmsk, outputnode, [("out_file", "rot_mask")])
    ])
    # fmt: on

    return workflow


def synthstrip_wf(name="synthstrip_wf", omp_nthreads=None):
    """Create a brain-extraction workflow using SynthStrip."""
    from nipype.interfaces.fsl import BET
    from nipype.interfaces.ants import N4BiasFieldCorrection
    from niworkflows.interfaces.nibabel import IntensityClip, ApplyMask
    from mriqc.interfaces.synthstrip import SynthStrip

    inputnode = pe.Node(niu.IdentityInterface(fields=["in_files", 'bspline']), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["out_corrected", "out_brain", "bias_image", 
                                      "out_mask", "out_skin_mask"]),
        name="outputnode",
    )

    # truncate target intensity for N4 correction
    pre_clip = pe.Node(IntensityClip(p_min=10, p_max=99.9), name="pre_clip")

    pre_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            num_threads=omp_nthreads,
            rescale_intensities=True,
            copy_header=True,
        ),
        name="pre_n4",
    )

    post_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            save_bias=True,
            num_threads=omp_nthreads,
            n_iterations=[50] * 4,
            copy_header=True,
        ),
        name="post_n4",
    )

    synthstrip = pe.Node(
        SynthStrip(num_threads=omp_nthreads),
        name="synthstrip",
        num_threads=omp_nthreads,
    )

    betted_skin = pe.Node(
        BET(frac =0.2),
        name="betted_skin",
        num_threads=omp_nthreads)
    betted_skin.inputs.args = '-A'
    betted_skin.inputs.surfaces = True
    final_masked = pe.Node(ApplyMask(), name="final_masked")
    final_inu = pe.Node(niu.Function(function=_apply_bias_correction), name="final_inu")

    workflow = pe.Workflow(name=name)
    # fmt: off
    workflow.connect([
        (inputnode, final_inu, [("in_files", "in_file")]),
        (inputnode, pre_n4, [('bspline', "bspline_fitting_distance")]),
        (inputnode, post_n4, [('bspline', "bspline_fitting_distance")]),
        (inputnode, pre_clip, [("in_files", "in_file")]),
        (pre_clip, pre_n4, [("out_file", "input_image")]),
        (pre_n4, synthstrip, [("output_image", "in_file")]),
        (synthstrip, post_n4, [("out_mask", "weight_image")]),
        (synthstrip, final_masked, [("out_mask", "in_mask")]),
        (pre_clip, post_n4, [("out_file", "input_image")]),
        (post_n4, final_inu, [("bias_image", "bias_image")]),
        (post_n4, final_masked, [("output_image", "in_file")]),
        (final_masked, outputnode, [("out_file", "out_brain")]),
        (post_n4, outputnode, [("bias_image", "bias_image")]),
        (synthstrip, outputnode, [("out_mask", "out_mask")]),
        (post_n4, outputnode, [("output_image", "out_corrected")]),
        (pre_n4, betted_skin,       [("output_image", "in_file")]),
        (betted_skin, outputnode,   [("outskin_mask_file", "out_skin_mask")]),
    ])
    # fmt: on
    return workflow


def _apply_bias_correction(in_file, bias_image, out_file=None):
    import os.path as op

    import numpy as np
    import nibabel as nb

    img = nb.load(in_file)
    data = np.clip(
        img.get_fdata() * nb.load(bias_image).get_fdata(),
        a_min=0,
        a_max=None,
    )
    out_img = img.__class__(
        data.astype(img.get_data_dtype()),
        img.affine,
        img.header,
    )

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_file))
        if ext == ".gz":
            fname, ext2 = op.splitext(fname)
            ext = ext2 + ext
        out_file = op.abspath("{}_inu{}".format(fname, ext))

    out_img.to_filename(out_file)
    return out_file


def _binarize(in_file, threshold=0.5, out_file=None):
    import os.path as op

    import nibabel as nb
    import numpy as np

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_file))
        if ext == ".gz":
            fname, ext2 = op.splitext(fname)
            ext = ext2 + ext
        out_file = op.abspath("{}_bin{}".format(fname, ext))

    nii = nb.load(in_file)
    data = nii.get_fdata() > threshold

    hdr = nii.header.copy()
    hdr.set_data_dtype(np.uint8)
    nb.Nifti1Image(data.astype(np.uint8), nii.affine, hdr).to_filename(out_file)
    return out_file


def _enhance(in_file, wm_tpm, modality, out_file=None):
    import numpy as np
    import nibabel as nb
    from mriqc.workflows.utils import generate_filename

    imnii = nb.load(in_file)
    data = imnii.get_fdata(dtype=np.float32)
    # Increase suppression of high-intensity WM voxels: improving FLAIR tissue segmentation
    thresh =  99.95 if modality.lower() == "flair" else 99.98
    range_max = np.percentile(data[data > 0], thresh)
    excess = data > range_max

    wm_prob = nb.load(wm_tpm).get_fdata()
    wm_prob[wm_prob < 0] = 0  # Ensure no negative values
    wm_prob[excess] = 0  # Ensure no outliers are considered

    # Calculate weighted mean and standard deviation
    wm_mu = np.average(data, weights=wm_prob)
    wm_sigma = np.sqrt(np.average((data - wm_mu) ** 2, weights=wm_prob))

    # Resample signal excess pixels
    data[excess] = np.random.normal(loc=wm_mu, scale=wm_sigma, size=excess.sum())

    out_file = out_file or str(generate_filename(in_file, suffix="enhanced").absolute())
    nb.Nifti1Image(data, imnii.affine, imnii.header).to_filename(out_file)
    return out_file


def image_gradient(in_file, brainmask, sigma=4.0, out_file=None):
    """Computes the magnitude gradient of an image using numpy"""
    import nibabel as nb
    import numpy as np
    from scipy.ndimage import gaussian_gradient_magnitude as gradient
    from mriqc.workflows.utils import generate_filename

    imnii = nb.load(in_file)
    mask = np.bool_(nb.load(brainmask).dataobj)
    data = imnii.get_fdata(dtype=np.float32)
    datamax = np.percentile(data.reshape(-1), 99.5)
    data *= 100 / datamax
    data[mask] = 100

    zooms = np.array(imnii.header.get_zooms()[:3])
    sigma_xyz = 2 - zooms / min(zooms)
    grad = gradient(data, sigma * sigma_xyz)
    gradmax = np.percentile(grad.reshape(-1), 99.5)
    grad *= 100.0
    grad /= gradmax
    grad[mask] = 100

    out_file = out_file or str(generate_filename(in_file, suffix="grad").absolute())
    nb.Nifti1Image(grad, imnii.affine, imnii.header).to_filename(out_file)
    return out_file


def gradient_threshold(in_file, brainmask, modality, thresh=15.0, out_file=None, aniso=False):
    """Compute a threshold from the histogram of the magnitude gradient image"""
    import nibabel as nb
    import numpy as np
    from scipy import ndimage as sim
    from mriqc.workflows.utils import generate_filename

    if not aniso:
        struct = sim.iterate_structure(sim.generate_binary_structure(3, 2), 2)
    else:
        # Generate an anisotropic binary structure, taking into account slice thickness
        img = nb.load(in_file)
        zooms = img.header.get_zooms()
        dist = max(zooms)
        dim = img.header["dim"][0]

        x = np.ones((5) * np.ones(dim, dtype=np.int8))
        np.put(x, x.size // 2, 0)
        dist_matrix = np.round(sim.distance_transform_edt(x, sampling=zooms), 5)
        struct = dist_matrix <= dist

    imnii = nb.load(in_file)

    hdr = imnii.header.copy()
    hdr.set_data_dtype(np.uint8)

    data = imnii.get_fdata(dtype=np.float32)

    mask = np.zeros_like(data, dtype=np.uint8)
    thresh = np.percentile(data[data != 0], 20) if modality == 'FLAIR' else thresh
    mask[data > thresh] = 1

    #Add padding
    pad = 100
    bmask = nb.load(brainmask).get_fdata()
    bmask = np.pad(bmask, pad, mode='constant', constant_values=0)
    mask = np.pad(mask, pad, mode='constant', constant_values=0)
    mask = sim.binary_closing(mask, struct, iterations=2).astype(np.uint8)
    mask = sim.binary_erosion(mask, sim.generate_binary_structure(3, 2)).astype(np.uint8)

    segdata = np.asanyarray(bmask) > 0
    segdata = sim.binary_dilation(segdata, struct, iterations=2, border_value=1).astype(np.uint8)
    mask[segdata] = 1

    # Remove small objects
    label_im, nb_labels = sim.label(mask)
    artmsk = np.zeros_like(mask)
    if nb_labels > 2:
        sizes = sim.sum(mask, label_im, list(range(nb_labels + 1)))
        ordered = list(reversed(sorted(zip(sizes, list(range(nb_labels + 1))))))
        for _, label in ordered[2:]:
            mask[label_im == label] = 0
            artmsk[label_im == label] = 1

    mask = sim.binary_fill_holes(mask, struct).astype(np.uint8)  # pylint: disable=no-member
    mask = mask[pad:-pad, pad:-pad, pad:-pad] #remove padding
    out_file = out_file or str(generate_filename(in_file, suffix="gradmask").absolute())
    nb.Nifti1Image(mask, imnii.affine, hdr).to_filename(out_file)
    return out_file


def _get_imgtype(in_file):
    from pathlib import Path

    return int(Path(in_file).name.rstrip(".gz").rstrip(".nii").split("_")[-1][1])


def _get_mod(in_file):
    from pathlib import Path

    return Path(in_file).name.rstrip(".gz").rstrip(".nii").split("_")[-1]


def _pop(inlist):
    if isinstance(inlist, (list, tuple)):
        return inlist[0]
    return inlist

def _get_info(in_file):
    from mriqc import config
    from templateflow.api import get as get_template
    from niworkflows.utils.misc import get_template_specs
    from mriqc.workflows.anatomical.base import _get_mod
    from mriqc.workflows.anatomical.flair_modules.normalisation import _get_custom_templates

    """
    Extracts image metadata and appropriate template resources for downstream anatomical processing

    Parameters
    ----------
    in_file : Path to a BIDS-compatible anatomical image file (e.g., T1w, T2w, FLAIR)

    Returns
    -------
    modality : Extracted modality from filename
    bspline : B-spline grid distance for N4 bias field correction (400 for FLAIR, 200 otherwise)
    tpl_target_path : Path to the selected MNI template image
    tpl_mask_path : Path to the binary brain mask of the selected MNI template.
    wm_tpl : Path to the white matter probabilistic segmentation map.
    tissue_tpls : List of paths to probabilistic segmentation maps for CSF, GM, and WM.
    likelihood_model : Tissue segmentation model type for Atropos (e.g., 'HistogramParzenWindows' for FLAIR, 
        'Gaussian' for others).
    """

    #1. Get modality
    modality = _get_mod(in_file)

    #2. Get bspline, atropos segmentation likelihood_model, MNI template ID  
    modality_params = {
        "FLAIR": {
            "bspline": 400, 
            "tpl_id": "GG853", 
            "likelihood_model": "HistogramParzenWindows"},
        "T1w":   {
            "bspline": 200, 
            "tpl_id": config.workflow.template_id, 
            "likelihood_model": "Linear"
            },
    }
    params = modality_params.get(modality, modality_params["T1w"])  # fallback to T1w
    #3. Get MNI template
    if modality == 'FLAIR':
       _get_custom_templates()
    tpl_id =params["tpl_id"]
    template_spec = {}
    template_spec["suffix"] = template_spec.get("suffix", modality)
    tpl_target_path, common_spec = get_template_specs(tpl_id, template_spec=template_spec, fallback=True,)

    #4. binary brain mask 
    tpl_mask_path =  get_template(tpl_id, desc="brain", suffix="mask", **common_spec)

    #5. Check species and configure tpls accordingly 
    tpl_target_path, tpl_mask_path = ((tpl_target_path, tpl_mask_path) 
                                      if config.workflow.species.lower() == "human" 
                                      else (str(get_template(tpl_id, suffix="T2w")), 
                                            str(get_template(tpl_id, desc="brain", suffix="mask")[0]))
                                            )
    #6. Get tissue templates
    tissue_tpls = [str(p) 
                   for p in get_template(tpl_id,
                                         suffix="probseg",
                                         resolution=(1 if config.workflow.species.lower() == "human" else None),
                                         label=["CSF", "GM", "WM"],
                                         )]
    #7. WM template                                         
    wm_tpl= next(tissue_tpl for tissue_tpl in tissue_tpls if 'label-WM' in tissue_tpl)

    message = f"""

    -----------------------------
    Anatomical preprocessing workflow parameters:
        * File: {in_file}
        * Modality:                                     {modality}
        * INU bspline fitting distance:                 {params["bspline"]}
        * Template:                                     {tpl_id}
            * Path:          {tpl_target_path}
    
            * Template specifications:              {common_spec}
        * tpl_mask_path:            {tpl_mask_path}
        * tissue_tpls:              {tissue_tpls}
        * wm_tpl:                   {wm_tpl}

    """
    config.loggers.workflow.info(message)

    return modality, params["bspline"], tpl_target_path, tpl_mask_path, wm_tpl, tissue_tpls, params["likelihood_model"]