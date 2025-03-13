from nipype import Node, Workflow, MapNode, Function
from nipype.interfaces import utility as niu
from nipype.interfaces.ants import (ImageMath, MultiplyImages,ThresholdImage)
from collections import OrderedDict

ATROPOS_MODELS = {
    'T1w': OrderedDict([('csf', 1), ('gm', 2), ('wm', 3)]),
    'T2w': OrderedDict([('csf', 3), ('gm', 2), ('wm', 1)]),
    'FLAIR': OrderedDict([('csf', 1), ('gm', 2), ('wm', 3)])
    }



def init_atropos_wf(name="atropos_wf", use_random_seed=True, omp_nthreads=None,mem_gb=3.0, padding = 10):
    """
    Create an ANTs' ATROPOS workflow for brain tissue segmentation.
    The workflow also executes steps 8 and 9 of the brain extraction workflow.
    """
    wf = Workflow(name)

    #inputnode = Node(niu.IdentityInterface(fields=["in_corrected", "in_mask", "modality"]), name="inputnode")

    inputnode = Node(niu.IdentityInterface(fields=["classified_image", "modality", 'posteriors']), name="inputnode")

    outputnode = Node(niu.IdentityInterface(fields=['out_segs_proc', 'out_segm', 'out_pvms', 'out_eroded']), name="outputnode")

    # Morphological dilation of brainmask, radius=2
    dil_bmsk = Node(ImageMath(operation="MD", op2="2", copy_header=True), name="bmask_dilate")
    # Get largest connected component
    bmsk_lcc = Node(ImageMath(operation="GetLargestComponent", copy_header=True),name="bmask_lcc")

    # Run atropos (core node)
    # atropos = Node(
    #     Atropos(
    #         convergence_threshold=0.0,
    #         dimension=3,
    #         initialization="KMeans",
    #         likelihood_model="Gaussian",
    #         mrf_radius=[1, 1, 1],
    #         mrf_smoothing_factor=0.1,
    #         n_iterations=3,
    #         number_of_tissue_classes=3,
    #         out_classified_image_name="segment.nii.gz",
    #         output_posteriors_name_template="segment_%02d.nii.gz",
    #         save_posteriors=True,
    #         use_random_seed=use_random_seed,
    #     ),
    #     name="atropos",
    #     n_procs=omp_nthreads,
    #     mem_gb=mem_gb,
    # )

    # Pad images with zeros before processing
    pad_segm = Node(
        ImageMath(operation="PadImage",
                  op2=f"{padding}"),
        name="segm_pad")

    # Split segmentation in binary masks
    sel_labels = Node(
            Function(function=_select_labels,
                         output_names=[ "out_csf", "out_gm","out_wm"]),
            name="segm_binarize")

    csf_wf = _tissue_wf(tissue='csf')
    wm_wf = _tissue_wf(tissue='wm')
    gm_wf = _tissue_wf(tissue='gm')

    seg_builder = _rebuild_segmentation()

    mask_pvms = MapNode(
        Function(
            input_names=["pvm", 'modality', 'seg_lst'],
            output_names=["out_pvm_msk", 'eroded'],
            function=_mask_pvms,
            execution={"keep_inputs": True, "remove_unnecessary_outputs": False},
        ),
        iterfield=["pvm"],
        name="mask_pvms",
    )
    wf.connect([
        #(inputnode, dil_bmsk,       [("in_mask", "op1")]),
        #(inputnode, atropos,        [("in_corrected", "intensity_images")]),
        #(dil_bmsk, bmsk_lcc,        [("output_image", "op1")]),
        #(bmsk_lcc, atropos,         [("output_image", "mask_image")]),
        (inputnode, pad_segm,       [("classified_image", "op1")]),
        (pad_segm, sel_labels,      [("output_image", "in_segm")]),
        (inputnode, sel_labels,     [("modality", "modality")]),
        (sel_labels, wm_wf,         [("out_wm", "inputnode.in_file")]),
        (sel_labels, gm_wf,         [("out_gm", "inputnode.in_file")]),
        (sel_labels, csf_wf,        [("out_csf", "inputnode.in_file")]),
        (wm_wf, seg_builder,        [("outputnode.out_tissue", "inputnode.wm")]),
        (csf_wf, seg_builder,       [("outputnode.out_tissue", "inputnode.csf")]),
        (gm_wf, seg_builder,        [("outputnode.out_tissue", "inputnode.gm")]),
        (inputnode, mask_pvms,      [("modality", "modality")]),
        (inputnode, mask_pvms,      [("posteriors", "pvm")]),
        (seg_builder, mask_pvms,    [("outputnode.proc_seg_lst", "seg_lst")]),
        (seg_builder, outputnode,   [("outputnode.sumd_tissues", "out_segm"),
                                     ("outputnode.proc_seg_lst", "out_segs_proc")]),
        # (inputnode, outputnode,     [("classified_image", "out_segm"),
        #                              ("posteriors", "out_posterior")]),
        (mask_pvms, outputnode,     [("out_pvm_msk", "out_pvms"),
                                     ("eroded", "out_eroded")]),
    ])
    return wf

def _rebuild_segmentation():
    """
    Reassemble segmentation after processing

    Outputs:
        1. Reassembled segmentation
        2. Segmentation list
    """
    inputnode = Node(niu.IdentityInterface(fields=['csf', 'gm', 'wm']), name = 'inputnode')
    outputnode = Node(niu.IdentityInterface(fields=['sumd_tissues', 'proc_seg_lst']), name = 'outputnode')

    def _merged(csf, wm, gm):
        return [csf, gm, wm]

    add_gm_csf = Node(ImageMath(operation="addtozero", output_image = '%s_addcsf.nii.gz'), name="add_gm_csf")
    add_all = Node(ImageMath(operation="addtozero", output_image = 'segment_processed.nii.gz'), name="add_all")
    merge_tpms = Node(Function(function=_merged, output_names=['proc_seg_lst']), name="merge_tpms")

    bin_gmcsf = Node(ThresholdImage(dimension=3, th_low=0.01, th_high=1e7, inside_value=0, outside_value=1), name="bin_gmcsf")
    #binMult_gmcsf = Node(ImageMath(operation='mul', op2=-1), name="binMult_gmcsf")
    #binMultAdd_gmcsf = Node(ImageMath(add='add', op2=1), name="binMultAdd_gmcsf")

    clean_wm = Node(MultiplyImages(dimension=3, output_product_image = f'seg-wm_multiply.nii.gz'), name = 'clean_wm')

    wf = Workflow(name ='segm_sum')
    wf.connect([
        (inputnode, add_gm_csf,             [('gm', 'op1'), ('csf', 'op2')]),
        (inputnode, add_all,                [('wm', 'op1')]),

        (add_gm_csf, bin_gmcsf,             [('output_image', 'input_image')]),
        # (bin_gmcsf, binMult_gmcsf,          [('output_image', 'op1')]),
        # (binMult_gmcsf, binMultAdd_gmcsf,   [('output_image', 'op1')]),
        # (binMultAdd_gmcsf, clean_wm,        [('output_image', 'first_input')]),

        (bin_gmcsf, clean_wm,               [('output_image', 'first_input')]),
        (inputnode, clean_wm,               [('wm', 'second_input')]),

        (add_gm_csf, add_all,               [('output_image', 'op2')]),
        (add_all, outputnode,               [('output_image', 'sumd_tissues')]),
        (inputnode, merge_tpms,             [('csf', 'csf'), ('gm', 'gm')]),
        (clean_wm, merge_tpms,              [('output_product_image', 'wm')]),

        (merge_tpms, outputnode,            [("proc_seg_lst", "proc_seg_lst")]),
        ])

    return wf


def _mask_pvms(pvm, modality, seg_lst):
    import nibabel as nb
    import re
    from nipype.utils.filemanip import fname_presuffix
    from os import getcwd
    from scipy.ndimage import binary_erosion, generate_binary_structure, label, binary_fill_holes
    from mriqc.config import ATROPOS_MODELS
    import numpy as np

    model = ATROPOS_MODELS[modality]
    idx = re.search(r'segment_(\d+).nii.gz', pvm).group(1)
    tissue = next(key for key, val in model.items() if f'0{val}' == idx)
    segmsk = next(file for file in seg_lst if re.search(f'seg-{tissue}', file))

    pvm_nii = nb.load(pvm)
    pvm_data = pvm_nii.get_fdata(dtype=np.float32)
    segmsk_data = nb.load(segmsk).get_fdata(dtype=np.float32)
    segmask_bin = (segmsk_data > 0).astype(np.int_) #boolean operation

    pvm_maskd = pvm_data * segmask_bin
    out_pvm_msk = fname_presuffix(pvm, suffix=f"_class-{idx}_tissue-{tissue.upper()}_mskseg", newpath=getcwd())
    nb.Nifti1Image(pvm_maskd.astype(np.float32), pvm_nii.affine, pvm_nii.header).to_filename(out_pvm_msk)
    eroded = out_pvm_msk

    if modality == 'T2w' and tissue == 'wm':
        erode_iter =2
        connect = 2

        # Erode WM to remove CSF overlap in T2w
        eroded_data = binary_erosion(segmask_bin, structure=generate_binary_structure(3, connectivity=connect), iterations=erode_iter).astype(np.uint8)
        arr_labels, features = label(eroded_data)
        count_labels = np.bincount(arr_labels.ravel())[1:]  #exclude background
        lcc = count_labels.argmax() +1
        lcc_mask = arr_labels == lcc
        fillh = binary_fill_holes(lcc_mask)
        out_t2_msk = fname_presuffix(pvm, suffix=f"_class-{idx}_tissue-{tissue.upper()}_mskseg_getLCC", newpath=getcwd())
        print('T2w')
        print(f'\t *{tissue}:')
        print(f'\t *{out_t2_msk}:')
        pvm_maskd_lcc = pvm_data * fillh
        nb.Nifti1Image(pvm_maskd_lcc.astype(np.float32), pvm_nii.affine, pvm_nii.header).to_filename(out_t2_msk)
        eroded = out_t2_msk


    return out_pvm_msk, eroded

def _select_labels(in_segm, modality):
    """
    Extract tissues from FSL segmentation and binarizes each tissue type

    Outputs:
        * List of binarized tissues
    """
    from os import getcwd
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix
    from mriqc.config import ATROPOS_MODELS

    model = ATROPOS_MODELS[modality]
    cwd = getcwd()
    nii = nb.load(in_segm)
    label_data = np.asanyarray(nii.dataobj).astype("uint8")

    for tissue, label in model.items():
        newnii = nii.__class__(np.uint8(label_data == label), nii.affine, nii.header)
        newnii.set_data_dtype("uint8")
        out_file = fname_presuffix(in_segm, suffix=f"_class-{label}_tissue-{tissue}", newpath=cwd)
        newnii.to_filename(out_file)
        #out_files.append(out_file)
        if tissue == 'csf':
            out_csf = out_file
        if tissue == 'gm':
            out_gm = out_file
        if tissue == 'wm':
            out_wm = out_file

        print(f'\t * Tissue: {tissue}')
        print(f'\t * label: {label}')
        print(f'\t * in_segm: {in_segm}')
        print(f'\t * out_file: {out_file}')

    return out_csf, out_gm, out_wm

def _tissue_wf(tissue, padding=10):
    """
    Perform processing for segmentations relative to each tissue type.

    Outputs:
        1. Processed segmentation tissue
        2. Eroded CSF
    """
    inputnode = Node(niu.IdentityInterface(fields=['in_file']), name = 'inputnode')
    outputnode = Node(niu.IdentityInterface(fields=['out_tissue', 'erode_csf']), name='outputnode')

    lcc = Node(ImageMath(operation="GetLargestComponent", copy_header=True), name ='getLCC')
    fill_holes = Node(ImageMath(operation="FillHoles", op2="2", copy_header=True), name="fill_holes")
    multiply_imgs = Node(MultiplyImages(dimension=3, output_product_image = f'seg-{tissue}_multiply-getLCCxfill_holes.nii.gz'), name = 'multiply_imgs')
    morph_erode = Node(ImageMath(operation="ME", op2="0.05", copy_header=True), name="morph_erode")
    reindex_tissue = Node(MultiplyImages(dimension=3, output_product_image = f'seg-{tissue}_processed.nii.gz'), name = 'reindex_tissue')
    depad = Node(ImageMath(operation="PadImage", op2="-%d" % padding), name="depad")

    workflow = Workflow(name = f'tissue-{tissue}_wf')

    def _get_idx(file):
        import re
        return int(re.search('class-(\d+)', file).group(1))

    if tissue == 'csf':
        workflow.connect([
            (inputnode, morph_erode,                [('in_file', 'op1')]),
            (inputnode, reindex_tissue,             [('in_file', 'first_input')]),
            (morph_erode, outputnode,               [('output_image', 'erode_csf')]),
            ])
    else:
        workflow.connect([
            (inputnode, lcc,                        [('in_file', 'op1')]),
            ])
    if tissue == 'wm':
        workflow.connect([
            (lcc, reindex_tissue,                   [('output_image', 'first_input')]),
            ])
    if tissue == 'gm':
        workflow.connect([
            (lcc, fill_holes,                       [("output_image", "op1")]),
            (lcc, multiply_imgs,                    [("output_image", "first_input")]),
            (fill_holes, multiply_imgs,             [("output_image", "second_input")]),
            (multiply_imgs, reindex_tissue,         [('output_product_image', 'first_input')]),
        ])
    workflow.connect([
            (inputnode, reindex_tissue,             [(('in_file', _get_idx), 'second_input')]),
            (reindex_tissue, depad,                 [('output_product_image', 'op1')]),
            (depad, outputnode,                     [('output_image', 'out_tissue')]),
        ])

    return workflow


# dir = '/home/unimelb.edu.au/mollyi/mnt/MASSIVE/vn36_scratch/IQMs/EMriqc/work/test/anatWorkflow/'
# mask = dir + 'skullstrip_wf/_in_file_..home..mollyi..vn36..IQMS..data..RADIOL_niis..sub-101423..ses-01..anat..sub-101423_FLAIR.nii.gz/synthstrip/clipped_corrected_desc-brain_mask.nii.gz'
# corrected = dir + 'HeadMaskWorkflow/_in_file_..home..mollyi..vn36..IQMS..data..RADIOL_niis..sub-101423..ses-01..anat..sub-101423_FLAIR.nii.gz/Enhance/clipped_corrected_enhanced.nii.gz'

# atropos_test = Workflow(name= 'atropos_test', base_dir = 'work_atropos/')

# fielders = ["in_corrected", 'modality', 'in_mask']
# inputnode = Node(niu.IdentityInterface(fields=fielders), name="inputnode")
# inputnode.inputs.modality = 'FLAIR'
# inputnode.inputs.out_denoised = corrected
# inputnode.inputs.in_mask = mask

# bts_kmeans = init_atropos_wf()
# atropos_test.connect([
#         (inputnode, bts_kmeans,         [('modality', "inputnode.modality")]),
#         (inputnode, bts_kmeans,         [("out_denoised", "inputnode.in_corrected")]),
#         (inputnode, bts_kmeans,         [("in_mask", "inputnode.in_mask")]),
# ])
# atropos_test.run()
