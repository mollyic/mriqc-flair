"""
Clean and refine brain tissue segmentations using using a simplified verions of ANTs
`antsBrainExtraction.sh`.

Takes outputted ATROPOS segmentation output and performs tissue-specific morphological 
cleaning steps to improve the quality of the segmentation masks.

Created by Molly Ireland

"""


from nipype import Node, Workflow, MapNode, Function
from nipype.interfaces import utility as niu
from nipype.interfaces.ants import (ImageMath, MultiplyImages,ThresholdImage)


def clean_tissue_segs(name="clean_segs", padding = 10):
    """
    These morphological changes are derived from the ANTS brain extraction workflow; antsBrainExtraction.sh
    """
    wf = Workflow(name)

    inputnode = Node(niu.IdentityInterface(fields=["segmentation", "modality", 'pvms']), name="inputnode")

    outputnode = Node(niu.IdentityInterface(fields=['out_segs_proc', 'out_segm', 'out_pvms']), name="outputnode")

    # Morphological dilation of brainmask, radius=2
    dil_bmsk = Node(ImageMath(operation="MD", op2="2", copy_header=True), name="bmask_dilate")
    # Get largest connected component
    bmsk_lcc = Node(ImageMath(operation="GetLargestComponent", copy_header=True),name="bmask_lcc")

    # Pad images with zeros before processing
    pad_segm = Node(
        ImageMath(operation="PadImage",
                  op2=f"{padding}"),
        name="segm_pad")

    # 1. Split segmentation in binary masks
    segm_binarise = Node(
            Function(function=_binarise_seg_labels,
                         output_names=[ "out_csf", "out_gm","out_wm"]),
            name="binarise_seg_labels")
    
    # 2.  Create individual workflows to clean and refine each tissue class 
    csf_wf = _tissue_wf(tissue='csf')
    wm_wf = _tissue_wf(tissue='wm')
    gm_wf = _tissue_wf(tissue='gm')

    # 3. Rebuild segmentation from cleaned tissue masks
    seg_builder = _rebuild_segmentation()

    # 4. Mask partial volume maps with refined segmentation masks 
    mask_pvms = MapNode(
        Function(
            input_names=["pvm", 'modality', 'seg_lst'],
            output_names=["out_pvm_msk"],
            function=_mask_pvms,
            execution={"keep_inputs": True, "remove_unnecessary_outputs": False},
        ),
        iterfield=["pvm"],
        name="mask_pvms",
    )
    wf.connect([
        (inputnode, pad_segm, [("segmentation", "op1")]),
        (pad_segm, segm_binarise, [("output_image", "in_segm")]),
        (inputnode, segm_binarise, [("modality", "modality")]),
        (segm_binarise, wm_wf, [("out_wm", "inputnode.in_file")]),
        (segm_binarise, gm_wf,         [("out_gm", "inputnode.in_file")]),
        (segm_binarise, csf_wf,        [("out_csf", "inputnode.in_file")]),
        (wm_wf, seg_builder,        [("outputnode.out_tissue", "inputnode.wm")]),
        (csf_wf, seg_builder,       [("outputnode.out_tissue", "inputnode.csf")]),
        (gm_wf, seg_builder,        [("outputnode.out_tissue", "inputnode.gm")]),
        (inputnode, mask_pvms,      [("modality", "modality")]),
        (inputnode, mask_pvms,      [("pvms", "pvm")]),
        (seg_builder, mask_pvms,    [("outputnode.proc_seg_lst", "seg_lst")]),
        (seg_builder, outputnode,   [("outputnode.sumd_tissues", "out_segm"),
                                     ("outputnode.proc_seg_lst", "out_segs_proc")]),
        (mask_pvms, outputnode,     [("out_pvm_msk", "out_pvms")]),
    ])
    return wf

def _binarise_seg_labels(in_segm, modality):
    """
    Extract binary CSF, GM, and WM masks from a labeled segmentation image.

    Parameters
    ----------
    in_segm : path to segmentation file
    modality : scan modality used to determine label-to-tissue mapping.

    Returns
    -------
    out_csf, out_gm, out_wm: Paths to binarized masks.

    """
    from os import getcwd
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix
    from mriqc.config import ATROPOS_MODELS

    model = ATROPOS_MODELS[modality]
    nii = nb.load(in_segm)
    label_data = np.asanyarray(nii.dataobj).astype(np.uint8)

    out_files = {}
    for tissue, label in model.items():
        mask_data = (label_data == label).astype(np.uint8)
        mask_data[mask_data > 0] = 1  # Ensure no values >1 sneak in
        out_mask = nii.__class__(mask_data, nii.affine, nii.header)
        out_mask.set_data_dtype(np.uint8)  # Enforce dtype in header

        out_file = fname_presuffix(
            in_segm,
            suffix=f"_class-{label}_tissue-{tissue}",
            newpath=getcwd()
        )
        out_mask.to_filename(out_file)
        out_files[tissue] = out_file

    return (
        out_files.get("csf"),
        out_files.get("gm"),
        out_files.get("wm")
    )

def _tissue_wf(tissue, padding=10):
    """
    Generate cleaned and eroded binary masks for a specific tissue class.

    Parameters
    ----------
    tissue : Tissue type to process 

    Returns
    -------
    workflow : nipype Workflow
        A workflow that outputs:
            - out_tissue: processed binary mask for the specified tissue.
            - erode_csf: eroded CSF mask (only for CSF tissue).
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


def _rebuild_segmentation():
    """
    Reconstruct cleaned segmentation image from individual CSF, GM, and WM masks.

    Combines processed tissue masks into a single labeled image and creates a list of
    binary masks for use in masking partial volume maps.

    Returns
    -------
    workflow : nipype Workflow
            - sumd_tissues: combined segmentation image (CSF + GM + WM).
            - proc_seg_lst: list of binarized tissue masks in CSF-GM-WM order.
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

    """
    Mask a partial volume map with a refined tissue segmentation mask.

    Parameters
    ----------
    pvm : Path to the partial volume map NIfTI file
    modality : Modality used to determine tissue labels 
    seg_lst : List of tissue segmentation mask file paths

    Returns
    -------
    out_pvm_msk : Path to the masked PVM NIfTI file.
    """
    import nibabel as nb
    import re
    from nipype.utils.filemanip import fname_presuffix
    from os import getcwd
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

    return out_pvm_msk
