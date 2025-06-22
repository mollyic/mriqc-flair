from niworkflows.interfaces.norm import (
    _SpatialNormalizationInputSpec,
    SpatialNormalization,
)

from nipype.interfaces.base import traits
from niworkflows.interfaces.reportlets.registration import SpatialNormalizationRPT,_SpatialNormalizationInputSpecRPT, _SpatialNormalizationOutputSpecRPT
from niworkflows.interfaces.reportlets import base as nrb
from pathlib import Path
import shutil
from typing import Optional, List
import urllib.request

from niworkflows.interfaces.norm import NIWORKFLOWS_LOG



class _WrapSpatialNormalizationInputSpec(_SpatialNormalizationInputSpec):
    reference = traits.Enum(
        "T1w",
        "T2w",
        "boldref",
        "PDw",
        "FLAIR",
        mandatory=True,
        usedefault=True,
        desc="set the reference modality for registration (extended)",
    )

class WrapSpatialNormalization(SpatialNormalization):
    input_spec = _WrapSpatialNormalizationInputSpec

class _WrapSpatialNormalizationInputSpecRPT(
    nrb._SVGReportCapableInputSpec, _WrapSpatialNormalizationInputSpec
):
    pass
class WrapSpatialNormalizationRPT(nrb.RegistrationRC, WrapSpatialNormalization):
    input_spec = _WrapSpatialNormalizationInputSpecRPT
    output_spec = _SpatialNormalizationOutputSpecRPT

    def _post_run_hook(self, runtime):
        # We need to dig into the internal ants.Registration interface
        self._fixed_image = self._get_ants_args()["fixed_image"]
        if isinstance(self._fixed_image, (list, tuple)):
            self._fixed_image = self._fixed_image[0]  # get first item if list

        if self._get_ants_args().get("fixed_image_mask") is not None:
            self._fixed_image_mask = self._get_ants_args().get("fixed_image_mask")
        self._moving_image = self.aggregate_outputs(runtime=runtime).warped_image
        NIWORKFLOWS_LOG.info(
            "Report - setting fixed (%s) and moving (%s) images",
            self._fixed_image,
            self._moving_image,
        )

        return super()._post_run_hook(runtime)


def _download_file(url: str, output_path: str) -> None:
    urllib.request.urlretrieve(url, output_path)
    print(f"\t * Outfile: {output_path}")
    return output_path

def _construct_mask(tpl_id:str, template_dir: Path, resolution: int = 1) -> None:

    """
    Note: excluding as too messy, just defaulting to standard MNI brain mask 
    Construct a brain mask from tissue probability maps for a given template.
    """
    import nibabel as nb
    import numpy as np
    import re
    print(f"Generating brain mask for {tpl_id} template from tissue probability maps.")

    template_dir = Path(template_dir)

    list_tissues = [fname for fname in template_dir.rglob("*") if re.match(rf'^tpl-{tpl_id}.*(CSF|GM|WM).*probseg\.nii(\.gz)?$', str(fname.name))]
    print(list_tissues)
    t1, t2, t3 = (nb.load(str(file)).get_fdata() for file in list_tissues)
    sum_tissues = t1+t2+t3

    bmask = (sum_tissues > 0.5).astype(np.uint8)

    in_nii = nb.load(list_tissues[0])
    out_file =Path(list_tissues[0]).parent / ('tpl-' + tpl_id + '_res-0' + str(resolution) +'_desc-brain_mask.nii.gz')
    hdr = in_nii.header.copy()
    hdr.set_data_dtype(np.uint8)
    nb.Nifti1Image(bmask.astype(np.uint8), in_nii.affine, hdr).to_filename(out_file)
 
    return out_file


def _get_custom_templates(modality: str, 
                          template_id: str = 'GG853', 
                          template_dir: Optional[Path] = None,
                          template_res: int = 1) -> List:
    """
    Retrieve or add a custom template to the users local TemplateFlow database.

    If the template is not already present in TemplateFlow, it will attempt to download if the template_id is 'GG853'
    or copy it from the specified directory.

    Parameters
    ----------
    modality : str
        Modality of the template (e.g., 'FLAIR').
    template_id : str
        Template identifier (default: 'GG853').
    template_dir : Optional[Path]
        Path to a local directory containing template files.
    template_res : int
        Resolution level of the template (default: 1).

    Returns
    -------
    int
        0 if the template is successfully registered or already exists.
        Raises FileNotFoundError if files are missing.
    """

    from pathlib import Path
    from templateflow import conf
    import re 

    template_dir = Path(template_dir) if template_dir else None
    tf_tpl_dir: Path = conf.TF_HOME / ('tpl-' + template_id)

    tpl_head_re = rf'{template_id}.*res-(0?{template_res}).*'
    tpl_patterns = [
            rf'{tpl_head_re}{suffix}'
            for suffix in [f'{modality}\.nii(\.gz)?$',
                           'desc-brain_mask\.nii(\.gz)?$',
                           'label-CSF.*probseg\.nii(\.gz)?$',
                           'label-WM.*probseg\.nii(\.gz)?$',
                           'label-GM.*probseg\.nii(\.gz)?$']]

    # Check if templateflow directory is present
    if tf_tpl_dir.exists():
        matches = [fname for fname in tf_tpl_dir.rglob("*") for ptn in tpl_patterns if re.match(ptn, str(fname.name))]
        print(f"Template {template_id} already exists in the templateflow database.")

        if len(matches) == len(tpl_patterns):
            print(f'Found {len(matches)}/5 matches')
            return 0
        else: 
            print(f"Insufficient files in local templateflow directory, attempting to source...")


    if template_id == 'GG853':
        # Download the template files if template_id is GG853
        tf_tpl_dir.mkdir(parents=True, exist_ok=True)
        dwnld_head = 'https://s3.us-east-2.amazonaws.com/brainder/software/flair/templates/'
        dwnld_lst = [{'url_name': 'GG-853-FLAIR-1.0mm.nii.gz', 'tf_name': f'tpl-{template_id}_res-0{template_res}_{modality}.nii.gz'},
                     {'url_name': 'GG-853-GM-1.0mm.nii.gz', 'tf_name': f'tpl-{template_id}_res-0{template_res}_label-GM_probseg.nii.gz'},
                     {'url_name': 'GG-853-WM-1.0mm.nii.gz', 'tf_name': f'tpl-{template_id}_res-0{template_res}_label-WM_probseg.nii.gz'},
                     {'url_name': 'GG-853-CSF-1.0mm.nii.gz', 'tf_name': f'tpl-{template_id}_res-0{template_res}_label-CSF_probseg.nii.gz'}]

        for file_dict in dwnld_lst:
            if not (tf_tpl_dir / file_dict['tf_name']).exists():
                _download_file(f"{dwnld_head}{file_dict['url_name']}", tf_tpl_dir / file_dict['tf_name'])

        # Create a brain mask GG853 from tissue probability maps
        bmask_files = [f for f in tf_tpl_dir.rglob("*")  if f'res-0{template_res}_desc-brain_mask' in f.name]
        if len(bmask_files) == 0:
            #bmask_file = _construct_mask(tpl_id= template_id, template_dir=tf_tpl_dir)
            bmask_file = _construct_mask(tpl_id= template_id, template_dir=tf_tpl_dir)

        else:
            bmask_file = bmask_files[0]
        
    else: 
        # Copy the template files from the specified directory
        if template_dir.exists():
            matches = [fname for fname in template_dir.rglob("*") for ptn in tpl_patterns if re.match(rf'{tpl_head_re}{ptn}', str(fname.name))]
            if len(matches) == len(tpl_patterns):
                for file in matches:
                    shutil.copy2(file, tf_tpl_dir / file.name)
                print(f"Template {template_id} added to the templateflow database.")
                return 0
            else: 
                raise FileNotFoundError(f"Template directory does not contain all required files: {template_dir}")
        else:
            raise FileNotFoundError(f"Template directory not found: {template_dir}")
    return [bmask_file, tf_tpl_dir]

