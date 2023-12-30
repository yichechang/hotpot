from pathlib import Path
from functools import partial

import numpy as np
import pandas as pd

from aicsimageio import AICSImage
from aicsimageio.writers import OmeTiffWriter


def main(
        roilist_path: str, 
        output_dir: str
    ) -> None:
    df = pd.read_csv(roilist_path)
    crop_func = partial(crop_rois_from_a_file, 
                        output_dir=output_dir)
    df.groupby('ImagesetFilepath').apply(crop_func)


def crop_rois_from_a_file(
        df: pd.DataFrame, *, 
        output_dir: str
    ) -> None:
    # open the image (all rois are from the same image)
    im_path = df.name
    im = AICSImage(im_path)

    roi_records = df.to_dict('records')

    for roi in roi_records:
        # TODO
        #   - Validate that we don't rely on actual spatial coordinates
        ri, rf, ci, cf = np.array([roi['ri'], roi['rf'], roi['ci'], roi['cf']]).astype(int)

        cropped = im.xarray_data.isel(Y=range(ri,rf+1), X=range(ci,cf+1))
        outpath = (
            (
                Path(output_dir) / 
                (Path(im_path).stem + '_roi' +  str(roi['RoiID']))
            )
            .with_suffix('.ome.tif')
        )
        OmeTiffWriter.save(cropped.data, str(outpath)) 


if __name__ == '__main__':
    if 'snakemake' in globals():
        main(snakemake.input[0], 
             snakemake.params['outdir'])