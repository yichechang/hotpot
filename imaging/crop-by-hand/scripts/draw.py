from typing import Dict, List, Optional, Any
from pathlib import Path
import shutil

import napari
from magicgui import magicgui
from magicgui.widgets import Container
from aicsimageio import AICSImage
import pandas as pd


# TODO:
#   Consider picking some defaults for chinfo_dict
def roi_drawer(
    imageset_records: List,
    output_dir: Path,
    chinfo_dict: Optional[Dict] = None
):
    """View and draw ROIs on each image set from the filepath list.

    Parameters
    ----------
    imageset_records: list of dict
        List of imageset records. Each record needs to contain minimally
        `ImagesetFilepath` (absolute path point to an image file).
    output_dir: pathlib.Path
        Path to the directry where compiled results (csv) to be saved.
    chinfo_dict: Dict
        Dictionary containing channel information, in the format of
        {key: [val_Ch0, val_Ch1, val_Ch2, ...]}. Keys required include
        `fluoro` and `colormap`.

    Notes
    -----
    Need to click 'next' or 'prev' to trigger saving. This is fine for
    all but the last imageset, as it is not something natural to do at
    the end.
    """
    # ====================
    # Start napari viewer using the first image set
    # ====================
    current_path = 0
    img = AICSImage(imageset_records[current_path]['ImagesetFilepath']).xarray_data

    if chinfo_dict is None:
        n_channels = img.coords['C'].size
        chinfo_dict = {
            'fluoro': ['channel:'+str(idx) for idx in range(n_channels)],
            'colormap': ['Greys_r' for _ in range(n_channels)]
        }

    viewer = napari.view_image(
        img,
        channel_axis=1,
        colormap=chinfo_dict['colormap'],
        name=chinfo_dict['fluoro'],
    )
    # Connect these two image layers to autoscale contrast when updated
    # From: https://forum.image.sc/t/function-to-autoadjust-lut-programatically/40396/2
    for l in viewer.layers:
        l.events.data.connect(lambda e: l.reset_contrast_limits())

    # Create ShapesLayer for drawing ROIs
    shapes_layer = viewer.add_shapes(
        shape_type='rectangle',
        edge_color='coral',
        face_color='royalblue',
        name='ROIs',
    )
    shapes_layer.mode = 'add_rectangle'

    #
    # Check and prepare output file directries
    #
    output_dir.mkdir(exist_ok=True)



    # ====================
    # GUI functions for processing
    # ====================
    @magicgui(call_button='next')
    def next_callback(viewer: napari.Viewer):
        save_current_rois(viewer)
        if current_path + 1 < len(imageset_records):
            set_up_new_file(viewer, 'next')


    @magicgui(call_button='prev')
    def prev_callback(viewer: napari.Viewer):
        save_current_rois(viewer)
        if current_path - 1 >= 0:
            set_up_new_file(viewer, 'prev')


    def save_current_rois(viewer):
        image_path = imageset_records[current_path]['ImagesetFilepath']
        outname = Path(image_path).stem + '.csv'
        outpath = output_dir / outname

        # Not sure how to access the dataframe describing ROI features.
        # But can save them asa csv and they will have the following
        # fields: index, shape-type, vertex-index, axis-0, axis-1.
        # => Since these are intermediate files, just set 'index' to be
        #    RoiID when combining the csv files. (See `compile_callback`)
        viewer.layers['ROIs'].save(str(outpath))

        # include original image filepath
        df = pd.read_csv(outpath)
        df['ImagesetFilepath'] = image_path
        df.to_csv(outpath)


    def set_up_new_file(viewer, direction):
        # Keep track of state of which image file to use
        nonlocal current_path
        if direction == 'next':
            current_path += 1
        elif direction == 'prev':
            current_path -= 1

        # Update layers for new images without adding/removing layers
        img = AICSImage(imageset_records[current_path]["ImagesetFilepath"]).xarray_data
        # viewer.layers['miRFP670'].data = img.isel(C=0).data
        # viewer.layers['GFP'].data = img.isel(C=1).data
        # TODO:
        #   Do not rely on the order (i.e., don't use index)
        for idx,fluoro in enumerate(chinfo_dict['fluoro']):
            viewer.layers[fluoro].data = img.isel(C=idx).data

        # Reset ROI layer for the new image
        viewer.layers['ROIs'].data = []
        viewer.layers.selection.active = viewer.layers['ROIs']
        viewer.layers['ROIs'].mode = 'add_rectangle'

    # Add widgets to napari viewer
    container = Container(
        widgets=[next_callback, prev_callback]
    )
    viewer.window.add_dock_widget(container)

    napari.run()


def compile_roilist(
    tmp_dir: str,
    output_path: str,
) -> None:
    """
    Compile a single roilist from all tmp roilist csv files.
    """
    #
    # Gather all csv files containing ROI coordinates from each imageset
    #
    files = list(Path(tmp_dir).glob('*.csv'))
    dfs_coords = [pd.read_csv(file) for file in files]
    df_coords = pd.concat(dfs_coords)
    # Add ID for each ROI coordinates. The Numeric ID already exists as
    # column named 'index' (note: not as the index of the dataframe).
    df_coords = df_coords.rename(columns={'index': 'RoiID'})
    df_coords['RoiID'] = df_coords['RoiID'].astype(str)

    #
    # Process ROI coordinates to ROI
    #
    colnames_roiID = ['ImagesetFilepath', 'RoiID']
    df_rois = pd.pivot(df_coords,
                       index=colnames_roiID,
                       columns=['vertex-index'],
                       values=['axis-0', 'axis-1'],
                      )
    df_rois['r1'] = df_rois[[('axis-0', 0)]]
    df_rois['r2'] = df_rois[[('axis-0', 2)]]
    df_rois['c1'] = df_rois[[('axis-1', 0)]]
    df_rois['c2'] = df_rois[[('axis-1', 2)]]
    df_rois['ri'] = df_rois[['r1', 'r2']].min(axis=1)
    df_rois['rf'] = df_rois[['r1', 'r2']].max(axis=1)
    df_rois['ci'] = df_rois[['c1', 'c2']].min(axis=1)
    df_rois['cf'] = df_rois[['c1', 'c2']].max(axis=1)

    df_rois = df_rois.reset_index()
    df_rois = df_rois[colnames_roiID + ['ri', 'rf', 'ci', 'cf']].droplevel('vertex-index', axis=1)


    #
    # save output to file
    #
    df_rois.to_csv(output_path, quoting=2, index=False, header=True)

# TODO:
#   - can be written using pandas.DataFrame.from_records()
def extract_channel_info(
    channel_info_list: Optional[List[Dict[str, Any]]]
) -> Optional[Dict[str, List[Any]]]:
    """
    From a pepfile, convert channel info from list of {k: v}'s to {key: list of v's}
    """
    channels = channel_info_list
    if channels is not None:
        keys = channels[0].keys()
        return {
            key: [ch[key] for ch in channels]
            for key in keys
        }
    else:
        return None

def main(
    imageset_paths: List[str],
    output_path: str,
    chinfo_dict: Optional[Dict] = None
) -> None:
    # TODO:
    #   If refactoring roi_drawer(), pick some sensible defaults.
    """
    Combined ROI information will be saved as file <output_path>. A tmp
    directory will be created in the same dir which <output_path> is in,
    to store individual csv files for each Imageset.
    """
    imageset_records = (
        pd.DataFrame({'ImagesetFilepath': imageset_paths})
        .to_dict('records')
    )

    output_pathobj = Path(output_path)

    # Create a new tmp dir for this specific drawing action. Note that
    # old 'tmp/roi_drawer/' will be removed first.
    tmp_dir = output_pathobj.parent / "tmp" / "roi_drawer"
    if tmp_dir.exists() and tmp_dir.is_dir():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(exist_ok=True, parents=True)
    roi_drawer(imageset_records, tmp_dir, chinfo_dict)

    compile_roilist(str(tmp_dir), output_path)


if __name__ == '__main__':
    if 'snakemake' in globals():
        chinfo_dict = extract_channel_info(
            snakemake.config['input']['channels'])

        main(snakemake.input.images[:],
             snakemake.output.csv,
             chinfo_dict=chinfo_dict)
