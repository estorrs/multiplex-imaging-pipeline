import anndata
import pandas as pd
import numpy as np
import tifffile
from skimage.measure import regionprops

from mip.utils import extract_ome_tiff


def get_spatial_features(label_fp, ome_fp):
    label_img = tifffile.imread(label_fp)
    cell_ids = np.unique(label_img)[1:]
    
    channel_to_img = extract_ome_tiff(ome_fp)
    channels = list(channel_to_img.keys())
    multichannel_img = np.stack([channel_to_img[c] for c in channels], axis=-1)
    
    props = regionprops(label_img, intensity_image=multichannel_img, )
    
    data = []
    columns = ['area', 'perimeter',
               'bbox_row_min', 'bbox_col_min', 'bbox_row_max', 'bbox_col_max',
               'centroid_row', 'centroid_col',
               'eccentricity']
    columns += [f'{c} intensity max' for c in channels]
    columns += [f'{c} intensity mean' for c in channels]
    columns += [f'{c} intensity min' for c in channels]

    for p in props:
        prop_data = [p.area, p.perimeter,
                     p.bbox[0], p.bbox[1], p.bbox[2], p.bbox[3],
                     p.centroid[0], p.centroid[1],
                     p.eccentricity]
        prop_data += p.intensity_max.tolist()
        prop_data += p.intensity_mean.tolist()
        prop_data += p.intensity_min.tolist()

        data.append(prop_data)

    df = pd.DataFrame(data=data, columns=columns, index=cell_ids)
    df.index.name = 'cell_id'
    
    return df


def save_spatial_features(label_fp, ome_fp, output_prefix):
    df = get_spatial_features(label_fp, ome_fp)
    
    adata = anndata.AnnData(
        X=df[[c for c in df.columns if 'intensity mean' in c]].values,
        obs=df[[c for c in df.columns if 'intensity mean' not in c]],
        var=pd.DataFrame(index=[c.replace(' intensity mean', '')
                                for c in df.columns if 'intensity mean' in c]))
    # for some reason needs this to make sure it saves correctly
    adata.var['marker'] = adata.var.index.to_list()
    
    df.to_csv(f'{output_prefix}.txt', sep='\t')
    adata.write_h5ad(f'{output_prefix}.h5ad')
