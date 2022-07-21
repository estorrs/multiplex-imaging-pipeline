from collections import Counter

import tifffile
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from skimage.transform import resize
from skimage.measure import label
from skimage.color import label2rgb
from skimage.segmentation import find_boundaries
from sklearn.neighbors import NearestNeighbors

from bokeh.io import output_notebook, show
from bokeh.plotting import figure
from bokeh.models import Div, RangeSlider, Spinner, Slider, ColumnDataSource, CustomJS
from bokeh.layouts import column, row, gridplot
output_notebook()

import mip.utils as utils


def convert_sparse(X):
    if 'sparse' in str(type(X)).lower():
        return X.toarray()
    return X


def get_ideal_window(adata, x_col='centroid_col', y_col='centroid_row',
        radius=200, cell_type=None, cell_type_col='cell_type',
	top_idx=1, return_filtered=True):
    if cell_type is not None and isinstance(cell_type, str):
        filtered = adata[adata.obs[cell_type_col]==cell_type]
    elif cell_type is not None:
        filtered = adata[[True if c in cell_type else False
                          for c in adata.obs[cell_type_col]]]
    else:
        filtered = adata

    X = filtered.obs[[x_col, y_col]]
    nbrs = NearestNeighbors(algorithm='ball_tree').fit(X)

    g = nbrs.radius_neighbors_graph(X, radius=radius).toarray()
    center_idx = np.argsort(g.sum(axis=1))[-top_idx]
    r, c = filtered.obs.iloc[center_idx][y_col], filtered.obs.iloc[center_idx][x_col]
    r, c = int(r), int(c)
    r1, r2 = max(0, r - radius), r + radius
    c1, c2 = max(0, c - radius), c + radius

    filtered = filtered[((filtered.obs[y_col]>r1)&(filtered.obs[y_col]<r2))]
    filtered = filtered[((filtered.obs[x_col]>c1)&(filtered.obs[x_col]<c2))]

    if return_filtered:
        return filtered, (r1, r2, c1, c2)
    return r1, r2, c1, c2


def gate_region(adata, channel, channel_img, boundary_img, x_col='centroid_col', y_col='centroid_row',
        radius=200, cell_type=None, cell_type_col='cell_type_macro', top_idx=1, default_value=None,
        use_raw=True):
    filtered, (r1, r2, c1, c2) = get_ideal_window(
        adata, x_col=x_col, y_col=y_col,
        radius=radius, cell_type=cell_type, cell_type_col=cell_type_col,
	top_idx=top_idx, return_filtered=True)

    df = filtered.obs[[y_col, x_col]].copy()
    if use_raw:
        df['marker'] = np.mean(convert_sparse(filtered.raw[:, channel].X), axis=1).flatten()
    else:
        df['marker'] = np.mean(convert_sparse(filtered[:, channel].X), axis=1).flatten()
    if default_value is None:
        default_value = min(df['marker'])
    df['color'] = ['blue' if m >= default_value else 'gray'
                   for m in df['marker']]
    source = ColumnDataSource(data=dict(
                            x=df[x_col],
                            y=df[y_col],
                            color=df['color'],
                            marker=df['marker']
    ))

    plot = figure(width=400, height=400)
    plot.xgrid.grid_line_color = None
    plot.ygrid.grid_line_color = None

    plot.image(image=[~boundary_img[r1:r2, c1:c2]],
        x=c1, y=r1, dw=radius * 2, dh=radius * 2, level="image")

    plot.image(image=[channel_img[r1:r2, c1:c2]],
        x=c1, y=r1, dw=radius * 2, dh=radius * 2, palette="Greens256", level="image", alpha=.92)


    plot.circle('x', 'y', source=source, size=5, color="color")

    callback = CustomJS(args=dict(source=source), code="""
    const data = source.data;
    const f = cb_obj.value;
    const x = data['x'];
    const y = data['y'];
    const color = data['color'];
    const marker = data['marker'];
    for (let i = 0; i < x.length; i++) {
        if (marker[i] < f) {
            color[i] = 'gray';
        } else {
            color[i] = 'blue';
        }
    }
    source.change.emit();""")

    slider = Slider(start=min(df['marker']), end=max(df['marker']), value=default_value, step=.01, title="marker threshold")
    slider.js_on_change('value', callback)

    layout = column(slider, plot)

    show(layout)
    
    return ((r1, r2, c1, c2),
            (r1 / channel_img.shape[0], r2 / channel_img.shape[0], c1 / channel_img.shape[1], c2 / channel_img.shape[1]))



