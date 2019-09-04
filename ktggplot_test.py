
import geoimagine.gis.mj_gis_v80 as mjgis

from plotnine import *
import pandas as pd
import numpy as np

import fiona
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import geopandas as gpd

#Tutorial site
#https://medium.com/@gscheithauer/data-visualization-in-python-like-in-rs-ggplot2-bc62f8debbf5

'''
g = ggplot(pd.DataFrame(data={"x": [-5, 5]}), aes(x="x"))\
    + stat_function(fun=lambda x: x**2+2.5)
print (g)
'''
#constants
c_remote_data ='https://raw.githubusercontent.com/nickhould/craft-beers-dataset/master/data/processed/beers.csv'
c_col = ["#2f4858", "#f6ae2d", "#f26419",
         "#33658a", "#55dde0", "#2f4858",
         "#2f4858", "#f6ae2d", "#f26419",
         "#33658a", "#55dde0", "#2f4858"]

#functions
def labels(from_, to_, step_):
    return pd.Series(np.arange(from_, to_ + step_, step_)).apply(lambda x: '{:,}'.format(x)).tolist()
def breaks(from_, to_, step_):
    return pd.Series(np.arange(from_, to_ + step_, step_)).tolist()

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return clr.LinearSegmentedColormap('CustomMap', cdict)

def testBeer():
    data = pd.read_csv(c_remote_data)
    data = (
        data.filter([
            'abv',
            'ibu',
            'id',
            'name',
            'style',
            'brewery_id',
            'ounces'
        ]).
        set_index('id')
    )
    '''
    
    fig = ggplot(data.dropna(subset = ['abv'])) +\
        geom_histogram(aes(x = 'abv'),bins=15)  
    
    
    fig = ggplot(data.dropna(subset = ['abv'])) +\
        geom_histogram(aes(x = 'abv'),fill = c_col[0], color = 'black',bins=15)
    '''
       
    fig = (
        ggplot(data.dropna(subset = ['abv', 'style'])) +
        geom_histogram(
            aes(x = 'abv',
            fill = 'style'),
            bins = 15
        ) +
        labs(
            title ='Distribution of The alcoholic content by volume (abv)',
            x = 'abv - The alcoholic content by volume',
            y = 'Count',
        )
    )
    
    print(fig)
    
    
def testImage01():
    srcFPN = '/Volumes/africa/export/trmm/region/rainfall-ixc-n/africasubsahara/199801-201712@M/obs-lag-best_3b43_africasubsahara_199801-201712@M_v7-f-m-30km-nadd-a1.tif'

    srcDS, rasterLayer = mjgis.RasterOpenGetFirstLayer(srcFPN)
    #source = rasterio.open('/Volumes/africa/export/trmm/region/rainfall-ixc-n/africasubsahara/199801-201712@M/obs-lag-best_3b43_africasubsahara_199801-201712@M_v7-f-m-30km-nadd-a1.tif', 'r')
    rasterLayer.ReadBand()
    
    bounds = (rasterLayer.bounds[0], rasterLayer.bounds[2], \
              rasterLayer.bounds[1], rasterLayer.bounds[3])
    
    srcDS.CloseDS()
    
    img = plt.imread(srcFPN)
    plt.imshow(img, extent=bounds, vmin=2, vmax=13)
    plt.colorbar(ticks=range(12), label='digit value')
    plt.clim(-0.5, 11.5)
    plt.show()

def testImage02():
    srcFPN = '/Volumes/africa/export/trmm/region/rainfall-ixc-n/africasubsahara/199801-201712@M/obs-lag-best_3b43_africasubsahara_199801-201712@M_v7-f-m-30km-nadd-a1.tif'
    
    
    colors=[(0.75, 0.15, 0.15), (1, 0.75, 0.15), (0.15, 0.75, 0.15)]
    
    nThresholds = 12
    cmap = clr.LinearSegmentedColormap.from_list(name='custom', colors = colors, N=nThresholds)
    
    
    srcDS, rasterLayer = mjgis.RasterOpenGetFirstLayer(srcFPN)
    #source = rasterio.open('/Volumes/africa/export/trmm/region/rainfall-ixc-n/africasubsahara/199801-201712@M/obs-lag-best_3b43_africasubsahara_199801-201712@M_v7-f-m-30km-nadd-a1.tif', 'r')
    rasterLayer.ReadBand()
    
    source = rasterLayer.NPBAND
    srcDS.CloseDS()
    bounds = (rasterLayer.bounds[0], rasterLayer.bounds[2], \
              rasterLayer.bounds[1], rasterLayer.bounds[3])
    
    
    f = plt.figure(figsize=(6, 6))
    
    ax = plt.imshow(source, extent=bounds, cmap=cmap , vmin=1, vmax=13)
    
    plt.show()
    
def testImage03():
    
    srcFPN = '/Volumes/africa/export/trmm/region/rainfall-ixc-n/africasubsahara/199801-201712@M/obs-lag-best_3b43_africasubsahara_199801-201712@M_v7-f-m-30km-nadd-a1.tif'

    cdict = {'red':   [(0.0,  0.0, 0.0),
                   (0.5,  1.0, 1.0),
                   (1.0,  1.0, 1.0)],

         'green': [(0.0,  0.0, 0.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 1.0, 1.0),
                   (1.0,  1.0, 1.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (0.5,  0.0, 0.0),
                   (1.0,  1.0, 1.0)]}
    
    cmap = clr.LinearSegmentedColormap('custom',cdict)
      
    srcDS, rasterLayer = mjgis.RasterOpenGetFirstLayer(srcFPN)
    #source = rasterio.open('/Volumes/africa/export/trmm/region/rainfall-ixc-n/africasubsahara/199801-201712@M/obs-lag-best_3b43_africasubsahara_199801-201712@M_v7-f-m-30km-nadd-a1.tif', 'r')
    rasterLayer.ReadBand()
    
    source = rasterLayer.NPBAND
    
    srcDS.CloseDS()
    
    bounds = (rasterLayer.bounds[0], rasterLayer.bounds[2], \
              rasterLayer.bounds[1], rasterLayer.bounds[3])
    
    f = plt.figure(figsize=(6, 6))
    
    ax = plt.imshow(source, extent=bounds, cmap=cmap , vmin=0, vmax=13)
    plt.colorbar(ticks=range(12), label='digit value')
    plt.clim(-0.5, 11.5)
    
    plt.show()

if __name__ == "__main__":
    testBeer()
    
    #testImage01()
    
    #testImage02()
    
    #testImage03()
    




