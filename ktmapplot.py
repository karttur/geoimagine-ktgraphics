
import geoimagine.gis.mj_gis_v80 as mjgis
from geoimagine.layout import SetDrawLegendText

from plotnine import *
import pandas as pd
import numpy as np

import fiona
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib.image as img
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import geopandas as gpd

class MapPlot():
    ''' For pyplot graphs
    '''
    def __init__(self, mapArr, process, colorRamp, legend, scaling, measure, srcComp, dstComp):
        ''' 
        '''
        self.mapArr = mapArr
        self.process = process
        self.params = process.params
        self.colorRamp = colorRamp
        self.srcComp = srcComp
        self.dstComp = dstComp
        self.legend = legend
        self.legendText = SetDrawLegendText(self.colorRamp, legend, scaling, measure)
        self._CraeteColorDict()
        self._RasterMatPlotLibContinuous()
        #self._RasterGGplotContinuous()
        
    def _CraeteColorDict(self):
        self.cdict = {'red': [], 'green': [], 'blue': []}
        for c in self.colorRamp:
            
            i = float(c[0])/255.0
            if i <= 1.0 and i >= 0.0:
                r = (i,float(c[1])/255.0, float(c[1])/255.0)
                g = (i,float(c[2])/255.0, float(c[2])/255.0)
                b = (i,float(c[3])/255.0, float(c[3])/255.0)
                a = (i,float(c[4])/255.0, float(c[4])/255.0)
                self.cdict['red'].append(r)
                self.cdict['green'].append(g)
                self.cdict['blue'].append(b)
                #self.cdict['alpha'].append(a)
        
    def _RasterMatPlotLibContinuous(self):
        '''
        '''
        #Set the colormap
        cmap = clr.LinearSegmentedColormap('custom',self.cdict)
        #Open the source raster layer - only to get the extent
        srcDS, rasterLayer = mjgis.RasterOpenGetFirstLayer(self.srcComp.FPN)
        #Read tbe first band, including the metadata
        rasterLayer.ReadBand()
        #Close the raster data
        srcDS.CloseDS()
        #Set the bouunds of the map
        bounds = (rasterLayer.bounds[0], rasterLayer.bounds[2], \
                  rasterLayer.bounds[1], rasterLayer.bounds[3])
                
        #Set the font
        '''
        font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}
        plt.rc('font', **font)
        
        #Set the labelsizes for the x and y axis
        plt.rc('xtick', labelsize=10) 
        plt.rc('ytick', labelsize=10)
         
        '''
        fig = plt.figure(figsize=(self.params.plotwidth, self.params.plotheight))
        #ax = fig.add_subplot(111)
        
        # I have no idea why I do like this
        ax = plt.subplot(1, 1, 1)
        
        #Set title
        if len(self.params.title) > 0:
            font = {}
            if len(self.params.titlefontcolor) > 0:
                font['color'] = self.params.titlefontcolor
            if len(self.params.titlefont) > 0:
                font['family'] = self.params.titlefont
            if self.params.titlefontsize:
                font['size'] = self.params.titlefontsize
            if 'bold' in self.params.titlefonteffect:
                font['weight'] = 'bold'
            ax.set_title(self.params.title, y=self.params.titleypos, fontdict=font)
  
        #Set the frame (spine) colors, style and width (colors can vary, style and width not)
        spinecolorL = self.params.spinecolors.split(',')
        spinelinestyle = self.params.spinelinestyle
        spinelinewidth = self.params.spinelinewidth

        ax.spines['top'].set_color( spinecolorL[0] )
        ax.spines['right'].set_color( spinecolorL[1] )
        ax.spines['bottom'].set_color( spinecolorL[2] )
        ax.spines['left'].set_color( spinecolorL[3] )
        
        ax.spines['top'].set_linestyle( spinelinestyle )
        ax.spines['right'].set_linestyle( spinelinestyle )
        ax.spines['bottom'].set_linestyle( spinelinestyle )
        ax.spines['left'].set_linestyle( spinelinestyle )
        if spinelinewidth: #if not it will be default
            ax.spines['top'].set_linewidth( spinelinewidth )
            ax.spines['right'].set_linewidth( spinelinewidth )
            ax.spines['bottom'].set_linewidth( spinelinewidth )
            ax.spines['left'].set_linewidth( spinelinewidth )
        
        #Set the axis major tick, these are always set but if given is set using locator
        if self.params.majorticks:
            majorLocator = MultipleLocator(self.params.majorticks)
            ax.xaxis.set_major_locator(majorLocator)
            ax.yaxis.set_major_locator(majorLocator)
        ax.tick_params(which='major', length=self.params.majorticklength, 
                    width=self.params.majortickwidth,direction=self.params.majortickdirection)
        
        if ',' in self.params.majortickcolors:
            majtickXcolor,majtickYcolor = self.params.majortickcolors.split(',')
            ax.tick_params(axis='x', colors=majtickXcolor)
            ax.tick_params(axis='y', colors=majtickYcolor)
            #Major tick mark labels
        ax.tick_params(labelbottom=True, labeltop=True, labelleft=True, labelright=True,
                     bottom=True, top=True, left=True, right=True)

        if self.params.axisfontsize:
            ax.tick_params(labelsize=self.params.axisfontsize)
        
        axisXfontcolors = axisYfontcolors = False
        if ',' in self.params.axisfontcolors: 
            axisXfontcolors,axisYfontcolors = self.params.axisfontcolors.split(',')[0:2]
        
        for tick in ax.get_xticklabels():
            if len(self.params.axisfont) > 0:
                tick.set_fontname(self.params.axisfont)
            if axisXfontcolors:
                tick.set_color(axisXfontcolors)
        for tick in ax.get_yticklabels():
            if len(self.params.axisfont) > 0:
                tick.set_fontname(self.params.axisfont)
            if axisYfontcolors:
                tick.set_color(axisYfontcolors)
      
        '''
        THE ALTERNATIVE BELOW DOES NOT WORK
        #Get the labels
        labels = [item.get_text() for item in ax.get_xticklabels()]
        #Rest the labels using fontdict

        print ('labels',labels)
        print ('font',font)
        ax.set_xticklabels(labels, fontdict=font)
        
        '''
        #minor ticks only set if defined
        if self.params.minorticks:
            if self.params.minorticks < 0:
                minorLocator = AutoMinorLocator()
            else:
                minorLocator = MultipleLocator(self.params.minorticks)
            mintickXcolor,mintickYcolor = self.params.minortickcolors.split(',')
            ax.tick_params(which='minor', length=self.params.minorticklength, 
                    width=self.params.minortickwidth,direction=self.params.minortickdirection)
            ax.tick_params(which='minor', axis='x', colors=mintickXcolor)
            ax.tick_params(which='minor', axis='y', colors=mintickYcolor)
            
            ax.xaxis.set_minor_locator(minorLocator)
            ax.yaxis.set_minor_locator(minorLocator)
            
            ax.tick_params(which='minor', bottom=True, top=True, left=True, right=True)
            
            #then also turn on, otherwise some are missing?? 
            ax.minorticks_on()
            #Never any labels for minor ticks
                
        #ax.xaxis.label.set_color('blue')
        #ax.tick_params(axis='x', colors='red')
                        
        if self.params.gridlinewidth: 
            gridXcolor,gridYcolor =  self.params.gridcolors.split(',')      
            ax.grid(linestyle=self.params.gridlinestyle, 
                    linewidth=self.params.gridlinewidth, alpha=self.params.gridalpha)
            ax.grid(axis='x', color=gridXcolor)
            ax.grid(axis='y', color=gridYcolor)
        
        '''
        with plt.rc_context({'axes.edgecolor':'orange', 'xtick.color':'red', 'ytick.color':'green', 'figure.facecolor':'white'}):
            # Temporary rc parameters in effect
            fig, (ax1, ax2) = plt.subplots(1,2)
            ax1.plot(range(10))
            ax2.plot(range(10))
        '''
            
        #Set north arrow
        if len(self.params.northarrow) > 6:
            im = img.imread(self.params.northarrow)
            imaspect = im.shape[0]/im.shape[1]

            xsize = self.params.northarrowsize*(bounds[1]-bounds[0])/self.params.plotwidth
            xmargin = self.params.northarrowmarginx*(bounds[1]-bounds[0])/100
            ymargin = self.params.northarrowmarginy*(bounds[3]-bounds[2])/100

            #fig, ax = plt.subplots()
            if self.params.northarrowanchor == 'coords':
                pass
            elif self.params.northarrowanchor == 'upperleft':

                left = bounds[0]+xmargin; top = bounds[3]-ymargin  
                right = left+xsize
                bottom = top-(xsize*imaspect)
                
            ax.imshow(im, aspect='auto', extent=(left, right, bottom, top), zorder=3)

        imgplot = plt.imshow(self.mapArr, extent=bounds, cmap=cmap,vmin=0, vmax=255)

        #extendfrac = 5/255
        #cbar = plt.colorbar(ticks=[0, 63, 125, 188,  250], orientation='horizontal',extend='max', extendfrac=extendfrac)
        ticks = []
        tickLabels = []
        for item in self.legendText:
            ticks.append(item['pos'])
            tickLabels.append(item['txt'])

        legsize = '%(ls)d%%' %{'ls':self.params.legendrelativesize}
        legwidth = '%(ls)d%%' %{'ls':self.params.legendrelativewidth}
        legheight = '%(ls)d%%' %{'ls':self.params.legendrelativeheight}
        
        
        if self.params.legendfitaxis:
            divider = make_axes_locatable(ax)
            
            cax = divider.append_axes(self.params.legendposition, size=legsize, pad=self.params.legendpad)
        
        elif self.params.legendinsetaxis:
            cax = inset_axes(ax, width=legwidth, height=legheight, loc=self.params.legendinsetposition, borderpad= self.params.legendborderpad) 
        #plt.colorbar(cax=cbaxes, ticks=[0.,1], orientation='horizontal')
        
        if self.params.legendposition in ['right','left']:
            cbar = plt.colorbar(ticks=ticks, orientation='vertical',cax=cax)
        
            #cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
            cbar.ax.set_yticklabels(tickLabels)  # horizontal colorbar
            # clim compresees the orignal data to the new range, not wh
        else:
            #cax = divider.append_axes('bottom', size='5%', pad=0.5)
            cbar = plt.colorbar(ticks=ticks, orientation='horizontal',cax=cax)
            #cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
            cbar.ax.set_xticklabels(tickLabels)  # horizontal colorbar
            # clim compresees the orignal data to the new range, not what I want
        cbar.set_label(self.legend.columnhead)
        #plt.clim(0, 250);
        if self.params.tightlayout:
            fig.tight_layout()
        plt.show()
        if self.process.dstpath.hdrfiletype == 'show':
            plt.show()
        else:
            if self.params.tightlayout:
                plt.savefig(self.dstComp.FPN, bbox_inches='tight')
            else:
                plt.savefig(self.dstComp.FPN)
           
    def _RasterGGplotContinuous(self):
        '''
        '''
        #Set the colormap
        cmap = clr.LinearSegmentedColormap('custom',self.cdict)
        #Open the source raster layer - only to get the extent
        srcDS, rasterLayer = mjgis.RasterOpenGetFirstLayer(self.srcComp.FPN)
        #Read tbe first band, including the metadata
        rasterLayer.ReadBand()
        #Close the raster data
        srcDS.CloseDS()
        #Set the bouunds of the map
        bounds = (rasterLayer.bounds[0], rasterLayer.bounds[2], \
                  rasterLayer.bounds[1], rasterLayer.bounds[3])
        
        
        #Convert the image array to a datarame
        
        #colourImg = Image.open("test.png")
        #colourPixels = colourImg.convert("RGB")
        #colourArray = np.array(colourPixels.getdata()).reshape(colourImg.size + (3,))

        #Flip the array along the lines (this is how the data must be given to p9
        mapArrF = np.flipud(self.mapArr)
        df = pd.DataFrame(mapArrF).reset_index().melt('index')
        df.columns = ['row', 'col', 'val']
        
        mapplot = (
                ggplot(df) +
                geom_tile(aes(x = 'col', y = 'row', fill = 'val')) 
            )

        print(mapplot)

        SNULEBULLE
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
        '''
                
        
        
        