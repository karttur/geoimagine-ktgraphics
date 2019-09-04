'''
Created on 14 Dec 2018

@author: thomasgumbricht
'''

import os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as pyplot
from math import exp

class GraphPlotV01():
    ''' For pyplot graphs
    '''
    def __init__(self,dstFP):
        ''' Set the dataframe, including the time index
        '''
        self.dstFP = dstFP
        
    def GraphDecomposeMKTSIni(self,aD,TS,origDF,trendDF,mk,slope,intercept,compinD, site = 'nosite'):
        #Calculate regression
        self.headerL = ["siteid","mann-kendall-Z","Theil-Sen-slope","lon" ,"lat"]
        self.valuesL = [site,mk.round(3),slope.round(3),aD['x'],aD['y']]
        if self.dstFP:
            FN = '%(site)s_%(p)s_data.csv' %{'site':site,'p':compinD.prefix}
            FPN = os.path.join(self.dstFP,FN)
            origDF.df.to_csv(FPN)

        tsr = []
        for i in range(trendDF.df['tendency'].values.shape[0]):
            tsr.append(intercept+slope*i)
        trendDF.df['regression'] = tsr
        titleD ={}
        figtitle = '%(site)s_decomp_%(p)s_%(decomp)s_%(trend)s' %{'site':site,'p':compinD.prefix,'decomp':aD['decompose'],'trend':aD['trend']}
        if aD['decompose'] == 'naive':
            if aD['trend'] == 'ma':
                smooth = 'moving average'
                titleD['seasons'] = 'Seasons (periods: %(k)s)' %{'k':TS.annualperiods}
            else:
                smooth = 'custom kernel'

        else:
            if aD['decompose'] =='trend':
                smooth = 'spline'
                titleD['seasons'] = 'Seasons (periods: %(k)s)' %{'k':TS.annualperiods}
            elif aD['decompose'] == 'cyclic':
                smooth = 'spline + kernel'
                titleD['seasons'] = 'Seasons (periods: %(k)s)' %{'k':TS.annualperiods}
                   
            elif aD['decompose'] ==  'detrend':
                smooth = 'spline + %(d)s' %{'d':TS.trend}
                titleD['detrended'] = 'Trend adjusted data (trend removed using %(y)d years spline)' %{'y':aD['periodfactor']}
                titleD['seasons'] = 'Seasons (periods: %(k)s)' %{'k':TS.annualperiods}
                         
        titleD['data'] = '%(p)s timeseries decomposition (smoothing: %(t)s)' %{'p': compinD.prefix,'t':smooth}
        if TS.trend == 'kernel':  
            titleD['tendency'] = 'Trend (filter size: %(n)d; MKtest Z: %(mk)1.2f; Theil-Sen slope: %(r)1.2f)' %{'n':TS.filtersize, 'mk':mk, 'r':slope}
    
        else:
            titleD['tendency'] = 'Trend (filter size: %(n)d; MKtest Z: %(mk)1.2f; Theil-Sen slope: %(r)1.2f)' %{'n':TS.periodfactor*TS.annualperiods, 'mk':mk, 'r':slope}

            
        titleD['residual'] = 'Residual (after removing trend and seasons)' 
        titleD['adjusted'] = 'Data adjusted for seasons (filter: spline; periods: %(p)s)' %{'p':TS.annualperiods}

        
        #Convert all to an ordered dict
        dataD = {}
        dataD['data'] = {'series':origDF.df.values, 'index':origDF.df.index, 'title':titleD['data']}
        for item in trendDF.df:
            
            if item == 'regression':
                continue
            if item not in titleD:
                continue
            dataD[item] = {'series':trendDF.df[item].values,'index':trendDF.df[item].index,'title':titleD[item]}
            if item == 'tendency':
                dataD[item]['regression'] = trendDF.df['regression'].values


        plotD = OrderedDict()
        self.plotOrder = ['data','adjusted','detrended','tendency','seasons','residual'] 
        for key in self.plotOrder:

            if key in dataD: 

                plotD[key] = dataD[key]

                
        return self.GraphDecomposeMKTSR(plotD,figtitle)
          
    def GraphDecomposeMKTSRtest(self,dataD,figtitle):
        fig, ax = pyplot.subplots(len(dataD), 1, figsize=(8,8), sharex=True)
        for a,key in zip(ax,dataD.keys() ):            
            a.plot(dataD[key]['series'].plot(),label=key)
            if 'regression' in dataD[key]:
                a.plot(dataD[key]['regression'])
            a.set_title(dataD[key]['title'])
            a.legend(loc=2) #2 is upper left
            # add labels/titles and such here 
        #fig.subplots_adjust(hspace=0.1)
        pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fig.canvas.set_window_title(figtitle)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.25)
        pyplot.show()
          
    def GraphDecomposeMKTSR(self,dataD,figtitle):
        figsizeY = min(len(dataD)*2,8)
        fig, ax = pyplot.subplots(len(dataD), 1, figsize=(8,figsizeY), sharex=True)
        #a = 0
      
        for a,key in zip(ax,dataD.keys() ):
            x = dataD[key]['index']
            y = dataD[key]['series']
            a.plot(x,y,label=key)
            if 'regression' in dataD[key]:
                x = dataD[key]['index']
                z = dataD[key]['regression']
                a.plot(x,z)
            #title = 'Decomposition method %(b)s (mk = %(mk)1.2f; ts = %(ts)1.2f)' %{'b':dataD[key]['band'],'mk':dataD[key]['mk'],'ts':dataD[key]['a']}
            a.set_title(dataD[key]['title'])
            a.legend(loc=2) #2 is upper left

        #fig.subplots_adjust(hspace=0.1)
        pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fig.canvas.set_window_title(figtitle)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.25)
        if self.dstFP:
            self.ResultsToFile(figtitle)
            '''
            FN = '%(f)s.png' %{'f':figtitle}
            FPN = os.path.join(self.dstFP, FN)
            pyplot.savefig(FPN, bbox_inches='tight')
            FN = '%(f)s.csv' %{'f':figtitle}
            FPN = os.path.join(self.dstFP, FN)
            CsvWrite(FPN,self.headerL,[self.valuesL])
            '''
        else:
            pyplot.show()
        pyplot.close()
        return dict(zip(self.headerL,self.valuesL))

    def GraphCrossCorrMultiTrends(self,compinD,aD,dataD,item, site='nosite'):
        self.headerL = ["siteid"]
        self.valuesL = [site]
        #figtitle = '%(site)s_crosscorr_%(p)s-%(i)s_%(decomp)s_%(trend)s' %{'site':aD['siteid'],'p':compinD.prefix,'i':item,'decomp':aD['decompose'],'trend':aD['trend']}
        figtitle = '%(site)s_%(p)s_%(decomp)s_%(trend)s' %{'site':site,'p':compinD.prefix,'decomp':aD['decompose'],'trend':aD['trend']}

        fig, ax = pyplot.subplots(len(dataD), 1, figsize=(8,8), sharex=True)
        for a,key in zip( ax,dataD.keys() ):
            head  = '%(item)s-%(i)s-lag' %{'item':item,'i': dataD[key]['i']}
            self.headerL.append(head)
            self.valuesL.append(dataD[key]['lag'])
            
            #head  = '%(i)s-corr' %{'i': dataD[key]['i']}
            head  = '%(item)s-%(i)s-corr' %{'item':item,'i': dataD[key]['i']}
            self.headerL.append(head)
            self.valuesL.append(dataD[key]['corr'].round(3))
            
            head  = '%(item)s-%(i)s-pearson' %{'item':item,'i': dataD[key]['i']}
            self.headerL.append(head)
            self.valuesL.append(dataD[key]['pearson'].round(3))
            
            #x = dataD[key]['index']
            y = dataD[key]['value']
            z = dataD[key]['item']
            #n = len(y)
            a.plot(dataD[key]['index'],y)
            a.plot(dataD[key]['crossindex'],z)
            title = 'crosscorr %(b)s vs %(i)s %(f)s (lag = %(l)d; corr = %(c)1.2f; pearson = %(p)1.2f)' %{'b':compinD.prefix,'i':dataD[key]['i'],'f':item,'l':dataD[key]['lag'],'c':dataD[key]['corr'], 'p':dataD[key]['pearson']}
            a.set_title(title)

        pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fig.canvas.set_window_title(figtitle)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.25)
        '''
        fig.subplots_adjust(hspace=0.25)
        pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        #pyplot.xlim((index[0],index[index.shape[0]-1]))
        fig.tight_layout()
        '''
        
        #fig.canvas.set_window_title(figtitle)
                
        if self.dstFP:
            figtitle = '%(site)s_crosscorr_%(p)s-%(i)s_%(decomp)s_%(trend)s' %{'site':site,'p':compinD.prefix,'i':item,'decomp':aD['decompose'],'trend':aD['trend']}
            self.ResultsToFile(figtitle)
            '''
            FN = '%(f)s.png' %{'f':figtitle}
            FPN = os.path.join(self.dstFP, FN)
            pyplot.savefig(FPN, bbox_inches='tight')
            
            FN = '%(f)s.csv' %{'f':figtitle}
            FPN = os.path.join(self.dstFP, FN)
            CsvWrite(FPN,self.headerL,[self.valuesL])
            '''
            
        else:
            pyplot.show()
        pyplot.close()
            
        return dict(zip(self.headerL,self.valuesL))
             
    def ResultsToFile(self,filename):
        filename = filename.replace(' ','-')
        FN = '%(f)s.png' %{'f':filename}
        FPN = os.path.join(self.dstFP, FN)
        pyplot.savefig(FPN, bbox_inches='tight')
        FN = '%(f)s.csv' %{'f':filename}
        FPN = os.path.join(self.dstFP, FN)
        CsvWrite(FPN,self.headerL,[self.valuesL])
        
class GraphPlotv02():
    def __init__(self):
        pyplot.ioff()
        pyplot.show()
        #Set the linestyles
        self._Linestyles()
            
    def _PdGraph(self,df):
        #pyplot.figure() 
        df.plot()
        pyplot.show()
    
    def _GraphSingleTimeLine(self,df,params):
        '''
        pyplot.figure()
        self.ts = ts
        title = figtitle
        pyplot.plot(self.ts, label='data')
        pyplot.suptitle(title)
        pyplot.legend(loc='upper left')
        pyplot.show()
        '''
        self.params = params
        #self.pyplot.plot(df.ts)
        #print ('df.ts',df.ts)
        #print ('df.date',df.date)
        form = self._SetLineFormat()
        #fig, ax = pyplot.subplots()
        if form:
            line, = pyplot.plot(df.date, df.ts, form)
        else:
            line, = pyplot.plot(df.date, df.ts, 
                        color=self.params.color, linewidth=self.params.linewidth, linestyle=self.params.linestyle )
        pyplot.xlabel(self.params.xlabel)
        pyplot.ylabel(self.params.ylabel)
        pyplot.title(self.params.title)
        pyplot.grid(params.grid)
        #plt.savefig("test.png")
        pyplot.show()
        
    def _GraphMultiTimeLine(self,df,paramsD,params):
        '''
        '''
        for item in paramsD:
            if paramsD[item]['obsformat'] == '0':
                line, = pyplot.plot(df.dates, df[item], color=paramsD[item]['obscolor'], 
                            linewidth=params.obslinewidth, linestyle=self.linestyleD[params.obslinestyle], label=item )
            else:
                line, = pyplot.plot(df.dates, df[item], paramsD[item]['obsformat'],
                            linewidth=params.obslinewidth, linestyle=self.linestyleD[params.obslinestyle], label=item )
        pyplot.ylabel(params.ylabel)
        pyplot.title(params.title)
        pyplot.grid(params.grid)
        if params.legend >= 0:
            pyplot.legend(loc=params.legend)
        #pyplot.savefig("test.png")
        pyplot.show()

    def _GraphMultiTimeLineSubPlots(self,df,paramsD,params):
        '''
        '''
        fig, axarr = pyplot.subplots(len(paramsD), sharex=True, figsize=(params.figwidth,params.figheight))
        fig.suptitle(params.title)
        for i , item in enumerate(paramsD):
            if paramsD[item]['format'] == '0':
                line, = axarr[i].plot(df.dates, df[item], 
                            color=paramsD[item]['color'], linewidth=paramsD[item]['linewidth'], linestyle=paramsD[item]['linestyle'] )
            else:
                line, = axarr[i].plot(df.dates, df[item], paramsD[item]['format'])

            axarr[i].set_ylabel(item)
   
        #only set the xlabel for the past subplot as sharex is set to True
        axarr[i].set_xlabel(params.xlabel)
        
        fig.subplots_adjust(hspace=0.1)

        #pyplot.savefig("test.png")
        pyplot.show()

    def _GraphMultiTimeLineDeCompose(self,dfD,regressD,paramsD,params):
        '''
        '''
        fig, axarr = pyplot.subplots(len(paramsD), sharex=True, figsize=(params.figwidth,params.figheight))
        fig.suptitle(params.title)

        for i , item in enumerate(paramsD):
            if params.plotoriginal:
                if paramsD[item]['obsformat'] == '0':
                    line, = axarr[i].plot(dfD[item]['dates'], dfD[item]['observed'], 
                                color=paramsD[item]['obscolor'], linewidth=params.obslinewidth, linestyle=self.linestyleD[params.obslinestyle], label= 'observed' )
                else:
                    line, = axarr[i].plot(dfD[item]['dates'], dfD[item]['observed'], 
                                paramsD[item]['obsformat'], linewidth=params.obslinewidth, linestyle=self.linestyleD[params.obslinestyle], label= 'observed' )
            
            if params.plottendency:
                if paramsD[item]['tendencyformat'] == '0':
                    line, = axarr[i].plot(dfD[item]['dates'], dfD[item]['tendency'], 
                                color=paramsD[item]['tendencycolor'], linewidth=params.tendencylinewidth, linestyle=self.linestyleD[params.tendencylinestyle], label= 'tendency')
                else:
                    line, = axarr[i].plot(dfD[item]['dates'], dfD[item]['tendency'], 
                                paramsD[item]['tendencyformat'], linewidth=params.tendencylinewidth, linestyle=self.linestyleD[params.tendencylinestyle], label= 'tendency'  )

            if params.plottrend:
                tsr = []
                for x in range(len(dfD[item]['dates'])):            
                    tsr.append(regressD[item]['trend']['tsintercept']+regressD[item]['trend']['tsslope']*x)
                regrA = np.array(tsr)
                               
                if paramsD[item]['trendformat'] == '0':
                    line, = axarr[i].plot(dfD[item]['dates'], regrA, 
                                color=paramsD[item]['trendcolor'], linewidth=params.trendlinewidth, linestyle=self.linestyleD[params.trendlinestyle], label= 'linear trend')
                else:
                    line, = axarr[i].plot(dfD[item]['dates'], regrA, 
                                paramsD[item]['trendformat'], linewidth=params.trendlinewidth, linestyle=self.linestyleD[params.trendlinestyle], label= 'linear trend')

            axarr[i].set_ylabel(item)
            if params.legend >= 0:
                axarr[i].legend(loc=params.legend)

        #only set the xlabel for the past subplot as sharex is set to True
        axarr[i].set_xlabel(params.xlabel)
        
        fig.subplots_adjust(hspace=0.1)

        #pyplot.savefig("test.png")
        pyplot.show()
  
    def _GraphMultiTimeLineDeComposeSingle(self,dfD,regressD,paramsD,params):
        '''
        '''
        fig, axarr = pyplot.subplots(1, sharex=True, figsize=(params.figwidth,params.figheight))
        fig.suptitle(params.title)
         
        for item in paramsD:
            
            if params.plotoriginal:
                if paramsD[item]['obsformat'] == '0':
                    line, = axarr.plot(dfD[item]['dates'], dfD[item]['observed'], 
                                color=paramsD[item]['obscolor'], linewidth=params.obslinewidth, linestyle=self.linestyleD[params.obslinestyle], label= '%s (observed)' %(item) )
                else:
                    line, = axarr.plot(dfD[item]['dates'], dfD[item]['observed'], 
                                paramsD[item]['obsformat'], linewidth=params.obslinewidth, linestyle=self.linestyleD[params.obslinestyle], label= '%s (observed)' %(item) )
            
            elif params.plottendency:
                if paramsD[item]['tendencyformat'] == '0':
                    line, = axarr.plot(dfD[item]['dates'], dfD[item]['tendency'], 
                                color=paramsD[item]['tendencycolor'], linewidth=params.tendencylinewidth, linestyle=self.linestyleD[params.tendencylinestyle], label= '%s (tendency)' %(item) )
                else:
                    line, = axarr.plot(dfD[item]['dates'], dfD[item]['tendency'], 
                                paramsD[item]['tendencyformat'], linewidth=params.tendencylinewidth, linestyle=self.linestyleD[params.tendencylinestyle], label= '%s (tendency)' %(item) )

            elif params.plottrend:
                tsr = []
                for x in range(len(dfD[item]['dates'])):            
                    tsr.append(regressD[item]['trend']['tsintercept']+regressD[item]['trend']['tsslope']*x)
                regrA = np.array(tsr)
                               
                if paramsD[item]['trendformat'] == '0':
                    line, = axarr.plot(dfD[item]['dates'], regrA, 
                                color=paramsD[item]['trendcolor'], linewidth=params.trendlinewidth, linestyle=self.linestyleD[params.trendlinestyle], label= '%s (linear trend)' %(item) )
                else:
                    line, = axarr[0].plot(dfD[item]['dates'], regrA, 
                                paramsD[item]['trendformat'], linewidth=params.trendlinewidth, linestyle=self.linestyleD[params.trendlinestyle], label= '%s (linear trend)' %(item) )
    
        axarr.set_ylabel(params.ylabel)
        if params.legend >= 0:
            axarr.legend(loc=params.legend)
        pyplot.grid(params.grid)
        #only set the xlabel for the past subplot as sharex is set to True
        axarr.set_xlabel(params.xlabel)
        #pyplot.savefig("test.png")
        pyplot.show()

    def _BarChartSingle(self, xpos, ydataD, paramsD, params):
        items = len(paramsD)
        width = params.width/items
        leftadjust = width*(items-1.0)/2.0
        for i,item in enumerate(paramsD):
            if params.errbars:
                pyplot.bar(xpos+width*i-leftadjust, ydataD[item]['acf'], width=width, color=paramsD[item]['color'], 
                           edgecolor=paramsD[item]['edgecolor'], alpha = params.alpha, yerr=ydataD[item]['err'], label = item)
            else:
                pyplot.bar(xpos+width*i-leftadjust, ydataD[item]['acf'], width=width, color=paramsD[item]['color'], 
                           edgecolor=paramsD[item]['edgecolor'], alpha = params.alpha, label = item)
        pyplot.xticks(xpos)
        pyplot.ylabel(params.ylabel)
        pyplot.xlabel(params.xlabel)
        pyplot.title(params.title)
        if params.legend >= 0:
            pyplot.legend(loc=params.legend)
        pyplot.grid(params.grid)
         
        #plt.xticks(index + bar_width / 2, ('A', 'B', 'C', 'D', 'E')) 
         
        pyplot.show()

    def _GraphDecomposedTrend(self, df, columns, regressD, layoutD, params):
        '''
        '''
        self.layoutD = layoutD
        self.params = params
        items = list(df.keys())
        fig, axarr = pyplot.subplots(len(columns), sharex=True, sharey=params.sharey, figsize=(params.figwidth, params.figheight))
        fig.suptitle(params.title)
        for a , key in enumerate(columns):
            for item in items:
                form = self._SetDecomposeLineFormat(key,item)
                if form != '0':
                    line, = axarr[a].plot(df[item].dates, df[item][key], form, label=item)
                else:
                    line, = axarr[a].plot(df[item].dates, df[item][key], 
                            color=self.params.color, linewidth=self.params.linewidth, 
                            linestyle=self.params.linestyle, label=item )
            title = '%(b)s;' %{'b':key}

            if key in regressD[item]:
                title = '%(title)s %(d)s ' %{'title':title,'d':self.params.trend}
                
                for item in items:
                    tsr = []
                    for i in range(len(df[item][key])):            
                        tsr.append(regressD[item][key]['tsintercept']+regressD[item][key]['tsslope']*i)
                        
                    regrA = np.array(tsr)

                    if not self.params.additive and key == 'tendency':
                        regrA = np.exp(regrA)
                        regressD[item][key]['tsslope'] = exp(regressD[item][key]['tsslope'])
                        #line, = axarr[a].plot(df[item].dates, df[item]['ts'], 'g-', label='ts')

                       
                    form = self._SetDecomposeLineFormat('regress',item)
                    label = 'trend %(item)s (mk = %(mk)1.2f; ts = %(ts)1.2f)' %{'item':item, 'mk':regressD[item][key]['mk'],'ts':regressD[item][key]['tsslope']}
 
                    if form != '0':
                        line, = axarr[a].plot(df[item].dates, regrA, form, label=label)
                    else:
                        line, = axarr[a].plot(df[item].dates, regrA, color=self.params.color, 
                                linewidth=self.params.linewidth, linestyle=self.params.linestyle, label=label ) 
                    '''   
                    if key == 'tendency':
                        tsr = []
                        for i in range(len(df[item][key])):            
                            tsr.append(306.0535077560561+0.12876836561228838*i)
                        regrA = np.array(tsr)
                        line, = axarr[a].plot(df[item].dates, regrA, 'b-', label=label)
                    '''

            else:
                if key == 'seasons':
                    title = '%(b)s; ' %{'b':key}
                    if self.params.prefilterseason:
                        title = '%(t)s (prefiltered)' %{'t':title}
                    if self.params.forceseason:
                        title = '%(t)s (forced)' %{'t':title}

            axarr[a].set_title(title)
            axarr[a].set_ylabel(self.params.ylabel)
            if params.legend >= 0:
                axarr[a].legend(loc=params.legend)

        #only set the xlabel for the past subplot as sharex is set to True
        axarr[a].set_xlabel(self.params.xlabel)
        
        fig.subplots_adjust(hspace=0.25)
        #pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        #fig.canvas.set_window_title(params.title)
        pyplot.show()
   
    def _GraphXCross(self, dfD, xcrossD, layoutD, indexLayoutD, params):
        self.layoutD = layoutD
        self.indexLayoutD = indexLayoutD
        self.params = params
        onlyone = False
        if len(dfD) == 1:
            onlyone = True
        fig, axarr = pyplot.subplots(len(dfD), sharex=True, sharey=params.sharey, figsize=(params.figwidth, params.figheight))
        
        title = params.title
        
        lagmax = params.maxlag 
        if params.mirrolag:
            lagmin = -params.maxlag   
        else:
            lagmin = 0
        title = '%(t)s; lag:%(min)d to %(max)d ' %{'t':title,'min':lagmin,'max':lagmax}
        
        if params.lagplot:
            title = '%(t)s; lag adj' %{'t':title}
            
        if params.abs:
            title = '%(t)s; abs' %{'t':title}
            
        if params.forceseason:
            title = '%(t)s; forced season' %{'t':title}
            
        if params.naive:
            title = '%(t)s; naive' %{'t':title}
        else:
            title = '%(t)s; %(x)s' %{'t':title,'x':params.trend}
            
            
        fig.suptitle(title)
        
        for a, xc in enumerate(dfD):
            #in each graph, do the index just once
            doIndex = True
            for xy in dfD[xc]:
                
                for i in dfD[xc][xy]: 
                    if doIndex:
                        iform =  self._SetXcrossIndexFormat(xc,i)
                        if iform != '0': 
                            if onlyone: 
                                line, = axarr.plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrIndex, iform, label=i)
                            else:
                                line, = axarr[a].plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrIndex, iform, label=i)
                        else:
                            if onlyone:
                                line, = axarr.plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrIndex, 
                                    color=self.params.color, linewidth=self.params.linewidth, 
                                    linestyle=self.linestyleD[self.params.linestyle], label=i )
                            else:
                                line, = axarr[a].plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrIndex, 
                                        color=self.params.color, linewidth=self.params.linewidth, 
                                        linestyle=self.linestyleD[self.params.linestyle], label=i )
                        doIndex = False
                    obsform = self._SetXcrossObsFormat(xc,xy)

                    label = '%(xy)s; lag:%(lag)s psonr:%(pearson)1.2f' %{'xy':xy, 
                                'lag':xcrossD[xc][xy][i]['lag'],'pearson':xcrossD[xc][xy][i]['pearson'] }
                    if obsform != '0': 
                        if onlyone: 
                            line, = axarr.plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrObs, obsform, label=label)
                        else:
                            line, = axarr[a].plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrObs, obsform, label=label)

                    else:
                        if onlyone:
                            line, = axarr.plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrObs, 
                                    color=self.params.color, linewidth=self.params.linewidth, 
                                    linestyle=self.linestyleD[self.params.linestyle], label=label )
                        else:
                            line, = axarr[a].plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrObs, 
                                    color=self.params.color, linewidth=self.params.linewidth, 
                                    linestyle=self.linestyleD[self.params.linestyle], label=label )
            if onlyone:       
                axarr.set_title(xc)
                axarr.set_ylabel(self.params.ylabel)
                if params.legend >= 0:
                    axarr.legend(loc=params.legend)
            else:
                axarr[a].set_title(xc)
                axarr[a].set_ylabel(self.params.ylabel)
                if params.legend >= 0:
                    axarr[a].legend(loc=params.legend)

        #only set the xlabel for the past subplot as sharex is set to True
        if onlyone:
            axarr.set_xlabel(self.params.xlabel)
        else:
            axarr[a].set_xlabel(self.params.xlabel)
            fig.subplots_adjust(hspace=0.25)
        #pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        #fig.canvas.set_window_title(params.title)
        pyplot.show()

        
    def _GraphXCrossMultiIndex(self, dfD, xcrossD, layoutD, indexLayoutD, params):
        self.layoutD = layoutD
        self.indexLayoutD = indexLayoutD
        self.params = params
        title = params.title
        
        for xc in dfD:
            title = '%(t)s; %(xc)s' %{'t':title,'xc':xc}
            #Create a list of the positions to plot (xy)
            xyL = list(dfD[xc].keys())
            for xy in dfD[xc]:
                
                fig, axarr = pyplot.subplots(len(dfD[xc][xy]), sharex=True, sharey=params.sharey, figsize=(params.figwidth, params.figheight))
                break
           
        lagmax = params.maxlag 
        if params.mirrorlag:
            lagmin = -params.maxlag   
        else:
            lagmin = 0
        title = '%(t)s; lag:%(min)d to %(max)d ' %{'t':title,'min':lagmin,'max':lagmax}
        
        if params.lagplot:
            title = '%(t)s; lag adj' %{'t':title}
            
        if params.abs:
            title = '%(t)s; abs' %{'t':title}
            
        if params.forceseason:
            title = '%(t)s; forced season' %{'t':title}
           
        if params.naive:
            if params.additive:
                title = '%(t)s; Additive' %{'t':title}
            else:
                title = '%(t)s; Multiplicative' %{'t':title}
            title = '%(t)s; naive' %{'t':title}
        else:
            if params.prefilterseason:
                title = '%(t)s; filtered season' %{'t':title}
            title = '%(t)s; %(x)s' %{'t':title,'x':params.trend}
            
        title = '%(t)s; %(y)d yrs' %{'t':title,'y':params.yearfac}
            
            
        fig.suptitle(title)
     
        #only one type of component
        for xc in dfD:
            
            #in each graph, do the index just once
            for a,i in enumerate(dfD[xc][xyL[0]]): 
                doIndex = True
                for xy in xyL:
                    if doIndex:
                        iform =  self._SetXcrossIndexFormat(xc,i)
                        if iform != '0':  
                            line, = axarr[a].plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrIndex, iform, label=i)
                        else:
                            line, = axarr[a].plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrIndex, 
                                    color=self.params.color, linewidth=self.params.linewidth, 
                                    linestyle=self.linestyleD[self.params.linestyle], label=i )
                        doIndex = False
                    obsform = self._SetXcrossObsFormat(xc,xy)

                    label = '%(xy)s; lag:%(lag)s psonr:%(pearson)1.2f' %{'xy':xy, 
                                'lag':xcrossD[xc][xy][i]['lag'],'pearson':xcrossD[xc][xy][i]['pearson'] }
                    if obsform != '0':  
                        line, = axarr[a].plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrObs, obsform, label=label)
                    else:
                        line, = axarr[a].plot(dfD[xc][xy][i].dates, dfD[xc][xy][i].corrObs, 
                                color=self.params.color, linewidth=self.params.linewidth, 
                                linestyle=self.linestyleD[self.params.linestyle], label=label )
                    
                axarr[a].set_ylabel(xc)
                if params.legend >= 0:
                    axarr[a].legend(loc=params.legend)

        #only set the xlabel for the past subplot as sharex is set to True
        axarr[a].set_xlabel(self.params.xlabel)
        fig.subplots_adjust(hspace=0.15)
        #pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        #fig.canvas.set_window_title(params.title)
        pyplot.show()
   
    def _SetLineFormat(self):
        #print ('self.params.format',self.params.format)
        if self.params.format == '0':
            #format is a digit and zero, that means that no format is set
            return False
        return self.params.format
    
    def _SetDecomposeLineFormat(self,component,item):
        #'observed', 'trend', 'seasons', 'residual'
        if component == 'observed':
            form = self.layoutD[item]['obsformat']; linewidth = self.params.obslinewidth;
            linestyle = self.params.obslinestyle; color  = self.layoutD[item]['obscolor'];
        elif component == 'regress':
            form = self.layoutD[item]['regressformat']; linewidth = self.params.regresslinewidth;
            linestyle = self.params.regresslinestyle; color  = self.layoutD[item]['regresscolor'];
        elif component == 'tendency':
            form = self.layoutD[item]['tendencyformat']; linewidth = self.params.tendencylinewidth;
            linestyle = self.params.tendencylinestyle; color  = self.layoutD[item]['tendencycolor'];
        elif component == 'seasons':
            form = self.layoutD[item]['seasonformat']; linewidth = self.params.seasonlinewidth;
            linestyle = self.params.seasonlinestyle; color  = self.layoutD[item]['seasoncolor'];
        elif component == 'residual':
            form = self.layoutD[item]['residualformat']; linewidth = self.params.residuallinewidth;
            linestyle = self.params.residuallinestyle; color  = self.layoutD[item]['residualcolor'];
        else:
            print (component)
            SNULLE
        #print ('self.params.format',self.params.format)
        if form == '0':
            self.params.linewidth = linewidth
            self.params.linestyle = linestyle
            self.params.color = color
            
            return form
        self.params.format = form
        return self.params.format

    def _SetXcrossObsFormat(self,component,item):
        #'observed', 'trend', 'seasons', 'residual'
        if component == 'observed':
            form = self.layoutD[item]['obsformat']; linewidth = self.params.obslinewidth;
            linestyle = self.params.obslinestyle; color  = self.layoutD[item]['obscolor'];
        elif component == 'regress':
            form = self.layoutD[item]['regressformat']; linewidth = self.params.regresslinewidth;
            linestyle = self.params.regresslinestyle; color  = self.layoutD[item]['regresscolor'];
        elif component == 'tendency':
            form = self.layoutD[item]['tendencyformat']; linewidth = self.params.tendencylinewidth;
            linestyle = self.params.tendencylinestyle; color  = self.layoutD[item]['tendencycolor'];
        elif component == 'seasons':
            form = self.layoutD[item]['seasonformat']; linewidth = self.params.seasonlinewidth;
            linestyle = self.params.seasonlinestyle; color  = self.layoutD[item]['seasoncolor'];
        elif component == 'residual':
            form = self.layoutD[item]['residualformat']; linewidth = self.params.residuallinewidth;
            linestyle = self.params.residuallinestyle; color  = self.layoutD[item]['residualcolor'];
        else:
            print (component)
            SNULLE
        #print ('self.params.format',self.params.format)
        if form == '0':
            self.params.linewidth = linewidth
            self.params.linestyle = linestyle
            self.params.color = color
            
            return form
        self.params.format = form
        return self.params.format

    def _SetXcrossIndexFormat(self,component,item):
        #'observed', 'trend', 'seasons', 'residual'
        if component == 'observed':
            form = self.indexLayoutD[item]['obsformat']; linewidth = self.params.obslinewidthindex;
            linestyle = self.params.obslinestyleindex; color  = self.indexLayoutD[item]['obscolor'];
        elif component == 'tendency':
            form = self.indexLayoutD[item]['tendencyformat']; linewidth = self.params.tendencylinewidthindex;
            linestyle = self.params.tendencylinestyleindex; color  = self.indexLayoutD[item]['tendencycolor'];
        elif component == 'seasons':
            form = self.indexLayoutD[item]['seasonformat']; linewidth = self.params.seasonlinewidthindex;
            linestyle = self.params.seasonlinestyleindex; color  = self.indexLayoutD[item]['seasoncolor'];
        elif component == 'residual':
            form = self.indexLayoutD[item]['residualformat']; linewidth = self.params.residuallinewidthindex;
            linestyle = self.params.residuallinestyleindex; color  = self.indexLayoutD[item]['residualcolor'];
        else:
            print (component)
            SNULLE
        #print ('self.params.format',self.params.format)
        if form == '0':
            self.params.linewidth = linewidth
            self.params.linestyle = linestyle
            self.params.color = color
            
            return form
        self.params.format = form
        return self.params.format
    
    def GraphDecomposeMKTSIniOld(self,TS,mk,slope,intercept,figtitle):
        self.ts = TS.ts
        tsr = []
        for i in range(len(TS.ts)):
            tsr.append(intercept+slope*i)
        regrA = np.array(tsr)
        title = 'Timeseries decomposition (type= %(d)s; period = %(p)d; smoothing = %(t)s)' %{'d':TS.decompose, 'p':TS.period, 't':TS.trend}
        if TS.decompose == 'trend': seastitle = 'Seasons (trend: %s)' %(TS.trend)
        elif TS.decompose == 'cyclic': seastitle = 'Seasons (trend: spline)'
        elif TS.decompose == 'naive': seastitle = 'Seasons (trend: naive)'
        elif TS.decompose == 'detrend': 
            if TS.splineseason:
                seastitle = 'Seasons (trend: spline)'
            else:
                seastitle = 'Seasons (trend: %s)' %(TS.trend)
        trendtitle = 'Trend (MKtest Z = %(mk)1.2f; regression = Theil Sen)' %{'i':TS.ts, 'mk':mk}

        self.GraphDecomposeMKTSR(TS,title,trendtitle,seastitle,regrA,TS.decompose,figtitle)
          
    def GraphSeasonalOld(self,TS, title, normalize=False, npA=False, labelA=False, npB=False, labelB=False, ts=False, tendency=False, seasons=False, adjusted=False, residual=False, longtendency=False , detrended=False): 
        pyplot.figure()
        if ts:
            if normalize:
                normA  = (TS.ts - mean(TS.ts)) /  std(TS.ts)
                pyplot.plot(normA, label='data')
            else:
                #plt.plot(data.index, data.amount)
                #pyplot.plot(TS.ts.index, TS.ts.data.amount)
                pyplot.plot(TS.ts, label='data')
        if tendency:
            if normalize:
                normA  = (self.tendency - np.average(self.tendency)) /  np.std(self.tendency)
                pyplot.plot(normA, label='trend')
            else:
                
                pyplot.plot(self.tendency, label='trend')
        if detrended:
            if normalize:
                normA  = (self.detrended - np.average(self.detrended)) /  np.std(self.tdetrended)
                pyplot.plot(normA, label='detrended')
            else:          
                pyplot.plot(self.detrended, label='detrended')
        if seasons:
            if normalize:
                normA  = (self.seasons - np.average(self.seasons)) /  np.std(self.seasons)
                pyplot.plot(normA, label='data')
            else:
                pyplot.plot(self.seasons, label='seasons')   
        if adjusted:
            pyplot.plot(self.adjusted, label='adjusted')
        if residual:
            pyplot.plot(self.residual, label='residual')
        if longtendency:
            if normalize:
                normA  = (TS.longtendency - np.average(TS.longtendency)) /  np.std(TS.longtendency)
                pyplot.plot(normA, label='data')
            else:
                pyplot.plot(TS.longtendency, label='long term tendency')
        if labelA:
            pyplot.plot(npA, label=labelA)
        if labelB:
            pyplot.plot(npB, label=labelB)
        if normalize:
            title = '%(t)s (normalized)' %{'t':title}
        pyplot.suptitle(title)
        pyplot.legend(loc='upper left')
        pyplot.show()
        
    def GraphPandasOld(self,TS, title, normalize=False, npA=False, labelA=False,npB=False, labelB=False, ts=False, tendency=False, seasons=False, adjusted=False, residual=False, longtendency=False , detrended=False): 
        pyplot.figure()
        if ts:
            if normalize:
                normA  = (TS.df - mean(TS.df)) /  std(TS.df)
                pyplot.plot(normA, label='data')
            else:
                pyplot.plot(TS.df.index, TS.df.values)
        if tendency:
            if normalize:
                normA  = (self.tendency - np.average(self.tendency)) /  np.std(self.tendency)
                pyplot.plot(normA, label='trend')
            else:
                
                pyplot.plot(TS.tendency, label='trend')
        if detrended:
            if normalize:
                normA  = (self.detrended - np.average(self.detrended)) /  np.std(self.tdetrended)
                pyplot.plot(normA, label='detrended')
            else:          
                pyplot.plot(self.detrended, label='detrended')
        if seasons:
            if normalize:
                normA  = (self.seasons - np.average(self.seasons)) /  np.std(self.seasons)
                pyplot.plot(normA, label='data')
            else:
                pyplot.plot(self.seasons, label='seasons')   
        if adjusted:
            pyplot.plot(self.adjusted, label='adjusted')
        if residual:
            pyplot.plot(self.residual, label='residual')
        if longtendency:
            if normalize:
                normA  = (TS.longtendency - np.average(TS.longtendency)) /  np.std(TS.longtendency)
                pyplot.plot(normA, label='data')
            else:
                pyplot.plot(TS.longtendency, label='long term tendency')
        if labelA:
            pyplot.plot(npA.df.index, npA.df.values)
            #pyplot.plot(npA, label=labelA)
        if labelB:
            pyplot.plot(npB.df.index, npB.df.values)
            #pyplot.plot(npB, label=labelB)
        if normalize:
            title = '%(t)s (normalized)' %{'t':title}
        pyplot.suptitle(title)
        pyplot.legend(loc='upper left')
        pyplot.show()
             
    def GraphDecomposeMKTSROld(self,TS,title,trendtitle,seastitle,tsr,decomp,figtitle):
        f, (ax1, ax2, ax3, ax4) = pyplot.subplots(4, sharex=True)
        if decomp == 'detrend':
            ax1.plot(TS.ts, label='detrended data')
        else:         
            ax1.plot(TS.ts, label='data')
        ax1.set_title(title)
        ax1.legend(loc='best')
        ax2.plot(TS.fullseasons, label='seasons')
        ax2.set_title(seastitle)
        ax2.legend(loc='best')
        ax3.plot(TS.residual, label='residual')
        ax3.set_title('residual')
        ax3.legend(loc='best')
        if decomp == 'detrend':
            ax4.plot(TS.longtendency, label='multi cycles trend')
        else:
            ax4.plot(TS.tendency, label='trend')
        ax4.plot(tsr)
        ax4.set_title(trendtitle)
        ax4.legend(loc='best')
        f.subplots_adjust(hspace=0.25)
        pyplot.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        f.canvas.set_window_title(figtitle)
        pyplot.show()
        
    def GraphDualDecomposeTrendsOld(self,periodfactor,title): 
        pyplot.figure() 
        pyplot.plot(self.trend, label='original trend')
        legstr = 'multiperiod trend (%s cycles)' %(periodfactor)
        pyplot.plot(self.longtrend, label=legstr)
        pyplot.plot(self.stationary, label='stationary data')
        pyplot.plot(self.stationarytrend, label='stationary trend')
        pyplot.legend(loc='best')
        pyplot.suptitle(title)
        pyplot.show()
        
    def GraphCrossCorrTrendsOld(self,y,x,title):
        pyplot.figure() 
        pyplot.plot(y, label='forcing')
        pyplot.plot(x, label='data')
        pyplot.legend(loc='best')
        pyplot.suptitle(title)
        pyplot.show()
        
    def GraphCrossCorrMultiTrendsOld(self, dataD, figtitle):
        fig, ax = pyplot.subplots(len(dataD), 1, figsize=(8,6))
        for a,key in zip(ax,dataD.keys() ):
            y = dataD[key]['index']
            z = dataD[key]['item']
            n = len(y)
            x = np.linspace(1,n,n)
            a.plot(x,y)
            a.plot(x,z)
            title = 'crosscorr %(b)s vs %(i)s (lag = %(l)d; corr = %(c)1.2f; pearson = %(p)1.2f)' %{'b':dataD[key]['band'],'i':dataD[key]['i'],'l':dataD[key]['lag'],'c':dataD[key]['corr'], 'p':dataD[key]['pearson']}
            a.set_title(title)
            # add labels/titles and such here 
        fig.subplots_adjust(hspace=0.25)
        pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fig.canvas.set_window_title(figtitle)
        pyplot.show()
        
    def GraphMultiDecompIniOld(self,dataD,figtitle,tendency=False):
        if tendency:
            for D in dataD:
                tsr = []
                for i in range(len(dataD[D]['ts'])):   
                    tsr.append(dataD[D]['b']+dataD[D]['a']*i)
                    regrA = np.array(tsr)
                dataD[D]['tsr'] = regrA
            self.GraphMultiTrends(dataD,figtitle)   
        else:
            self.GraphMultiDecomp(dataD,figtitle)
            #self.GraphMultiDecompAllInOne(dataD,figtitle)
        
    def GraphMultiDecompOld(self,dataD,figtitle):
        fig, ax = pyplot.subplots(len(dataD), 1, figsize=(8,8), sharex=True)

        for a,key in zip(ax,dataD.keys() ):
            x = dataD[key]['index'].index
            y = dataD[key]['index'].values
            #z = dataD[key]['item']
            #n = len(y)

            #x = np.linspace(1,n,n)
            a.plot(x,y)
            #a.plot(x,z)
            title = 'Decomposition method %(b)s' %{'b':dataD[key]['band']}
            a.set_title(title)
            # add labels/titles and such here 
        
        pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        #pyplot.setp(a.index, a.values, visible=False)
        #pyplot.plot(npB.df.index, npB.df.values)
        fig.canvas.set_window_title(figtitle)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.25)
        pyplot.show()
        
    def GraphMultiDecompAllInOneOld(self,dataD,figtitle):
        fig, ax = pyplot.subplots(len(dataD), 1, figsize=(8,8), sharex=True)

        for a,key in zip(ax,dataD.keys() ):
            x = dataD[key]['index'].index
            y = dataD[key]['index'].values
            #z = dataD[key]['item']
            #n = len(y)

            #x = np.linspace(1,n,n)
            #a.plot(x,y)
            #a.plot(x,z)
            #title = 'Decomposition method %(b)s' %{'b':dataD[key]['band']}
            #a.set_title(title)
            # add labels/titles and such here 
            pyplot.plot(x, y)
        #pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        #pyplot.setp(a.index, a.values, visible=False)
        #pyplot.plot(npB.df.index, npB.df.values)
        fig.canvas.set_window_title(figtitle)
        fig.tight_layout()
        #fig.subplots_adjust(hspace=0.25)
        pyplot.show()
        
    def GraphMultiTrendsOld(self,dataD,figtitle):
        fig, ax = pyplot.subplots(len(dataD), 1, figsize=(8,8), sharex=True)
        for a,key in zip(ax,dataD.keys() ):
            #y = dataD[key]['index']
            x = dataD[key]['index'].index
            y = dataD[key]['index'].values
            z = dataD[key]['tsr']

            a.plot(x,y)
            a.plot(x,z)
            title = 'Decomposition method %(b)s (mk = %(mk)1.2f; ts = %(ts)1.2f)' %{'b':dataD[key]['band'],'mk':dataD[key]['mk'],'ts':dataD[key]['a']}
            a.set_title(title)
            # add labels/titles and such here 
        #fig.subplots_adjust(hspace=0.1)
        pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fig.canvas.set_window_title(figtitle)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.25)
        pyplot.show()
        
    def GraphCellSeasonIniOld(self, origdataD, seasonD, residualD, tendencyD, figtit):

        dataD = {'origdata':origdataD,'season':seasonD, 'residual':residualD, 'trend':tendencyD}
        self.GraphCellSeason(dataD,figtit)
        
    def GraphCellSeasonOld(self,dataD,figtitle):
        fig, ax = pyplot.subplots(4, 1, figsize=(8,8), sharex=True)

        for a,key in zip(ax,dataD.keys() ):
            x = dataD[key]['index'].index
            y = dataD[key]['index'].values
            a.plot(x,y)
            #a.plot(x,z)
            if dataD[key]['comp'] == 'trend':
                title = 'Component: %(b)s (mk = %(mk)1.2f; ts = %(ts)1.2f)' %{'b':dataD[key]['band'],'mk':dataD[key]['mk'],'ts':dataD[key]['a']}

                tsr = []
                for i in range(len(dataD[key]['ts'])):   
                    tsr.append(dataD[key]['b']+dataD[key]['a']*i)
                    regrA = np.array(tsr)
                #dataD[D]['tsr'] = regrA
                a.plot(x,regrA)
            else:
                title = 'Component: %(b)s' %{'b':dataD[key]['comp']}
            a.set_title(title)
            # add labels/titles and such here 
        
        pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        #pyplot.setp(a.index, a.values, visible=False)
        #pyplot.plot(npB.df.index, npB.df.values)
        fig.canvas.set_window_title(figtitle)
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.25)
        pyplot.show()
        
    def _Linestyles(self): 
        self.linestyleD = OrderedDict(
        [('solid',               (0, ())),
         ('loosely dotted',      (0, (1, 10))),
         ('dotted',              (0, (1, 5))),
         ('densely dotted',      (0, (1, 1))),
    
         ('loosely dashed',      (0, (5, 10))),
         ('dashed',              (0, (5, 5))),
         ('densely dashed',      (0, (5, 1))),
    
         ('loosely dashdotted',  (0, (3, 10, 1, 10))),
         ('dashdotted',          (0, (3, 5, 1, 5))),
         ('densely dashdotted',  (0, (3, 1, 1, 1))),
    
         ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
         ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
         ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
        
if __name__ == "__main__":
    
    graph = GraphPlotv02()
    npA = [296, 170, 145, 21,  8,  9 , 16 , 20 , 13 , 59, 104, 224, 316, 665, 210, 36 , 24 , 11,
           16 , 24 ,  0 , 21 , 72 , 46, 152, 369, 416 , 64 , 59 , 69 , 34 , 22 , 13 ,  9, 121, 182,
           70, 478, 212 , 41 , 14 , 31 , 13 ,  3 ,  8 ,  7 , 91, 230 , 22 , 37 , 54 , 31 ,  4 , 30,
           89 ,  4 , 29, 219 , 30 , 27 , 92 , 51, 144 , 31 , 11, 133 , 10 ,  3 , 22 , 81 , 66 , 16,
           107 , 84, 228 , 89 , 15 , 17 ,  9 ,  0 ,  3 , 81 , 26, 164 , 68 , 66 , 65 , 14 , 26 ,  7,
           18 ,  1 , 58 ,  4 , 24, 393, 110, 100, 332 , 69 , 20, 105 ,  3 ,  3 ,  2 ,  2 , 11, 132,
           81, 249 , 58, 290 ,  5 , 38 , 12 ,  1 ,  0 ,  9 , 78, 246 , 92 , 28 , 84 , 16 , 11 , 19,
           6 ,  2 ,  3 ,  6 , 12, 293, 116, 334 , 75 , 14 , 13 ,  1 , 32 ,  5 ,  3 , 11 , 86 , 70,
           102, 106 , 28 , 66 , 26 ,  6 , 29 ,  2 ,  1 ,  7 , 38, 138, 161 , 73 , 67, 111 ,  7 ,  0,
           62 ,  4 ,  1 , 52 , 46 , 47, 142 , 35 , 81 , 88 ,  6 ,  5 ,  2 ,  1 , 15 , 25 , 20 , 60,
           172, 143 , 81 , 25 ,  6 , 15 , 24 ,  9 , 14 , 79 , 16, 182, 343, 108 , 35 , 82 ,  2 ,  2,
           3 , 51 ,  7 , 15 , 11, 202 , 47 , 17 , 37 , 24 ,  2 ,  4 ,  3 ,  1 ,  1 , 42 , 37, 159,
           55 , 37 , 39 , 21 , 24 ,  9 , 32 ,  4 ,  0 ,  4, 149, 122, 324, 246 , 44, 101 ,  4 , 15,
           16 ,  3 ,  0 , 33 , 20, 120]
    TS = np.array(npA)
    figtitle = 'snullebulle'
    graph._GraphSingleTimeLine(TS,figtitle)
    