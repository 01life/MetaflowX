#!/usr/bin/env python

import markdown,sys,plotly,codecs,dash_bio,base64,os
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering

from math import pi
import math
import pandas as pd
import numpy as np
import argparse as ap
from jinja2 import Environment, FileSystemLoader
import time
 
from bokeh.models import ColumnDataSource, CustomJS, TextInput,HoverTool,Button,Div,CategoricalColorMapper
from bokeh.models import ColorBar,FactorRange, LinearAxis, Range1d,BasicTicker, PrintfTickFormatter
from bokeh.layouts import column,row,gridplot
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.embed import components
from bokeh.plotting import figure
from bokeh.transform import linear_cmap


def ParsReceiver():
    p = ap.ArgumentParser(
    formatter_class=ap.RawDescriptionHelpFormatter,	
    description='''
-------------------
    Report the whole pipeline result\n
    Version:
        V1	2023-07-11
    Author:
        Liang Lifeng
------------------a'''
)
    ars = p.add_argument
    ars('-r',dest="result",action="store",help='report result list',required=True)
    ars('-c',dest="templatehtml",action="store",help='template html',required=True)
    ars('-j',dest="js",action="store",help='resport html javascript file',required=False,default='NA')
    ars('-t',dest="txt",action="store",help='resport html txt file',required=True)
    ars('-i',dest="plotOrder",action="store",help='resport html result plot order file',required=True)
    ars('-m',dest="img",action="store",help='image path, default:shellPath/img/',required=False,default=os.path.split(os.path.realpath(__file__))[0]+'/img/')
    ars('-p',dest="prefix",action="store",help='output file prefix ',required=False,default='Test')
    #ars('-M',dest="mode",action="store",help='the whole pipeline mode',required=False,default='0')
    ars('-o',dest="outPath",action="store",help='outpathh, default:./',required=False,default=os.getcwd())

    return vars(p.parse_args())

#Function
def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def file2Dir(onefile, oneDir,valueindex,keyOrder):
    with open(onefile,'r') as oneF:
        for a in oneF:
            al = a.strip().split('\t')
            oneDir[al[0]] = al[valueindex]
            keyOrder.append(al[0])

def get_file_line_count(file_path):
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        return len(lines)

def check_data_frame_more_than_three(dataFrame):
    # check row and col number more than 2
    num_rows, num_cols = dataFrame.shape
    if num_rows <= 2 or num_cols <= 2:
        return('''<p style="font-size: 12pt; color: red; font-style: italic;" >Oops! Number of Smaple less than 2, nor item less than 2.  Socouldn't show this figure.<br></p>''')
    else:
        return('Pass')


def checkResultCompletely(resultDir):
    completeResultItem = {"0":["readstat","contigstat","mpaspecies","krakenspecies","krakenspeciesT","genesetLenStat","genesetSampleStat","genesetAbundance","genesetCOG","genesetGO","genesetKEGG","genesetCAZY","rawbinQS","binInfo","eachBinFunction","gtdb","binDepth","binRPKM","binTPM"],
                            }
    failed =[]
    for r in resultDir:
        if float(resultDir[r]) < 2:
            failed.append(r)
    return(failed)

def read_first_line(filepath):
    with open(filepath, 'r') as f:
        first_line = f.readline().strip().split("\t")
        newTitle = {}
        newcolumn = []
        for f in first_line:
            ff = f.strip().split("(")[0]
            newTitle[ff] = f
            newcolumn.append(ff)
    return(newcolumn,newTitle)

def plotlyPCA(oneProfileF, pointcolor,topTitle,returnType,inputType):
    if inputType =="F":
        oneProfile=pd.read_table(oneProfileF,sep="\t",low_memory=False,index_col=0)
        geneT = oneProfile.T
        data_frame_check = check_data_frame_more_than_three(geneT)
        if data_frame_check != 'Pass':
            return(data_frame_check,None,None,None,None)
            
    else:
        geneT = oneProfileF.T
        data_frame_check = check_data_frame_more_than_three(geneT)
        if data_frame_check != 'Pass':
            return(data_frame_check,None,None,None,None)

    if len(geneT.index) >= 2 :
        if len(geneT.index) == 2:
            pca = PCA(n_components=2)
            components = pca.fit_transform(geneT)
            total_var = pca.explained_variance_ratio_.sum() * 100
            pc1Label = f'PC1 :{pca.explained_variance_ratio_[0]*100:.2f}%'
            pc2Label = f'PC2 :{pca.explained_variance_ratio_[1]*100:.2f}%'
            pc3Label = 'PC3 :0%'
            componentsPD = pd.DataFrame(components, columns =['pc1', 'pc2'])
            componentsPD['Sample'] = geneT.index
            #PCA
            onePCA = px.scatter(componentsPD, x="pc1", y="pc2",
                hover_data = 'Sample' ,
                title=f'Total Explained Variance: {total_var:.2f}%\n{topTitle}',
                labels={'pc1': pc1Label, 'pc2': pc2Label})
            onePCA.update_traces(marker_size=20, marker_color = pointcolor,opacity=0.7)
            onePCA.update_layout(margin=dict(l=500, r=500, b=100, t=100), height = 600, )
        

        else:
            pca = PCA(n_components=3)
            components = pca.fit_transform(geneT)
            total_var = pca.explained_variance_ratio_.sum() * 100
            pc1Label = f'PC1 :{pca.explained_variance_ratio_[0]*100:.2f}%'
            pc2Label = f'PC2 :{pca.explained_variance_ratio_[1]*100:.2f}%'
            pc3Label = f'PC3 :{pca.explained_variance_ratio_[2]*100:.2f}%'
            componentsPD = pd.DataFrame(components, columns =['pc1', 'pc2','pc3'])
            componentsPD['Sample'] = geneT.index
            pc3 = '%.4f'%float(pca.explained_variance_ratio_[2]*100)
            if float(pc3) != 0:
                    onePCA = px.scatter_3d(
                        componentsPD, x="pc1", y="pc2", z="pc3", 
                        hover_data='Sample',
                        title=f'Total Explained Variance: {total_var:.2f}%\n{topTitle}',
                        labels={'pc1': pc1Label, 'pc2': pc2Label, 'pc3': pc3Label},)
                    onePCA.update_traces(marker_size=5, marker_color = pointcolor,opacity=0.7)
                    onePCA.update_layout(margin=dict(l=1, r=1, b=100, t=100), height = 600,)

            
            else:
                onePCA = px.scatter(componentsPD, x="pc1", y="pc2",
                    hover_data = 'Sample' ,
                    title=f'Total Explained Variance: {total_var:.2f}%\n{topTitle}',
                    labels={'pc1': pc1Label, 'pc2': pc2Label})
                onePCA.update_traces(marker_size=20, marker_color = pointcolor,opacity=0.7)
                onePCA.update_layout(margin=dict(l=500, r=500, b=100, t=100), height = 600, )


        #onePCA.show()
        if returnType == "px":
            return(onePCA.to_html(full_html=False, include_plotlyjs=False))

        else:
            return(componentsPD,total_var,pc1Label,pc2Label,pc3Label)

    else:
        if returnType == "px":
            return('''<p style="font-size: 12pt; color: red; font-style: italic;" >Oops! Number of Smaple less than 2, couldn't show as 3D PCA.<br></p>''')
        else:
            return(None,None,None,None,None)


def bokehPCA(dataframe,pc1Label,pc2Label,pointcolor,figTitle):
    data_frame_check = check_data_frame_more_than_three(dataframe)
    if data_frame_check != 'Pass':
        return(data_frame_check)

    source = ColumnDataSource(dataframe)
    p = figure(width=200, height=200, tools="pan,wheel_zoom,box_select,reset", active_drag="box_select",title=figTitle,sizing_mode="stretch_width")
    showPoint = p.circle(x="pc1", y="pc2", fill_color=pointcolor, size=10, alpha=0.6, source=source)
    p.xaxis.axis_label = pc1Label
    p.yaxis.axis_label = pc2Label

    tooltips = [("SampleID", "@Sample"),]
    showPoint_hover_tool = HoverTool(renderers=[showPoint],tooltips=tooltips)
    p.add_tools(showPoint_hover_tool)

    return(p)


def multiPCA(filelistPath,pointcolor):
    plot_col_dir = {1:2,
                    2:2,
                    3:3,
                    4:2,
                    5:3,
                    6:3,
                    7:4,
                    8:4,
                    9:3,
                    10:2}
    fig_list = []
    for onefile in filelistPath.strip().split(","):
        coverM_method = onefile.strip().split("bin_")[1].split(".xls")[0]
        if coverM_method not in ["coverage_histogram","length"]:
            componentsPD,total_var,pc1Label,pc2Label,pc3Label =  plotlyPCA(onefile,pointcolor,coverM_method,"pd","F")
            if componentsPD is not None and not isinstance(componentsPD, str):
                onepca = bokehPCA(componentsPD,pc1Label,pc2Label,pointcolor,coverM_method)
                fig_list.append(onepca)
    if fig_list :
        colsNumber = plot_col_dir[len(fig_list)]
        rows = math.ceil(len(fig_list) / colsNumber) # 计算行数
        grid = gridplot(fig_list, ncols=colsNumber,height=rows*150,sizing_mode="stretch_width",merge_tools=False,toolbar_location=None)
        grid.margin = (50, 50, 50, 50)
        script, div = components(grid)
        multigridtxt= div+script
        return(multigridtxt)
    else:
        return('''<p style="font-size: 12pt; color: red; font-style: italic;" >Oops! Number of Smaple less than 2, couldn't show as 3D PCA.<br></p>''')

##Geneset Function
def summaryDB2Plot(oneCountDF,databaseName,pngHeight):
    oneCountDF_sorted = pd.read_csv(oneCountDF,sep="\t",low_memory=False,index_col = 0)
    oneCountDFBar = go.Figure(go.Bar(
        x=oneCountDF_sorted['log2geneNumber'],
        y=oneCountDF_sorted.index,
        # width=[0.8]*len(list(oneCountDF_sorted.index)),
        width= 0.7,
        marker=dict(
            color='rgba(50, 171, 96, 0.6)',
            line=dict(color='rgba(50, 171, 96, 1.0)',width=0.5),),
            orientation='h',)
            )

    oneCountDFBar.update_layout(
        legend=dict(orientation="h",yanchor="bottom",y=-0.2,xanchor="right",x=1,),
        xaxis=dict(title='log2(gene numbers)'),
        yaxis=dict(title=str(databaseName)+'|Level1 Category'),
        height=pngHeight
        )
    return(oneCountDFBar.to_html(full_html=False, include_plotlyjs=False))

def figTable(tableFile, valueList):
    tf = pd.read_csv(tableFile,sep="\t",low_memory=False)
    cellsValue=[]
    for f in valueList:
        cellsValue.append(tf[f])
    figtable = go.Figure(data=[go.Table(
        header=dict(values=list(tf.columns),fill_color='rgba(50, 171, 96, 1)',align='left'),
        cells=dict(values=cellsValue,fill_color='rgba(50, 171, 96, 0.1)',align='left'))]
        )

    figtable.update_layout(margin=dict(b=50, t=50), height = 20*len(list(tf.index))+150,)
    return(figtable.to_html(full_html=False, include_plotlyjs=False))

def contigFigure(contigFile):
    contig=pd.read_csv(contigFile,sep="\t",low_memory=False)

    contig['Sample'] = contig['Sample'].astype('str')

    #F1 contig's length histogram

    # set histogram borderline and labels
    bins = list(np.arange(2000, 21000, 500))
    bins.append(np.inf)
    labels = [f"{i}-{i+1000}" for i in range(2000, 21000, 500)]
    labels[-1] = ">20K"
    hist, edges = np.histogram(contig['Length'], bins=bins)
    contigLengthHistogram = px.histogram(x=labels,y=hist,color_discrete_sequence=['rgba(50, 171, 96, 0.6)'])
    contigLengthHistogram.update_layout(
        title_text='Total Contigs Number :' + str(len(contig['Length'])), # title of plot
        xaxis_title_text='Contigs Length', # xaxis label
        yaxis_title_text='Count', # yaxis label
        #bargap=0.05, # gap between bars of adjacent location coordinates
        #bargroupgap=0, # gap between bars of the same location coordinates
        margin=dict(l=50, r=50, b=50, t=50),
    )

    #F2 contig's GC histogram
    contigGCHistogram = px.histogram(contig, x="GC",nbins = 100,color_discrete_sequence=['rgba(50, 171, 96, 0.6)'])
    contigGCHistogram.update_layout(height=500,margin=dict(l=50, r=50, b=50, t=50),)

    #F3 each sample contig numbers
    #sample contig number and contig length
    sampleContigCount = contig['Sample'].value_counts()
    sampleNum = len(list(sampleContigCount.index))

    if sampleNum <= 100:
        cotigNumber = make_subplots(rows=1, cols=2,shared_yaxes=True, specs = [[{}, {}]],  horizontal_spacing = 0.01,column_widths=[0.7, 0.3])
        cotigNumber.add_box(x=np.log2(contig["Length"]), y=contig["Sample"],
            name = "log2(Contig Length)",
            row=1, col=1,
            width = 0.8, 
            orientation = "h",
            marker=dict(color="rgba(50, 171, 96, 0.6)"))
        cotigNumber.add_bar(x=list(sampleContigCount.values), y=list(sampleContigCount.index),
                name = "Contig Number",
                row=1, col=2, 
                width = 0.8,
                orientation = "h",
                marker=dict(color="rgba(50, 171, 96, 0.6)"))
        cotigNumber.update_layout(
            showlegend=False,
            height = 10*sampleNum+200,
            margin=dict(l=50, r=50, b=50, t=50),
        )
        cotigNumber.update_xaxes(title_text="log2(Contig Length)", row=1, col=1)
        cotigNumber.update_xaxes(title_text="Contig Number", row=1, col=2)
        cotigNumber_html = cotigNumber.to_html(full_html=False, include_plotlyjs=False)

    else:
        sampleContig = contig.groupby("Sample").agg({"Length": "mean", "GC": "mean", "ContigID": "count"}).reset_index()
        sampleContig.columns = ["Sample", "avg_Length", "avg_GC", "Contig_Count"]
        source = ColumnDataSource(sampleContig)
        # Define the categories and data
        samples = sampleContig['Sample'].tolist()
        contig_numbers = sampleContig['Contig_Count'].tolist()
        contig_average_lengths = sampleContig['avg_Length'].tolist()

        # Create the figure
        fill_color =  "rgba(50, 171, 96, 0.6)"
        line_color = "rgba(50, 171, 96, 1.0)"
        toolslist = "hover,pan,box_zoom,wheel_zoom,undo,redo,reset,save"
        p = figure(x_range=FactorRange(*samples), height=600,
                sizing_mode="stretch_width",
                y_axis_label="Contig Metrics", 
                tools=toolslist,
                tooltips=[('sample', '@Sample'),('ContigNumber', '@Contig_Count'),('ContigAverageLength', '@avg_Length')]
                )

        # Plot the gene numbers
        p.vbar(x='Sample', top='Contig_Count', width=0.8,
            source=source,
            fill_color=fill_color, fill_alpha=0.8, line_color=line_color, line_width=1.2,
            legend_label="Total Contig Number of Sample")

        # Add the second y-axis for gene average length
        p.extra_y_ranges = {"average_contig_length": Range1d(start=0, end=max(contig_average_lengths))}
        p.add_layout(LinearAxis(y_range_name="average_contig_length", axis_label="Contig Average Length"), 'right')

        # Plot the gene average lengths
        p.line(x='Sample', y='avg_Length', source=source, color=line_color, line_width=2,y_range_name="average_contig_length", legend_label="Contig Average Length")
        p.circle(x='Sample', y='avg_Length', source=source,size=8,line_color=line_color, fill_color="white", line_width=1.5,y_range_name="average_contig_length")

        #Create a button.
        xbutton = Button(label="Toggle X Axis", button_type="success")

        # Define a JavaScript callback function to toggle the display state of Xaxis labels.
        callback = CustomJS(args=dict(xaxis=p.xaxis[0]), code="""
            if (xaxis.visible) {
                xaxis.visible = false;
            } else {
                xaxis.visible = true;
            }
        """)

        # Bind the JavaScript callback to the button.
        xbutton.js_on_click(callback)

        # Customize layout
        p.legend.location = "top_left"
        p.legend.orientation = "horizontal"
        p.legend.title = ""
        p.xaxis.axis_label_text_font_size = "10pt"
        p.yaxis.axis_label_text_font_size = "10pt"
        p.legend.label_text_font_size = "10pt"
        p.y_range.start = 0
        p.x_range.range_padding = 0.1
        p.xaxis.major_label_orientation = 1
        p.xgrid.grid_line_color = None
        p.legend.click_policy="hide"
        # Create the final layout
        layout = column(p,xbutton)
        layout.sizing_mode = "stretch_width"
        layout.margin = (50, 50, 50, 50)  # (top, right, bottom, left)
        script, div = components(layout)
        cotigNumber_html= div+script

    contigLengthHistogram_html = contigLengthHistogram.to_html(full_html=False, include_plotlyjs=False)
    contigGCHistogram_html = contigGCHistogram.to_html(full_html=False, include_plotlyjs=False)
    return(contigLengthHistogram_html,contigGCHistogram_html,cotigNumber_html)


#Gene set
def uniqueGeneSetNumber(genelenF):
    geneset=pd.read_csv(genelenF,sep="\t",low_memory=False)
    gbins = list(np.arange(200, 2100, 100))
    gbins.append(np.inf)
    glabels = [f"{i}-{i+100}" for i in range(200, 2100, 100)]
    glabels[-1] = ">2K"
    ghist, gedges = np.histogram(geneset['length'], bins=gbins)
    genesetLengthHistogram = px.histogram(x=glabels,y=ghist,color_discrete_sequence=['rgba(50, 171, 96, 0.6)'])
    genesetLengthHistogram.update_layout(
        title_text='Total Gene Number :' + str(len(geneset['length'])), # title of plot
        xaxis_title_text='Genes Length', # xaxis label
        yaxis_title_text='Count', # yaxis label
        bargap=0.2, # gap between bars of adjacent location coordinates
        bargroupgap=0, # gap between bars of the same location coordinates
        height = 500,
        margin=dict(l=50, r=50, b=50, t=50),
    )
    return(genesetLengthHistogram.to_html(full_html=False, include_plotlyjs=False))

#each sample gene stat
def sampleGene(sampleGenestatF):
    gene=pd.read_csv(sampleGenestatF,sep="\t",low_memory=False)

    gene['Sample'] = gene['Sample'].astype('str')

        #Sample  geneNumber      geneAverageLength
    
    if len(list(gene['Sample'])) <= 100:

        geneNumner = go.Figure(
            data=go.Bar(
                y=gene['Sample'],
                x=gene['geneNumber'],
                name="Gene Number",
                width=0.8,
                marker=dict(color="rgba(50, 171, 96, 0.6)"),
                orientation = "h",))
        geneNumner.add_trace(
            go.Scatter(
                y=gene['Sample'],
                x=gene['geneAverageLength'],
                xaxis="x2",
                name="Gene Length",
                marker=dict(color="#109618"),
                orientation = "h",))
        geneNumner.update_layout(
            margin=dict(l=50, r=50, b=50, t=50),
            # legend=dict(orientation="h",yanchor="bottom",y=-0.5,xanchor="right",x=1,title=''),
            xaxis=dict(title=dict(text="Total Gene Number of sample"),side="top",),
            xaxis2=dict(title=dict(text="Gene Average Length"),overlaying="x",tickmode="sync",),
            height = 10*len(list(gene['Sample']))+200,)
        # geneNumner.update_layout(legend=dict(y = -1))
        return(geneNumner.to_html(full_html=False, include_plotlyjs=False))
    else:
        source = ColumnDataSource(gene)
        # Define the categories and data
        samples = gene['Sample'].tolist()
        gene_numbers = gene['geneNumber'].tolist()
        gene_average_lengths = gene['geneAverageLength'].tolist()

        # Create the figure
        fill_color =  "rgba(50, 171, 96, 0.6)"
        line_color = "rgba(50, 171, 96, 1.0)"
        toolslist = "hover,pan,box_zoom,wheel_zoom,undo,redo,reset,save"
        p = figure(x_range=FactorRange(*samples), height=600,
                sizing_mode="stretch_width",
                y_axis_label="Gene Metrics", 
                tools=toolslist,
                tooltips=[('sample', '@Sample'),('geneNumber', '@geneNumber'),('geneAverageLength', '@geneAverageLength')]
                )

        # Plot the gene numbers
        p.vbar(x='Sample', top='geneNumber', width=0.8,
            source=source,
            fill_color=fill_color, fill_alpha=0.8, line_color=line_color, line_width=1.2,
            legend_label="Total Gene Number of Sample")

        # Add the second y-axis for gene average length
        p.extra_y_ranges = {"gene_length": Range1d(start=0, end=max(gene_average_lengths))}
        p.add_layout(LinearAxis(y_range_name="gene_length", axis_label="Gene Average Length"), 'right')

        # Plot the gene average lengths
        p.line(x='Sample', y='geneAverageLength', source=source, color=line_color, line_width=2,y_range_name="gene_length", legend_label="Gene Average Length")
        p.circle(x='Sample', y='geneAverageLength', source=source,size=8,line_color=line_color, fill_color="white", line_width=1.5,y_range_name="gene_length")

        #Create a button.
        xbutton = Button(label="Toggle X Axis", button_type="success")

        # Define a JavaScript callback function to toggle the display state of Xaxis labels.
        callback = CustomJS(args=dict(xaxis=p.xaxis[0]), code="""
            if (xaxis.visible) {
                xaxis.visible = false;
            } else {
                xaxis.visible = true;
            }
        """)

        # Bind the JavaScript callback to the button.
        xbutton.js_on_click(callback)

        # Customize layout
        p.legend.location = "top_left"
        p.legend.orientation = "horizontal"
        p.legend.title = ""
        p.xaxis.axis_label_text_font_size = "10pt"
        p.yaxis.axis_label_text_font_size = "10pt"
        p.legend.label_text_font_size = "10pt"
        p.y_range.start = 0
        p.x_range.range_padding = 0.1
        p.xaxis.major_label_orientation = 1
        p.xgrid.grid_line_color = None
        p.legend.click_policy="hide"
        # Create the final layout
        layout = column(p,xbutton)
        layout.sizing_mode = "stretch_width"
        layout.margin = (50, 50, 50, 50)  # (top, right, bottom, left)
        script, div = components(layout)
        bartxt= div+script
        return(bartxt)


#gene abundance PCA
#gene profile
def top100Profile(geneProfile):
    df=pd.read_csv(geneProfile,sep="\t",low_memory=False,index_col=0)
    gene_occurrence = (df > 0).sum(axis=1)
    filtered_genes = gene_occurrence[gene_occurrence == df.shape[1]].index
    gene_abundance = df.loc[filtered_genes].sum(axis=1)
    top_100_genes = gene_abundance.nlargest(100)
    top100 = df.loc[top_100_genes.index,:]
    return(top100)

#Gene Function
def cogfig(cogF):
    cogLevelCount_sorted=pd.read_csv(cogF,sep="\t",low_memory=False,index_col=0)
    cogLevelCountBar = px.bar(cogLevelCount_sorted, x='loggeneNumber', y='label', color='Cog_CATEGORIE', 
                color_discrete_map={'POORLY CHARACTERIZED': '#fc8d62', 'METABOLISM': '#8da0cb', 'INFORMATION STORAGE AND PROCESSING': '#e78ac3',"CELLULAR PROCESSES AND SIGNALING":"#a6d854"},
                orientation='h')
    cogLevelCountBar.update_layout(
        legend=dict(orientation="h",yanchor="bottom",y=-0.2,xanchor="right",x=1,title='Level1 Category'),
        xaxis=dict(title='log2(gene numbers)'),
        yaxis=dict(title='COG|DataBase Category'),
        height = 30*len(list(cogLevelCount_sorted.index))+50,
        margin=dict(l=50, r=50, b=50, t=50),
        )
    return(cogLevelCountBar.to_html(full_html=False, include_plotlyjs=False))

def binQCFig(binQCF):
    binQC = pd.read_csv(binQCF,sep="\t",low_memory=False)
    binQC["QC"] = binQC[["Contamination", "Completeness"]].apply(lambda x: "HQ" if x["Completeness"] > 90 and x["Contamination"] < 10 else "LQ", axis=1)
    binQS = px.scatter(binQC, x="Contamination", y="Completeness",  color="QC",
        size_max=30,
        hover_name="Name",
        color_discrete_map={"HQ": "rgba(50, 171, 96, 0.6)", "LQ": "grey"},
        category_orders={"QC": ["HQ", "LQ"]},
        # color_discrete_sequence=["rgba(50, 171, 96, 0.6)","grey"],
    )
    binQS.add_vline(x=10, line_width=2, line_dash="dash", line_color="grey")
    binQS.add_hline(y=90, line_width=2, line_dash="dash", line_color="grey")
    binQS.update_layout(
        xaxis=dict(title='Contamination'),
        yaxis=dict(title='Completeness'),
        margin=dict(l=400, r=400, b=50, t=50), height = 600,
    )
    return(binQS.to_html(full_html=False, include_plotlyjs=False))

def binSizeFig(binsF):
    bins = pd.read_csv(binsF,sep="\t",low_memory=False)
    bins["bublesize"] = bins[["contigNumber","GenomeSize"]].apply(lambda x:x["GenomeSize"]/x["contigNumber"],axis=1)
    binsContigNumberSize = px.scatter(bins, x="contigNumber", y="GenomeSize",size = "bublesize", size_max=40,
    hover_name="BinID",
    color_discrete_sequence=['rgba(50, 171, 96, 0.6)'])
    binsContigNumberSize.update_layout(
        xaxis=dict(title='contig Number'),
        yaxis=dict(title='Total Contig Length'),
        margin=dict(l=400, r=400, b=50, t=50), height = 600,
    )
    return(binsContigNumberSize.to_html(full_html=False, include_plotlyjs=False))

def pasteKraken2Species(speciesF):
    allLevelPD = pd.read_csv(speciesF,sep="\t",low_memory=False, index_col = 0)
    filtered_df = allLevelPD[allLevelPD.index.str.contains('s_')]
    filtered_df.index = filtered_df.index.str.extract('s_(.*)')[0]
    return(filtered_df)


#taxonomay
def taxonomayBarFig(speciesProfile,ftype):
    if ftype =="pd" : 
        df = speciesProfile
    else:
        df =pd.read_csv(speciesProfile,sep="\t",index_col=0,low_memory=False)
    df_relative = df.div(df.sum()) * 100
    bin_sums = df_relative.sum(axis=1)
    top20_bins = bin_sums.nlargest(20).index
    df_top20 = df_relative.loc[top20_bins]
    df_others = df_relative.loc[~df_relative.index.isin(top20_bins)].sum().to_frame().T
    df_top20.loc['Others'] = df_others.values[0]
    df_sorted = df_top20.T.sort_values(by = [top20_bins[0],top20_bins[1],top20_bins[2],top20_bins[3]], ascending = [False,True, False,True], na_position = 'last')
    # df_sorted['item'] = df_sorted.index
    # df_melt = df_sorted.melt(id_vars=['item'],var_name='sample', value_name='abundance')
    
    speecies_list = df_sorted.columns.to_list()
    sample_list = df_sorted.index.to_list()
    df_sorted['sample'] = df_sorted.index
    data = {key: df_sorted[key].tolist() for key in df_sorted.columns}

    
    # colorlist = list(px.colors.qualitative.Pastel[:-1])+list(px.colors.qualitative.Bold[:-1]) +list(px.colors.qualitative.Set2)
    colorlist = list(px.colors.qualitative.Pastel[:-1])+list(px.colors.qualitative.Bold[:-1]) +list(px.colors.qualitative.Set2)
    colorlist = colorlist[:20]
    colorlist.append('#d9d9d9')

    if len(sample_list) <= 60:

        #go.bar
        taxonomay_bar = go.Figure()
        for n,species in enumerate(speecies_list):
            hover_text = [f'{sample}|{species}:{abundance}' for sample,abundance in zip(sample_list, data[species])]
            taxonomay_bar.add_trace(go.Bar(
                y=sample_list,
                x=data[species],
                orientation='h',
                name= species,
                hovertext=hover_text,  # 指定悬停标签
                marker_color=colorlist[n],
                width=0.95
            )
        )
        taxonomay_bar.update_layout(
            # title='Microbiome Abundance',
            yaxis=dict(title='Samples'),
            xaxis=dict(title='Relative Abundance',side="top"),
            barmode='stack',
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-1,
                xanchor="right",
                x=1,
                title='Species'
            ),
            height = 10*len(sample_list)+300,
            margin=dict(l=50, r=50, t=50, b=50),
        )

        ##px.bar
        # taxonomay_bar = px.bar(df_melt, y="sample", x="abundance", color="item",
        # 				width=0.5,
        # 				color_discrete_sequence= colorlist,
        # 				category_orders={"sample": df_sorted.columns.to_list(),
        # 							"item": df_sorted.index.to_list(),})
        # taxonomay_bar.update_layout(legend=dict(
        # 	orientation="h",
        # 	yanchor="bottom",
        # 	y=-2.0,
        # 	xanchor="right",
        # 	x=1,
        # 	title='Species'),
        # 	height = 50*len(list(df_sorted.columns))+200,
        # 	title='Microbiome Abundance',
        # 	yaxis=dict(title='Samples'),
        # 	xaxis=dict(title='Relative Abundance'),
        # 	barmode='stack' 
        # 	)
        return(taxonomay_bar.to_html(full_html=False, include_plotlyjs=False))
    else:
        toolslist = "hover,pan,box_zoom,wheel_zoom,undo,redo,reset,save"

        p = figure(
            x_range=sample_list,
            y_range=(0, 100),
            title="Microbiome Abundance",
            x_axis_label="Samples",
            y_axis_label="Relative Abundance %",
            height=600,
            sizing_mode="stretch_width",
            tools=toolslist, 
            tooltips=[('sample', '@sample'),('Species', '$name'),('Abundance', '@$name')],
        )

        p.vbar_stack(speecies_list, x='sample', width=0.9, color=colorlist, source=data)

        p.y_range.start = 0
        p.x_range.range_padding = 0.1
        p.xgrid.grid_line_color = None
        p.axis.minor_tick_line_color = None
        p.xaxis.major_label_orientation = 1.2
        p.outline_line_color = None
        p.toolbar.autohide = True

        cmap = CategoricalColorMapper(palette=colorlist[:21], factors=speecies_list)
        color_bar = ColorBar(color_mapper=cmap, 
                    major_label_text_font_size="10px",
                    orientation = "vertical", 
                    height = 400, 
                    width=20,
                    location="right",
                    major_label_text_line_height=2) 
        p.add_layout(color_bar, 'right')

        #Create a button.
        xbutton = Button(label="Toggle X Axis", button_type="success")

        # Define a JavaScript callback function to toggle the display state of Xaxis labels.
        callback = CustomJS(args=dict(xaxis=p.xaxis[0]), code="""
            if (xaxis.visible) {
                xaxis.visible = false;
            } else {
                xaxis.visible = true;
            }
        """)

        # Bind the JavaScript callback to the button.
        xbutton.js_on_click(callback)

        #hide legend
        # Create a button.
        legend_button = Button(label="Toggle Legend label", button_type="success")

        # Define a JavaScript callback function to toggle the display state of the ColorBar.
        legend_callback = CustomJS(args=dict(color_bar=color_bar), code="""
            if (color_bar.visible) {
                color_bar.visible = false;
            } else {
                color_bar.visible = true;
            }
        """)

        # Bind the JavaScript callback function to the button.
        legend_button.js_on_click(legend_callback)

        # Combine the graph and the button together.
        layout = column(p, row(xbutton,legend_button))
        layout.margin = (50, 50, 50, 50)  # (top, right, bottom, left)
        layout.sizing_mode = "stretch_width"
        script, div = components(layout)
        bartxt= div+script
        return(bartxt)




#taxonomy sunburst
def taxonomySunFig(sunburstFile):
    hierarchy = pd.read_csv(sunburstFile,sep="\t",low_memory=False)

    hierarchyPD = dict(
        character=hierarchy["child"].to_list(),
        parent=hierarchy["parent"].to_list(),
        abundance=hierarchy["abundance"].to_list(),
        grandsonNumber=hierarchy["grandsonNumber"].to_list()
        )
    taxonomy_sunburst = px.sunburst(hierarchyPD, names='character',parents='parent',values='grandsonNumber',
                    color='abundance', 
                    color_continuous_scale='RdBu',
                    color_continuous_midpoint=np.average(hierarchy['abundance'], weights=hierarchy['grandsonNumber'])
                    )
    taxonomy_sunburst.update_layout(height = 800,margin=dict(l=50, r=50, b=50, t=50),)
    return(taxonomy_sunburst.to_html(full_html=False, include_plotlyjs=False))


#bin gtdb 
#taxonomy sunburst
def treemapFig(treemapFile):
    gtdbPd = pd.read_csv(treemapFile,sep="\t",low_memory=False)
    taxonomy_treemap = px.treemap(gtdbPd, names='child',parents='parent',values='value',
                    color='value', 
                    color_continuous_scale='RdBu',
                    color_continuous_midpoint=np.average(gtdbPd['value'], weights=gtdbPd['value'])
                    )
    taxonomy_treemap.update_layout(height = 800,margin=dict(l=50, r=50, b=50, t=50),)
    return(taxonomy_treemap.to_html(full_html=False, include_plotlyjs=False))


def binFunctionAnimation(binFunctionDatabasePD,batabaseName):
    keggPD = binFunctionDatabasePD
    ymin = min(keggPD["count"])- (max(keggPD["count"] - min(keggPD["count"])))*0.2
    ymax = max(keggPD["count"]) + (max(keggPD["count"] - min(keggPD["count"])))*0.2
    eachBinFunction_KEGG_A = px.scatter(keggPD, x="class", y="count", animation_frame="binID", animation_group="class", size="count", color="class", hover_name="class",  size_max=45, range_y=[ymin,ymax],)
    eachBinFunction_KEGG_A.update_xaxes(showticklabels=False)
    eachBinFunction_KEGG_A.update_layout(legend=dict( orientation="h",yanchor="bottom", y=-0.6,xanchor="right",x=1, ),
        height = 700,
        #title="Each bin's "+ batabaseName+" Function statistics",
        yaxis=dict(title='Count'),
        xaxis=dict(title=''),
        margin=dict(l=50, r=50, b=50, t=50),
        )
    return(eachBinFunction_KEGG_A.to_html(full_html=False, include_plotlyjs=False))





# =============================== #
# ====== add bokeh module  ====== #
# =============================== #

#@1 circle link a table 
def circle_taable_link(profile,columnsNameList,newTitleDir,xaxisID,yaxisID,xaxislabel,yaxislabel):
    data=pd.read_csv(profile,sep="\t",low_memory=False)
    data.columns=columnsNameList
    source = ColumnDataSource(data)

    table_height = 400
    if len(data.index) <= 10:
        table_height = 25*len(data.index) + 50


    columns = []
    tooltips = []
    for onecolumn in columnsNameList:
        columns.append(TableColumn(field=onecolumn, title=newTitleDir[onecolumn]))
        tooltips.append((f"{newTitleDir[onecolumn]}", f"@{onecolumn}"))
            
    data_table = DataTable(source=source, columns=columns,index_position=-1, index_header="row index", index_width=20,sizing_mode="stretch_width",height=table_height)

    p = figure(tools="pan,wheel_zoom,xbox_select,reset", active_drag="xbox_select",sizing_mode="stretch_width",height=400)

    showPoint = p.circle(x=xaxisID, y=yaxisID, fill_color="rgba(50, 171, 96, 1)", size=10, alpha=0.6, source=source)
    
    p.background_fill_color = "rgba(229,236,246,1)"
    p.grid.grid_line_color = "rgba(255, 255, 255, 1)"
    # p.background_fill_alpha = 0.9
    p.xaxis.axis_label = xaxislabel
    p.yaxis.axis_label = yaxislabel

    
    showPoint_hover_tool = HoverTool(renderers=[showPoint],tooltips=tooltips)

    p.add_tools(showPoint_hover_tool)

    layout = column(p,data_table)
    layout.margin = (50, 50, 50, 50)  # (top, right, bottom, left)
    layout.sizing_mode = "stretch_width"

    # # Convert the Bokeh plot object to an HTML string.
    script, div = components(layout)
    linktxt= div+script
    return(linktxt)
#1=================================================

#@2 A spreadsheet that can filter multiple columns.


#2=================================================

def createCustomJSCode(columnsNameList, filter_column):
    codetxt = 'const data = source.data;\n'
    colTxt = ''
    pushTxt = ''
    for i in columnsNameList:
        colTxt += f'{i}:[] ,'
        pushTxt += f'filtered_data.{i}.push(data.{i}[i]);\n'

    codetxt += 'const filtered_data = { '+colTxt.strip(",")+' };\n'
    codetxt +='const filter_value = filter_input.value.toLowerCase();\n'
    codetxt +='''for (let i = 0; i < data.'''+filter_column+'''.length; i++) {
    if (data.'''+filter_column+'''[i].toLowerCase().includes(filter_value)) {
        '''+pushTxt+'''
        }
    }
    filtered_source.data = filtered_data;'''
    return(codetxt)


def multiConditionCustomJSCode(columnsNameList, multi_filter_columns):
    codetxt = 'const data = source.data;\n'
    colTxt = ''
    pushTxt = ''
    for i in columnsNameList:
        colTxt += f'{i}:[] ,'
        pushTxt += f'filtered_data.{i}.push(data.{i}[i]);\n'

    #ceate a new empty dataframe
    codetxt += 'const filtered_data = { '+colTxt.strip(",")+' };\n'
    
    # define all filter item
    filter_condition = []
    for n,columnName in enumerate(multi_filter_columns):
        codetxt +=f'const filter_item_{columnName} = multi_input[{n}].value.toLowerCase();\n'
        filter_condition.append(f'data.{columnName}[i].toLowerCase().includes(filter_item_{columnName})')

    filter_condition_txt = ' && '.join(filter_condition)

    codetxt +='''for (let i = 0; i < data.'''+multi_filter_columns[0]+'''.length; i++) {
    if ( ''' + filter_condition_txt +''' ) {
        '''+pushTxt+'''
        }
    }
    filtered_source.data = filtered_data;'''
    return(codetxt)





def tableFilter(profile, columnsNameList, filter_columns, newTitleDir):
    """
    Create a DataTable with filtering capability

    Parameters:
        profile (str): Path to the data file.
        columnsNameList (list): List of column names to be displayed.
        filter_columns (list): List of column names used for filtering.
        newTitleDir (dict): Dictionary mapping original column names to new titles.

    Returns:
        str: HTML code containing the DataTable and filter inputs.
    """

    # Read data and set column names
    data = pd.read_csv(profile, sep="\t", low_memory=False)

    if len(data.columns) > len(columnsNameList):
        new_data = data[columnsNameList]
        data = new_data
    else:
        data.columns = columnsNameList

    table_height = 600
    if len(data.index) <= 30:
        table_height = 25*len(data.index) + 100


    # Create original data source
    original_source = ColumnDataSource(data)

    # Create filtered data source
    filtered_source = ColumnDataSource(data)

    columns = [TableColumn(field=col, title=newTitleDir.get(col, col)) for col in columnsNameList]

    data_table = DataTable(
        source=filtered_source,
        columns=columns,
        index_position=-1,
        index_header="row index",
        index_width=40,
        height=table_height,
        editable = True ,
        sizing_mode="stretch_width"
    )

    # Create TextInputs for filtering
    filter_inputs = []
    for col in filter_columns:
        filter_input = TextInput(title=f"Filter by {col}:", value="")
        filter_inputs.append(filter_input)

    # Generate custom JS code with updated filtering logic
    custom_js_codes = [createCustomJSCode(columnsNameList, col) for col in filter_columns]


    # Create and bind custom JS callbacks for each filter input
    for i, (filter_input, custom_js_code) in enumerate(zip(filter_inputs, custom_js_codes)):
        js_data = original_source if i == 0 else filtered_source
        filter_callback = CustomJS(
            args=dict(
                source=js_data,
                filtered_source=filtered_source,
                filter_input=filter_input,
            ),
            code=custom_js_code
        )
        filter_input.js_on_change('value', filter_callback)

    #create search button
    multi_js = multiConditionCustomJSCode(columnsNameList, filter_columns)
    search_button = Button(label="Search", button_type="success")

    search_callback = CustomJS(
        args=dict(
            source=original_source,
            filtered_source=filtered_source,
            multi_input = filter_inputs,
            ),
        code=multi_js
        )

    search_button.js_on_click(search_callback)




    # Create a button to clear all filters primary
    clear_button = Button(label="Clear All Filters", button_type="primary")

    clear_callback = CustomJS(
        args=dict(
            source=original_source,
            filtered_source=filtered_source,
            filter_inputs=filter_inputs,
        ),
        code="""
        // Clear the values of all filter inputs
        for (let i = 0; i < filter_inputs.length; i++) {
            filter_inputs[i].value = "";
        }
        // Restore the original data source
        filtered_source.data = source.data;
        """
    )
    clear_button.js_on_click(clear_callback)

    # Combine all filter inputs and the clear button into a row layout
    filter_inputs_row = row(*filter_inputs )

    button_row = row(search_button,clear_button)
    # Combine the filter inputs row and the data table into a column layout
    layout = column(filter_inputs_row,button_row,data_table, sizing_mode="stretch_width")
    layout.margin = (50, 50, 50, 50)

    # Convert the Bokeh plot object to an HTML string
    script, div = components(layout)
    tabletxt = div + script
    return(tabletxt)


#@3 heatmap.

def bokehHeatmapFig(dataFrame, titletxt, removeIndexList=[]):
    # input data
    rawdata = pd.read_csv(dataFrame,sep="\t",index_col=0,low_memory=False)

    if removeIndexList != []:
        rawdata = rawdata[~rawdata.index.isin(["UNMAPPED", "UNINTEGRATED"])]

    # check row and col number more than 2
    # num_rows, num_cols = rawdata.shape
    # if num_rows < 2 or num_cols < 2:
    # 	return('''<p style="font-size: 12pt; color: red; font-style: italic;" >Oops! Number of Smaple less than 2, nor MetaCyc pathway less than 2.  Socouldn't show as Heatmap.<br></p>''')
    
    data_frame_check = check_data_frame_more_than_three(rawdata)
    if data_frame_check != 'Pass':
        return(data_frame_check)

    else:
        # Find the minimum value greater than zero.
        min_positive = np.min(rawdata[rawdata > 0].min())
        # Set the offset.
        offset = min_positive * 0.01
        # Take the logarithm of the data.
        data = np.log(rawdata + offset)
        # Perform hierarchical clustering on rows and columns.
        row_clusters = AgglomerativeClustering().fit_predict(data.values)
        col_clusters = AgglomerativeClustering().fit_predict(data.values.T)
        # Reorder the matrix data based on the clustering results.
        sorted_data = data.iloc[np.argsort(row_clusters)]
        sorted_data = sorted_data.iloc[:, np.argsort(col_clusters)]
        order_pathway = list(sorted_data.index)
        order_pathway_simple = [element.split(":")[0] for element in order_pathway]
        order_sample = list(sorted_data.columns)
        # Convert the DataFrame to long format.
        long_data = pd.melt(sorted_data.reset_index(), id_vars='index', var_name='sample', value_name='logAbundance')
        long_data['index'] = long_data['index'].astype(str)
        long_data['sample'] = long_data['sample'].astype(str)
        long_data['simpleLabel'] = long_data['index'].str.split(':').str.get(0)
        orgLongData =  pd.melt(rawdata.reset_index(), id_vars='index', var_name='sample', value_name='Abundance')
        merged_data = pd.merge(long_data, orgLongData, on=['index', 'sample'])

        # this is the colormap from the original NYTimes plot
        colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]

        #TOOLS = "hover,save,pan,box_zoom,undo,redo,reset,save,wheel_zoom"
        TOOLS = "hover,pan,box_zoom,wheel_zoom,undo,redo,reset,save"

        p = figure(title=f"{titletxt}",
            x_range=order_sample, 
            y_range=order_pathway_simple,
            x_axis_location="below", 
            y_axis_location="right", 
            height=800,
            sizing_mode="stretch_width",
            tools=TOOLS, 
            toolbar_location='above',
            tooltips=[('Pathway', '@index'),('sample', '@sample'), ('abundance', '@Abundance')])

        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = "10px"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = pi / 3

        r = p.rect(x="sample", y="simpleLabel", width=1, height=1, source=merged_data,
                    fill_color=linear_cmap("logAbundance", colors, low=merged_data.logAbundance.min(), high=merged_data.logAbundance.max()),
                    line_color=None)

        p.add_layout(r.construct_color_bar(
            major_label_text_font_size="7px",
            ticker=BasicTicker(desired_num_ticks=len(colors)),
            formatter=PrintfTickFormatter(format="%6.3f"),
            label_standoff=6,
            border_line_color=None,
            padding=5,
        ), 'left')

        # Set autohide to true to only show the toolbar when mouse is over plot
        p.toolbar.autohide = True

        #Create a button.
        ybutton = Button(label="Toggle Y Axis Label", button_type="success")

        # Define a JavaScript callback function to toggle the display state of Y-axis labels.
        ycallback = CustomJS(args=dict(yaxis=p.yaxis[0]), code="""
            if (yaxis.visible) {
                yaxis.visible = false;
            } else {
                yaxis.visible = true;
            }
        """)

        # Bind the JavaScript callback to the button.
        ybutton.js_on_click(ycallback)

        xbutton = Button(label="Toggle X Axis Label", button_type="success")

        # Define a JavaScript callback function to toggle the display state of Y-axis labels.
        xcallback = CustomJS(args=dict(xaxis=p.xaxis[0]), code="""
            if (xaxis.visible) {
                xaxis.visible = false;
            } else {
                xaxis.visible = true;
            }
        """)

        # Bind the JavaScript callback to the button.
        xbutton.js_on_click(xcallback)

        # Combine the graph and the button together.
        layout = column(p, row(xbutton,ybutton))
        layout.sizing_mode = "stretch_width"
        layout.margin = (50, 50, 50, 50)

        # Convert the Bokeh plot object to an HTML string.
        script, div = components(layout)
        heatmaptxt= div+script
        return(heatmaptxt)
#3=================================================

# =============================== #
# ====== add bokeh module  ====== #
# =============================== #





def plotAllResult(reportFile,templateHTML):
    runTime =  time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) 
    timeTable = go.Figure(data=[go.Table(header=dict(values=['Run times'],
                line_color='white',
                fill_color='white',
                align='center',
                height=20,),
            cells=dict(values=[runTime],
                line_color='white',
                fill_color='white',
                align='center',
                height=20,
                ))])
    timeTable.update_layout(margin=dict(l=1, r=1, b=1, t=1),height=70)
    timeTable_html = timeTable.to_html(full_html=False, include_plotlyjs='cdn')
    
    plotData = [
        {'type': 'home', 'content': 'Nice! Workflow execution completed successfully!'},
        {'type': 'plot', 'content': timeTable_html},]


    plotFunDir = {
        #"readstat":  lambda a : figTable(a,["ID","raw_reads","raw_bases","clean_reads","clean_bases","effective_rate(%)","GC(%)","lq_reads","host_reads"]),
        "readstat":  lambda a : circle_taable_link(a,read_first_line(a)[0],read_first_line(a)[1],"effective_rate","GC","Effective Rate %","GC %"),
        "mpaspecies" : lambda a : taxonomySunFig(a),
        "mpaspeciesT" : lambda a : taxonomayBarFig(a,"file"),
        "mpaspeciesTPCA" : lambda a : plotlyPCA(a,"#24B064","MetaPhlAn Abundance PCA","px","F"),
        "krakenspecies" : lambda a : taxonomySunFig(a),
        "krakenspeciesT" : lambda a : taxonomayBarFig(pasteKraken2Species(a),"pd"),
        "krakenspeciesTPCA" : lambda a : plotlyPCA(a,"#24B064","Kraken2 Abundance PCA","px","F"),
        "genesetLenStat" : lambda a : uniqueGeneSetNumber(a),
        "genesetSampleStat": lambda a : sampleGene(a),
        "genesetAbundance" : lambda a : plotlyPCA(top100Profile(a),"#24B064","Gene PCA","px","df"),
        "genesetCOG": lambda a : cogfig(a),
        "genesetGO" : lambda a : summaryDB2Plot(a,"GO",280),
        "genesetKEGG" : lambda a : summaryDB2Plot(a,"KEGG",350),
        "genesetCAZY" : lambda a : summaryDB2Plot(a,"CAzY",310),
        "antismash" : lambda a : plotlyPCA(a,"#24B064","Geneset antisMASH function abundance PCA","px","F"),
        "VFDB" : lambda a : plotlyPCA(a,"#24B064","Geneset VFDB function abundance PCA","px","F"),
        "CARD" : lambda a : plotlyPCA(a,"#24B064","Geneset CARD function abundance PCA","px","F"),
        "customPR" : lambda a : plotlyPCA(a,"#24B064","Geneset custom protein database function abundance PCA","px","F"),
        "customNT" : lambda a : plotlyPCA(a,"#24B064","Geneset gene custom nucleic acid database function abundance PCA","px","F"),
        "bintable" : lambda a : tableFilter(a, read_first_line(a)[0], ['BinID','originalID','GTDB_taxonomy','NCBI_taxonomy'],read_first_line(a)[1]),
        "rawbinQS" : lambda a : binQCFig(a),
        "binInfo" : lambda a : binSizeFig(a),
        "gtdb" : lambda a : treemapFig(a),
        "binabundance" : lambda a : multiPCA(a,"#24B064"),
        "binDepth" : lambda a : plotlyPCA(a,"#24B064","Bin's Depth Abundance PCA","px","F"),
        "binRPKM": lambda a : plotlyPCA(a,"#24B064","Bin's RPKM Abundance PCA","px","F"),
        "binTPM" : lambda a : plotlyPCA(a,"#24B064","Bin's TPM Abundance PCA","px","F"),
        # "bin_res_stat" : lambda a : tableFilter(a, ['BinID','Orignal_GTDB_taxonomy','Orignal_Completeness','ReAssembly_Completeness','ReAssembly_Deepurify_Completeness','Orignal_Contamination','ReAssembly_Contamination','ReAssembly_Deepurify_Contamination'],['BinID','Orignal_GTDB_taxonomy'],read_first_line(a)[1]),
        "bin_res_stat" : lambda a : tableFilter(a,read_first_line(a)[0], ['BinID','Completeness','Contamination','Level','Staus','Stage'],read_first_line(a)[1]),
        "metacyc" : lambda a : bokehHeatmapFig(a,"HUMAnN MetaCyc log realtive abundance",["UNMAPPED", "UNINTEGRATED"]),

    }

    H1TitleDir = {
        "readstat" : "Reads Quality",
        "mpaspeciesT" : "Taxonomy Classification | MetaPhlAn",
        "metacyc" : "Metabolic Analysis Network | HUMAnN",
        "krakenspeciesT" : "Taxonomy Classification | Kraken2",
        "genesetLenStat" : "Gene Set",
        "bintable" : "Binning",
        # "rawbinQS" : "Binning",
    }

    plotContent = {}

    with open(reportFile,'r') as reportF:
        for r in reportF:
            rl = r.strip().split(' ')
            if rl[0] in plotFunDir:
                if rl[0] == 'binabundance':
                    resultFilelineDir[rl[0]] = get_file_line_count(rl[1].strip().split(",")[0])
                else:
                    resultFilelineDir[rl[0]] = get_file_line_count(rl[1])
                    if int(get_file_line_count(rl[1])) < 2:
                        continue

            if rl[0] == "contigstat":
                contigStatList = contigFigure(rl[1])
                oneplotDir = [
                    {'type': 'H1', 'content': "Assembly Contigs"},
                    {'type': 'text', 'content': resultPlottxtDir[rl[0]+"1"]},
                    {'type': 'plot', 'content': contigStatList[0]},
                    {'type': 'text', 'content': resultPlottxtDir[rl[0]+"2"]},
                    {'type': 'plot', 'content': contigStatList[1]},
                    {'type': 'text', 'content': resultPlottxtDir[rl[0]+"3"]},
                    {'type': 'plot', 'content': contigStatList[2]},
                ]
                plotContent[rl[0]] = oneplotDir

            if rl[0] == "eachBinFunction":
                eachBinfp = pd.read_csv(rl[1],sep="\t",low_memory=False)
                cazyPD = eachBinfp[eachBinfp['database'] == 'CAzY']
                cogPD = eachBinfp[eachBinfp['database'] == 'COG']
                goPD = eachBinfp[eachBinfp['database'] == 'GO']
                keggPD = eachBinfp[eachBinfp['database'] == 'KEGG']
                oneplotDir = [
                    {'type': 'text', 'content': resultPlottxtDir[rl[0]+"KEGG"]},
                    {'type': 'plot', 'content': binFunctionAnimation(keggPD,'KEGG')},
                    {'type': 'text', 'content': resultPlottxtDir[rl[0]+"CAzY"]},
                    {'type': 'plot', 'content': binFunctionAnimation(cazyPD,'CAzY')},
                    {'type': 'text', 'content': resultPlottxtDir[rl[0]+"COG"]},
                    {'type': 'plot', 'content': binFunctionAnimation(cogPD,'COG')},
                    {'type': 'text', 'content': resultPlottxtDir[rl[0]+"GO"]},
                    {'type': 'plot', 'content': binFunctionAnimation(goPD,'GO')},
                ]
                plotContent[rl[0]] = oneplotDir

            else:
                if rl[0] in plotFunDir :
                    onehtml = plotFunDir[rl[0]](rl[1])
                    if rl[0] in H1TitleDir:
                        H1 = H1TitleDir[rl[0]]
                        oneplotDir = [{'type': 'H1', 'content': H1}]
                    else:
                        oneplotDir = []
                    oneplotDir.append({'type': 'text', 'content': resultPlottxtDir[rl[0]]})
                    oneplotDir.append({'type': 'plot', 'content': onehtml})
                    plotContent[rl[0]] = oneplotDir


    for figindex in  resultPlotOrderKeyList:
        lebel = resultPlotOrderDir[figindex]
        if lebel in plotContent:
            for h in plotContent[lebel]:
                plotData.append(h)
            
    reference = resultPlottxtDir['reference']

    with open(imgPath+"/nf-core-logo-square.png", "rb") as img_file:
        logoB64 = base64.b64encode(img_file.read()).decode()

    with open(imgPath+"/MetaFlowX.png", "rb") as img_file2:
        pipelinePNG64 = base64.b64encode(img_file2.read()).decode()
    
    template_dir = os.path.dirname(templateHTML)
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template(os.path.basename(templateHTML))
    
    html = template.render(logo=logoB64,
                            pipelinePNG = pipelinePNG64,
                            plotdata = plotData,
                            ReferencesTxt = reference,
                            username = "Hi!")
    nowTime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    errList = checkResultCompletely(resultFilelineDir)
    if len(errList) == 0:
        with open(outpath+'/'+str(filePrefix)+'_Report_'+str(nowTime)+'.html', 'w') as f:
            f.write(html)
        print("Good Luck ~ Report successfully !\n "+nowTime)

    else:
        with open("Report_Reuslt_Failed.txt",'w') as errF:
            errList.append("\n")
            errF.write("\tFailed\n".join(errList))
        print("Oops~ Something is wrong, please check your result report list !")
        sys.exit()



def main():
    if len(sys.argv) == 1:
        sys.argv.append('-h')
        ParsReceiver()
        sys.exit()
    global pars,outpath,shellPath,resultPlotOrderDir,resultPlottxtDir,resultPlotOrderKeyList,imgPath,resultPlottxtKeyList,filePrefix,resultFilelineDir
    pars = ParsReceiver()
    makedir(pars['outPath'])
    outpath = pars['outPath']
    filePrefix = pars['prefix']
    imgPath = pars['img']
    shellPath = os.path.split(os.path.realpath(__file__))[0]
    resultFilelineDir = {}
    resultPlotOrderDir,resultPlotOrderKeyList = {},[]
    resultPlottxtDir,resultPlottxtKeyList = {},[]
    file2Dir(pars['plotOrder'],resultPlotOrderDir,1,resultPlotOrderKeyList)
    file2Dir(pars['txt'],resultPlottxtDir,1,resultPlottxtKeyList)
    plotAllResult(pars['result'],pars['templatehtml'])


if __name__ == '__main__':
    main()
