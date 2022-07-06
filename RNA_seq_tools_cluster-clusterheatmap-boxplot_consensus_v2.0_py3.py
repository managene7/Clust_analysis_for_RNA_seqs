# -*- coding: utf-8 -*-
"""
Created on Sat May 16 00:56:30 2020

@author: minkj
"""


#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-out':"",'-dot_size':"2.5"}
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print("""
_________________________________________________________________________________

Usage;

This code is for analysis of RNA-seq and data visualization
-def    cl: generate individual files for cluster in teh clust result.
            -tsv: file name of cluster info by clust analysis
            -dge: file name of digital gene expression
        hmamp: generate clustered heat map
            -consensus  consensus of names for multiple running
        hmamp_fc: generate clustered heat map by fold change
            -consensus  consensus of names for multiple running
        bp: generate box plot
            -consensus  consensus of names for multiple running
            -dot_size   size of dots in the plot (default is 2.5)
        vp: generate violin plot
            -consensus  consensus of names for multiple running
            -dot_size   size of dots in the plot (default is 2.5)
_________________________________________________________________________________
""")
                quit()


        
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt

import os
import csv 
dot_size=float(option_dict['-dot_size'])
#__________convert clust result into individual cluster files_____________________
if option_dict['-def']=='cl':
    tsv_name=option_dict['-tsv']
    tsv=open(tsv_name,'r')
    data_name=option_dict['-dge']
    cluster_dic={}
    tsv_lines=tsv.readlines()
    init=1
    for row in tsv_lines:
        row=row.strip()
        row=row.split('\t')
        if init==1:
            cl_names=row
            for cl in row:
                cluster_dic[cl]=[]
            init=2
        elif init==2:
            init=3
            pass
        else:
            for i in range(len(cl_names)):
                #print row
                try:
                    if row[i]!="":
                        cluster_dic[cl_names[i]].append(row[i])
                except:
                    pass



    data_file=open(data_name,'r')
    line=data_file.readline().split()
    new_first_line=(',').join(line)
    
    data_dic={}
    while 1:
        line=data_file.readline().strip()
        if line=="":
            break
        line_list=line.split()
        new_line=(",").join(line_list)
        data_dic[line_list[0]]=new_line
    
    for cl, g_names in list(cluster_dic.items()):
        cl_list=cl.split()
        cl_concatenated="_".join(cl_list)
        out_file=open(cl_concatenated+".csv",'w')
        out_file.write(new_first_line+'\n')
        for g in g_names:
            #print g
            if g in data_dic:
                out_file.write(data_dic[g]+'\n')
            else:
                print(g)
    print("Completed!!")
            

#___________draw boxplot________________________________    
elif option_dict['-def']=='bp':
    
    name=option_dict['-consensus']
    file_list_row=os.listdir('.')
    file_list=[]
    for f in file_list_row:
        if name in f:
            file_list.append(f)
    
    for name in file_list:
    
        df = pd.read_csv(name)
        df = df.set_index(df.columns[0])
        df.head()
        bplot=sns.boxplot(data=df, width=0.7)
        bplot=sns.swarmplot(data=df, color="gray", size=dot_size)
        bplot.axes.set_title(name[:-4], fontsize=16)
        plt.savefig(name[:-4]+'_box_plot.jpg',dpi=150, figsize=(4,4))
        plt.close()


elif option_dict['-def']=='vp':
    
    name=option_dict['-consensus']
    file_list_row=os.listdir('.')
    file_list=[]
    for f in file_list_row:
        if name in f:
            file_list.append(f)
    
    for name in file_list:
        
        df = pd.read_csv(name)
        df = df.set_index(df.columns[0])
        df.head()
        bplot=sns.violinplot(data=df, width=0.7, inner=None)
        bplot=sns.swarmplot(data=df, color="white", size=dot_size, edgecolor="black")
        bplot.axes.set_title(name[:-4], fontsize=16)
        plt.savefig(name[:-4]+'_violin_plot.jpg',dpi=150, figsize=(4,4))
        plt.close()



#_____________draw clusterd heat map______________________        
elif option_dict['-def']=='hmamp' or option_dict['-def']=='hmamp_fc':

    name=option_dict['-consensus']
    file_list_row=os.listdir('.')
    file_list=[]
    for f in file_list_row:
        if name in f:
            file_list.append(f)
    
    for name in file_list:
        
        df = pd.read_csv(name)
        df = df.set_index(df.columns[0])
        df.head()
        
        if option_dict['-def']=='hmamp':
            sns.clustermap(df, col_cluster=False,cmap=plt.cm.gist_heat_r, figsize=(4,8), method='average')
        elif option_dict['-def']=='hmamp_fc':
            sns.clustermap(df, col_cluster=False,cmap=plt.cm.seismic,figsize=(4,8), method='average')
        plt.savefig(name[:-4]+'_heat_map.jpg', dpi=150, figsize=(4,4))
        plt.close()
