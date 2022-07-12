# -*- coding: utf-8 -*-
"""
Created on Sat May 16 00:56:30 2020

@author: minkju3003@gmail.com (Minkyu Park)
"""


#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-out':"",'-dot_size':"1.0",'-rep':"",'-log2':'no','-fig_w':"10", "-fig_h":"20",'-dpi':"300",'-clust':"Clusters_Objects.tsv", '-min':'0'}
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print("""
___________________________________________________________________________________________

Usage;

This code is for analysis of RNA-seq and data visualization
-def    cl: generate individual files for cluster in teh clust result.
            -clust      file name of cluster info by clust analysis (default is 'Clusters_Objects.tsv')
            -dge        file name of digital gene expression
            -rep        file name of replicate info if applicable (option)
            -log2       yes => convert dge into log2 value, no => do not convert (default is 'no')
        hmap: generate clustered heat map
            -consensus  consensus of names for multiple running
            -fig_w      width of the figure (default is 10)
            -fig_h      height of the figure (default is 20)
            -dpi        resolution of figure (default is 300)
        hmap_fc: generate clustered heat map by fold change
            -consensus  consensus of names for multiple running
            -fig_w      width of the figure (default is 10)
            -fig_h      height of the figure (default is 20)
            -dpi        resolution of figure (default is 300)
        bp: generate box plot
            -consensus  consensus of names for multiple running
            -fig_w      width of the figure (default is 10)
            -fig_h      height of the figure (default is 20)
            -dot_size   size of dots in the plot (default is 1.0)
            -dpi        resolution of figure (default is 300)
        vp: generate violin plot
            -consensus  consensus of names for multiple running
            -fig_w      width of the figure (default is 10)
            -fig_h      height of the figure (default is 20)
            -dot_size   size of dots in the plot (default is 1.0)
            -dpi        resolution of figure (default is 300)

        convert: convert replicates into single values by average (and log2 value). 
            -dge        file name of digital gene expression
            -rep        file name of replicate info (optional)
            -log2       yes => convert dge into log2 value, no => do not convert (default is 'no')
            -min        minimum value threshold of dges in each gene (default is 0)
___________________________________________________________________________________________
""")
                quit()


        
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import math
import os
import csv 
#dot_size=float(option_dict['-dot_size'])
#__________convert clust result into individual cluster files_____________________

def clust_divider(cl_file):
    print ("Dividing clust into separate files..\n")
    cl_name=cl_file
    cl=open(cl_name,'r')
    
    #__clust info parsing__
    clust_dic={}
    cl_lines=cl.readlines()
    init=1
    for row in cl_lines:
        row=row.strip()
        row=row.split('\t')
        if init==1:
            cl_names=row
            for cl in row:
                clust_dic[cl]=[]
            init=2
        elif init==2:
            init=3
            pass
        else:
            for i in range(len(cl_names)):
                #print row
                try:
                    if row[i]!="":
                        clust_dic[cl_names[i]].append(row[i])
                except:
                    pass
    return clust_dic

def Parsing_dge_data(dge_file, log2, cl_write):
    data_dic={}
    name_list=[]
    name_list.append("head")
    #minimum_value=float(min)
    #__parsing dge data__
    data_file=open(dge_file,'r')
    line=data_file.readline().split()
    new_first_line=(',').join(line) # generate head line
    data_dic['head']=new_first_line
    #__dge info parsing__
    while 1:
        line=data_file.readline().strip()
        if line=="":
            break
        line_list=line.split()
        #__filtering_with_min_value__
        if log2=="yes":
            if cl_write=="yes":
                new_line_list=[line_list[0]]
                temp=list(map(lambda x: str(round(math.log(float(x)+1,2),2)), line_list[1:]))
                line_list=new_line_list+temp

        new_line=(",").join(line_list)
        key=line_list[0]
        data_dic[key]=new_line
        name_list.append(key)
    return data_dic, name_list

def process_dge_data(rep_file, data_dic, log2, min, cl_write):#__replicate info parsing__
    minimum_value=float(min)
    print ("Calculate average dge of replicates..\n")
    if rep_file!="":
        rep_open=open(rep_file, 'r')
        rep_read=rep_open.readlines()
        sample_id_list=data_dic['head'].split(",")

        rep_cont=[]
        for line in rep_read:
            line=line.strip().split('\t')
            rep_list=list(map(lambda x: x.strip(), line[2].split(",")))
            index=[]
            for rep in rep_list:
                index.append(sample_id_list.index(rep))
            if log2!="yes":
                rep_cont.append([rep_list[0]+"_avg", index])
            else:
                rep_cont.append([rep_list[0]+"_log2_avg", index])
        
    #__convert dge info to average value based on replicates__
    new_data_dic={}
    for key, cont in data_dic.items():
        new_first_line=data_dic['head']
        new_first_line=[new_first_line.split(",")[0]] # generate new head line
        if key !="head":
            cont=cont.split(",")
            new_cont=[cont[0]]
            if rep_file!="":
                for rep in rep_cont:
                    new_first_line.append(rep[0])
                    value=0
                    #print (cont)
                    for ind in rep[1]:
                        value=value+float(cont[ind])
                    value=value/float(len(rep[1]))
                    new_cont.append(value)
                
                new_first_line=",".join(new_first_line)
                new_data_dic['head']=new_first_line

            else:
                new_data_dic["head"]=data_dic['head']
                new_cont=[cont[0]]+list(map(lambda x: float(x), cont[1:]))
            
            if max(new_cont[1:])>=minimum_value:
                if log2=="yes" and cl_write!="yes":
                    new_cont=[new_cont[0]]+list(map(lambda x: str(round(math.log(x+1,2),2)), new_cont[1:]))
                else:
                    new_cont=[new_cont[0]]+list(map(lambda x: str(round(x,2)), new_cont[1:]))
                new_cont=",".join(new_cont)
                new_data_dic[key]=new_cont
            
    
    data_dic=new_data_dic
    return data_dic

def generate_clust_files(clust_dic, data_dic, log2):#__generate clust files with dge data__
    for cl, g_names in list(clust_dic.items()):
        cl_list=cl.split()
        cl_concatenated="_".join(cl_list)
        if log2=="yes":
            out_file=open(cl_concatenated+"_log2.csv",'w')
        else:
            out_file=open(cl_concatenated+".csv",'w')
        new_first_line=data_dic['head']
        out_file.write(new_first_line+'\n')
        for g in g_names:
            #print g
            if g in data_dic:
                out_file.write(data_dic[g]+'\n')
            else:
                print(g)
    print("Completed!!")
            
def write_converted_data(dge_file_name, data_dic_converted, gene_list, rep):
    if rep !="":
        out_file=open(dge_file_name+"_converted_avg.txt",'w')
    else:
        out_file=open(dge_file_name+"_converted.txt",'w')
    for gene in gene_list:
        if gene in data_dic_converted:
            out_file.write("\t".join(data_dic_converted[gene].split(','))+"\n")

def boxplotting(consensus, w, h, dot, dpi):#__draw boxplot__ 
    dpi_value=int(dpi)
    dot_size=float(dot)
    width=float(w)
    height=float(h)
    name=consensus
    file_list_row=os.listdir('.')
    file_list=[]
    for f in file_list_row:
        if name in f:
            file_list.append(f)
    
    for name in file_list:
    
        df = pd.read_csv(name)
        df = df.set_index(df.columns[0])
        df.head()
        bplot=sns.boxplot(data=df, width=0.7, fliersize=1)
        bplot=sns.swarmplot(data=df, color="gray", size=dot_size)
        bplot.axes.set_title(name[:-4], fontsize=16)
        plt.savefig(name[:-4]+'_box_plot.jpg',dpi=dpi_value, figsize=(width,height))
        plt.close()

def violin_plotting(consensus, w, h, dot, dpi):#__draw violin plot__ 
    dpi_value=int(dpi)
    dot_size=float(dot)
    width=float(w)
    height=float(h)
    name=consensus
    file_list_row=os.listdir('.')
    file_list=[]
    for f in file_list_row:
        if name in f:
            file_list.append(f)
    
    for name in file_list:
        
        df = pd.read_csv(name)
        df = df.set_index(df.columns[0])
        df.head()
        bplot=sns.violinplot(data=df, width=0.7, inner=None, fliersize=1)
        bplot=sns.swarmplot(data=df, color="white", size=dot_size, edgecolor="black")
        bplot.axes.set_title(name[:-4], fontsize=16)
        plt.savefig(name[:-4]+'_violin_plot.jpg',dpi=dpi_value, figsize=(width,height))
        plt.close()

def heatmap_plotting(def_option, w, h, consensus, dpi):#__draw heatmap__ 
    dpi_value=int(dpi)
    width=float(w)
    height=float(h)

    name=consensus
    file_list_row=os.listdir('.')
    file_list=[]
    for f in file_list_row:
        if name in f:
            file_list.append(f)
    
    for name in file_list:
        
        df = pd.read_csv(name)
        df = df.set_index(df.columns[0])
        df.head()
        
        if def_option=='hmap':
            sns.clustermap(df, col_cluster=False,cmap=plt.cm.gist_heat_r, method='average', figsize=(width,height))
        elif def_option=='hmap_fc':
            sns.clustermap(df, col_cluster=False,cmap=plt.cm.seismic, method='average',figsize=(width,height))
        plt.savefig(name[:-4]+'_heat_map.jpg', dpi=dpi_value, figsize=(width,height))
        plt.close()

def main():
    if option_dict['-def']=="cl":
        clust_dic=clust_divider(option_dict['-clust'])
        data_dic=Parsing_dge_data(option_dict['-dge'], option_dict['-log2'], "yes")[0]
        if option_dict['-rep']!="":
            data_dic=process_dge_data(option_dict['-rep'], data_dic, option_dict['-log2'],0,"yes")
        generate_clust_files(clust_dic, data_dic, option_dict['-log2'])

    elif option_dict['-def']=="convert":
        data_dic_tuple=Parsing_dge_data(option_dict['-dge'], option_dict['-log2'], "no")
        data_dic=data_dic_tuple[0]
        gene_list=data_dic_tuple[1]
        data_dic_processed=process_dge_data(option_dict['-rep'], data_dic, option_dict['-log2'], option_dict['-min'],"no")
        write_converted_data(option_dict['-dge'], data_dic_processed, gene_list, option_dict["-rep"])

    elif option_dict['-def']=="hmap" or option_dict['-def']=="hmap_fc" :
        heatmap_plotting(option_dict['-def'], option_dict['-fig_w'], option_dict['-fig_h'], option_dict['-consensus'], option_dict['-dpi'])
    elif option_dict['-def']=="bp":
        boxplotting(option_dict['-consensus'], option_dict['-fig_w'], option_dict['-fig_h'], option_dict['-dot_size'], option_dict['-dpi'])
    elif option_dict['-def']=="vp":
        violin_plotting(option_dict['-consensus'], option_dict['-fig_w'], option_dict['-fig_h'], option_dict['-dot_size'], option_dict['-dpi'])

if __name__=='__main__':
    main()
