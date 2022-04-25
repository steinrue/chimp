import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import cm
from matplotlib.ticker import PercentFormatter
import sys
from pandas import DataFrame, read_csv
import pandas as pd
import re
import csv        







def getRunTimeStats(filename):
    cpu_time= "NA       "
    wall_use= "NA       "
    wall_max = "NA       "
    mem_use = "NA       "
    vmem_use = "NA       "
    EM_use = "NA       "
    EM_complete = "False    "
    LLconverge = "False    "
    Q_eval = 0

    with open(filename, "rt") as ofile:
        line = ofile.readline()
        while(line):

            #cpu_time
            s_s = "cput="
            e_s = ",vmem="
            search = re.search(f'{s_s}(.*){e_s}', line)
            if search:
                cpu_time = search.group(1)
                
            #wall_use
            s_s = ",walltime="
            e_s = ",mem="
            search = re.search(f'{s_s}(.*){e_s}', line)
            if search:
                wall_use = search.group(1)

            #wall_max
            s_s = "walltime="
            e_s = ",nodes"
            search = re.search(f'{s_s}(.*){e_s}', line)
            if search:
                wall_max = search.group(1)


            #mem_use
            s_s = ",mem="
            e_s = ",energy"
            search = re.search(f'{s_s}(.*){e_s}', line)
            if search:
                mem_use = search.group(1)

            #vmem_use
            s_s = ",vmem="
            e_s = ",walltime"
            search = re.search(f'{s_s}(.*){e_s}', line)
            if search:
                vmem_use = search.group(1)

            #EM_use
            s_s = "after "
            e_s = " iteration, ll"
            search = re.search(f'{s_s}(.*){e_s}', line)
            if search:
                EM_use = search.group(1) + "       "

            #EM_completed
            search = re.search('Done!!', line)
            if search:
                EM_complete = "True    "
                
            #LLconverge
            search = re.search("CONVERGENCE COMPLETE!", line)
            if search:
                LLconverge = "True    "

            #LLconverge
            search = re.search("Q_eval:", line)
            if search:
                Q_eval+=1

            # last line of while loop, move to next line
            line = ofile.readline()

    
    return [cpu_time,wall_use,wall_max,mem_use,vmem_use, EM_use,EM_complete, LLconverge, Q_eval]

   
if __name__=="__main__":
    
    # name of folder where we want to read files and create consolidated result
    
    folder = sys.argv[1]
    prefix = sys.argv[2]
    replicates = int(sys.argv[3])
    time_catalog = False
    if(len(sys.argv) >4):
        time_catalog = bool(sys.argv[4])


    # CREATE CONSOLIDATED DEMINF RESULTS FILE

    # get x_column from file 1
    iset_df =  pd.read_csv(f"{folder}/{prefix}16/deminf_results.csv", header = 0)
    x_ind = -2
    for head in iset_df.columns.tolist():
        if(head.strip() == 'x'):
            x_ind = head
    x_col = iset_df.get(x_ind)
        
    iset_df=pd.DataFrame()
    iset_df.insert(0,'x',x_col)


    # create dataframe for the runtimes
    times_df = pd.DataFrame()
    time_fields = ["cpu_time ","wall_use","wall_max","mem_use ","vmem_use","EM_use  ","EM_complete","LLconverge", "Q_eval     "]
    if(time_catalog):
        times_df.insert(0,"fields  ",time_fields)

    # cycle through each one individually

    for i in range(replicates):


        print("----")

        # GET THE DEMHIST from deminf_results

        try:
            t_df =  pd.read_csv(f"{folder}/{prefix}{i+1}/deminf_results.csv", header = 0)
        except FileNotFoundError:
            print(f"!!! DID NOT FIND FILE !!!")
            continue

        # see if finished properly, otherwise raise alarm
        try:
            outFilename = f"{folder}/{prefix}{i+1}/{prefix}{i+1}.out"
            ifs =  open (outFilename, "r")
            properFinish = False
            exitCode = None
            for line in ifs:
                fields = line.strip().split(":")
                if(fields[0].startswith('Exit Code')):
                    exitCode = int(fields[1])
                    if (exitCode != 0):
                        break
                    else:
                        properFinish = True
            if (not properFinish): 
                    print (f"!!! RUN DID NOT FINISH (YET) !!!")
                    print (outFilename)
                    if (exitCode is None):
                        print ("!!! NO EXIT CODE GIVEN !!!")
                    else:
                        print (f"!!! EXIT CODE: {exitCode} !!!")
                    continue
        except FileNotFoundError:
            print(f"!!! NO OUTPUT FILE FROM CLUSTER !!!")
            continue


        y_head = ''
    
        # isolate the headers, find column with y vals
        for header in  t_df.columns.tolist():
            if('y' == header.strip()):
                y_head = header

        # write the correct row over into the consolidated structure
        vals = t_df.to_numpy()
        print (i+1)
        iset_df.insert(len(iset_df.columns), 'y'+str(i+1), np.array(t_df[y_head]) )

        # print the column we just copied
        print( np.array(t_df[y_head]) )


        # PERFORM THE TIME CATALOG FOR CHIMP OUT FILES
        if(time_catalog):
            out_file_stats = getRunTimeStats(f"{folder}/{prefix}{i+1}/{prefix}{i+1}.out")
            times_df.insert(len(times_df.columns), "set"+str(i+1)+"    ", out_file_stats )


    
    iset_df.to_csv(path_or_buf= f"{folder}/{prefix}_CONSOLIDATED.csv", index = False)

    # write tsv for the runtimes
    with open(f"{folder}/{prefix}_RUNTIMES.tsv", 'w', newline='') as times_out:
        tsv_output = csv.writer(times_out, delimiter='\t')
        df_data = times_df.to_numpy()
        tsv_output.writerow(np.insert(df_data[:,0],0,"set     "))
        for i in range(1,df_data.shape[1]):
            tsv_output.writerow(np.insert(df_data[:,i],0,times_df.columns[i]))





