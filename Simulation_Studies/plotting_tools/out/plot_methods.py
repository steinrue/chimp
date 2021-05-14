import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import cm
from matplotlib.ticker import PercentFormatter
from collections import deque
from scipy.interpolate import interp1d
from matplotlib.ticker import *#(FixedLocator, FixedFormatter, NullLocator)

from pandas import DataFrame, read_csv
import pandas as pd 


def extract_times(label_files):

    def string_to_sec(t_string):
        hms= t_string.strip().split(':',2)
        return 3600. * int(hms[0]) + 60. * int(hms[1]) + 1.0 * int(hms[2])

    def sec_to_string(seconds):
        seconds = int(seconds)
        hrs = seconds / 3600.0
        return f"{hrs:.02f}"
    
    output = dict()
    labels = list(label_files)
    for lab in labels:

        times_file = label_files[lab].rsplit("CONSOLIDATED.csv",1)[0] + "RUNTIMES.tsv"
        times_df = pd.read_csv(times_file ,sep='\t', header=0)
        #rtimes = times_df.loc[:,'wall_use'].to_numpy()
        rtimes = times_df.loc[:,'cpu_time '].to_numpy()

        tot_sec = 0.0
        num_ent = 0
        for time in rtimes:
            tot_sec = tot_sec + string_to_sec(time)
            num_ent = num_ent + 1

        avg_time = tot_sec / num_ent 

        lab = lab.rsplit("*",1)[0]
        
        output[lab] = sec_to_string(avg_time)

    return output


# Perform smoothing on dif_avg, each value is avg of 10 adjacent values.
def smooth_list(func_vals, w_size):

        temp = func_vals.copy()
        
        def win_val(www): # function to get value from sliding window
            return (np.sum(www) - np.amin(www) - np.amax(www))/ (np.size(www)-2)
            #return np.mean(www)

        window = deque(temp[0:w_size+1],maxlen=(2 * w_size + 1)) # initialize with window centered on first index (half is empty)
        temp[0] = win_val(window) # get first val

        for i in range(w_size+1,len( temp)): # add values till reach end of main list
            window.append( temp[i])
            temp[i - w_size] = win_val(window)

        for i in range(w_size): # shorten and find vals for remaining half window
            window.popleft()
            temp[i - w_size] = win_val(window)

        return temp
        


def plotting_method(out_file, truth_file, trajectories, xlims, ylims,int_lims, title='', colors=['blue', 'red', 'green', 'purple', 'orange', 'cyan'], show=False, **kwargs):

    ####
    # kwargs can contain: 
    #'lloc' (legend location)
    #'ncol' (number of legend columns)
    #'min_xpoints' and min_xlabels (lists of minor x tick points and their respective labels (to add in ADDITION to what's there))
    
    i_colors = colors

    # SETUP TRUTH DATA
    df_tr = pd.read_csv( truth_file, header = 0).values
    df_tr= np.insert(df_tr, 0, df_tr[0,:],0)
    df_tr = np.row_stack((df_tr,df_tr[-1,:]))
    df_tr[0,0] = min(xlims[0],df_tr[1,0])
    df_tr[-1,0] = max(xlims[1], df_tr[-2,0])

    # FUNCTION FOR EVALUATING TRUTH

    def truth_eval(domain):
        if(isinstance(domain,np.ndarray)):
            out = np.copy(domain)
            for i in range(out.size):
                out[i] = truth_eval(out[i])
            return out

        temp = np.where((df_tr[:,0] > domain))[0]
        if temp.size == 0:
            return df_tr[:,1][-1]
        else:
            return df_tr[:,1][temp[0]]

 
    
    fig, [ax1,ax2] = plt.subplots(2,1,figsize=(10,8), sharex=True, gridspec_kw={'height_ratios':[3.2,1]}) # use if legends are on plot

    ###################################
    # FIRST PLOT
    ax1.loglog()
    ax1.set_xlim(xlims) 
    ax1.set_ylim(ylims)

    nan_switch = False

    ax1.plot(df_tr[:,0], df_tr[:,1],  "black", linewidth = 4.5 , label='Truth' )
    
    all_labels = list(trajectories.keys())
    for iii in range(len(all_labels)):
        i_label = all_labels[iii]
        iset = trajectories[i_label]
    
        # setup data frames and values for the datasets
        iset_df = pd.read_csv(iset , header=0)
        iset_values = iset_df.values
        iset_values = np.insert(iset_values, 0, iset_values[0,:],0)
        iset_values = np.row_stack((iset_values,iset_values[-1,:]))
        iset_values[0,0] = min(xlims[0],iset_values[1,0])
        iset_values[-1,0] = max(xlims[1], iset_values[-2,0])

    
        # setup and plot the mean of the inferred trajectories
        iset_avg = np.copy(iset_values[:,0])
        iset_ln_std = np.copy(iset_avg)

        for i in range(len(iset_avg)):
            ln_pop =  np.nan_to_num( np.log(iset_values[i,1::]), posinf = np.nan )
            if(np.any(np.isnan(ln_pop))):
                nan_switch=True
            ln_pop = ln_pop[np.logical_not(np.isnan(ln_pop))]

            
            iset_avg[i] = np.exp(np.mean(ln_pop))
            iset_ln_std[i] = np.exp(np.std(ln_pop))
        
        iset_stds = np.array([iset_avg * (1. - 1.0/iset_ln_std) ,iset_avg * (iset_ln_std - 1)])
        ax1.plot( iset_values[:,0],  iset_avg, color = i_colors[iii], linewidth = 2., label=i_label)
        ax1.fill_between(iset_values[:,0], iset_avg - iset_stds[0], iset_avg+iset_stds[1], alpha=.15, color = i_colors[iii])

        print(len(iset_values[0,1::]))
        print(nan_switch)

    #############################
    # SECOND PLOT
    
    epsilon = .001
    x_evengrid = np.exp(np.arange(np.log(xlims[0] + epsilon), np.log(xlims[1] - epsilon), np.log(xlims[1]/xlims[0])/1000.  ))

    
    ax2.plot(df_tr[:,0], [0]* df_tr[:,0].size,  "black", linewidth = 4.5 , label='Truth' )


    for iii in range(len(all_labels)):
        i_label = all_labels[iii]
        iset = trajectories[i_label]
    
        # setup data frames and values for the datasets
        iset_df = pd.read_csv(iset , header=0)
        iset_values = iset_df.values
        iset_values = np.insert(iset_values, 0, iset_values[0,:],0)
        iset_values = np.row_stack((iset_values,iset_values[-1,:]))
        iset_values[0,0] = min(xlims[0],iset_values[1,0])
        iset_values[-1,0] = max(xlims[1], iset_values[-2,0])
        
        # setup and plot the mean of the inferred trajectories
        dif_avg = np.copy(iset_values[:,0])
        dif_std = np.copy(dif_avg)

        for i in range(len(dif_avg)):
            true_v =  truth_eval(iset_values[i,0])

            # define distance measure
            distance_plot = np.log(iset_values[i,1::]/true_v)
            #distance_plot = np.square(distance_plot) ##comment for mean distance, uncomment for mean square distance

            # clean to remove positive infinities
            distance_plot = np.nan_to_num(distance_plot, posinf = np.nan)
            distance_plot = distance_plot[np.logical_not(np.isnan(distance_plot))]

            dif_avg[i] = np.mean(np.abs(distance_plot)) # abs val here stops under and over trajectories from averaging out
            #dif_std[i] = np.exp(np.std(ln_pop))

        #######

        # get single distance measure if limits provided

        phi_val_string = ""

        if(int_lims != None):
            # get single distance measure
            delx = np.log(iset_values[:,0]) # log of the domain
            delx = np.nan_to_num(delx - np.roll(delx,1), nan=0., posinf = 0., neginf= 0.) # avoid issues with log of zeros 
            delx[0] = 0 # first one is meaningless since end rolled over to beginning

            #distance measures
            distances = np.abs(dif_avg * delx) # distance along entire domain
            distances = distances[ np.where(np.logical_and( iset_values[:,0] >= int_lims[0], iset_values[:,0] <= int_lims[1])) ] # distances in the xlims for plot
            int_dist = np.sum(distances) # integrated distance measure

            phi_val_string = fr" ($\phi$={int_dist:.1f})"

        ###################

        # use uniform grid for all the methods based on x_evengrid
        # this ensures that smoothing treats them all the same way
        interp_f = interp1d(iset_values[:,0], dif_avg)
        y_evengrid = interp_f(x_evengrid)
        
        ax2.plot(x_evengrid, smooth_list(y_evengrid,40), color = i_colors[iii], linewidth = 2., label=i_label+phi_val_string)
        

    ##
    # WRAP UP FIGURE DETAILS

    textSize=22
    ax2.set_xlabel('Generations Before Present', fontsize=textSize)
    ax1.set_ylabel('Population Size', fontsize=textSize)
    vv = r"$\Delta(k)$"
    ax2.set_ylabel(f"Mean Signed\nError {vv}", wrap=True, fontsize=textSize)
    ax1.set_title(title, fontsize=textSize)

    
    lloc = "upper right"
    if 'lloc' in kwargs.keys():
        lloc = kwargs['lloc']

    col_ns = 2
    if 'ncol' in kwargs.keys():
        col_ns = kwargs['ncol']

    ax1.legend(*ax2.get_legend_handles_labels(),ncol = col_ns,fontsize = textSize-2, loc=lloc)

    if('min_ypoints' in  kwargs.keys()):
        ylabs = ax1.yaxis.get_minorticklocs()
        l_i = 0
        print(ylabs)
        print(kwargs['min_ypoints'])
        for iii in range(len(ylabs)):
            if not(ylabs[iii] in kwargs['min_ypoints']):
                ylabs[iii]=""
            else:
                print("Bingo")
                ylabs[iii]=kwargs['min_ylabels'][l_i]
                l_i = l_i + 1

        print(ylabs)
        ax1.yaxis.set_minor_formatter(FixedFormatter(ylabs))
        #ax1.yaxis.set_minor_formatter(FixedFormatter(kwargs['min_ylabels']))
        

    if ('maj_xpoints' in kwargs.keys()):

        ax1.xaxis.set_major_locator( matplotlib.ticker.FixedLocator(kwargs['maj_xpoints']))
        ax1.xaxis.set_major_formatter(FixedFormatter(kwargs['maj_xlabels']))


        
    ax1.tick_params(axis='both', which='major', labelsize=textSize-2, length=7, width = .9)
    ax2.tick_params(axis='both', which='major', labelsize=textSize-2, length=7, width = .9)

    ax2.tick_params(axis='x', which = 'minor', labelsize=textSize-2)
    ax1.tick_params(axis='y', which = 'minor', labelsize=textSize-8)




    plt.tight_layout()
    if show:
        plt.show()
    fig.savefig(out_file)
    plt.close(fig)













    
if __name__=="__main__":

    print("X")
  
