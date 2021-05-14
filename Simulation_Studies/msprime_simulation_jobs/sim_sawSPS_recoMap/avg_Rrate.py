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

    

    

## Will Extract Avg Per Site Recombination Rate

def avg_rrate(file):

    t_df = pd.read_csv(file, sep='\t', header = 0)
    locs = t_df.loc[:,'Position(bp)'].to_numpy()
    rates = t_df.loc[:,'Rate(cM/Mb)'].to_numpy()

    weights =  np.roll(locs,-1) - locs
    weights[-1] = 0
    rates[-1] = 0.0

    return (np.sum(weights * rates) / np.sum(weights))
    

if __name__=="__main__":
    
    arrate = avg_rrate("./final_sims/msprime_simulation_jobs/sim_sawSPS_recoMap/chr3map_msp.txt")
    print(f"average recombination rate is: {arrate}")




    
