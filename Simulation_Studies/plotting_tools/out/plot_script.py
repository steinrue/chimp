import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import shutil
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import cm
from matplotlib.ticker import PercentFormatter

from pandas import DataFrame, read_csv
import pandas as pd 

import os

from plot_methods import *


runtimes = pd.DataFrame(columns=["CHIMP-$\mathcal{T}_{10}$", "CHIMP-$\mathcal{T}_{2,5,10}$","CHIMP-$\mathcal{L}_{10}$", "Relate", "MSMC2"])

#########
# directories that we need
#########

# if you want to overwrite a previous analysis, you have to manually comment this out
os.makedirs ("./truth_trajectories/")
os.makedirs ("./Fplots/")


######################
# BOTTLE DEMO
###################################################################################
print ("---------- BOTTLE DEMO ----------")

dem_type = "bottle" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottle_3p_10n.pdf'

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/3p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/3p_10n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/3p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/3p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/3p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [400,5000],[3000,32000], [500,4000], maj_xpoints=[1000,2000],maj_xlabels=[r"$10^3$",r"$2\times10^3$"], min_ypoints=[5000,20000],min_ylabels=[r"$5\times10^3$",r"$2\times10^4$"])

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Bottleneck: n=10"))

###################################################################################

dem_type = "bottle" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottle_3p_200n.pdf'

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/3p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/3p_200n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/3p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/3p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/3p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [400,5000],[3000,32000], [500,4000], maj_xpoints=[1000,2000], maj_xlabels=[r"$10^3$",r"$2\times10^3$"], min_ypoints=[5000,20000], min_ylabels=[r"$5\times10^3$",r"$2\times10^4$"])

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Bottleneck: n=200"))


######################
# EXP GROWTH DEMO
#########################################################################################################
print ("---------- EXP GROWTH DEMO ----------")


thisAnalysis = ""

dem_type = "exp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = f'./Fplots/exp_3p_10n{thisAnalysis}.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/3p_10n_TH_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/3p_10n_TH_n2510_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/3p_10n_TL_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/3p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/3p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [400,5000],[6000,2000000], [500,4000], maj_xpoints=[1000,2000], maj_xlabels=[r"$10^3$",r"$2\times10^3$"])

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Piecewise Growth: n=10"))

###################################################################################

dem_type = "exp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/exp_3p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/3p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/3p_200n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/3p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/3p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/3p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [400,5000],[6000,2000000],[500,4000], maj_xpoints=[1000,2000], maj_xlabels=[r"$10^3$",r"$2\times10^3$"])

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Piecewise Growth: n=200"))


######################
# PWCSAWX DEMO
#########################################################################################################
print ("---------- PWCSAWX DEMO ----------")


thisAnalysis = ""

dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = f'./Fplots/pwcsawX_15p_10n{thisAnalysis}.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/15p_10n_TH_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/15p_10n_TH_n2510_no_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/15p_10n_TL_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/15p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/15p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000],  [28,895484.])

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Piecewise Sawtooth: n=10"))

###################################################################################

dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/pwcsawX_15p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/15p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/15p_200n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/15p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/15p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/15p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000], [28,895484.])

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Piecewise Sawtooth: n=200"))


######################
# SAWSPS DEMO
#########################################################################################################
print ("---------- SAWSPS DEMO ----------")


dem_type = "sawSPS" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawSPS_20p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/20p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories,[13,120000],[500,300000],[27,58712])

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Sawtooth: n=10"))

###################################################################################

dem_type = "sawSPS" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawSPS_20p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/20p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/20p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories,[13,120000],[500,300000],[27,58712])

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Sawtooth: n=200"))


######################
# BOTTLE+EXP DEMO
#####################################################################################################
print ("---------- BOTTLE+EXP DEMO ----------")


thisAnalysis = ""

dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = f'./Fplots/bottleExp_20p_10n{thisAnalysis}.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_n2510_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/20p_10n_TL_{dem_type}_dataset{thisAnalysis}/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [100,40000],[200,80000],[155,25831], lloc = 'lower right' )

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Bottleneck+Growth: n=10"))

###################################################################################

dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottleExp_20p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/20p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/20p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [100,40000],[200,80000],[155,25831], lloc = 'lower right')

runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Bottleneck+Growth: n=200"))


#########################################################################################################
# PWCSAWX, PSEUDOHAPLOID COMPARISON
#######################################
print ("---------- PWCSAWX, PSEUDOHAPLOID COMPARISON ----------")

dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/pseudo_pwcsawX_15p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$-PH"]= f"../../analysis_CHIMP/out/15p_10n_TH_n2510_PH_{dem_type}_pseudo_dataset/{dem_type}_pseudo_dataset_CONSOLIDATED.csv"
trajectories[r"Relate"]= f"../../analysis_relate/out/15p_10n_Reg_{dem_type}_pseudo_dataset/{dem_type}_pseudo_dataset_CONSOLIDATED.csv"
trajectories[r"MSMC2"]= f"../../analysis_msmc2/out/15p_10n_Reg_{dem_type}_pseudo_dataset/{dem_type}_pseudo_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000], None, ncol=2,colors=['red', 'purple', 'orange'])


############################################
# BOTTLEEXP HEURISTIC PLOTS
##########################################################################################################################
print ("---------- BOTTLEEXP HEURISTIC PLOTS ----------")


dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottleExp_00p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/00p_10n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"Relate"]= f"../../analysis_relate/out/00p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"MSMC2"]= f"../../analysis_msmc2/out/00p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [40,50000],[200,80000],None, lloc = 'lower right',colors=['red', 'purple', 'orange'] )

###################################################################################

dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottleExp_00p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/00p_200n_TH_n2510_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"Relate"]= f"../../analysis_relate/out/00p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"MSMC2*"]= f"../../analysis_msmc2/out/00p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [40,50000],[200,80000],None, lloc = 'lower right',colors=['red', 'purple', 'orange'] )


############################################
# SAWSPS with RECOMBINATION MAP   PLOTS
#####################################################################################################
print ("---------- SAWSPS with RECOMBINATION MAP ----------")


dem_type = "sawSPS" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawSPS_RM_10n.pdf'
truth_file = f"./truth_trajectories/sawSPS_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_n2510_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_10n_RM_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/20p_10n_RMrf_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories,[13,120000],[500,300000],[27,58712],colors=['red', 'purple', 'orange'])

###################################################################################

dem_type = "sawSPS" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawSPS_RM_200n.pdf'
truth_file = f"./truth_trajectories/sawSPS_truth.csv"

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_n2510_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_200n_RM_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/20p_100n_RMrf_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories,[13,120000],[500,300000],[27,58712],colors=['red', 'purple', 'orange'])


######################
# PWCSAWX DEMO -- different single ns
#########################################################################################################
print ("---------- PWCSAWX DEMO -- different single ns ----------")


dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/pwcsawX_15p_10n_TH_single_n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{2}$-o"]= f"../../analysis_CHIMP/out/15p_10n_TH_n2_o_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2}$"]= f"../../analysis_CHIMP/out/15p_10n_TH_n2_no_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{5}$"]= f"../../analysis_CHIMP/out/15p_10n_TH_n5_no_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{10}$"]= f"../../analysis_CHIMP/out/15p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/15p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000],  [28,895484.])

###################################################################################

dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/pwcsawX_15p_10n_TL_single_n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{L}_{2}$-o"]= f"../../analysis_CHIMP/out/15p_10n_TL_n2_o_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{2}$"]= f"../../analysis_CHIMP/out/15p_10n_TL_n2_no_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{5}$"]= f"../../analysis_CHIMP/out/15p_10n_TL_n5_no_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{10}$"]= f"../../analysis_CHIMP/out/15p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/15p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000],  [28,895484.])


######################
# PWCSAWX DEMO -- mutliple ns
#########################################################################################################
print ("---------- PWCSAWX DEMO -- mutliple ns ----------")


dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/pwcsawX_15p_10n_multi_n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}$"]= f"../../analysis_CHIMP/out/15p_10n_TH_n2510_no_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}_{2,5,10}$"]= f"../../analysis_CHIMP/out/15p_10n_TL_n2510_no_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/15p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000],  [28,895484.],colors=['blue', 'red', 'orange'])


######################
# BOTTLE+EXP REGULARIZATION
#####################################################################################################
print ("---------- BOTTLE+EXP REGULARIZATION ----------")


dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = f'./Fplots/bottleExp_20p_10n_reg.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

shutil.copyfile (f"../../data/{dem_type}/{dem_type}_dataset_truth.csv", f"./truth_trajectories/{dem_type}_truth.csv")

trajectories = dict()

trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}^{\:\:(-7)}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_n2510_l-7_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}^{\:\:(-6)}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_n2510_l-6_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}^{\:\:(-5)}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_n2510_l-5_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}^{\:\:(-4)}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_n2510_l-4_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{T}_{2,5,10}^{\:\:(-3)}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_n2510_l-3_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [100,40000],[200,80000],[155,25831], lloc = 'lower right' )


################################################################
# SAVE FINAL RUNTIMES
#########################################################################################################

runtimes.to_csv("./Fplots/avg_run_times.csv",index=True,header=True)
