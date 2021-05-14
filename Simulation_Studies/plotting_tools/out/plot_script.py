import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import cm
from matplotlib.ticker import PercentFormatter

from pandas import DataFrame, read_csv
import pandas as pd 
from plot_methods import *



runtimes = pd.DataFrame(columns=["CHIMP-$T_{MRCA}$","CHIMP-$\mathcal{L}$","Relate","MSMC2"])#,"SMC++"])



######################
# BOTTLE DEMO
###################################################################################

dem_type = "bottle" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottle_3p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/3p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/3p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/3p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/3p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [400,5000],[3000,32000], [500,4000], maj_xpoints=[1000,2000],maj_xlabels=[r"$10^3$",r"$2\times10^3$"], min_ypoints=[5000,20000],min_ylabels=[r"$5\times10^3$",r"$2\times10^4$"])
extract_times(trajectories)
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Bottleneck: n=10"))



###################################################################################

dem_type = "bottle" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottle_3p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/3p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/3p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/3p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/3p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [400,5000],[3000,32000], [500,4000], maj_xpoints=[1000,2000], maj_xlabels=[r"$10^3$",r"$2\times10^3$"], min_ypoints=[5000,20000], min_ylabels=[r"$5\times10^3$",r"$2\times10^4$"])
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Bottleneck: n=200"))








######################
# EXP GROWTH DEMO
#########################################################################################################

dem_type = "exp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/exp_3p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/3p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/3p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/3p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/3p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [400,5000],[6000,2000000], [500,4000], maj_xpoints=[1000,2000], maj_xlabels=[r"$10^3$",r"$2\times10^3$"])
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Growth: n=10"))




###################################################################################

dem_type = "exp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/exp_3p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/3p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/3p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/3p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/3p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [400,5000],[6000,2000000],[500,4000], maj_xpoints=[1000,2000], maj_xlabels=[r"$10^3$",r"$2\times10^3$"])
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Growth: n=200"))












######################
# PWCSAWX DEMO
#########################################################################################################

dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/pwcsawX_15p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/15p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/15p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/15p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/15p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000],  [28,895484.])
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Piecewise Saw: n=10"))

###################################################################################

dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/pwcsawX_15p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/15p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/15p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/15p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/15p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000], [28,895484.])
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Piecewise Saw: n=200"))



#########################################################################################################
# PWCSAWX, PSEUDOHAPLOID COMPARISON
#######################################

dem_type = "pwcsawX" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/pseduo_pwcsawX_15p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$-PH"]= f"../../analysis_CHIMP/out/15p_10n_PH_{dem_type}_pseudo_dataset/{dem_type}_pseudo_dataset_CONSOLIDATED.csv"
#trajectories[r"CHIMP-$T_{MRCA}$-Reg"]= f"../../analysis_CHIMP/out/15p_10n_Reg_{dem_type}_pseudo_dataset/{dem_type}_pseudo_dataset_CONSOLIDATED.csv"
trajectories[r"Relate"]= f"../../analysis_relate/out/15p_10n_Reg_{dem_type}_pseudo_dataset/{dem_type}_pseudo_dataset_CONSOLIDATED.csv"
trajectories[r"MSMC2*"]= f"../../analysis_msmc2/out/15p_10n_Reg_{dem_type}_pseudo_dataset/{dem_type}_pseudo_dataset_CONSOLIDATED.csv"

plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,2500000],[1000,1000000],  [28,895484.], ncol=2,colors=['blue', 'green', 'purple'])








'''
######################
# SAWTOOTH DEMO
#########################################################################################################

dem_type = "sawtooth" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawtooth_20p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/20p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
#trajectories["SMC++"]= f"../analysis_smcpp/out/20pwp_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,100000],[2000,2000000],[34,73390])
#runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Saw: n=10"))


###################################################################################

dem_type = "sawtooth" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawtooth_20p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/20p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/20p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
#trajectories["SMC++"]= f"../analysis_smcpp/out/20pwp_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [20,100000],[2000,2000000],[34,73390])
#runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Saw: n=200"))


'''




######################
#SAWSPS DEMO
#########################################################################################################

dem_type = "sawSPS" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawSPS_20p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/20p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories,[13,120000],[500,300000],[27,58712])
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Saw: n=10"))


###################################################################################

dem_type = "sawSPS" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawSPS_20p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/20p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/20p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories,[13,120000],[500,300000],[27,58712])
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Saw: n=200"))















######################
# BOTTLE+EXP DEMO
#####################################################################################################

dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottleExp_20p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/20p_10n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/20p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
#trajectories["SMC++"]= f"../analysis_smcpp/out/20pwp_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [100,40000],[200,80000],[155,25831], lloc = 'lower right' )
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Collapse+Growth: n=10"))


###################################################################################

dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottleExp_20p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"CHIMP-$\mathcal{L}$"]= f"../../analysis_CHIMP/out/20p_200n_TL_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/20p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
#trajectories["SMC++"]= f"../analysis_smcpp/out/20pwp_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [100,40000],[200,80000],[155,25831], lloc = 'lower right')
runtimes = runtimes.append(pd.Series(extract_times(trajectories), name="Collapse+Growth: n=200"))















############################################
# SAWSPS with RECOMBINATION MAP   PLOTS
#####################################################################################################


dem_type = "sawSPS" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawSPS_RM_10n.pdf'
truth_file = f"./truth_trajectories/sawSPS_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/20p_10n_TH_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_10n_RM_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"
trajectories["MSMC2"]= f"../../analysis_msmc2/out/20p_10n_RMrf_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories,[13,120000],[500,300000],[27,58712],colors=['blue', 'green', 'purple'])





dem_type = "sawSPS" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/sawSPS_RM_200n.pdf'
truth_file = f"./truth_trajectories/sawSPS_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/20p_200n_TH_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"
trajectories["Relate"]= f"../../analysis_relate/out/20p_200n_RM_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"
trajectories["MSMC2*"]= f"../../analysis_msmc2/out/20p_100n_RMrf_sawSPS_RM_dataset/sawSPS_RM_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories,[13,120000],[500,300000],[27,58712],colors=['blue', 'green', 'purple'])















############################################
# BOTTLEEXP HEURISTIC PLOTS
##########################################################################################################################

dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottleExp_00p_10n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/00p_10n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"Relate"]= f"../../analysis_relate/out/00p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"MSMC2"]= f"../../analysis_msmc2/out/00p_10n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [40,50000],[200,80000],None, lloc = 'lower right',colors=['blue', 'green', 'purple'] )


###################################################################################

dem_type = "bottleExp" #'sawtooth', 'bottle', 'bottleExp', 'exp', 'pwcsaw'
write_file = './Fplots/bottleExp_00p_200n.pdf'
truth_file = f"./truth_trajectories/{dem_type}_truth.csv"

isLabels = []
isets = []
trajectories = dict()

trajectories[r"CHIMP-$T_{MRCA}$"]= f"../../analysis_CHIMP/out/00p_200n_TH_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"Relate"]= f"../../analysis_relate/out/00p_200n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"
trajectories[r"MSMC2*"]= f"../../analysis_msmc2/out/00p_100n_{dem_type}_dataset/{dem_type}_dataset_CONSOLIDATED.csv"


plotting_method(write_file,f"./truth_trajectories/{dem_type}_truth.csv", trajectories, [40,50000],[200,80000],None, lloc = 'lower right',colors=['blue', 'green', 'purple'] )










################################################################

runtimes.to_csv("./Fplots/avg_run_times.csv",index=True,header=True)
