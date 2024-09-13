#!/usr/bin/env python3

# Runs spin permutation for within and between network dispersion hormonal contrasts with lmer (on data excluding day 26)


### Define parameters

surface_name='fsa5'
parcellation_name='schaefer_400'  # I don't specify that it's Schaefer 400 17 networks but I think it's ok because here what matters is the number and shape of parcels for spin rotations
n_rot=1000  # number of permutations


### Load required packages
import os
import sys
import numpy as np
import pandas as pd

from enigmatoolbox.permutation_testing import centroid_extraction_sphere
from enigmatoolbox.permutation_testing import rotate_parcellation
import statsmodels.regression.mixed_linear_model as sm
from enigmatoolbox.datasets import load_fsa5


### Define directories

datadir = '/data/p_02667/28Me/data/'
resdir = '/data/p_02667/28Me/results/'





### Load the required data

print(f'\nLoad the required data')


# aligned gradient
array_aligned_fc_G1_excl26 = np.genfromtxt(resdir+'array_aligned_fc_G1_excl26.csv', delimiter=',')

# hormone data
hormones_28Me_excl26 = pd.read_csv(resdir+'hormones_28Me_excl26.csv', index_col = 0)


# Yeo network array (Sofie's 17 to 7 network conversion for Schaefer 400 (17 network) parcel order)
# labels: 1=visual, 2=sensory motor, 3=dorsal attention, 4=ventral attention, 5=limbic, 6=fronto parietal, 7= DMN
yeo_17_2_7 = np.genfromtxt(datadir+'yeo_17_2_7.csv', delimiter=',')


# root_pth = os.path.dirname(__file__)  

# root_pth via os.pathdirname(__file__) did not work so me so manually directed to path to downlowaded enigmatoolbox 
# ie via 3 commands (see https://enigma-toolbox.readthedocs.io/en/latest/pages/01.install/index.html): 
## git clone https://github.com/MICA-MNI/ENIGMA.git
## cd ENIGMA
## python setup.py install
root_pth = '/data/p_02667/sex_diff_gradients/scripts/ENIGMA/enigmatoolbox/permutation_testing/'



# load fsa5 surface
sphere_lh, sphere_rh = load_fsa5(as_sphere=True)

# get sphere coordinates of parcels
annotfile_lh = os.path.join(root_pth, 'annot', surface_name + '_lh_' + parcellation_name + '.annot')
annotfile_rh = os.path.join(root_pth, 'annot', surface_name + '_rh_' + parcellation_name + '.annot')

lh_centroid = centroid_extraction_sphere(sphere_lh.Points, annotfile_lh)
rh_centroid = centroid_extraction_sphere(sphere_rh.Points, annotfile_rh)




### Generate permutation maps (rotated centroids n_rot times)

print(f'\nGenerate permutation maps')

perm_id = rotate_parcellation(lh_centroid, rh_centroid, n_rot)

nroi = perm_id.shape[0]   # number of regions
nperm = perm_id.shape[1]  # number of permutations




### Permutation of yeo7 network labels 

print(f'\nPermutate Yeo network labels')

yeo7_networks_array_perm = np.empty((0, nroi))

for rr in range(nperm):
    
    print(f'---- Permuation No. {rr+1} ----')
    
    yeo7_networks_array_perm2 = []
    
    for ii in range(nroi):
        yeo7_networks_array_perm2 = np.append(yeo7_networks_array_perm2, yeo_17_2_7[int(perm_id[ii, rr])])

    yeo7_networks_array_perm = np.vstack((yeo7_networks_array_perm, yeo7_networks_array_perm2))




### Compute within and between network dispersion & hormonal contrast -- same code as for original, only differences: 
# WN: don't take the real yeo7_networks_array labels, but iterate over the permuted ones in yeo7_networks_array_perm
# BN: how data is saved

print(f'\nCompute within and between network dispersion - as well as hormonal contrast - for {n_rot} permutations')

# dictionaries of lists that will contain the permutation t-values (null distribution of t-values) for WN and BN dispersion hormonal contrasts
WN_dispersion_perm_tval_hormonal_contrast = {'estradiol_t_val': [], 'progesterone_t_val': []}
BN_dispersion_perm_tval_hormonal_contrast = {'estradiol_t_val': [], 'progesterone_t_val': []}



# iterate over the i permutations of the yeo network labels
for i in range(len(yeo7_networks_array_perm)):
    
    print(f'\n---- Permuation No. {i+1} ----')
    
    
    print(f'---- Computed WN dispersion ----')
    
    yeo_cog_perm = []  # center of gravity (median) for each network (i.e., network centroid position) (7 x 29)
    WN_dispersion_perm = []  # Within network dispersion: sum squared Euclidean distance of network nodes to the network centroid at experimental day level (7 x 29)


    # gradient values
    g1 = array_aligned_fc_G1_excl26.T  # transpose to obtain shape (400 x 29) in order to access/index the relevant network nodes


    # dictionary of lists that will contain the current permuation's t values (len of lists is 7 because will contain all 7 networks per permutation iteration)
    WN_dispersion_perm_tval_hormonal_contrast_temp = {'estradiol_t_val': [], 'progesterone_t_val': []}

    
    # iterate over the 7 Yeo networks
    for n in range(7):
        
        print(f'Network: {n+1}')

        # identify the nodes of given Yeo network
        netNodes = np.where(yeo7_networks_array_perm[i] == (n + 1))  # here we take the ith permutation of the yeo network array labels !!
        netNodes = np.squeeze(netNodes)

        # get the gradient loadings of the nodes of the given Yeo network, for each experimental day (shape: number of nodes in network x N)
        G1_net = g1[netNodes]


        ### identify the centroid / center of gravity (= median) of the given Yeo network for each experumental day (shape: N)
        yeo_cog_perm_net = np.median(G1_net, axis=0)  
        yeo_cog_perm.append(yeo_cog_perm_net)


        ### within network dispersion: 1 within network dispersion value per experimental day (per network)

        # compute (per subject) the Eucledian distance between each gradient loading (in Yeo network) and that network's centroid
        dist_nodes_to_cog = G1_net - yeo_cog_perm_net  # shape: number of nodes in network x N

        # take the sum of squares of this Eucledian distance 
        sum_of_squares = np.sum((dist_nodes_to_cog**2), axis = 0)  # shape: N

        # append to list
        WN_dispersion_perm.append(sum_of_squares)


        
        ### compute linear model for current network (WN dispersion)
        print(f'lmer in progress...')
        
        
        # define model
        model = sm.MixedLM(endog = sum_of_squares, exog = pd.DataFrame({'estradiol': hormones_28Me_excl26.z_estradiol, 'progesterone': hormones_28Me_excl26.z_progesterone}), groups = hormones_28Me_excl26.Experiment, exog_re=hormones_28Me_excl26.ExperimentDay, exog_vc=None)
			
        # fit model
        results = model.fit()
        #results.summary()

        # save results
        WN_dispersion_perm_tval_hormonal_contrast_temp['estradiol_t_val'].append(results.tvalues['estradiol'])
        WN_dispersion_perm_tval_hormonal_contrast_temp['progesterone_t_val'].append(results.tvalues['progesterone'])
        
    
    # append to permutation results
    WN_dispersion_perm_tval_hormonal_contrast['estradiol_t_val'].append(WN_dispersion_perm_tval_hormonal_contrast_temp['estradiol_t_val'])
    WN_dispersion_perm_tval_hormonal_contrast['progesterone_t_val'].append(WN_dispersion_perm_tval_hormonal_contrast_temp['progesterone_t_val'])

    
    
    
    ### Compute between network dispersion (Eucledian distance between two network centroids) only once for each pariwise network comparison, whilst keeping track what networks are being compared
    
    print(f'---- Computed BN dispersion ----')
    
    # dictionary of lists that will contain the current permuation's t values (len 21 because will contain all 21 comparisons of networks per permutation iteration)
    BN_dispersion_perm_tval_hormonal_contrast_temp = {'estradiol_t_val': [], 'progesterone_t_val': []}
    
    
    # to keep track the order in which the pairwise comparisons are computed
    pairwise_comparison_order = []

    # iterate over 7 Yeo networks
    for n1 in range(7):
    
        for n2 in range(7):

            current_pairwise_comparison = [n1+1, n2+1]


            if n1!=n2 and [n1+1, n2+1] not in pairwise_comparison_order and list(reversed([n1+1, n2+1])) not in pairwise_comparison_order:
                
                print(f'Networks: {current_pairwise_comparison}')

                # append the pairwise comparison to dict
                pairwise_comparison_order.append(current_pairwise_comparison)

                # compute the distance between centroid of network 1 and centroid of network 2 and append it to dict
                BN_dispersion = yeo_cog_perm[n1] - yeo_cog_perm[n2]
                


                ### compute linear model for current pairwise between network dispersion (distance between centroids)
                print(f'lmer in progress...')
                
                # define model
                model = sm.MixedLM(endog = BN_dispersion, exog = pd.DataFrame({'estradiol': hormones_28Me_excl26.z_estradiol, 'progesterone': hormones_28Me_excl26.z_progesterone}), groups = hormones_28Me_excl26.Experiment, exog_re=hormones_28Me_excl26.ExperimentDay, exog_vc=None)

                # fit model
                results = model.fit()
                #results.summary()

                # save results
                BN_dispersion_perm_tval_hormonal_contrast_temp['estradiol_t_val'].append(results.tvalues['estradiol'])
                BN_dispersion_perm_tval_hormonal_contrast_temp['progesterone_t_val'].append(results.tvalues['progesterone'])
    
    
    # append to permutation results
    BN_dispersion_perm_tval_hormonal_contrast['estradiol_t_val'].append(BN_dispersion_perm_tval_hormonal_contrast_temp['estradiol_t_val'])
    BN_dispersion_perm_tval_hormonal_contrast['progesterone_t_val'].append(BN_dispersion_perm_tval_hormonal_contrast_temp['progesterone_t_val'])


    
    
    
### Clean results for export

print(f'\n---- Export results at /data/p_02667/28Me/results/WN and BN_dispersion_perm_tval_hormonal_contrast_null_distr.csv ----')

# pack into arrays (estradiol and progesterone separately)
WN_dispersion_perm_tval_estradiol_contrast = np.array(WN_dispersion_perm_tval_hormonal_contrast['estradiol_t_val'])
WN_dispersion_perm_tval_progesterone_contrast = np.array(WN_dispersion_perm_tval_hormonal_contrast['progesterone_t_val'])

BN_dispersion_perm_tval_estradiol_contrast = np.array(BN_dispersion_perm_tval_hormonal_contrast['estradiol_t_val'])
BN_dispersion_perm_tval_progesterone_contrast = np.array(BN_dispersion_perm_tval_hormonal_contrast['progesterone_t_val'])


# contains the null distribution of t values for the hormonal contrasts (estradiol and progesterone separately) on the Within Network dispersion per network (7)
WN_dispersion_perm_tval_estradiol_contrast_null_distr = {'visual': WN_dispersion_perm_tval_estradiol_contrast.T[0],
							 'sensory motor': WN_dispersion_perm_tval_estradiol_contrast.T[1], 
                                                   	 'dorsal attention': WN_dispersion_perm_tval_estradiol_contrast.T[2], 
                                                 	 'ventral attention': WN_dispersion_perm_tval_estradiol_contrast.T[3], 
                                                   	 'limbic': WN_dispersion_perm_tval_estradiol_contrast.T[4], 
                                                 	 'fronto parietal': WN_dispersion_perm_tval_estradiol_contrast.T[5], 
                                                  	 'DMN': WN_dispersion_perm_tval_estradiol_contrast.T[6]}

WN_dispersion_perm_tval_progesterone_contrast_null_distr = {'visual': WN_dispersion_perm_tval_progesterone_contrast.T[0],
							    'sensory motor': WN_dispersion_perm_tval_progesterone_contrast.T[1], 
                                                   	    'dorsal attention': WN_dispersion_perm_tval_progesterone_contrast.T[2], 
                                                 	    'ventral attention': WN_dispersion_perm_tval_progesterone_contrast.T[3], 
                                                      	    'limbic': WN_dispersion_perm_tval_progesterone_contrast.T[4], 
                                                 	    'fronto parietal': WN_dispersion_perm_tval_progesterone_contrast.T[5], 
                                                  	    'DMN': WN_dispersion_perm_tval_progesterone_contrast.T[6]}


# contains the null distribution of t values for the sex contrast on the Between Network pairwise comparisons (21)
# column names is the pairwise comparison labels (turned into strings)
BN_dispersion_perm_tval_estradiol_contrast_null_distr = pd.DataFrame(BN_dispersion_perm_tval_estradiol_contrast, columns = [str(pair) for pair in pairwise_comparison_order])  
BN_dispersion_perm_tval_progesterone_contrast_null_distr = pd.DataFrame(BN_dispersion_perm_tval_progesterone_contrast, columns = [str(pair) for pair in pairwise_comparison_order])  



# export 
pd.DataFrame(WN_dispersion_perm_tval_estradiol_contrast_null_distr).to_csv(resdir+'WN_dispersion_perm_tval_estradiol_contrast_null_distr_Me.csv', header = True, index = False)
pd.DataFrame(WN_dispersion_perm_tval_progesterone_contrast_null_distr).to_csv(resdir+'WN_dispersion_perm_tval_progesterone_contrast_null_distr_Me.csv', header = True, index = False)

BN_dispersion_perm_tval_estradiol_contrast_null_distr.to_csv(resdir+'BN_dispersion_perm_tval_estradiol_contrast_null_distr_Me.csv', header = True, index = False)
BN_dispersion_perm_tval_progesterone_contrast_null_distr.to_csv(resdir+'BN_dispersion_perm_tval_progesterone_contrast_null_distr_Me.csv', header = True, index = False)


    
    
