# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 15:48:23 2022

@author: Matteo
"""
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
import numpy as np
import pandas as pd
import scipy.stats as st


# Number of monitored channels
n_channels = 8
# Which dof monitored
which_channels = [1,2,3,4,5,6,7,8] 
# Number of classes in the dataset
n_class = 8
# Specify if you are using accelerations =1 or displacements =0
accelerations = 0
# Specify the length of the data time series 
seq_len = 200
# Specify the amount of noise you want to add to the data
addedd_SNR = 100

# Path containing the data to be loaded
path_data = '.\\istantest_7_1'
# Path containing the models to be loaded
path_classifier_total = '.\\Classifier_total_7_1'
path_regressors = ['.\\Level_Regressor_7_' + str(i1+1) for i1 in range(n_class-1)]
    
# Load the statistics for each channel, computed for the different models
signals_means_classifier = np.load(path_classifier_total+'\\Corrupted_signals_means.npy')
signals_stds_classifier = np.load(path_classifier_total+'\\Corrupted_signals_stds.npy')
signals_means_regressors = [np.load(path_regressors[i1]+'\\Corrupted_signals_means.npy') for i1 in range(n_class-1)]
signals_stds_regressors = [np.load(path_regressors[i1]+'\\Corrupted_signals_stds.npy') for i1 in range(n_class-1)]
# statistics to normalize the data for the ML models
statistics = [signals_means_classifier,signals_stds_classifier,signals_means_regressors,signals_stds_regressors]
# Load the statistics used to normalize the damage levels for the different models
level_means = [np.load(path_regressors[i1]+'\\Damage_level_mean.npy') for i1 in range(n_class-1)]
level_stds  = [np.load(path_regressors[i1]+'\\Damage_level_std.npy') for i1 in range(n_class-1)]
statistics_level = [level_means, level_stds]
    
def rms(vect):
    """
    Compute the root mean squared value of vect
    """
    return np.sqrt(np.mean(np.square(vect)))

# Load the labels relative to the observations
label_path_class = path_data + '\\Damage_class.csv'                                
labels_class     = np.genfromtxt(label_path_class)
labels_class     = labels_class.astype('int')
label_path_level = path_data + '\\Damage_level.csv'                                
labels_level     = np.genfromtxt(label_path_level)
label_path_amplitude = path_data + '\\Amplitude.csv'                                
labels_amplitude     = np.genfromtxt(label_path_amplitude)
label_path_frequency = path_data + '\\Frequency.csv'                                
labels_frequency     = np.genfromtxt(label_path_frequency)

N_ist = len(labels_level)  
if accelerations == 1:
    path_rec  = path_data + '\\U2_concat_'
elif accelerations == 0:
    path_rec  = path_data + '\\U_concat_'
# Create the dataset output structure for clean and corrupted data
X       = np.zeros((n_channels, seq_len, N_ist))
X_noise = np.zeros((n_channels, seq_len, N_ist))
# Create the dataset output structure for different models
X_classifier = np.zeros((N_ist, seq_len, n_channels))
X_regressors = np.array([np.zeros((N_ist, seq_len, n_channels)) for _ in range(n_class-1)])
# Fill the dataset output structures 
for i1 in range(n_channels):
    path_data_var  = path_rec + str(which_channels[i1]) + '.csv'
    X_singledof = pd.read_csv(path_data_var, header=None).to_numpy()[:,0]
    for i2 in range(N_ist):
        # Reshape the signals as a dataset of multivariate time series
        X[i1,:,i2]= X_singledof[1 + i2*(1+seq_len) : (i2+1)*(1+seq_len)]   
        # Corrupt the signals with noise to obtain "addedd_SNR"
        rms_signal = rms(X[i1,:,i2]);
        dev_std    = rms_signal / np.sqrt(addedd_SNR);
        # Define the Gaussian distribution to sample the noise
        sample = st.norm(0, dev_std)
        X_noise[i1,:,i2] = X[i1,:,i2] + sample.rvs(size=seq_len)
    # Normalize the data accoring to the training of different models
    X_classifier[:,:,i1] = (np.transpose(X_noise[i1,:,:]) - statistics[0][i1])/statistics[1][i1]
    for i3 in range(n_class-1):
        X_regressors[i3,:,:,i1] = (np.transpose(X_noise[i1,:,:]) - statistics[2][i3][i1])/statistics[3][i3][i1]

np.save('.\\recordings_classifier', X_classifier)
np.save('.\\recordings_regressors_1of4', X_regressors[:,0:1000])
np.save('.\\recordings_regressors_2of4', X_regressors[:,1000:2000])
np.save('.\\recordings_regressors_3of4', X_regressors[:,2000:3000])
np.save('.\\recordings_regressors_4of4', X_regressors[:,3000:4000])
np.save('.\\labels_class', labels_class)
np.save('.\\labels_level', labels_level)
np.save('.\\statistics_level', statistics_level)
np.save('.\\labels_amplitude', labels_amplitude)
np.save('.\\labels_frequency', labels_frequency)