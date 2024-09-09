%% PLOT CONFUSION MATRIX

%% READ CSV confusion matrix
clc 
clear all
close all

start='/Users/matteotorzoni/Desktop/Corigliano'
dim_char=20;
%Prompt user for filename
[fname, pname] = uigetfile('*.csv','Choose input file',start);  
%Create fully-formed filename as a string
filename = fullfile(pname, fname);
%Read in the data, skipping the first row 
data = csvread(filename,0,0);
plotConfMat(data',{'\Omega_0','\Omega_1','\Omega_2','\Omega_3','\Omega_4','\Omega_5','\Omega_6'},dim_char)