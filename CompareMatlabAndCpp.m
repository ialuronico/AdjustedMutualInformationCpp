close all;
clear all;
clc;

% Script to compare the Adjusted Mutual Information computed in 
% Matlab and the one computed in Cpp
% 
% I am computing AMI_max (normalized by the maximum entropy)

% Number of clusters 
k = 50;
% Number of points 
N = 100;
% Generate one clustering
A = randi(k,1,N);
% Generate the other clustering
B = randi(k,1,N);

tic;

valAMImatlab = ami(A,B);

tMatlab = toc;

tic; 

valAMIcpp = AMIcpp(A,B);

tCpp = toc;

disp(['AMImatlab = ' num2str(valAMImatlab) ' in ' num2str(tMatlab) ' seconds']);
disp(['AMIcpp = ' num2str(valAMIcpp) ' in ' num2str(tCpp) ' seconds']);
disp(['Speedup = ' num2str(tMatlab/tCpp)]);