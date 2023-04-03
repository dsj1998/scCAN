% This is a demo showing how to running CAN-Impute using the uso
clear;clc
 addpath(genpath('./'))

%----------------------load mat processed data---------------------

load('data.mat')
data = Brain';
data = log10(1+data);
%----------------------------parameters-----------------------------
K = 50;
r = 20;
lambda1 = 0.2;
lambda2 = 0.4;
paras = [K,r,lambda1,lambda2];

%----------------------------find missing positions-----------------------------
cn = 3;
P = findMP(data,cn,threshold); % cn is the number of clusters 

%--------------------------Initialization---------------------------
[W0,H0,S0] = Initialize(data, paras,P, 100); % data;K;max_neighbors;lambda1;MaxIter
Init = [{data},{W0},{H0},{S0},P];
%----------------------------  run-----------------------------------


[data_full] = adaptive_nmf(Init, paras);% return the imputated data matrix;D is the imputed data
X=data_full';
Y= max(10.^X-1,0);

fileDir_csv= 'result/';
filename_mat='data.mat';
savePath_csv = strcat(fileDir_csv,filename_csv); 
csvwrite(savePath_csv, Y);


