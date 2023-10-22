disp('--------------------start--------------------')
tic
%----------------------load mat processed data---------------------

data_sc = csvread('/data/liver.csv',1,1); 
data_raw = data_sc';
data = log10(1+data_raw);       

cn = 7;                
K = 7;
r = 15;
lambda1 = 2;
lambda2 = 4;

threshold = 0.1;
paras = [K,r,lambda1,lambda2];
                   
                    
%----------------------------find missing positions-----------------------------
 P = findMP(data,cn,threshold);% cn is the number of clusters
                    
                    
 %--------------------------Initialization---------------------------
 [W0,H0,S0] = Initialize(data, paras,P,100); % data;K;max_neighbors;lambda1;MaxIter
                     
 Init = [{data},{W0},{H0},{S0},P];
                    
%----------------------------  run-----------------------------------
[data_full,~,~,~] = adaptive_nmf(Init, paras);

                   
X=data_full';
Y= max(10.^X-1,0);
                   
fileDir_csv = 'result/'; 
filename_csv="liver_result.csv";
savePath_csv = strcat(fileDir_csv,filename_csv); 
csvwrite(savePath_csv, Y);

disp('--------------------end--------------------')
addpath(genpath('./'))
