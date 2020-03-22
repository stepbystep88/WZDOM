% Show results  of all optimization methods
%
% show_all_results_FINAL1.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axislimits_f=[0 2000  1e-16 1e4]; axislimits_snr=[0 2000 -2 17];
load reports_ComprSens_Sparse_734127.4818_smooth_1e-3.mat;ind=[1:6]; name='ComprSens_Sparse_FISTA';plot_all;   

%%%%%%%   LATEST RESULTS FOR THE PAPER  %%%%%%%%%%%%%%

axislimits_f=[0 990  1e-7 1e4]; axislimits_snr=[0 1000 -2 12];
load reports_ComprSens_Sparse_733988.4394_smooth_1e-3.mat; ind=[4:9]; name='ComprSens_Sparse1';plot_all;            
load reports_ComprSens_Sparse_733988.4394_smooth_1e-3.mat; ind=[1 3 6 10 11 12];  name='ComprSens_Sparse2'  ;plot_all;  

axislimits_f=[0 1000  1e-5 1e3]; axislimits_snr=[0 1000 10 45];
load reports_ComprSens_Sparse_733988.4228_concave.mat; ind=[4,5,7,8];name='ComprSens_Sparse_concave' ;plot_all;   

axislimits_f=[0 55  1e-3 1e3]; axislimits_snr=[0 60 0 14];
load reports_Tomogr_Phantom_734108.7581_smooth_1e-3.mat
%load reports_Tomogr_Phantom_733988.4817_smooth_1e-3.mat; 
ind=4:9;    name='Tomogr_Phantom1'   ;plot_all;    
ind=[1 3 6 10 11 12];   name='Tomogr_Phantom2'   ;plot_all;  

axislimits_f=[0 90  1e-3 1e3]; axislimits_snr=[0 100 2 12];
load reports_Deblur_Phantom_733989.3502_smooth_1e-3.mat; ind=4:9;  name='Deblur_Phantom1'   ;plot_all;  
load reports_Deblur_Phantom_733989.3502_smooth_1e-3.mat; ind=[1 3 6 10 11];  name='Deblur_Phantom2'  ;plot_all;  


