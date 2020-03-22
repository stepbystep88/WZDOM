% Show results  of all optimization methods
%
% show_all_results_FINAL1.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load reports_Matrix_Sparse_734131.4333_smooth_1e-3_weight_1e-6.mat
%load reports_Matrix_Sparse_734131.4424_smooth_1e-3.mat
%load reports_Matrix_Sparse_geogaussmatrix_734131.5007_smooth_1e-3_weight_1e-6.mat;
%load reports_Matrix_Sparse_tonytaperwavreal_734131.4770_smooth_1e-3_weight_1e-6.mat

%load reports_Tomogr_Phantom32_734134.6232_smooth_1e-3.mat
%load reports_Tomogr_Phantom64_734134.6264_smooth_1e-3.mat
%load reports_Tomogr_Phantom128_734134.6511_smooth_1e-3.mat
%load reports_Tomogr_Phantom256_734134.6749_smooth_1e-3.mat

%load reports_Matrixgeogaussmatrix_Sparse128_734135.7210_smooth_1e-3_weight_1e-6.mat
%load reports_Matrixgeogaussmatrix_Sparse128_734135.7587_smooth_1e-3.mat
%load reports_Matrixbadgaussmatrix_Sparse128_734135.6777_smooth_1e-3_weight_1e-6.mat
%load reports_Matrixbadgaussmatrix_Sparse128_734135.6777_smooth_1e-3.mat
%load reports_Matrixgaussmatrix_Sparse128_734135.7963_smooth_1e-3_weight_1e-6.mat
%load reports_Matrixgaussmatrix_Sparse128_734135.7798_smooth_1e-3.mat



axislimits_f=[0 600  5e-17 1e3]; axislimits_snr=[0 600 -20 56];tmax=50;
load reports_Matrixgaussmatrix_Sparse128_734135.7798_smooth_1e-3.mat
i_run=numel(reports); ind=1:i_run; name='Gaussmatrix_lam_1e-3';plot_all;   


axislimits_f=[0 1500  1e-8 1e0]; axislimits_snr=[0 1500 -10 10];tmax=120;
load reports_Matrixgaussmatrix_Sparse128_734135.7963_smooth_1e-3_weight_1e-6.mat
i_run=numel(reports); ind=1:i_run; name='Gaussmatrix_lam_1e-6';plot_all;   



axislimits_f=[0 300  1e-9 1e3]; axislimits_snr=[0 300 -3 2]; tmax=30;
load reports_Matrixbadgaussmatrix_Sparse128_734135.6777_smooth_1e-3.mat
i_run=numel(reports); ind=1:i_run; name='BadGaussmatrix_lam_1e-3';plot_all;   


axislimits_f=[0 300  3e-8 1e1]; axislimits_snr=[0 300 -3 2]; tmax=30;
load reports_Matrixbadgaussmatrix_Sparse128_734135.6853_smooth_1e-3_weight_1e-6.mat  
i_run=numel(reports); ind=1:i_run; name='BadGaussmatrix_lam_1e-6';plot_all;   


axislimits_f=[0 600  1e-12 1e3]; axislimits_snr=[0 600 -3 -1.6]; tmax=60;
load reports_Matrixgeogaussmatrix_Sparse128_734135.7587_smooth_1e-3.mat
i_run=numel(reports); ind=1:i_run; name='GeoGaussmatrix_lam_1e-3';plot_all;   


axislimits_f=[0 3000  1e-5 1e1]; axislimits_snr=[0 3000 -3 1]; tmax=250;
load reports_Matrixgeogaussmatrix_Sparse128_734135.7210_smooth_1e-3_weight_1e-6.mat
i_run=numel(reports); ind=1:i_run; name='GeoGaussmatrix_lam_1e-6';plot_all;   






axislimits_f=[0 1000  1e-8 1e2]; axislimits_snr=[0 1000 -5 16]; tmax=12;
load reports_ComprSens_Sparse_all_128_734135.5321_smooth_1e-3.mat
i_run=numel(reports); %ind=1:i_run; name='ComprSens_Sparse';plot_all;   
ind=[4:9]; name='ComprSens_Sparse1';plot_all;  
ind=[1 6 10 11 12];  name='ComprSens_Sparse2'  ;plot_all; 

axislimits_f=[0 120  1e-4 1e2]; axislimits_snr=[0 120 0 15]; tmax=60;
load reports_Deblur_Phantom128_all_734135.5212_smooth_1e-3.mat
%i_run=numel(reports); ind=1:i_run; name='Deblur_Phantom128new';plot_all;   
ind=[4:9]; name='Deblur_Phantom128-1';plot_all;  
ind=[1 6 10 11 12];  name='Deblur_Phantom128-2'  ;plot_all; 

axislimits_f=[0 120  1e-4 1e2]; axislimits_snr=[0 120 0 15]; tmax=60;
load reports_Tomogr_Phantom128_all_734135.3694_smooth_1e-3.mat
%i_run=numel(reports); ind=1:i_run; name='Tomogr_Phantom128new';plot_all;   
ind=[4:9]; name='Tomogr_Phantom128-1';plot_all;  
ind=[1 6 10 11 12];   name='Tomogr_Phantom128-2';plot_all; 




axislimits_f=[0 120  1e-4 1e2]; axislimits_snr=[0 120 0 15]; tmax=10;
load reports_Tomogr_Phantom32_734134.6232_smooth_1e-3.mat
i_run=numel(reports); ind=1:i_run; name='Tomogr_Phantom32';plot_all;   

axislimits_f=[0 120  1e-4 1e2]; axislimits_snr=[0 120 0 15]; tmax=60;
load reports_Tomogr_Phantom128_734134.6511_smooth_1e-3.mat
i_run=numel(reports); ind=1:i_run; name='Tomogr_Phantom128';plot_all;   

axislimits_f=[0 120  1e-4 1e2]; axislimits_snr=[0 120 0 15]; tmax=500;
load reports_Tomogr_Phantom256_734134.6749_smooth_1e-3.mat; 
i_run=numel(reports); ind=1:i_run; name='Tomogr_Phantom256';plot_all;   

 

% axislimits_f=[0 2000  1e-16 1e2]; axislimits_snr=[0 2000 -2 17]; tmax=30;
% load reports_ComprSens_Sparse_734127.4818_smooth_1e-3.mat;ind=[1:6]; 
% name='ComprSens_Sparse_FISTA';plot_all;   


% axislimits_f=[0 120  1e-4 1e2]; axislimits_snr=[0 120 0 15]; tmax=30;
% load reports_Matrix_Sparse_734131.4424_smooth_1e-3.mat
% i_run=numel(reports); ind=1:i_run; name='Matrix_Sparse';plot_all;   



% %%%%%%%   AUGUST 2009 RESULTS FOR THE PAPER  %%%%%%%%%%%%%%
% 
% axislimits_f=[0 990  1e-7 1e2]; axislimits_snr=[0 1000 -2 12]; tmax=30;
% load reports_ComprSens_Sparse_733988.4394_smooth_1e-3.mat; ind=[4:9]; name='ComprSens_Sparse1';plot_all;            
% load reports_ComprSens_Sparse_733988.4394_smooth_1e-3.mat; ind=[1 3 6 10 11];  name='ComprSens_Sparse2'  ;plot_all;  
% 
 axislimits_f=[0 1000  1e-5 1e2]; axislimits_snr=[0 1000 0 45]; tmax=12;
 load reports_ComprSens_Sparse_733988.4228_concave.mat; ind=[4,5,7,8];name='ComprSens_Sparse_concave' ;plot_all;   
% 
% axislimits_f=[0 55  1e-3 1e2]; axislimits_snr=[0 60 0 14]; tmax=30;
% load reports_Tomogr_Phantom_734108.7581_smooth_1e-3.mat
% %load reports_Tomogr_Phantom_733988.4817_smooth_1e-3.mat; 
% ind=4:9;    name='Tomogr_Phantom1'   ;plot_all;    
% ind=[1 3 6 10 11 12];   name='Tomogr_Phantom2'   ;plot_all;  
% 
% axislimits_f=[0 90  1e-3 1e2]; axislimits_snr=[0 100 2 12]; tmax=30;
% load reports_Deblur_Phantom_733989.3502_smooth_1e-3.mat; ind=4:9;  name='Deblur_Phantom1'   ;plot_all;  
% load reports_Deblur_Phantom_733989.3502_smooth_1e-3.mat; ind=[1 3 6 10 11];  name='Deblur_Phantom2'  ;plot_all;  
% %load reports_Deblur_Phantom_733989.3502_smooth_1e-3.mat; ind=[1:11];  name='Deblur_Phantom_all'  ;plot_all;  
% 
% 
