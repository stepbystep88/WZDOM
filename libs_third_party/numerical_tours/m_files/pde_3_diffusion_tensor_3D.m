% close all;

% n1 = 60;
% n2 = 60;
% n3 = 60;
n = 100;

% f = randn(n3, n2, n1);
segyInfo = bsCreateGSegyInfo();
segyInfo.t0 = 0;
segyInfo.isPosZero = 0;
segyInfo.inlineId = 9;
segyInfo.crosslineId = 21;
basePath = 'E:/HRS projects/GEO_POST/';

fileName = sprintf('%s/seismic/seismic.sgy', basePath);
shiftfileName = sprintf('%s/seismic/phase_shift_90.sgy', basePath);

[trData, inIds, crossIds] = bsReadAllTraces(fileName, segyInfo);

test = bsAddNoise(trData(1001:1500, :), 1, 10, [], [], [], []);

volumeData = bsReshapeDataAs3D(test, 142, 110);

% inputData = volumeData(1:n, 1:n, 1:n);
inputData = volumeData;
% inputData = permute(inputData, [2 3 1]);

% name = 'hibiscus';
% f0 = load_image(name,n);
% f0 = rescale( sum(f0,3) );
% inputData = repmat(f0, [1, 1, n]);

% inputData = permute(inputData, [3 1 2]);


volumeData = inputData;

f = inputData;



GST_options = bsCreateGSTParam(3);
GST_options.sigma = 4;
GST_options.iterNum = 30;
GST_options.lambda = 1e-4;
GST_options.show_mid_results = false;
GST_options.sub = 4;
GST_options.filter_fcn = 1;
GST_options.ttv = 0.0001;
GST_options.tau = 0.2;

f2 = bsSmoothByGST3D(f, [],  GST_options);


figure;
set(gcf, 'position', [154         199        1457         779]);
subplot(1, 3, 1);
 range = [prctile(inputData(:), 10), prctile(inputData(:), 90)];
bsShow3DVolume(inputData, 1, range, 1, 1, 1, [], 'colormap', bsGetColormap('seismic'));

subplot(1, 3, 2);
% outputData = permute(f2, [3 1 2]);
outputData = f2;
bsShow3DVolume(outputData, 1, range, 1, 1, 1, [], 'colormap', bsGetColormap('seismic'));

subplot(1, 3, 3);
% outputData = permute(f2, [3 1 2]);
outputData = imgaussfilt3(f2, 0.1);
bsShow3DVolume(outputData, 1, range, 1, 1, 1, [], 'colormap', bsGetColormap('seismic'));

