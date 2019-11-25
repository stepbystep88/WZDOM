function [wavelet, GInvParam] = bsExtractWavelet(GInvParam, timeLine, wellLogs, type)

    if ~exist('type', 'var')
        type = 'ricker';
    end
    
    [GInvParam, postSeisData, horizon] = bsPrepareForExtractingWavelet(...
        GInvParam, ...
        timeLine{GInvParam.usedTimeLineId}, ...
        wellLogs, ...
        type);

    wellNum = length(wellLogs);      % the number of wells
    % get synthetic data
    scaleFactors = zeros(wellNum, 1);
    similarities = zeros(wellNum, 1);
    
    for i = 1 : wellNum
        synData = bsGetSynthetic(GInvParam, wellLogs{i}, horizon(i));

%         scaleFactors(i) = bsCalcScaleFactor(postSeisData(:, i), synData);
        scaleFactors(i) = bsComputeGain(postSeisData(:, i), synData);
        correlation = corrcoef(synData, postSeisData(:, i));
        similarities(i) = correlation(1, 2);
    end
    
    [~, index] = bsMaxK(similarities, ceil(0.3*wellNum));
    meanScaleFactor = mean(scaleFactors(index));

    wavelet = GInvParam.wavelet * meanScaleFactor;
    
%     for i = 1 : wellNum
%         synData = bsGetSynthetic(GInvParam, wellLogs{i}, horizon(i));
%         
%         figure;
%         plot(1:length(synData), synData*meanScaleFactor, 'r', 'linewidth', 2); hold on;
%         plot(1:length(synData), postSeisData(:, i), 'k', 'linewidth', 2);
%         legend('Synthetic', 'Real seismic data');
%         
%     end
    
end


function synseis = bsGetSynthetic(GInvParam, wellInfo, horizon)
    sampNum = GInvParam.upNum + GInvParam.downNum; 
    
    if isfield(GInvParam.indexInWellData, 'Ip')
        Ip = bsExtractWellDataByHorizon(...
                    wellInfo.wellLog, ...
                    horizon, ...
                    GInvParam.indexInWellData.Ip, ...
                    GInvParam.indexInWellData.time, ...
                    GInvParam.upNum, ...
                    GInvParam.downNum, ...
                    1);
    else
        data = bsExtractWellDataByHorizon(...
                    wellInfo.wellLog, ...
                    horizon, ...
                    [GInvParam.indexInWellData.vp, GInvParam.indexInWellData.rho], ...
                    GInvParam.indexInWellData.time, ...
                    GInvParam.upNum, ...
                    GInvParam.downNum, ...
                    1);
        Ip = data(:, 1) .* data(:, 2);
    end
    
    G = bsPostGenGMatrix(GInvParam.wavelet, sampNum);
    
    synseis = G * log(Ip);
end

function scaleFactor = bsCalcScaleFactor(v1, v2)
    scaleFactor = (v2' * v1) / (v2' * v2);
end



    
