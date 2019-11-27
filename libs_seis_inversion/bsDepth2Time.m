function [GInvParam, outputWelllogs, wavelet] = bsDepth2Time(GInvParam, timeLine, inputWelllogs, waveletType)
%% get time information of welllog data
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
     
    if ~exist('waveletType', 'var')
        waveletType = 'waveletType';
    end
    
    sampNum = GInvParam.upNum + GInvParam.downNum;
    [GInvParam, postSeisData, horizon] = bsPrepareForExtractingWavelet(...
       GInvParam, ...
       timeLine{GInvParam.usedTimeLineId}, ...
       inputWelllogs, ...
       waveletType);
   
    [~, wellNum] = size(postSeisData);                    
    
    % generate wavelet convolution matrix
    G = bsPostGenGMatrix(GInvParam.wavelet, sampNum);
                  
    scaleFactors = zeros(wellNum, 1);
    similarities = zeros(wellNum, 1);
    synSeisData = zeros(sampNum-1, wellNum);
    
    expandNum = GInvParam.depth2time.expandNum;
    indexInWell = GInvParam.indexInWellData;
    outputWelllogs = inputWelllogs;
    searchOffsetNum = GInvParam.depth2time.searchOffsetNum;
    
    for i = 1 : wellNum
        wellInfo = inputWelllogs{i};
        
        wellData = wellInfo.wellLog;                       % welllog data including vp, vs, rho
        layerNum = size(wellData, 1);
        
        % expand welllog data
        if isfield(indexInWell, 'Ip')
            wellData = bsExpandPostWellData(wellData, expandNum);
            dist = abs(wellData(:, indexInWell.time) - horizon);
        else
            wellData = bsExpandPreWellData(wellData, ...
                expandNum, ...
                GInvParam.dt, ...
                indexInWell.depth, ...
                indexInWell.vp);
            targetDep = wellInfo.targetDepth;               % the depth at target region
            dist = abs(wellData(:, indexInWell.depth) - targetDep);  
        end
        [~, targetIndex] = min(dist);
        


        %% find the segment of welllog data that matches the seismic data best
        maxCorrelation = -1;
        
        for j = -searchOffsetNum : searchOffsetNum
            index = targetIndex + j;
            upIndex = index - GInvParam.upNum;
            downIndex = index + GInvParam.downNum - 1;
            
            if isfield(indexInWell, 'Ip')
                Ip = wellData(upIndex : downIndex, indexInWell.Ip);
            else
                vp = wellData(upIndex : downIndex, indexInWell.vp);
                rho = wellData(upIndex : downIndex, indexInWell.rho);
                Ip =  vp .* rho;
            end
            
            synthSeis = G * log(Ip);                                   

            % only compare the valid part of welllog data
            if( upIndex <= expandNum)
                num = expandNum - upIndex + 2;
                correlation = corrcoef(synthSeis(num:end), postSeisData(num:end, i)) ;
            elseif ( downIndex > layerNum + expandNum )
                num = layerNum + expandNum - upIndex + 1;
                correlation = corrcoef(synthSeis(1:num), postSeisData(1:num, i)) ;
            else
                correlation = corrcoef(synthSeis, postSeisData(:, i));
            end
            
            if correlation(1,2) > maxCorrelation
                maxCorrelation = correlation(1,2);
                bestIndex = index;
                synSeisData(:, i) = synthSeis;
            end
        end
            
        upIndex = bestIndex - GInvParam.upNum;
        downIndex = bestIndex + GInvParam.downNum - 1;
            
        bestWelllog = wellData(upIndex : downIndex, :);
        
        % add time information
        timeData = (0:(sampNum-1))'*GInvParam.dt;
        timeData = timeData - timeData(GInvParam.upNum) + horizon(i);
        
        outputWelllogs{i}.wellLog = [bestWelllog, timeData];
        outputWelllogs{i}.t0 = timeData(1);
        
        similarities(i) = maxCorrelation;
        scaleFactors(i) = bsComputeGain(postSeisData(:, i), synSeisData(:, i));
        
    end
    
    % re-scale wavelet
    [~, index] = bsMaxK(similarities, ceil(0.3*wellNum));
    meanScaleFactor = mean(scaleFactors(index));

    wavelet = GInvParam.wavelet * meanScaleFactor;
    GInvParam.wavelet = wavelet;
    GInvParam.indexInWellData.time = size(inputWelllogs{1}.wellLog, 2) + 1;
    
    if GInvParam.depth2time.isShowCompare
        
        for i = 1 : min(wellNum, GInvParam.depth2time.showCompareNum)
            figure;
            plot(1:sampNum-1, postSeisData(:, i), 'k', 'linewidth', 2); hold on;
            plot(1:sampNum-1, synSeisData(:, i)*meanScaleFactor, 'r', 'linewidth', 2);
            legend('real seismic data', 'synthetic seismic data');
            bsSetDefaultPlotSet(bsGetDefaultPlotSet());
        end
        
    end

end

function wellData = bsExpandPostWellData(wellData, expandNum)
    wellData = [...
        repmat(wellData(1, :), expandNum, 1); ...
        wellData; ...
        repmat(wellData(end, :), expandNum, 1)];
end

function wellData = bsExpandPreWellData(wellData, expandNum, dt, depthIndex, vpIndex)
    [layerNum, propertyNum] = size(wellData);          

    otherIndex = setdiff(1:propertyNum, depthIndex);
    
    dz = dt * wellData(1, vpIndex) * 0.001 * 0.5;
    tmp1 = zeros(expandNum, propertyNum);
    for i = 1 : expandNum
        tmp1(expandNum-i+1, 1) = wellData(1, depthIndex) - i * dz;
        tmp1(expandNum-i+1, 2 : propertyNum) = wellData(1, otherIndex);
    end
    
    dz = dt * wellData(layerNum, vpIndex) * 0.001 * 0.5;    
    tmp2 = zeros(expandNum, propertyNum);
    for i = 1 : expandNum
        tmp2(i, 1) = wellData(layerNum, depthIndex) + i * dz;
        tmp2(i, 2 : propertyNum) = wellData(layerNum, otherIndex);
    end

    wellData = [tmp1; wellData; tmp2];
end