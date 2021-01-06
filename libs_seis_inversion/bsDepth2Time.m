function [GInvParam, outputWelllogs, wavelet] ...
    = bsDepth2Time(GInvParam, timeLine, inputWelllogs, waveletType)
%% get time information of welllog data
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
     
    if ~exist('waveletType', 'var')
        waveletType = 'waveletType';
    end
    
    bsCheckIsValidWellData(GInvParam, inputWelllogs);
    
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
    saveOffsetNum = GInvParam.depth2time.saveOffsetNum;
    
    for i = 1 : wellNum
        wellInfo = inputWelllogs{i};
        
        wellData = wellInfo.wellLog;                       % welllog data including vp, vs, rho
        layerNum = size(wellData, 1);
        
        % expand welllog data
        if isfield(indexInWell, 'ip') && indexInWell.ip > 0
            wellData = bsExpandPostWellData(wellData, expandNum);
            dist = abs(wellData(:, indexInWell.time) - horizon(i));
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
            upIndex = index - GInvParam.upNum + 1;
            downIndex = index + GInvParam.downNum;
            
            try
                if isfield(indexInWell, 'ip') && indexInWell.ip > 0
                    ip = wellData(upIndex : downIndex, indexInWell.ip);
                else
                    vp = wellData(upIndex : downIndex, indexInWell.vp);
                    rho = wellData(upIndex : downIndex, indexInWell.rho);
                    ip =  vp .* rho;
                end
            catch
                error('Data of well %s is not enough.', inputWelllogs{i}.name);
            end
            
            
            synthSeis = G * log(ip);  
            % only compare the valid part of welllog data
%             if( upIndex <= expandNum)
%                 num = expandNum - upIndex + 2;
%                 correlation = corrcoef(synthSeis(num:end), postSeisData(num:end, i)) ;
%             elseif ( downIndex > layerNum + expandNum )
%                 num = layerNum + expandNum - upIndex + 1;
%                 correlation = corrcoef(synthSeis(1:num), postSeisData(1:num, i)) ;
%             else
            correlation = corrcoef(synthSeis, postSeisData(:, i));
%             end
            
            try
                if correlation(1,2) > maxCorrelation
                    maxCorrelation = correlation(1,2);
                    bestIndex = index;
                    synSeisData(:, i) = synthSeis;
                end
            catch
                error('Data of well %s is not enough.', inputWelllogs{i}.name);
            end
        end
            
        
        targetDep = wellData(bestIndex, indexInWell.depth);

        % 长度不够填充方式， 注释后长度不够不填充
        wellData = wellData(expandNum+1:end-expandNum, :);
        bestIndex = bestIndex - expandNum;
        
        upIndex = bestIndex - GInvParam.upNum - saveOffsetNum;
        downIndex = bestIndex + GInvParam.downNum - 1 + saveOffsetNum;
        
        if upIndex < 1 
            upIndex = 1;
        end
        
        if downIndex > size(wellData, 1)
            downIndex = size(wellData, 1);
        end
        
        bestWelllog = wellData(upIndex : downIndex, :);            
        dist = abs(bestWelllog(:, indexInWell.depth) - targetDep);  
        [~, targetIndex] = min(dist);
        
        % add time information
        timeData = (0:(size(bestWelllog, 1)-1))'*GInvParam.dt;
        timeData = timeData - timeData(targetIndex) + horizon(i);
        
        outputWelllogs{i}.wellLog = [bestWelllog, timeData];
        outputWelllogs{i}.t0 = timeData(1);
        outputWelllogs{i}.targetDepth = targetDep;
        
        similarities(i) = maxCorrelation;
        scaleFactors(i) = bsComputeGain(postSeisData(:, i), synSeisData(:, i));
        
    end
    
    % re-scale wavelet
    [~, index] = bsMaxK(similarities, ceil(0.8*wellNum));
    meanScaleFactor = mean(scaleFactors(index));
    bestIndex = index(1);
    bestSeisData = postSeisData(:, bestIndex);
    bestScale = scaleFactors(bestIndex);
    
    
%     meanScaleFactor = 1;
    wavelet = GInvParam.wavelet * meanScaleFactor;
    GInvParam.wavelet = wavelet;
    
    GInvParam.indexInWellData.time = size(inputWelllogs{1}.wellLog, 2) + 1;
    
    if GInvParam.isScale 
        GInvParam.bestPostSeisData = bestSeisData / bestScale * meanScaleFactor;
    else
        GInvParam.bestPostSeisData = bestSeisData;
    end
    
    
    if GInvParam.depth2time.isShowCompare
        
        figure;
        bsSetPosition(0.78, 0.56);
        
        nSubFigure = min(wellNum, GInvParam.depth2time.showCompareNum);
        for i = 1 : nSubFigure
            
            switch nSubFigure
                case {1, 2, 3, 4, 5}
                    subplot(1, nSubFigure, i);
                otherwise
                    subplot(2, ceil(nSubFigure/2), i);
            end
            
            plot(postSeisData(:, i), 1:sampNum-1,  'k', 'linewidth', 2); hold on;
            plot(synSeisData(:, i)*meanScaleFactor, 1:sampNum-1,  'r', 'linewidth', 2);
            set(gca, 'ydir', 'reverse');
            
            if i == 1
                legend('real seismic data', 'synthetic seismic data');
                ylabel('Sample number');
            end
            
            if isfield(inputWelllogs{i}, 'name')
                title(inputWelllogs{i}.name);
            elseif isfield(inputWelllogs{i}, 'wellName')
                title(inputWelllogs{i}.wellName);
            end
            
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