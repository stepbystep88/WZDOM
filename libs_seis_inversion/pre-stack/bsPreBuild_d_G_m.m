function [d, G, initX, lsdCoef, angleData] = bsPreBuild_d_G_m(GPreInvParam, inline, crossline, startTime, initLog, G)
    sampNum = GPreInvParam.upNum + GPreInvParam.downNum; 
    
    % load prestack seismic data
    preDataInfo = GPreInvParam.preSeisData;
    switch lower(preDataInfo.mode)
        case 'angle_separate_files'
            separates = preDataInfo.separates;
            angleSeisData = bsReadMultiSegyFiles(separates, inline, crossline, ...
                startTime, sampNum-1, GPreInvParam.dt);
            angleData = GPreInvParam.angleData;
            
        case 'angle_one_file'
            gather = bsReadGathersByIds(preDataInfo.fileName, preDataInfo.segyInfo, ...
                inline, crossline, startTime, sampNum-1, GPreInvParam.dt);
            angleSeisData = gather{1}.data;
            if ~isempty(GPreInvParam.angleData)
                angleData = GPreInvParam.angleData;
            else
                angleData = gather{1}.offsets;
            end
            
            if angleData(end) > 10
                angleData = angleData / 180 * pi;
            end
            
        case 'offset_one_file'
            gather = bsReadGathersByIds(preDataInfo.fileName, preDataInfo.segyInfo, ...
                inline, crossline, startTime, sampNum, GPreInvParam.dt);
            preData = gather{1}.data;
            offsets = gather{1}.offsets;
            
            [angleSeisData, angleData, ~] = bsOffsetData2AngleData(GPreInvParam, preData, offsets, ...
                initLog(:, 1), initLog(:, 2), initLog(:, 3), initLog(:, 4));
            
        case 'function'
            angleData = GPreInvParam.angleData;
            angleSeisData = preDataInfo.fcn(inline, crossline, startTime);
            if angleData(end) > 10
                angleData = angleData / 180 * pi;
            end
            
        otherwise
            validatestring(GPreInvParam.preSeisData.mode, ...
                'angle_separate_files', 'angle_one_file', 'offset_one_file');
    end

% -------------------------------------------------------------------------
    % build model parameter
    [initX, lsdCoef] = bsPreBuildModelParam(initLog, GPreInvParam.mode, GPreInvParam.lsdCoef);
    
    if GPreInvParam.isInitDeltaZero
        initX(sampNum+1:end) = 0;
    end
    
    % build forward matrix G
    if nargin < 6 || isempty(G)
        G = bsPreBuildGMatrix(...
                    GPreInvParam.mode, ...
                    initLog(:, 2), ...
                    initLog(:, 3), ...
                    angleData, ...
                    GPreInvParam.wavelet, ...
                    lsdCoef);
    end
            
    % reshape angle seismic data as a vector
    d = reshape(angleSeisData, GPreInvParam.angleTrNum*(sampNum-1), 1);
    
end


