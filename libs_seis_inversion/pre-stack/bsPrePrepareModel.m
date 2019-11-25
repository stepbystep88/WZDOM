function model = bsPrePrepareModel(GPreInvParam, inline, crossline, horizonTime, trueLog, model)
%% create model package for prestack inversion which involves d, G, m, etc.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    
    
    sampNum = GPreInvParam.upNum + GPreInvParam.downNum; 
    
    % load model
    initModel = GPreInvParam.initModel;
    switch initModel.mode
        % the source of initial model
        % return initLog: 2D matrix, first-fourth columns are depth, vp,
        % vs, rho, respectively.
        case 'segy' % get initial model from segy file
            % start location of the inverted time interval
            initLog = bsReadMultiSegyFiles(GPreInvParam, ...
                [initModel.depth, initModel.vp, initModel.vs, initModel.rho], ...
                inId, crossId, horizonTime, sampNum);
 
            initLog = bsFiltWelllog(initLog, initModel.filtCoef);
            
        case 'filter_from_true_log' % get initial model by filtering the true model
            initLog = bsFiltWelllog(trueLog, initModel.filtCoef);
            
        case 'function' % get initial model by calling a function
            
            if isempty(initModel.fcn)
                error('When GPreInvParam.initModel.mode is function, the GPreInvParam.initModel.fcn could not be empty!\n');
            end
            initLog = initModel.fcn(GPreInvParam, inline, crossline);
            
        otherwise
            validatestring(GPreInvParam.initModel.mode, ['segy', 'filter_from_true_log', 'function']);
    end
    
    % load prestack seismic data
    switch GPreInvParam.preSeisData.mode
        case 'angle_separate_files'
            separates = GPreInvParam.preSeisData.separates;
            angleSeisData = bsReadMultiSegyFiles(GPreInvParam, separates, inId, crossId, horizonTime, sampNum-1);
            angleData = GPreInvParam.angleData;
            
        case 'angle_one_file'
            pos = bsCalcT0Pos(GPreInvParam, GPreInvParam.preSeisData.segyInfo, horizonTime);
            gather = bsReadGathersByIds(GPreInvParam.preSeisData.fileName, inId, crossId, pos, sampNum-1);
            angleSeisData = gather{1}.data;
            angleData = GPreInvParam.angleData;
            
        case 'offset_one_file'
            pos = bsCalcT0Pos(GPreInvParam, GPreInvParam.preSeisData.segyInfo, horizonTime);
            gather = bsReadGathersByIds(GPreInvParam.preSeisData.fileName, inId, crossId, pos, sampNum);
            preData = gather{1}.data;
            offsets = gather{1}.offsets;
            
            [angleSeisData, angleData, ~] = bsOffsetData2AngleData(GPreInvParam, preData, offsets, depth, initVp, initVs, initRho);
   
        otherwise
            validatestring(GPreInvParam.preSeisData.mode, ...
                'angle_separate_files', 'angle_one_file', 'offset_one_file');
    end

% -------------------------------------------------------------------------
    % build forward matrix G
    model.G = bsPreBuildGMatrix(...
                GPreInvParam.mode, ...
                initLog(:, 2), ...
                initLog(:, 3), ...
                angleData, ...
                GPreInvParam.wavelet, ...
                GPreInvParam.lsdCoff);
            
    % build model parameter
    [model.initX, model.lsdCoef] = bsPreBuildModelParam(initLog, GPreInvParam.mode, GPreInvParam.lsdCoff);
    GPreInvParam.lsdCoef = model.lsdCoef;
    
    % reshape angle seismic data as a vector
    model.d = reshape(angleSeisData, GInvParam.angleTrNum*(sampNum-1), 1);
    
% -------------------------------------------------------------------------            
    
    % start time of the inverted time interval
    model.t0 = round(horizonTime / GPreInvParam.dt) * GPreInvParam.dt - GPreInvParam.upNum * GPreInvParam.dt;
    model.inId = inline;
    model.crossId = crossline;
    model.initLog = initLog;
    
    if exist('trueLog', 'var') && ~isempty(trueLog)
        model.trueLog = trueLog;
        
        [model.trueX, ~] = bsPreBuildModelParam(trueLog, GPreInvParam.mode, model.lsdCoef);
    end
    
    % set boundary information
    [model.Lb, model.Ub] = bsGetBound(GPreInvParam, initLog);
    
    % normalize
    if GPreInvParam.isNormal
        model.maxAbsD = max(abs(model.d));
        model.d = model.d / model.maxAbsD;
        model.G = model.G / model.maxAbsD;    % we have to use the original G to normalize
    end
end

function seisData = bsReadMultiSegyFiles(GPreInvParam, separates, inId, crossId, horizonTime, sampNum)
    nFile = length(separates);
    seisData = size(sampNum, nFile);
    for i = 1 : nFile
        separate = separates(i);
        pos = bsCalcT0Pos(GPreInvParam, separate.segyInfo, horizonTime);
        seisData(:, i) = stpReadTracesByIds(separate.fileName, inId, crossId, pos, sampNum);
    end
end

function [Lb, Ub] = bsGetBound(GPreInvParam, initLog)
    bound = GPreInvParam.bound;
    sampNum = size(initLog, 1);
    
    switch bound.mode
        case 'off'
            Lb = [];
            Ub = [];
        case 'fixed'
            if length(bound.vp.Lb) == 1
                OneVector = ones(sampNum, 1);
                lvp = OneVector * bound.vp.Lb;
                uvp = OneVector * bound.vp.Ub;
                lvs = OneVector * bound.vs.Lb;
                uvs = OneVector * bound.vs.Ub;
                lrho = OneVector * bound.rho.Lb;
                urho = OneVector * bound.rho.Ub;
                
                Lb = bsPreBuildModelParam(...
                        [initLog(:, 1), lvp, lvs, lrho], ...
                        GPreInvParam.mode, ...
                        GPreInvParam.lsdCoef);
                Ub = bsPreBuildModelParam(...
                        [initLog(:, 1), uvp, uvs, urho], ...
                        GPreInvParam.mode, ...
                        GPreInvParam.lsdCoef);
            else
                Lb = bsPreBuildModelParam(...
                        [initLog(:, 1), bound.vp.Lb, bound.vs.Lb, bound.rho.Lb], ...
                        GPreInvParam.mode, ...
                        GPreInvParam.lsdCoef);
                Ub = bsPreBuildModelParam(...
                        [initLog(:, 1), bound.vp.Ub, bound.vs.Ub, bound.rho.Ub], ...
                        GPreInvParam.mode, ...
                        GPreInvParam.lsdCoef);
            end
        case 'based_on_init'
            offset = repmat([0, bound.vp.offset_init, bound.vs.offset_init, bound.rho.offset_init], sampNum, 1);
            
            Lb = bsPreBuildModelParam(...
                    initLog - offset, ...
                    GPreInvParam.mode, ...
                    GPreInvParam.lsdCoef);
            Ub = bsPreBuildModelParam(...
                    initLog + offset, ...
                    GPreInvParam.mode, ...
                    GPreInvParam.lsdCoef);
        otherwise
            validatestring(GPreInvParam.bound.mode, ['off', 'fixed', 'based_on_init']);
    end
end
