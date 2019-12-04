function bsPreGetOtherAttributesByFile(GPreInvParam, timeLine, loadInfo)
    nFile = length(loadInfo.fileName);
    
    % horion of the whole volume
    usedTimeLine = timeLine{GPreInvParam.usedTimeLineId};
    sampNum = GPreInvParam.upNum + GPreInvParam.downNum;
    dt = GPreInvParam.dt;
    
    vp = [];
    vs = [];
    rho = [];
    
    for iFile = 1 : nFile
        if length(loadInfo.segyInfo) == 1
            % share the common segy basic info
            iSegyInfo = loadInfo.segyInfo;
        else
            % use different segy info for different files
            iSegyInfo = loadInfo.segyInfo(iFile);
        end
        
        fileName = loadInfo.fileName{iFile};
        if iFile == 1
            [rangeInline, rangeCrossline] = bsGetWorkAreaRange(iSegyInfo, fileName);
            [inIds, crossIds] = bsGetCDPsByRange(rangeInline, rangeCrossline);
            
            % horizon of given traces
            horizonTimes = bsCalcHorizonTime(usedTimeLine, inIds, crossIds, ...
                    GPreInvParam.isParallel, GPreInvParam.numWorkers);

            startTimes = horizonTimes - dt * GPreInvParam.upNum;
    
        end
        
        [data, outInIds, outCrossIds] = bsReadAllTraces(fileName, iSegyInfo, startTimes, sampNum, dt);

        switch lower(loadInfo.type{iFile})
            case 'vp'
                vp = data;
            case 'vs'
                vs = data;
            case {'rho', 'density'}
                rho = data;
        end
    end
    
    vp_vs = bsGetVp_Vs(vp, vs);
    possion = bsGetPossion(vp, vs);
    
    profile.inIds = inIds;
    profile.crossIds = crossIds;
    profile.horizon = horizonTimes;
    profile.upNum = GPreInvParam.upNum;
    profile.dt = GPreInvParam.dt;
    
    [filepath, name, ~] = fileparts(fileName);
    names = split(name, '_');
    
    dstFileName = sprintf('%s/vp_vs_%s.sgy', filepath, bsJointStr(names(2:end), '_'));
    bsWriteInvResultIntoSegyFile(profile, vp_vs, fileName, iSegyInfo, dstFileName);
    
    dstFileName = sprintf('%s/possion_%s.sgy', filepath, bsJointStr(names(2:end), '_'));
    bsWriteInvResultIntoSegyFile(profile, possion, fileName, iSegyInfo, dstFileName);
    
end

function vp_vs = bsGetVp_Vs(vp, vs)
    vp_vs = vp ./ vs;
end

function possion = bsGetPossion(vp, vs)
    vp_vs = (vp ./ vs).^2;
    possion = (0.5*vp_vs - 1)./(vp_vs - 1);
end

function str = bsJointStr(strs, dchar)
    str = '';
    for i = 1 : length(strs)
        str = [str, dchar, strs{i}];
    end
end