function [angleSeisData, angleData, superTrData] = bsOffsetData2AngleData(GPreInvParam, preData, offsets, depth, initVp, initVs, initRho)
    
    % get super gather from pre-stack CDP gather
    [superTrData, offsetMax, offsetMin] = bsGetSuperGather(preData, offsets, GPreInvParam.oldSuperTrNum);
    
    % the last super traces is dropped
    superTrData = superTrData(:, 1:GPreInvParam.newSuperTrNum);
    
    % reset the range of offset
    GPreInvParam.offsetMax = (offsetMax - offsetMin)/(GPreInvParam.oldSuperTrNum - 1) * (GPreInvParam.newSuperTrNum - 1) + offsetMin;
    GPreInvParam.offsetMin = offsetMin;
   
    % calculate the angle information 
    [meanTheta, angleData] = bsGetAngleByRayTracking(GPreInvParam, depth, initVp, initVs, initRho);
    [angleIndex, angleCoef] = bsGetAngleIndexAndCoef(meanTheta, angleData);
    [angleSeisData] = bsOffset2Angle(angleIndex, angleCoef, superTrData, GPreInvParam.angleTrNum);

end

function [superTrData, offsetMax, offsetMin] = bsGetSuperGather(preData, offsets, superTrNum)
    offsetMax = max(offsets);
    offsetMin = min(offsets);
    
    [sampNum, nTrace] = size(preData);
    superTrData = zeros(sampNum, superTrNum);
    superNum = zeros(1, superTrNum);
    invOffset = (offsetMax-offsetMin) / superTrNum;
    
    for i = 1 : nTrace
        % get the group id which the current trace should be 
        superTraceId = (offsets(i) - offsetMin) / invOffset;
        superTraceId = floor(superTraceId) + 1;
        if(superTraceId > superTrNum)
            superTraceId = superTrNum;
        elseif(superTraceId < 1)
            superTraceId = 1;
        end
        
        % add the current trace onto the corresponding group
        superTrData(:, superTraceId) = superTrData(:, superTraceId) + preData(:, i); 
        superNum(superTraceId) = superNum(superTraceId) + 1;
    end
    
    superNum( superNum == 0 ) = 1;
    for i = 1 : superTrNum
        superTrData(:, i) = superTrData(:, i) / superNum(i);
    end
end

function [angleIndex, angleCoef] = bsGetAngleIndexAndCoef(angle, angleData)

    sampNum = size(angle, 1);
    angleTrNum = length(angleData);
    
    bb = cell(sampNum, angleTrNum);
    tt = cell(sampNum, angleTrNum);
    
    for j = 1  : sampNum
        for i = 1 : angleTrNum
            index0 = find(angle(j, : ) < angleData(i));
            index1 = find(angle(j, : ) > angleData(i));

            if isempty(index0)
                bb{j}{i} =[index1(1)];
                tt{j}{i} =1;
            end
            if isempty(index1)
                bb{j}{i} =[index0(length(index0))];
                tt{j}{i} =1;
            end
            if length(index0) == 1 || length(index1) == 1
                bb{j}{i} = [index0(length(index0)),index1(1)];
                tt{j}{i} = [angle(j,index0(length(index0))),angle(j,index1(1))]/sum([angle(j,index0(length(index0))),angle(j,index1(1))]);
            end
            if length(index0)>=2 && length(index1)>=2
                bb{j}{i} = [index0(length(index0)-1 : length(index0)),index1(1 : 2)];
                tt{j}{i} = [angle(j,index0(length(index0)-1 : length(index0))),angle(j,index1(1 : 2))]/sum([angle(j,index0(length(index0)-1 : length(index0))),angle(j,index1(1 : 2))]);
            end
        end
    end
    angleIndex = bb;
    angleCoef = tt;
end

function [angleTrace] = bsOffset2Angle(angleIndex, angleCoef, superTrData, angleTrNum)

    sampNum = length(angleIndex);
    angleTrace = zeros(sampNum, angleTrNum);

    for i = 1 : sampNum 
        for j = 1 : angleTrNum
            if angleIndex{i}{j}  ==  0
                angleTrace(i,j) = 0;
            else    
                angleTrace(i,j) = superTrData(i, angleIndex{i}{j})*angleCoef{i}{j}';
            end
        end
    end
end