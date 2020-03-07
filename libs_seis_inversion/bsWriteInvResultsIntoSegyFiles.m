function bsWriteInvResultsIntoSegyFiles(GInvParam, invResults, title)
    
    if isempty(invResults)
        return;
    end
    
    switch GInvParam.flag
        case {'prestack', 'pre-stack'}
            if isempty(GInvParam.postSeisData.fileName)
                fileName = GInvParam.preSeisData.fileName;
                GSegyInfo = GInvParam.preSeisData.segyInfo;
            else
                fileName = GInvParam.postSeisData.fileName;
                GSegyInfo = GInvParam.postSeisData.segyInfo;
            end

        case {'poststack', 'post-stack'}
            fileName = GInvParam.postSeisData.fileName;
            GSegyInfo = GInvParam.postSeisData.segyInfo;
    end
            
    for i = 1 : length(invResults)
        data = invResults{i}.data;
        
        
        if ~iscell(data)
            dstFileName = getDstFileName(GInvParam, invResults{i}.name, invResults{i}.type, title);
            bsWriteInvResultIntoSegyFile(invResults{i}, data, fileName, GSegyInfo, dstFileName);
                
        else
            for j = 1 : length(data)
                dstFileName = getDstFileName(GInvParam, invResults{i}.name, invResults{i}.type{j}, title);
                bsWriteInvResultIntoSegyFile(invResults{i}, data{j}, fileName, GSegyInfo, dstFileName);
            end
        end
        
    end
end

function str = getDstFileName(GInvParam, name, type, title)
    warning('off');
    mkdir(sprintf('%s/sgy_results', GInvParam.modelSavePath));
    str = sprintf('%s/sgy_results/%s-%s-%s.sgy', GInvParam.modelSavePath, type, name, title);
    warning('on');
end
