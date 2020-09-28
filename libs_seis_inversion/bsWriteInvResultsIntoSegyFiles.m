function bsWriteInvResultsIntoSegyFiles(GInvParam, invResults, title, isSort, fileName, GSegyInfo)
    
    if isempty(invResults)
        return;
    end
    
    if nargin <= 3
        isSort = 1;
    end
    
    if nargin < 6
        switch GInvParam.flag
            case {'prestack', 'pre-stack'}
                if ~isempty(GInvParam.initModel.vp.fileName)
                    fileName = GInvParam.initModel.vp.fileName;
                    GSegyInfo = GInvParam.initModel.vp.segyInfo;
                elseif ~isempty(GInvParam.postSeisData.fileName)
                    fileName = GInvParam.postSeisData.fileName;
                    GSegyInfo = GInvParam.postSeisData.segyInfo;
                else
                    fileName = GInvParam.preSeisData.fileName;
                    GSegyInfo = GInvParam.preSeisData.segyInfo;
                end

            case {'poststack', 'post-stack'}
                fileName = GInvParam.postSeisData.fileName;
                GSegyInfo = GInvParam.postSeisData.segyInfo;
        end
    end
            
    for i = 1 : length(invResults)
        data = invResults{i}.data;
        
        
        if ~iscell(data)
            dstFileName = getDstFileName(GInvParam, invResults{i}.name, invResults{i}.type, title);
            bsWriteInvResultIntoSegyFile(invResults{i}, data, fileName, GSegyInfo, dstFileName, isSort);
                
        else
            for j = 1 : length(data)
                dstFileName = getDstFileName(GInvParam, invResults{i}.name, invResults{i}.type{j}, title);
                bsWriteInvResultIntoSegyFile(invResults{i}, data{j}, fileName, GSegyInfo, dstFileName, isSort);
            end
        end
        
    end
end

function str = getDstFileName(GInvParam, name, type, title)
    warning('off');
    mkdir(sprintf('%s/%s/sgy_results', GInvParam.modelSavePath, name));
    str = sprintf('%s/%s/sgy_results/%s-%s-%s.sgy', GInvParam.modelSavePath, name, type, name, title);
    warning('on');
end
