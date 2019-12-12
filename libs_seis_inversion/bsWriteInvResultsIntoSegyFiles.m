function bsWriteInvResultsIntoSegyFiles(GInvParam, invResults, title)
    
    if isempty(invResults)
        return;
    end
    
    for i = 1 : length(invResults)
        data = invResults{i}.data;
        
        
        if ~iscell(data)
            dstFileName = getDstFileName(GInvParam, invResults{i}.name, invResults{i}.type, title);
            
            switch GInvParam.flag
                case {'prestack', 'pre-stack'}
                bsWriteInvResultIntoSegyFile(invResults{i}, data, ...
                    GInvParam.preSeisData.fileName, ...
                    GInvParam.preSeisData.segyInfo, dstFileName);
                case {'poststack', 'post-stack'}
                bsWriteInvResultIntoSegyFile(invResults{i}, data, ...
                    GInvParam.postSeisData.fileName, ...
                    GInvParam.postSeisData.segyInfo, dstFileName);
            end
            
        else
            for j = 1 : length(data)
                dstFileName = getDstFileName(GInvParam, invResults{i}.name, invResults{i}.type{j}, title);
                switch GInvParam.flag
                    case {'prestack', 'pre-stack'}
                        bsWriteInvResultIntoSegyFile(invResults{i}, data{j}, ...
                            GInvParam.preSeisData.fileName, ...
                            GInvParam.preSeisData.segyInfo, dstFileName);
                    case {'poststack', 'post-stack'}
                        bsWriteInvResultIntoSegyFile(invResults{i}, data{j}, ...
                            GInvParam.postSeisData.fileName, ...
                            GInvParam.posteisData.segyInfo, dstFileName);
                end
                
            end
        end
        
    end
end

function str = getDstFileName(GInvParam, name, type, title)
    str = sprintf('%s/sgy_results/%s-%s-%s.sgy', GInvParam.modelSavePath, type, name, title);
end
