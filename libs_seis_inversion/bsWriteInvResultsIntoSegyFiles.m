function bsWriteInvResultsIntoSegyFiles(GPreInvParam, invResults, title)
    
    if isempty(invResults)
        return;
    end
    
    for i = 1 : length(invResults)
        data = invResults{i}.data;
        
        
        if ~iscell(data)
            dstFileName = getDstFileName(GPreInvParam, invResults{i}.name, invResults{i}.type, title);
            
            bsWriteInvResultIntoSegyFile(invResults{i}, data, ...
                GPreInvParam.preSeisData.segyFileName, ...
                GPreInvParam.preSeisData.segyInfo, dstFileName);
            
        else
            for j = 1 : length(data)
                dstFileName = getDstFileName(GPreInvParam, invResults{i}.name, invResults{i}.type{j}, title);
                
                bsWriteInvResultIntoSegyFile(invResults{i}, data{j}, ...
                    GPreInvParam.preSeisData.fileName, ...
                    GPreInvParam.preSeisData.segyInfo, dstFileName);
            end
        end
        
    end
end

function str = getDstFileName(GPreInvParam, name, type, title)
    str = sprintf('%s/sgy_results/%s-%s-%s.sgy', GPreInvParam.modelSavePath, type, name, title);
end
