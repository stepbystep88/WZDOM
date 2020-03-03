function wellData = bsResampleWellData(data, dataIndex, vpIndex, depthIndex, timeIndex, dt)
    
    ct = 0;
    sampNum = size(data, 1);
    
    wellData = [];
    sumData = data(1, dataIndex);
    num = 1;
    
    for i = 2 : sampNum
        depth = data(i, depthIndex);
        
        
        dz = depth - data(i-1, depthIndex);
        it = dz * 2 / data(i, vpIndex) * 1000;
        ct = ct + it;
        sumData = sumData + data(i, dataIndex);
        
        num = num + 1;
        if ct > dt
            if isempty(timeIndex) || timeIndex<0
                wellData = [wellData; [depth, sumData/num]];
            else
                
%                 if isempty(wellData)
%                     time = data(i, timeIndex) * 1000;
%                   	wellData = [wellData; [depth, sumData/num, time]];
%                 else
%                     time = wellData(end, end) + dt;
%                     wellData = [wellData; [depth, sumData/num, time]];
%                 end
                time = data(i, timeIndex);
                wellData = [wellData; [depth, sumData/num, time]];
                
            end
            num = 0;
            sumData = zeros(1, length(dataIndex));
            ct = ct - dt;
        end
    end
end