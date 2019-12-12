function bsShowErrorOfWellogs(GInvParam, GShowProfileParam, timeLine, wellLogs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the seismic error distribution of welllogs
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------

    data = cell(1, length(wellLogs));
    
    for i = 1 : length(wellLogs)
        wellInfo = wellLogs{i};
        inline = wellInfo.inline;
        crossline = wellInfo.crossline;
        
        [~, ~, horizonTime] = bsGetWellBaseInfo(timeLine{GInvParam.usedTimeLineId}, ...
            inline, crossline, 1, 2, 1, 2, 3);
        startTime = horizonTime - GInvParam.upNum * GInvParam.dt;
        
        switch GInvParam.flag
            case {'post-stack', 'poststack'}
                % obtain target welllog data based on horizon information
                trueLog = bsExtractWellDataByHorizon(...
                                wellInfo.wellLog, ...
                                horizonTime, ...
                                GInvParam.indexInWellData.ip, ...
                                GInvParam.indexInWellData.time, ...
                                GInvParam.upNum, ...
                                GInvParam.downNum, ...
                                1);

                [d, G, m] = bsPostBuild_d_G_m(GInvParam, ...
                    inline, crossline, startTime, trueLog, []);
        
            case {'pre-stack', 'prestack'}
            
                trueLog = bsExtractWellDataByHorizon(...
                    wellInfo.wellLog, ...
                    horizonTime, ...
                    [   indexInWellData.depth, ...
                        indexInWellData.vp, ...
                        indexInWellData.vs, ...
                        indexInWellData.rho], ...
                    indexInWellData.time, ...
                    GPreInvParam.upNum, ...
                    GPreInvParam.downNum, ...
                    1);

                [d, G, m] = bsPreBuild_d_G_m(GInvParam, ...
                    inline, crossline, startTime, trueLog, []);
            
            otherwise
                warning('No valid type of GInvParam.flag is found');
        end
        
        synData = G * m;
        realData = d;
        data{i} = [realData, synData, realData - synData];
        
    end
    
    bsShowDistribution(GShowProfileParam, wellLogs, data);
end


function bsShowDistribution(GShowProfileParam, wellLogs, data)
    nWell = length(wellLogs);
    figure;
    bsSetPosition(0.78, 0.59);
        
    nSubFigure = min(nWell, 15);
    nRow = ceil(nSubFigure / 5);
    
    for i = 1 : nSubFigure
        
        subplot(nRow, 5, i);
        
        histogram(data{i}(:, 3), 20, 'facecolor', 'b');
    end
    
    bsSetDefaultPlotSet(GShowProfileParam.plotParam);
end