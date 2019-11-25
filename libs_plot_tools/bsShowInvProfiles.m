function bsShowInvProfiles(GPostInvParam, GShowProfileParam, profiles, wellLogs)
%% Show the inversion results
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    nProfile = length(profiles);
 
    figure;
    % set the screen size
    switch nProfile
        case 1
            bsSetPosition(0.4, 0.25);
            nRow = 1;
            nCol = 1;
            loc = [0.9, 0.9, 0.01, 0.08, 0.07, 0.01];
        case 2
            bsSetPosition(0.4, 0.5);
            nRow = 2;
            nCol = 1;
            loc = [0.9, 0.95, 0.01, 0.08, 0.07, 0.02];
        case 3
            bsSetPosition(0.4, 0.75);
            nRow = 3;
            nCol = 1;
            loc = [0.9, 0.96, 0.01, 0.06, 0.07, 0.01];
        case 4
            bsSetPosition(0.8, 0.5);
            nRow = 2;
            nCol = 2;
            loc = [0.98, 0.95, 0.05, 0.08, 0.04, 0.02];
        case {5, 6}
            bsSetPosition(0.8, 0.75);
            nRow = 3;
            nCol = 2;
            loc = [0.98, 0.96, 0.05, 0.06, 0.04, 0.01];
        case {7, 8, 9}
            bsSetPosition(0.95, 0.75);
            nRow = 3;
            nCol = 3;
            loc = [0.98, 0.96, 0.03, 0.06, 0.04, 0.01];
    end
    
    % show profiles
    for iProfile = 1 : nProfile
        
        profile = profiles{iProfile};
        
        bsSubPlotFit(nRow, nCol, iProfile, loc(1), loc(2), loc(3), loc(4), loc(5), loc(6));
        
        % change the name
        profile.name = sprintf('(%s) %s', 'a'+iProfile-1, profile.name);
        bsShowOneProfile(profile);
    end
    
    % show one profile
    function bsShowOneProfile(profile)
        
        profileData = profile.data;
        
        if profile.inIds(1) == profile.inIds(2)
            traceIds = profile.crossIds;
        else
            traceIds = profile.inIds;
        end
        
        % data preprocessing base on the type of profile
        switch profile.type
            case 'IP'
                
                profileData(profileData<=0) = nan;
                
                if(max(max(profileData)) > 1000)
                    profileData = profileData / 1000;
                end
                
                range = GShowProfileParam.rangeIP;
                if range(1) > 1000
                    range = range / 1000;
                end
                
                attName = 'Impedance (g/cm^3\cdotkm/s)';
                [wellPos, wellData] = bsFindWellLocation(GPostInvParam, ...
                    wellLogs, ...
                    profile.inIds, ...
                    profile.crossIds, ...
                    profile.horizon, ...
                    1, ...
                    2, ...
                    GShowProfileParam.showWellFiltCoef);
                if ~isempty(wellData)
                    wellData = wellData / 1000;
                end
                
            case 'Seismic'
                attName = 'Seismic (Amplitude)';
                wellPos = [];
                wellData = [];
                range = GShowProfileParam.rangeSeismic;
                
        end
        
        profile.horizon = bsButtLowPassFilter(profile.horizon, 0.1);
        
        % filter data along with horizon
        profileData = bsFilterData(profileData, profile.showFiltCoef);
        
        % replace the traces that are near by well location as welllog data
        profileData = bsReplaceWellLocationData(GShowProfileParam, profileData, wellPos, wellData);
        
        % fill data by using horion information
        [newProfileData, minTime] = bsHorizonRestoreData(profileData, profile.horizon, GPostInvParam.upNum, GPostInvParam.dt);

        [newProfileData, newDt, newTraceIds, newHorizon] = bsReScaleData(GShowProfileParam.scaleFactor, ...
            newProfileData, GPostInvParam.dt, traceIds, profile.horizon);
        
        % show the filled data by using imagesc
        bsShowHorizonedData(GShowProfileParam, ...
            newProfileData, ...
            newHorizon, minTime, newDt, profile.name, newTraceIds, range, GShowProfileParam.dataColorTbl);
        
        % set attribute name of colorbar
        ylabel(colorbar(), attName, ...
            'fontsize', GShowProfileParam.plotParam.fontsize, ...
            'fontweight', 'bold', ...
            'fontname', GShowProfileParam.plotParam.fontname);
        bsSetDefaultPlotSet(GShowProfileParam.plotParam);
        
    end
end

function [wellPos, wellData] = bsFindWellLocation(GPostInvParam, wellLogs, inIds, crossIds, horizon, dataIndex, timeIndex, showWellFiltCoef)
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    
    wellPos = [];
    wellData = [];
    for i = 1 : length(inIds)
        for j = 1 : length(wellInIds)
            if wellInIds(j) == inIds(i) && wellCrossIds(j) == crossIds(i)
                wellPos = [wellPos, i];
                
                subWellData = bsExtractWellDataByHorizon(...
                    wellLogs{j}.wellLog, horizon(i), dataIndex, timeIndex, ...
                    GPostInvParam.upNum, GPostInvParam.downNum, showWellFiltCoef);
                
                wellData = [wellData, subWellData];
            end
        end
    end
end

function profileData = bsReplaceWellLocationData(GShowProfileParam, profileData, wellPos, wellData)
    [~, trNum] = size(profileData);
    
    if ~isempty(wellPos)
        % replace the data at well location by wellData
        for i = 1 : length(wellPos)
            s = wellPos(i) - GShowProfileParam.showWellOffset;
            if s < 1
                s = 1;
            end
            
            e = wellPos(i) + GShowProfileParam.showWellOffset;
            if e > trNum
                e = trNum;
            end
            
            % replace serveral traces neary by the current well as the
            % corresponding welllog data
            profileData(:, s:e) = repmat(wellData(:, i), [1, e-s+1]);
        end
    end
    
end

function profileData = bsFilterData(profileData, showFiltCoef)
    [sampNum, ~] = size(profileData);
    
    % filter data along with horizon
    if showFiltCoef > 0 && showFiltCoef < 1
        try
            for i = 1 : sampNum
                profileData(i, :) = bsButtLowPassFilter(profileData(i, :), showFiltCoef);
            end
        catch
        end
    end
end



function [newProfileData, newDt, newTraceIds, newHorizon] = bsReScaleData(scaleFactor, profileData, dt, traceIds, horizon)
    if exist('scaleFactor', 'var') && scaleFactor > 1
        scaleFactor = round(scaleFactor);
        [sampNum, traceNum] = size(profileData);
        
        newDt = dt / scaleFactor;
        newTraceIds = linspace(traceIds(1), traceIds(end), traceNum * scaleFactor);
        [X, Y] = meshgrid(traceIds, linspace(1, sampNum, sampNum));
        [Xq,Yq] = meshgrid(newTraceIds, linspace(1, sampNum, sampNum * scaleFactor));
        
        newProfileData = interp2(X, Y, profileData, Xq, Yq,'cubic');
        newHorizon = interp1(traceIds, horizon, newTraceIds, 'spline');
    else
        newProfileData = profileData;
        newDt = dt;
        newTraceIds = traceIds;
    end
end

function bsShowHorizonedData(GShowProfileParam, profileData, horizon, minTime, dt, name, traceIds, range, colorTbl)
    
    [sampNum, traceNum] = size(profileData);
%     min_val = min(min(profileData));
%     min_val = min_val - abs(min_val)*10;
%     profileData(isnan(profileData )) = min_val;
    profileData(isinf(profileData )) = nan;
    
    h = imagesc(profileData); hold on;
    
    if ~isempty(range)
        set(gca, 'clim', range);
    end
    
    if GShowProfileParam.isColorReverse
        colorTbl = flipud(colorTbl);
        set(gcf, 'colormap', colorTbl);
    else
        set(gcf, 'colormap', colorTbl);
    end
    
    set(h, 'AlphaData', ~isnan(profileData));
    colorbar;
    
    % set labels of x and y axises
    [label, data] = bsGenLabel(minTime, minTime+sampNum*dt, sampNum, GShowProfileParam.yLabelNum);
    data = floor(data / 10) / 100;
    set(gca,'Ytick', label, 'YtickLabel', data);
    
    [label, data] = bsGenLabel(traceIds(1), traceIds(end), traceNum, GShowProfileParam.xLabelNum);
    set(gca,'Xtick', label, 'XtickLabel', data);
        
    % set title
    xlabel('');
    if(~isempty(name))
        title(name); 
    end
    ylabel('Time (s)');
    
%     set(gca, 'ydir', 'reverse');
%     set(gca, 'xlim', [traceIds(1), traceIds(end)]);
    
    % show horizon
    if( GShowProfileParam.isShowHorizon)
        y = 1 + round((horizon - minTime) / dt);
%         y = horizon / 1000;
        plot(1:traceNum, y, 'k-','LineWidth', GShowProfileParam.plotParam.linewidth); hold on;
    end
end