function bsShowInvProfiles(GInvParam, GShowProfileParam, profiles, wellLogs)
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
        switch lower(profile.type)
            case 'ip'
                
                [profileData, range, wellPos, wellData, wellTime] = bsGetRangeByType(...
                    profile, ...
                    wellLogs, ...
                    GShowProfileParam.rangeIP, ...
                    1000, ...
                    GInvParam.indexInWellData.Ip, ...
                    GInvParam.indexInWellData.time, ...
                    GShowProfileParam.showWellFiltCoef);

                attName = 'Impedance (g/cm^3\cdotkm/s)';
                
            case 'vp'
                
                [profileData, range, wellPos, wellData, wellTime] = bsGetRangeByType(...
                    profile, ...
                    wellLogs, ...
                    GShowProfileParam.rangeVP, ...
                    1000, ...
                    GInvParam.indexInWellData.vp, ...
                    GInvParam.indexInWellData.time, ...
                    GShowProfileParam.showWellFiltCoef);

                attName = 'V_P km/s';
                
            case 'vs'
                
                [profileData, range, wellPos, wellData, wellTime] = bsGetRangeByType(...
                    profile, ...
                    wellLogs, ...
                    GShowProfileParam.rangeVS, ...
                    1000, ...
                    GInvParam.indexInWellData.vs, ...
                    GInvParam.indexInWellData.time, ...
                    GShowProfileParam.showWellFiltCoef);

                attName = 'Density g/cm^3';
                
            case 'rho'
                
                [profileData, range, wellPos, wellData, wellTime] = bsGetRangeByType(...
                    profile, ...
                    wellLogs, ...
                    GShowProfileParam.rangeRho, ...
                    1, ...
                    GInvParam.indexInWellData.rho, ...
                    GInvParam.indexInWellData.time, ...
                    GShowProfileParam.showWellFiltCoef);

                attName = ' km/s';
                
            case 'seismic'
                attName = 'Seismic (Amplitude)';
                wellPos = [];
                wellData = [];
                range = GShowProfileParam.rangeSeismic;
                
        end
        
        profile.horizon = bsButtLowPassFilter(profile.horizon, 0.1);
        
        % filter data along with horizon
        profileData = bsFilterData(profileData, profile.showFiltCoef);
        
        % fill data by using horion information
%         [newProfileData, minTime] = bsHorizonRestoreData(profileData, profile.horizon, GInvParam.upNum, GInvParam.dt);

%         [newProfileData, newDt, newTraceIds, newHorizon] = bsReScaleData(GShowProfileParam.scaleFactor, ...
%             newProfileData, GInvParam.dt, traceIds, profile.horizon);
        
        [newProfileData, newDt, newTraceIds, newHorizon, minTime] = bsReScaleAndRestoreData(...
                GShowProfileParam.scaleFactor, ...
                profileData, ...
                profile.horizon, ...
                traceIds, ...
                GInvParam.upNum, ...
                GInvParam.dt);

                % replace the traces that are near by well location as welllog data
        newProfileData = bsReplaceWellLocationData(GShowProfileParam, ...
            newProfileData, ...
            traceIds, ...
            newTraceIds, ...
            wellPos, ...
            wellData, ...
            wellTime, ...
            minTime, ...
            newDt);
        
        
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

function [profileData, range, wellPos, wellData, wellTime] ...
    = bsGetRangeByType(profile, wellLogs, range, scale, dataIndex, timeIndex, showWellFiltCoef)
    
    profileData = profile.data;
    profileData(profileData<=0) = nan;    
    profileData = profileData / scale;

    range = range / scale;

    [wellPos, wellData, wellTime] = bsFindWellLocation(...
        wellLogs, ...
        profile.inIds, ...
        profile.crossIds, ...
        dataIndex, ...
        timeIndex, ...
        showWellFiltCoef);
    
    if ~isempty(wellData)
        wellData = wellData / scale;
    end
end

function [wellPos, wellData, wellTime] = bsFindWellLocation(wellLogs, inIds, crossIds, dataIndex, timeIndex, showWellFiltCoef)
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    
    wellPos = [];
    wellData = [];
    wellTime = [];
    
    for i = 1 : length(inIds)
        for j = 1 : length(wellInIds)
            if wellInIds(j) == inIds(i) && wellCrossIds(j) == crossIds(i)
                wellPos = [wellPos, i];
                
                data = wellLogs{j}.wellLog;
                for k = dataIndex
                    data(:, k) = bsButtLowPassFilter(data(:, k), showWellFiltCoef);
                end
    
%                 subWellData = bsExtractWellDataByHorizon(...
%                     wellLogs{j}.wellLog, horizon(i), dataIndex, timeIndex, ...
%                     GInvParam.upNum, GInvParam.downNum, showWellFiltCoef);
                
                wellData = [wellData, data(:, dataIndex)];
                wellTime = [wellTime, data(:, timeIndex)];
            end
        end
    end
end

function profileData = bsReplaceWellLocationData(GShowProfileParam, ...
    profileData, traceIds, newTraceIds, ...
    wellPos, wellData, time, t0, dt)
    [newSampNum, trNum] = size(profileData);
    
    offsetNum = GShowProfileParam.showWellOffset * GShowProfileParam.scaleFactor;
    if ~isempty(wellPos)
        % replace the data at well location by wellData
        for i = 1 : length(wellPos)
            
            traceId = traceIds(wellPos(i));
            [~, index] = min(abs(traceId - newTraceIds));
            s = index - offsetNum;
            if s < 1
                s = 1;
            end
            
            e = index + offsetNum;
            if e > trNum
                e = trNum;
            end
            
            z = zeros(newSampNum, 1);
            
            for j = 1 : newSampNum
                tj = t0 + dt * (j-1);
                
                if tj < time(1, i)
                    z(j) = wellData(1, i);
                elseif tj > time(end, i)
                    z(j) = wellData(end, i);
                else
                    z(j) = bsCalVal(tj, time(:, i), wellData(:, i));
                end
            end
            
            
            
            % replace serveral traces neary by the current well as the
            % corresponding welllog data
            profileData(:, s:e) = repmat(z, [1, e-s+1]);
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
        newTraceIds = linspace(traceIds(1), traceIds(end), traceNum * 1);
        [X, Y] = meshgrid(traceIds, linspace(1, sampNum, sampNum));
        [Xq,Yq] = meshgrid(newTraceIds, linspace(1, sampNum, sampNum * scaleFactor));
        
        newProfileData = interp2(X, Y, profileData, Xq, Yq,'cubic');
        newHorizon = interp1(traceIds, horizon, newTraceIds, 'spline');
    else
        newProfileData = profileData;
        newDt = dt;
        newTraceIds = traceIds;
        newHorizon = horizon;
    end
end

function [newProfileData, newDt, newTraceIds, newHorizon, minTime] = bsReScaleAndRestoreData(...
    scaleFactor, profileData, horizon, traceIds, ...
    upNum, dt, minTime)

    scaleFactor = round(scaleFactor);
    [sampNum, traceNum] = size(profileData);

    newDt = dt / scaleFactor;
    newTraceIds = linspace(traceIds(1), traceIds(end), traceNum * scaleFactor);
    
    time0 = horizon - upNum * dt;
    if ~exist('minTime', 'var')
        minTime = min(time0) - 5 * dt;
    end
    maxTime = max(time0) + sampNum * dt + 5 * dt;
%     t = minTime : dt : maxTime;
    newT = minTime : newDt : maxTime;
    
    [X, Y] = meshgrid(traceIds, newT);
    Z = inf(size(X));
    [Xq,Yq] = meshgrid(newTraceIds, newT);
    
    sequence = linspace(0, dt*(sampNum-1), sampNum);
    time = repmat(sequence', 1, traceNum)...
        + repmat(time0, sampNum, 1);
    newSampNum = length(newT);
    
    for j = 1 : traceNum
        for i = 1 : newSampNum
            ti = newT(i);

            if ti >= time(1, j) && ti<= time(sampNum, j)
                Z(i, j) = bsCalVal(ti, time(:, j), profileData(:, j));
            end
        end
    end
    
    newProfileData = interp2(X, Y, Z, Xq, Yq,'cubic');
    newHorizon = interp1(traceIds, horizon, newTraceIds, 'spline');
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