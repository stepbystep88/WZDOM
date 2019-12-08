function bsShowInvProfiles(GInvParam, GShowProfileParam, profiles, wellLogs, timeLine)
%% Show the inversion results
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    nProfile = length(profiles);
    GPlotParam = GShowProfileParam.plotParam;
    isSameType = bsCheckIsSameType(profiles);
    nHorizon = length(timeLine);
    
    if nProfile < 1
        return;
    end
    
    % get horizon information
    horizons = bsGetHorizons();
    newHorizons = [];
    
    figure;
    % set the screen size
    if isSameType
        [nRow, nCol, loc, colorbar_pos] = setShareFigureSize(nProfile);
    else
        [nRow, nCol, loc] = setFigureSize(nProfile);
    end
    
    % show profiles
    for iProfile = 1 : nProfile
        
        profile = profiles{iProfile};
        
        bsSubPlotFit(nRow, nCol, iProfile, loc(1), loc(2), loc(3), loc(4), loc(5), loc(6));
        
        % change the name
        profile.name = sprintf('(%s) %s', 'a'+iProfile-1, profile.name);
        [attName, range] = bsShowOneProfile(profile, isSameType);
        
        if mod(iProfile, nCol) ~= 1 && nCol ~= 1
            set(gca, 'ylabel', []);
            set(gca, 'ytick', [], 'yticklabel', []);
        end
        
        if ceil(iProfile/nCol) == nRow || iProfile+nCol>nProfile
%             set(gca, 'ylabel', []);
        else
            set(gca, 'xtick', [], 'xticklabel', []);
            set(gca, 'xlabel', []);
        end
    end
    
    if isSameType
        bsShowShareColorbar();
    end

    function horizons = bsGetHorizons()
        if isfield(GShowProfileParam, 'isShowHorizon') ...
                && GShowProfileParam.isShowHorizon

            inIds = profiles{1}.inIds;
            crossIds = profiles{1}.crossIds;
            horizons = zeros(nHorizon, length(inIds));

            for i = 1 : nHorizon
                horizons(i, :) = bsGetHorizonTime(timeLine{i}, ...
                    inIds, crossIds, GInvParam.isParallel, GInvParam.numWorkers);
            end
        else
            horizons = [];
        end
        
    end
    
    function bsShowShareColorbar()
        subplot('Position', colorbar_pos);
        axis off;
        
        % set colormap
        colorTbl = GShowProfileParam.colormap.allTheSame;
        if GShowProfileParam.isColorReverse
            colorTbl = flipud(colorTbl);
            set(gca, 'colormap', colorTbl);
        else
            set(gca, 'colormap', colorTbl);
        end
    
        set(gca, 'clim', range);

        hc = colorbar( ...
            'fontsize', GPlotParam.fontsize,...
            'fontweight', 'bold', ...
            'fontname', GPlotParam.fontname);
        hc.Position = hc.Position + [0 0 0.01 0];

        ylabel(hc, attName, ...
            'fontsize', GShowProfileParam.plotParam.fontsize, ...
            'fontweight', 'bold', ...
            'fontname', GShowProfileParam.plotParam.fontname);
    end

%         gtext(attName, 'fontsize', GPlotParam.fontsize,'fontweight', 'bold', 'fontname', GPlotParam.fontname);
        
    % show one profile
    function [attName, range] = bsShowOneProfile(profile, isSameColorbar)
        
        fprintf('Showing the %s profile of %s...\n', upper(profile.type), profile.name);
        dataIndex = [];
        scale = 1;
        dataColorTbl = [];
        
        % data preprocessing base on the type of profile
        switch lower(profile.type)
            case 'ip'
                
                range = GShowProfileParam.range.ip;
                scale = 1000;
                dataIndex = GInvParam.indexInWellData.ip;
                dataColorTbl = GShowProfileParam.colormap.ip;
                attName = 'Impedance (g/cm^3\cdotkm/s)';
                
            case 'vp'
                range = GShowProfileParam.range.vp;
                scale = 1000;
                dataIndex = GInvParam.indexInWellData.vp;
                dataColorTbl = GShowProfileParam.colormap.vp;
                attName = 'V_P km/s';
                
            case 'vs'
                range = GShowProfileParam.range.vs;
                scale = 1000;
                dataIndex = GInvParam.indexInWellData.vs;
                dataColorTbl = GShowProfileParam.colormap.vs;
                attName = 'V_S km/s';
                
            case {'rho', 'density'}
                range = GShowProfileParam.range.rho;
                dataIndex = GInvParam.indexInWellData.rho;
                dataColorTbl = GShowProfileParam.colormap.rho;
                attName = 'Density g/cm^3';
            
            case {'vp_vs', 'vpvs_ratio'}
                attName = 'Vp/Vs ratio';
                range = GShowProfileParam.range.vpvs_ratio;
                dataIndex = GInvParam.indexInWellData.vpvs_ratio;
                dataColorTbl = GShowProfileParam.colormap.vpvs_ratio;
                
            case {'possion'}
                attName = 'Possion ratio';
                range = GShowProfileParam.range.possion;
                dataIndex = GInvParam.indexInWellData.possion;
                dataColorTbl = GShowProfileParam.colormap.possion;
            case 'seismic'
                attName = 'Seismic (Amplitude)';
%                 [wellPos, ~] = bsFindWellLocation(wellLogs, profile.inIds, profile.crossIds);
                range = GShowProfileParam.range.seismic;
                dataColorTbl = GShowProfileParam.colormap.seismic;
                
            otherwise
                validatestring({'seismic', 'ip', ...
                    'vp', 'vs', 'rho', 'density', ...
                    'vp_vs', 'vpvs_ratio', 'possion'})
        end
        
        if ~isempty(GShowProfileParam.colormap.allTheSame) && ~strcmpi(profile.type, 'seismic')
            dataColorTbl = GShowProfileParam.colormap.allTheSame;
        end
        
        [profile.data, range, wellPos, wellData, wellTime, wellNames] = bsGetRangeByType(...
                    profile, ...
                    wellLogs, ...
                    range, ...
                    scale, ...
                    dataIndex, ...
                    GInvParam.indexInWellData.time, ...
                    GShowProfileParam.showWellFiltCoef);
                
        % smooth the horizon
        profile.horizon = bsButtLowPassFilter(profile.horizon, 0.1);
        
        % extract part of the whole data to show based on 
        % GShowProfileParam.showLeftTrNumByWells and 
        % GShowProfileParam.showRightTrNumByWells
        [profile, traceIds, wellPos] = bsExtractShowData(GShowProfileParam, profile, wellPos);

        % filter data along with horizon
        profile.data = bsFilterProfileData(profile.data, profile.showFiltCoef);
        
        % fill data by using horion information
%         [newProfileData, minTime] = bsHorizonRestoreData(profileData, profile.horizon, GInvParam.upNum, GInvParam.dt);

%         [newProfileData, newDt, newTraceIds, newHorizon] = bsReScaleData(GShowProfileParam.scaleFactor, ...
%             newProfileData, GInvParam.dt, traceIds, profile.horizon);
        
        % fill data by using horizon information, it also does some
        % interpolatation for visual effects
        [newProfileData, newDt, newTraceIds, newHorizon, minTime] = bsReScaleAndRestoreData(...
                GShowProfileParam.scaleFactor, ...
                profile.data, ...
                profile.horizon, ...
                traceIds, ...
                GInvParam.upNum, ...
                GInvParam.dt, ...
                [], ...
                GInvParam.isParallel);
            
        
        [newProfileData, newHorizons, minTime] = bsShowVerticalPartData(GShowProfileParam, horizons, ...
                traceIds, newProfileData, newTraceIds, newHorizon, newHorizons, newDt, minTime);    

        % replace the traces that are near by well location as welllog data
        newProfileData = bsReplaceWellLocationData(GShowProfileParam, ...
            newProfileData, ...
            traceIds, ...
            newTraceIds, ...
            newHorizon, ...
            wellPos, ...
            wellData, ...
            wellTime, ...
            minTime, ...
            newDt);
        
        
        % show the filled data by using imagesc
        bsShowHorizonedData(GShowProfileParam, ...
            newProfileData, ...
            newHorizons, minTime, newDt, profile.name, ...
            newTraceIds, range, dataColorTbl, wellPos, wellNames);
        
        % set attribute name of colorbar
        if ~isSameColorbar
            ylabel(colorbar(), attName, ...
                'fontsize', GShowProfileParam.plotParam.fontsize, ...
                'fontweight', 'bold', ...
                'fontname', GShowProfileParam.plotParam.fontname);
        end
        bsSetDefaultPlotSet(GShowProfileParam.plotParam);
        
    end
end

function [newProfileData, newHorizons, minTime] = bsShowVerticalPartData(GShowProfileParam, horizons, ...
    traceIds, newProfileData, newTraceIds, newHorizon, newHorizons, newDt, minTime)

    mode = lower(GShowProfileParam.showPartVert.mode);
    [sampNum, nTrace] = size(newProfileData);
    
    if (GShowProfileParam.isShowHorizon ...
            || strcmp(mode, 'in_2_horizons') ...
            || strcmp(mode, 'in2horizons')) && isempty(newHorizons)
        nHorizon = size(horizons, 1);
        newHorizons = zeros(nHorizon, nTrace);
        
        for i = 1 : nHorizon
            newHorizons(i, :) = interp1(traceIds, horizons(i, :), newTraceIds, 'spline');
        end
    end
    
    switch mode
        case {'in_2_horizons', 'in2horizons'}
            horizonIds = GShowProfileParam.showPartVert.horizonIds;
            
            if isempty(horizonIds)
                error("When GShowProfileParam.showPartVert.mode is 'in_2_horizons' or 'in2horizons', ...GShowProfileParam.showPartVert.horizonIds can not be empty.");
            end
            
            time = repmat(((0:sampNum-1)*newDt)', 1, nTrace) + minTime;

            % get horizons
            time1 = repmat(newHorizons(horizonIds(1), :), sampNum, 1);
            time2 = repmat(newHorizons(horizonIds(2), :), sampNum, 1);

            index = time<time1 | time>time2;
            newProfileData(index) = nan;
            
            index = sum(index, 2) / nTrace;
            [newProfileData, minTime] = bsSetNanDataAsEmpty(newProfileData, index, minTime, newDt);
        case {'up_down_time', 'updowntime'}
            [sampNum, nTrace] = size(newProfileData);
            time = repmat(((0:sampNum-1)*newDt)', 1, nTrace) + minTime;
            
            tmp = repmat(newHorizon, sampNum, 1);
            time1 = tmp - GShowProfileParam.showPartVert.upTime;
            time2 = tmp + GShowProfileParam.showPartVert.downTime;
            index = time < time1 | time > time2;
            newProfileData(index) = nan;
            [newProfileData, minTime] = bsSetNanDataAsEmpty(newProfileData, index, minTime, newDt);
        case {'off'}
         
        otherwise
            validatestring(mode, ...
                {'in_2_horizons', 'in2horizons', 'off', 'up_down_time', 'updowntime'})
    end
    
end

function [newProfileData, minTime] = bsSetNanDataAsEmpty(newProfileData, index, minTime, newDt)
    index = mean(index, 2);
    sampNum = length(index);
    sPos = 0;
    ePos = sampNum+1;
    
    for i = 1 : sampNum
        if index(i) == 1
            sPos = i;
        else
            break;
        end
    end
    
    for i = 1 : sampNum
        if index(sampNum - i + 1) == 1
            ePos = sampNum - i + 1;
        else
            break;
        end
    end
    
    sPos = sPos - 5;
    ePos = ePos + 5;
    
    minTime = minTime + newDt * sPos;
    newProfileData(ePos+1:end, :) = [];
    newProfileData(1:sPos, :) = [];
    
end
        
function [profile, traceIds, wellPos] = bsExtractShowData(GShowProfileParam, profile, wellPos)

    if length(unique(profile.inIds)) < length(unique(profile.crossIds))
        traceIds = profile.crossIds;
    else
        traceIds = profile.inIds;
    end
    
    if isempty(wellPos)
        return;
    end
    
    left = min(wellPos) - GShowProfileParam.showLeftTrNumByWells;
    right = max(wellPos) + GShowProfileParam.showRightTrNumByWells;
    
    if left < 1
        left = 1;
    else
        wellPos = wellPos - min(wellPos) + 1 + GShowProfileParam.showLeftTrNumByWells;
    end
    
    if right > length(profile.inIds)
        right = length(profile.inIds);
    end
    
    profile.horizon = profile.horizon(left:right);
    profile.data = profile.data(:, left:right);
    profile.inIds = profile.inIds(:, left:right);
    profile.crossIds = profile.crossIds(:, left:right);
    traceIds = traceIds(left:right);
end

function [nRow, nCol, loc, colorbar_pos] = setShareFigureSize(nProfile)
    switch nProfile
        case 1
            bsSetPosition(0.3, 0.2);
            nRow = 1;
            nCol = 1;
            loc = [0.8, 0.89, 0.01, 0.08, 0.08, -0.01];
            colorbar_pos = [0.87 0.06 0.01 0.88];
        case 2
            bsSetPosition(0.3, 0.42);
            nRow = 2;
            nCol = 1;
            loc = [0.8, 0.94, 0.01, 0.07, 0.08, 0.02];
            colorbar_pos = [0.87 0.06 0.01 0.9];
        case 3
            bsSetPosition(0.3, 0.65);
            nRow = 3;
            nCol = 1;
            loc = [0.8, 0.95, 0.01, 0.04, 0.08, 0.01];
            colorbar_pos = [0.87 0.05 0.01 0.92];
        case 4
            bsSetPosition(0.6, 0.42);
            nRow = 2;
            nCol = 2;
            loc = [0.9, 0.94, 0.02, 0.07, 0.04, 0.02];
            colorbar_pos = [0.93 0.05 0.01 0.9];
        case {5, 6}
            bsSetPosition(0.6, 0.65);
            nRow = 3;
            nCol = 2;
            loc = [0.9, 0.95, 0.02, 0.04, 0.04, 0.01];
            colorbar_pos = [0.93 0.05 0.01 0.92];
        case {7, 8, 9}
            bsSetPosition(0.8, 0.65);
            nRow = 3;
            nCol = 3;
            loc = [0.9, 0.95, 0.01, 0.04, 0.04, 0.01];
            colorbar_pos = [0.94 0.05 0.01 0.92];
    end
end

function [nRow, nCol, loc] = setFigureSize(nProfile)
    switch nProfile
        case 1
            bsSetPosition(0.3, 0.2);
            nRow = 1;
            nCol = 1;
            loc = [0.89, 0.84, 0.01, 0.08, 0.08, 0.00];
        case 2
            bsSetPosition(0.3, 0.42);
            nRow = 2;
            nCol = 1;
            loc = [0.89, 0.93, 0.01, 0.07, 0.08, 0.02];
        case 3
            bsSetPosition(0.3, 0.65);
            nRow = 3;
            nCol = 1;
            loc = [0.89, 0.94, 0.01, 0.04, 0.08, 0.00];
        case 4
            bsSetPosition(0.6, 0.42);
            nRow = 2;
            nCol = 2;
            loc = [0.95, 0.94, 0.02, 0.07, 0.04, 0.02];
        case {5, 6}
            bsSetPosition(0.6, 0.65);
            nRow = 3;
            nCol = 2;
            loc = [0.95, 0.95, 0.02, 0.04, 0.04, 0.01];
        case {7, 8, 9}
            bsSetPosition(0.8, 0.65);
            nRow = 3;
            nCol = 3;
            loc = [0.95, 0.95, 0.01, 0.04, 0.04, 0.01];
    end
end

function isSameType = bsCheckIsSameType(profiles)
    nProfiles = length(profiles);
    isSameType = 1;
    if nProfiles == 0
        error('There profile data to show is empty.');
    elseif nProfiles == 1
        isSameType = 0;
        return;
    end
    
    type = profiles{1}.type;
    for i = 2 : nProfiles
        itype = profiles{i}.type;
        
        if ~strcmpi(type, itype)
            isSameType = 0;
            break;
        end
    end
end

function [profileData, range, wellPos, wellData, wellTime, wellNames] ...
    = bsGetRangeByType(profile, wellLogs, range, scale, dataIndex, timeIndex, showWellFiltCoef)
    
    
    profileData = profile.data;
%     profileData(profileData<=0) = nan;    
    profileData = profileData / scale;

    if ~isempty(range)
        range = range / scale;
    end
    
    if isempty(dataIndex)
        wellData = [];
        wellTime = [];
        wellPos = [];
        wellNames = [];
        return;
    end
    
    [wellPos, wellIndex, wellNames] = bsFindWellLocation(wellLogs, profile.inIds, profile.crossIds);
    
    [wellData, wellTime] = bsGetWellData(wellLogs, wellIndex, dataIndex, timeIndex, showWellFiltCoef);
    
    if ~isempty(wellData)
        for i = 1 : length(wellData)
            wellData{i} = wellData{i} / scale;
        end
    end
end

function [wellPos, wellIndex, wellNames] = bsFindWellLocation(wellLogs, inIds, crossIds)

    wellPos = [];
    wellIndex = [];
    wellNames = {};
    
    if isempty(wellLogs)
        return;
    end
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];
    
    for i = 1 : length(inIds)
        for j = 1 : length(wellInIds)
            if wellInIds(j) == inIds(i) && wellCrossIds(j) == crossIds(i)
                wellPos = [wellPos, i];
                wellIndex = [wellIndex, j];
                wellNames = [wellNames, wellLogs{j}.name];
            end
        end
    end
end

function [wellData, wellTime] = bsGetWellData(wellLogs, wellIndex, dataIndex, timeIndex, showWellFiltCoef)

    
    wellData = cell(0);
    wellTime = cell(0);

    for j = wellIndex
        data = wellLogs{j}.wellLog;
        for k = dataIndex
            data(:, k) = bsButtLowPassFilter(data(:, k), showWellFiltCoef);
        end

        wellData = [wellData, data(:, dataIndex)];
        wellTime = [wellTime, data(:, timeIndex)];
    end
end

function profileData = bsReplaceWellLocationData(GShowProfileParam, ...
    profileData, traceIds, newTraceIds, newHorizons, ...
    wellPos, wellData, wellTime, t0, dt)

    if GShowProfileParam.isShowColorWells == 0
        return;
    end
    
    [newSampNum, trNum] = size(profileData);
    
    offsetNum = GShowProfileParam.showWellOffset * GShowProfileParam.scaleFactor;
    if ~isempty(wellData) && ~isempty(wellPos)
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
            
            for k = s : e
                z = zeros(newSampNum, 1);
                ctime = wellTime{i} - newHorizons(index) + newHorizons(k);
                for j = 1 : newSampNum
                    tj = t0 + dt * (j-1);

                    if tj < ctime(1)
                        z(j) = wellData{i}(1);
                    elseif tj > ctime(end)
                        z(j) = wellData{i}(end);
                    else
                        z(j) = bsCalVal(tj, ctime, wellData{i});
                    end
                end
                
                % replace serveral traces neary by the current well as the
                % corresponding welllog data
                profileData(:, k) = z;
            end
        end
    end
    
end

function [newProfileData, newDt, newTraceIds, newHorizon, minTime] = bsReScaleAndRestoreData(...
    scaleFactor, profileData, horizon, traceIds, ...
    upNum, dt, minTime, isParallel)

    scaleFactor = round(scaleFactor);
    [sampNum, traceNum] = size(profileData);

    newDt = dt / scaleFactor;
    newTraceIds = linspace(traceIds(1), traceIds(end), traceNum * scaleFactor);
    
    time0 = horizon - upNum * dt;
    if ~exist('minTime', 'var') || isempty(minTime)
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
    
    if isParallel
        parfor j = 1 : traceNum
            for i = 1 : newSampNum
                ti = newT(i);

                if ti >= time(1, j) && ti<= time(sampNum, j)
                    Z(i, j) = bsCalVal(ti, time(:, j), profileData(:, j));
                end
            end
        end
    else
        for j = 1 : traceNum
            for i = 1 : newSampNum
                ti = newT(i);

                if ti >= time(1, j) && ti<= time(sampNum, j)
                    Z(i, j) = bsCalVal(ti, time(:, j), profileData(:, j));
                end
            end
        end
    end
    
    newProfileData = interp2(X, Y, Z, Xq, Yq,'cubic');
    newHorizon = interp1(traceIds, horizon, newTraceIds, 'spline');
end


function bsShowHorizonedData(GShowProfileParam, profileData, horizons, minTime, dt, ...
    name, traceIds, range, colorTbl, wellPos, wellNames)
    
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
        set(gca, 'colormap', colorTbl);
    else
        set(gca, 'colormap', colorTbl);
    end
    
    set(h, 'AlphaData', ~isnan(profileData));
%     colorbar;
    
    % set labels of x and y axises
    [label, data] = bsGenLabel(minTime, minTime+sampNum*dt, sampNum, GShowProfileParam.yLabelNum);
    data = floor(data / 10) / 100;
    set(gca,'Ytick', label, 'YtickLabel', data);
    
    [label, data] = bsGenLabel(traceIds(1), traceIds(end), traceNum, GShowProfileParam.xLabelNum);
    set(gca,'Xtick', label, 'XtickLabel', data);
        
    % set title
    xlabel('Trace numbers');
    if(~isempty(name))
        title(name, 'fontweight', 'bold'); 
    end
    ylabel('Time (s)');
    
%     set(gca, 'ydir', 'reverse');
%     set(gca, 'xlim', [traceIds(1), traceIds(end)]);
    
    % show horizon
    if( GShowProfileParam.isShowHorizon)
        for i = 1 : size(horizons, 1)
            y = 1 + round((horizons(i, :) - minTime) / dt);
    %         y = horizon / 1000;
            plot(1:traceNum, y, 'k-', 'LineWidth', GShowProfileParam.plotParam.linewidth); hold on;
        end
    end
    
    % show the name of wells
    if ~isempty(wellPos)
        for i = 1 : length(wellPos)
            ipos = GShowProfileParam.scaleFactor * wellPos(i) - 0.05*traceNum;
            text(ipos, 20, wellNames{i}, ...
                'fontsize', GShowProfileParam.plotParam.fontsize,...
                'fontweight', 'bold', ...
                'fontname', GShowProfileParam.plotParam.fontname, ...
                'Interpreter','none');
        end
    end
end



