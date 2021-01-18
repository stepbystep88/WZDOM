function bsShowInvProfiles(GInvParam, GShowProfileParam, profiles, wellLogs, timeLine)
%% Show the inversion results
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Nov 2019
% -------------------------------------------------------------------------
    
    nProfile = length(profiles);
    GPlotParam = GShowProfileParam.plotParam;
    [isSameType, type] = bsCheckIsSameType(profiles);
    
    if nProfile < 1
        return;
    end
    
    
    
    figure;
    if strcmp(GShowProfileParam.layout.mode, 'auto')
        % set the screen size
        if isSameType
            [nRow, nCol, loc, colorbar_pos] = setShareFigureSize(nProfile);
        else
            [nRow, nCol, loc] = setFigureSize(nProfile);
        end
    else
        nRow = GShowProfileParam.layout.nRow;
        nCol = GShowProfileParam.layout.nCol;
        loc = GShowProfileParam.layout.loc;
        colorbar_pos = GShowProfileParam.layout.colorbar_pos;
    end
    
    
    % show profiles
    for iProfile = 1 : nProfile
        [basicInfo] = bsInitBasicInfoForShowingProfile(GShowProfileParam, GInvParam, wellLogs, timeLine, profiles{iProfile});
        
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

    
    
    function bsShowShareColorbar()
        subplot('Position', colorbar_pos);
        axis off;
        
        % set colormap
        if strcmpi(type, 'seismic')
            colorTbl = bsGetColormap('seismic');
        else
            colorTbl = GShowProfileParam.colormap.allTheSame;
        end
        if GShowProfileParam.isColorReverse
            colorTbl = flipud(colorTbl);
            set(gca, 'colormap', colorTbl);
        else
            set(gca, 'colormap', colorTbl);
        end
    
        set(gca, 'clim', range);

        hc = colorbar( ...
            'fontsize', GPlotParam.fontsize,...
            'fontweight', GShowProfileParam.plotParam.fontweight, ...
            'fontname', GPlotParam.fontname);
        hc.Position = hc.Position + [0 0 0.01 0];

        ylabel(hc, attName, ...
            'fontsize', GShowProfileParam.plotParam.fontsize, ...
            'fontweight', GShowProfileParam.plotParam.fontweight, ...
            'fontname', GShowProfileParam.plotParam.fontname);
    end

%         gtext(attName, 'fontsize', GPlotParam.fontsize,'fontweight', 'bold', 'fontname', GPlotParam.fontname);
        
    % show one profile
    function [attName, range] = bsShowOneProfile(profile, isSameColorbar)
        
        fprintf('Showing the %s profile of %s...\n', upper(profile.type), profile.name);
        
        % data preprocessing base on the type of profile
        [range, scale, dataIndex, dataColorTbl, attName] = ...
            bsGetInfoByType(GShowProfileParam, GInvParam, profile.type);
        
        if ~isempty(GShowProfileParam.colormap.allTheSame) && ~strcmpi(profile.type, 'seismic')
            dataColorTbl = GShowProfileParam.colormap.allTheSame;
        end
        
        [data, range, wellData, wellTime] = bsGetRangeByType(...
                    profile, ...
                    wellLogs, ...
                    range, ...
                    scale, ...
                    basicInfo.wellIndex, ...
                    dataIndex, ...
                    GInvParam.indexInWellData.time, ...
                    GShowProfileParam.showWellFiltCoef);
                
        
        % extract part of the whole data to show based on 
        % GShowProfileParam.showLeftTrNumByWells and 
        % GShowProfileParam.showRightTrNumByWells
%         [profile, traceIds, wellPos, horizons] = bsExtractShowData(GShowProfileParam, profile, wellPos, horizons);
        data = data(:, basicInfo.left:basicInfo.right);
        
        % filter data along with horizon
        data = bsFilterProfileData(data, profile.showFiltCoef);
        
        % fill data by using horion information
%         [newProfileData, minTime] = bsHorizonRestoreData(profileData, profile.horizon, GInvParam.upNum, GInvParam.dt);

%         [newProfileData, newDt, newTraceIds, newHorizon] = bsReScaleData(GShowProfileParam.scaleFactor, ...
%             newProfileData, GInvParam.dt, traceIds, profile.horizon);
        
        % fill data by using horizon information, it also does some
        % interpolatation for visual effects
        [newProfileData] = bsReScaleAndRestoreData(basicInfo, data, GInvParam.isParallel);
            
        
        [newProfileData, minTime] = bsShowVerticalPartData(GShowProfileParam, basicInfo, newProfileData);    

        if ~isempty(wellLogs) && strcmpi(GShowProfileParam.showWellMode, 'color')
            % replace the traces that are near by well location as welllog data
            newProfileData = bsReplaceWellLocationData(GShowProfileParam, basicInfo, newProfileData, wellData, wellTime, minTime);
        end
        
        % show the filled data by using imagesc
        bsShowHorizonedData(GShowProfileParam, basicInfo, newProfileData, minTime, ...
            profile.name, range, dataColorTbl, wellData)

        % set attribute name of colorbar
        if ~isSameColorbar
            ylabel(colorbar(), attName, ...
                'fontsize', GShowProfileParam.plotParam.fontsize, ...
                'fontweight', GShowProfileParam.plotParam.fontweight, ...
                'fontname', GShowProfileParam.plotParam.fontname);
        end
        bsSetDefaultPlotSet(GShowProfileParam.plotParam);
        
    end
end


function [nRow, nCol, loc, colorbar_pos] = setShareFigureSize(nProfile)
    switch nProfile
        case 1
            bsSetPosition(0.3, 0.2);
            nRow = 1;
            nCol = 1;
            loc = [0.79, 0.89, 0.01, 0.08, 0.09, -0.01];
            colorbar_pos = [0.87 0.06 0.01 0.88];
        case 2
            bsSetPosition(0.3, 0.42);
            nRow = 2;
            nCol = 1;
            loc = [0.87, 0.83, 0.025, 0.08, 0.05, -0.005];
            colorbar_pos = [0.91 0.05 0.01 0.9];
%             loc = [0.79, 0.94, 0.01, 0.07, 0.09, 0.02];
            
%             loc = [0.84, 0.91, 0.01, 0.07, 0.06, 0.02];
%             colorbar_pos = [0.89 0.05 0.01 0.92];
        case 3
            bsSetPosition(0.3, 0.65);
            nRow = 3;
            nCol = 1;
            loc = [0.79, 0.93, 0.01, 0.06, 0.09, 0.01];
            colorbar_pos = [0.87 0.05 0.01 0.92];
%             loc = [0.84, 0.94, 0.02, 0.05, 0.07, 0.01];
%             colorbar_pos = [0.89 0.05 0.01 0.92];
        case 4
            bsSetPosition(0.6, 0.42);
            nRow = 2;
            nCol = 2;
            loc = [0.83, 0.9, 0.03, 0.07, 0.08, 0.01];
            colorbar_pos = [0.89 0.05 0.01 0.9];
        case {5, 6}
            bsSetPosition(0.6, 0.65);
            nRow = 3;
            nCol = 2;
            loc = [0.89, 0.95, 0.02, 0.06, 0.05, 0.01];
            colorbar_pos = [0.93 0.05 0.01 0.92];
        case {7, 8}
            bsSetPosition(0.6, 0.7);
            nRow = 4;
            nCol = 2;
%             loc = [0.83, 0.9, 0.03, 0.07, 0.08, 0.01];
%             colorbar_pos = [0.89 0.05 0.01 0.9];
            loc = [0.83, 0.95, 0.03, 0.05, 0.08, 0.01];
            colorbar_pos = [0.89 0.05 0.01 0.92];
            
        case {9}
            bsSetPosition(0.8, 0.65);
            nRow = 3;
            nCol = 3;
            loc = [0.9, 0.92, 0.01, 0.045, 0.04, -0.01];
            colorbar_pos = [0.94 0.05 0.01 0.92];
        otherwise
            bsSetPosition(0.9, 0.7);
            nRow = 3;
            nCol = 4;
            loc = [0.9, 0.92, 0.015, 0.045, 0.04, -0.01];
            colorbar_pos = [0.94 0.05 0.01 0.92];
    end
end

function [nRow, nCol, loc] = setFigureSize(nProfile)
    switch nProfile
        case 1
            bsSetPosition(0.3, 0.2);
            nRow = 1;
            nCol = 1;
            loc = [0.81, 0.82, 0.01, 0.08, 0.12, 0.00];
        case 2
            bsSetPosition(0.3, 0.42);
%             nRow = 2;
%             nCol = 1;
%             loc = [0.87, 0.91, 0.01, 0.07, 0.09, 0.02];
            nRow = 1;
            nCol = 2;
            loc = [0.93, 0.82, 0.02, 0.08, 0.05, 0.00];
        case 3
            bsSetPosition(0.3, 0.65);
            nRow = 3;
            nCol = 1;
            loc = [0.87, 0.92, 0.01, 0.04, 0.09, 0.00];
        case 4
            bsSetPosition(0.6, 0.42);
            nRow = 2;
            nCol = 2;
            loc = [0.93, 0.92, 0.02, 0.07, 0.05, 0.02];
        case {5, 6}
            bsSetPosition(0.6, 0.65);
            nRow = 3;
            nCol = 2;
            loc = [0.93, 0.93, 0.02, 0.04, 0.05, 0.00];
        case {7, 8, 9}
            bsSetPosition(0.8, 0.65);
            nRow = 3;
            nCol = 3;
            loc = [0.93, 0.93, 0.02, 0.05, 0.05, 0.00];
    end
end

function [isSameType, type] = bsCheckIsSameType(profiles)
    nProfiles = length(profiles);
    type = profiles{1}.type;
    
    isSameType = 1;
    if nProfiles == 0
        error('There profile data to show is empty.');
    elseif nProfiles == 1
        isSameType = 0;
        return;
    end
    
    
    for i = 2 : nProfiles
        itype = profiles{i}.type;
        
        if ~strcmpi(type, itype)
            isSameType = 0;
            break;
        end
    end
end

function [profileData, range, wellData, wellTime] ...
    = bsGetRangeByType(profile, wellLogs, range, scale, wellIndex, dataIndex, timeIndex, showWellFiltCoef)
    
    
    profileData = profile.data;
%     profileData(profileData<=0) = nan;    
    profileData = profileData / scale;

    if ~isempty(range)
        range = range / scale;
    else
        range = [prctile(profileData(:), 5), prctile(profileData(:), 95)];
    end
    
    if isempty(dataIndex)
        wellData = [];
        wellTime = [];
        return;
    end
    
    [wellData, wellTime] = bsGetWellData(wellLogs, wellIndex, dataIndex, timeIndex, showWellFiltCoef);
    
    if ~isempty(wellData)
        for i = 1 : length(wellData)
            wellData{i} = wellData{i} / scale;
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

function profileData = bsReplaceWellLocationData(GShowProfileParam, basicInfo, profileData, wellData, wellTime, minTime)

    if GShowProfileParam.isShowColorWells == 0
        return;
    end
    
    [newSampNum, trNum] = size(profileData);
    
    offsetNum = GShowProfileParam.showWellOffset * GShowProfileParam.scaleFactor + GShowProfileParam.emptyOffset;
    wellPos = basicInfo.wellPos;
    
    if ~isempty(wellData) && ~isempty(wellPos)
        % replace the data at well location by wellData
        for i = 1 : length(wellPos)
            traceId = basicInfo.traceIds(wellPos(i));
            [~, index] = min(abs(traceId - basicInfo.newTraceIds));
            s = index - offsetNum;
            if s < 1
                s = 1;
            end
            
            e = index + offsetNum;
            if e > trNum
                e = trNum;
            end
            
%             profileData(:, s : s + GShowProfileParam.emptyOffset-1) = inf;
%             profileData(:, e - GShowProfileParam.emptyOffset + 1 : e) = inf;
            
            for k = s : e
                z = profileData(:, k);
                ctime = wellTime{i} - basicInfo.newHorizon(index) + basicInfo.newHorizon(k);
                inf_set = 0;
                
                for j = 1 : newSampNum
                    tj = minTime + basicInfo.newDt * (j-1);
                    
                    tk = basicInfo.newHorizon(k) +  GShowProfileParam.scaleFactor * basicInfo.downNum * basicInfo.newDt;
                    
                    if tj < ctime(1)
%                         z(j) = wellData{i}(1);
                        z(j) = inf;
                    elseif tj > ctime(end) || tj > tk
%                         z(j) = wellData{i}(end);
                        inf_set = inf_set + 1;
                        z(j) = inf;
                        
                        if inf_set > GShowProfileParam.emptyOffset*2
                            break;
                        end
                    else
                        if k < s + GShowProfileParam.emptyOffset || k > e - GShowProfileParam.emptyOffset
                            z(j) = inf;
                        else
                            z(j) = bsCalVal(tj, ctime, wellData{i});
                        end
                    end
                end
                
                % replace serveral traces neary by the current well as the
                % corresponding welllog data
                profileData(:, k) = z;
            end
        end
    end
    
end


function [newProfileData, minTime] = bsShowVerticalPartData(GShowProfileParam, basicInfo, newProfileData)

    mode = lower(GShowProfileParam.showPartVert.mode);
    [sampNum, nTrace] = size(newProfileData);
    minTime = basicInfo.minTime;
    
    switch mode
        case {'in_2_horizons', 'in2horizons'}
            horizonIds = GShowProfileParam.showPartVert.horizonIds;
            
            if isempty(horizonIds)
                error("When GShowProfileParam.showPartVert.mode is 'in_2_horizons' or 'in2horizons', ...GShowProfileParam.showPartVert.horizonIds can not be empty.");
            end
            
            time = repmat(((0:sampNum-1)*basicInfo.newDt)', 1, nTrace) + basicInfo.minTime;

            % get horizons
            time1 = repmat(basicInfo.newHorizons(horizonIds(1), :), sampNum, 1);
            time2 = repmat(basicInfo.newHorizons(horizonIds(2), :), sampNum, 1);

            index = time<time1 | time>time2;
            newProfileData(index) = nan;
            
            index = sum(index, 2) / nTrace;
            [newProfileData, minTime] = bsSetNanDataAsEmpty(newProfileData, index, basicInfo.minTime, basicInfo.newDt);
        case {'up_down_time', 'updowntime'}
            [sampNum, nTrace] = size(newProfileData);
            time = repmat(((0:sampNum-1)*basicInfo.newDt)', 1, nTrace) + basicInfo.minTime;
            
            tmp = repmat(basicInfo.newHorizon, sampNum, 1);
            time1 = tmp - GShowProfileParam.showPartVert.upTime;
            time2 = tmp + GShowProfileParam.showPartVert.downTime;
            index = time < time1 | time > time2;
            newProfileData(index) = nan;
            [newProfileData, minTime] = bsSetNanDataAsEmpty(newProfileData, index, basicInfo.minTime, basicInfo.newDt);
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


function bsShowHorizonedData(GShowProfileParam, basicInfo, profileData, minTime, ...
    name, range, colorTbl, wellData)
    
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
    [label, data] = bsGenLabel(minTime, minTime+sampNum*basicInfo.newDt, sampNum, GShowProfileParam.yLabelNum);
    data = floor(data / 10) / 100;
    set(gca,'Ytick', label, 'YtickLabel', data);
    
    [label, data] = bsGenLabel(basicInfo.traceIds(1), basicInfo.traceIds(end), traceNum, GShowProfileParam.xLabelNum);
    set(gca,'Xtick', label, 'XtickLabel', data);
        
    % set title
    if strcmpi(GShowProfileParam.language, 'en')
        xlabel('Trace numbers');
    else
        xlabel('µÀºÅ');
    end
    
    if(~isempty(name))
        title(name, 'fontweight', GShowProfileParam.plotParam.fontweight); 
    end
    
    if strcmpi(GShowProfileParam.language, 'en')
        ylabel('Time (s)');
    else
        ylabel('Ê±¼ä \fontname{Times New Roman}(s)');
    end
    
%     set(gca, 'ydir', 'reverse');
%     set(gca, 'xlim', [traceIds(1), traceIds(end)]);
    
    % show horizon
    if( GShowProfileParam.isShowHorizon)
        for i = 1 : size(basicInfo.newHorizons, 1)
            y = 1 + round((basicInfo.newHorizons(i, :) - minTime) / basicInfo.newDt);
    %         y = horizon / 1000;
            plot(1:traceNum, y, 'k-', 'LineWidth', GShowProfileParam.plotParam.linewidth); hold on;
        end
    end
    
    nWell = length(basicInfo.wellPos);
    % show the curves of wells
    if ~isempty(wellData) && strcmpi(GShowProfileParam.showWellMode, 'curve')
        for i = 1 : nWell
            
            ipos = GShowProfileParam.scaleFactor * basicInfo.wellPos(i);
            offset = GShowProfileParam.scaleFactor*GShowProfileParam.showWellOffset;
            x = normalize(wellData{i}, 'range', [ipos-offset, ipos+offset]); 
            y = (-basicInfo.upNum : basicInfo.downNum-1)*basicInfo.dt + basicInfo.horizon(ipos);
            
            y = 1 + (y - minTime) / basicInfo.newDt;
            
            plot(x, y, 'k-', 'linewidth', 2);
        end
    end
    
    % show the names of wells
    if ~isempty(wellData) && GShowProfileParam.isShowWellNames
        for i = 1 : nWell
            ipos = GShowProfileParam.scaleFactor * basicInfo.wellPos(i);
            x = ipos - 5*GShowProfileParam.scaleFactor;
            y = round((basicInfo.newHorizon(ipos) - minTime) / basicInfo.newDt) +  GShowProfileParam.scaleFactor * (basicInfo.downNum + 5);
            
            
            text(x, y, basicInfo.wellNames{i}, ...
                'color', basicInfo.wellColors{i}, ...
                'fontsize', GShowProfileParam.plotParam.fontsize-2,...
                'fontweight', GShowProfileParam.plotParam.fontweight, ...
                'fontname', GShowProfileParam.plotParam.fontname, ...
                'Interpreter','none');
        end
    end
end

