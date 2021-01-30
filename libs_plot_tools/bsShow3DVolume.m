function [houts, horizonSlices] = bsShow3DVolume(data, dt, clim, xslices, yslices, zslices, horizons, varargin)
%% This function is used to display a 3D volume
%
% Inputs
%
% data : a 3d volume of size (sampNum*nInline*nCrossline), where sampNum is
% the number of sample points for each trace, nInline is the number of
% inline sections, nCrossline is the number of traces for each inline
% section. If the data is not organized as (sampNum*nInline*nCrossline),
% you can adjust the order by using function "permute".
%
% dt : the sample interval in microsecond.
%
% clim : which is used to set the colormap limits for the current axes. It
% should be a vector of size 1*2.
%
% xslices : the inline sections will be shown. The value should be a vector
% representing the inline numbers.
%
% yslices : the crossline sections will be shown. The value should be a vector
% representing the crossline numbers.
%
% zslices : the slices with the same time will be shown. The value should be a vector
% representing a series of time nodes.
%
% horizons : the horizons indicate the slices which not has the same time,
% so it need to indicate the time information of each horizon. It should be
% a array of structures, for each structure, it has .horizon of size
% (nInline*nCrossline) and .shift representing the shift time in
% microsecond.
%
% startTime : the time of the first sample point of this volume in
% microsecond. When startTime a scalar, it means each trace has the same
% starting time. Otherwise, it should be the size of (nInline*nCrossline).
%
% firstInlineId : the minimum inline number of the volume, the default
% value is 1.
%
% firstCrosslineId : the minimum crossline number of each inline section,
% the default value is 1.
%
% colormap : the colormap for showing the data
%
% attributeName : the attribute name of the data. For instance, it would be
% "Amplitude" for post-stack seismic data, or "Velocity (m/s)" for velocity data. 
%
% fontName : defualt is 'Times new roman'
%
% fontsize : default is 11.
% 
% fontweight : default is 'bold'.
%
% welllogs : the welllogs of the voume. It is a structure array. Each
% structrue records the information of one well, including .inline,
% .crossline, .data, .name, .color, .startTime
%
% showsize : the size of each shown well in the figure.
%
% filtCoef : the cutoff angle frequency of well-log data. The range of it is (0, 1]. 
% The default value is 1 which means no filtering is performed.
%
% isShading : whether the background of is white. The defaut value is
% false.
%
% Written by Bin She, University of Electronic Science and Technology of
% China, 01-31-2019 
% Email: bin.stepbystep@gmail.com
%

    [sampNum, nInline, nCrossline] = size(data);
    
    if isempty(clim)
        clim = [prctile(data(:), 5), prctile(data(:), 95)];
    end
    
    %% check inputs
    p = inputParser;
    
    addRequired(p, 'dt', @(x)(x>0));
    addRequired(p, 'clim', @(x) validateattributes(x,{'numeric'},{'size', [1, 2]}));
    
    addParameter(p, 'startTime', 0, @(x) (isscalar(x)&&(x>0)) ||  (size(x, 1) == nInline && size(x, 2) == nCrossline) );
    addParameter(p, 'firstInlineId', 1, @(x) (isscalar(x)&&(x>0)) );
    addParameter(p, 'firstCrosslineId', 1, @(x) (isscalar(x)&&(x>0)) );
    addParameter(p, 'colormap', bsGetColormap('velocity') );
    addParameter(p, 'fontname', 'Times new roman' );
    addParameter(p, 'fontsize', 11 );
    addParameter(p, 'fontweight', 'bold' );
    addParameter(p, 'attributeName', 'AttributeName' );
    addParameter(p, 'filtCoef', 1, @(x) ((isscalar(x)&&(x>0)&&(x<=1))));
    addParameter(p, 'welllogs', []);
    addParameter(p, 'dataIndex', 1);
    addParameter(p, 'showsize', 1, @(x) (isscalar(x)&&(x>0)));
    addParameter(p, 'isShading', 1, @(x) (islogical(x)));
    
    addParameter(p, 'isSmooth', 0);
    addParameter(p, 'smoothFiltCoef', 0.3);
    
    addParameter(p, 'scaleFactor', 1);
    addParameter(p, 'view', [-49.3414   37.6342]);
    
    addParameter(p, 'zlabel_str', 'Time (s)');
    
    p.parse(dt, clim, varargin{:});  
    params = p.Results;
    
    if params.isShading
        params.colormap = [params.colormap];
    end

    %% calculate the wells which will be shown
    welllogs = params.welllogs;
%     wellData = [];
    wellPos = [];
    wellIndex = [];
    for i = 1 : length(welllogs)
        iwell = welllogs{i};
        
        if ismember(iwell.crossline, yslices) || ismember(iwell.inline, xslices)
            tmp = bsButtLowPassFilter(iwell.wellLog(:, params.dataIndex), params.filtCoef);
            
            if( length(iwell.wellLog(:, params.dataIndex)) ~= sampNum)
               fprintf('The length of the well-log data must be the same as the data volume!!!\n');
               return;
            end
            
            ix = iwell.inline - params.firstInlineId + 1;
            iy = iwell.crossline - params.firstCrosslineId + 1;
            
%             wellData = [wellData, tmp];
            wellPos = [wellPos; [ix, iy]];
            wellIndex = [wellIndex; i];
            
            showsize = params.showsize;
            
            for m = ix-showsize:ix+showsize-1
                for n = iy-showsize:iy+showsize-1
                    if m<1 || n<1 || m>nInline || n>nCrossline
                        continue;
                    end
                    
                    data(:, m, n) = tmp;
                end
            end
        end
    end
    
    % smooth the slice
    if params.isSmooth
        for i = xslices
            profileData = reshape(data(:, i, :), sampNum, []);
            data(:, i, :) = bsFilterProfileData(profileData, params.smoothFiltCoef, 1);
        end    
        for i = yslices
            profileData = reshape(data(:, :, i), sampNum, []);
            data(:, :, i) = bsFilterProfileData(profileData, params.smoothFiltCoef, 1);
        end
    end
    
    firstInline = params.firstInlineId;
    firstCrossline = params.firstCrosslineId;
    endInline = firstInline + nInline - 1;
    endCrossline = firstCrossline + nCrossline - 1;
    
    
    %% prepare basical information
    if isscalar(params.startTime)
        minTime = params.startTime;
        maxTime = minTime + (sampNum-1) * dt;
        
        t = minTime + (0:sampNum-1)*dt;
        t = t / 1000;
        [x3, y3, z3] = meshgrid(firstInline : endInline, firstCrossline : endCrossline, t);
        
        V = data;
        
    else
%         minTime = min(min(startTime));
        % create horizons
        interval = 1 / params.scaleFactor;
        newDt = dt * interval;
        
        [x3, y3, z3] = meshgrid(firstInline : interval : endInline, 1 : interval : sampNum, firstCrossline : interval : endCrossline);
        [xx,yy,zz] = meshgrid(firstInline : endInline, 1 : sampNum, firstCrossline : endCrossline);
        newData = interp3(xx,yy,zz,data,x3,y3,z3);
        
        [xx,yy] = meshgrid(firstInline : endInline, firstCrossline : endCrossline);
        [x2,y2] = meshgrid(firstInline : interval : endInline, firstCrossline : interval : endCrossline);
        newStartTime = interp2(xx,yy,params.startTime',x2,y2);
        
        [V, minTime, maxTime, newSampNum] = bsHorizontalData(newData, newStartTime', newDt);
        
        
        t = minTime + (0:newSampNum-1)*newDt;
        t = t / 1000;
        [x3, y3, z3] = meshgrid(firstInline : interval : endInline, firstCrossline : interval : endCrossline, t);
    end
    
    
    %% calculate the horizons
    nHorizon = length(horizons);
    zHorizons = cell(1, nHorizon);
    
    for iHorizon = 1 : nHorizon
        zh.xd = zeros(nInline, nCrossline);
        zh.yd = zeros(nInline, nCrossline);
        zh.zd = zeros(nInline, nCrossline);
        
        for i = 1 : nInline
                for j = 1 : nCrossline
                    zh.xd(i, j) = i + firstInline - 1;
                    zh.yd(i, j) = j + firstCrossline - 1;
                    zh.zd(i, j) = (horizons(iHorizon).horizon(i, j) ...
                            + horizons(iHorizon).shift) / 1000;
                end
        end
            
        residual = horizons(iHorizon).horizon - params.startTime;
        if max(residual(:)) == min(residual(:))
            % 层位与起始时间一致
            zh.horizonSameAsStartime = 1;
            k = round((residual(1) + horizons(iHorizon).shift) / params.dt);
            zh.data = squeeze(data(k, :, :));
        else
            % 层位与起始时间不一致
            zh.horizonSameAsStartime = 0;
        end
        
        zHorizons{iHorizon} = zh;
    end
    
    %% start to plot the data
    p_V = permute(V, [3 2 1]);
    
%     newp_V = interp3(xx,yy,zz,p_V,x3,y3,z3);
    houts = slice(x3, y3, z3, p_V, xslices, yslices, zslices/1000, 'cubic'); hold on;
    horizonSlices = cell(1, nHorizon);
    
    for iHorizon = 1 : nHorizon
        xd = interp2(xx,yy,(zHorizons{iHorizon}.xd)',x2,y2);
        yd = interp2(xx,yy,(zHorizons{iHorizon}.yd)',x2,y2);
        zd = interp2(xx,yy,(zHorizons{iHorizon}.zd'),x2,y2);
        
        if zHorizons{iHorizon}.horizonSameAsStartime
%             surf(xd, yd, zd);
            tmp = repmat(zHorizons{iHorizon}.data', 1, 1, newSampNum);
            horizonSlices{iHorizon} = slice(x3, y3, z3, tmp, xd, yd, zd, 'cubic');
        else
            horizonSlices{iHorizon} = slice(x3, y3, z3, p_V, xd, yd, zd, 'cubic');
        end
        
    end
    
    %% show the titles of all crossed well
    for i = 1 : size(wellPos, 1)
        ix = wellPos(i, 1);
        iy = wellPos(i, 2);

        if ~isfield(welllogs{wellIndex(i)}, 'color')
        	welllogs{wellIndex(i)}.color = 'black';
        end
        
        if isscalar(params.startTime)
            text(ix, iy, (minTime - 0.02*(maxTime-minTime))/1000, welllogs{wellIndex(i)}.name, ...
            'color', welllogs{wellIndex(i)}.color, 'fontsize', params.fontsize,...
            'fontweight','bold', 'fontname', params.fontname);
        else
            text(ix, iy, params.startTime(ix, iy)/1000 - 0.02*(maxTime-minTime)/1000, welllogs{wellIndex(i)}.name, ...
            'color', welllogs{wellIndex(i)}.color, 'fontsize', params.fontsize,...
            'fontweight','bold', 'fontname', params.fontname);
        end
    end
    
    
    annotation(gcf,'textbox',[0.1 0.95 0.8 0.05],...
        'String',{params.attributeName},...
        'LineStyle','none',...
        'HorizontalAlignment', 'center', ...
        'FontWeight', params.fontweight,...
        'FontSize',  params.fontsize,...
        'FontName', params.fontname,...
        'FitBoxToText','off', 'EdgeColor',[0.94 0.94 0.94]);
    
%     title(params.attributeName);
    
    if params.isShading
        shading interp;
    else
        shading flat;
    end
    
%     for i = 1 : length(houts)
%         bsChangeColorData(houts(i), params.isShading, clim, params.colormap);
%     end
%     
%     for iHorizon = 1 : nHorizon
%         bsChangeColorData(horizonSlices{iHorizon}, params.isShading, clim, params.colormap);
%     end
    
%     shading interp;
    
    set(gca, 'clim', clim);
    set(gca, 'colormap', params.colormap);
%     colormap(params.colormap);
    
    
    colorbar('location', 'southoutside', 'position', [0.092 0.06 0.88 0.02]);
    
    set(gca,'zdir','reverse');
    xlabel('Inline', 'fontsize', params.fontsize,'fontweight', params.fontweight, 'fontname', 'Times New Roman');
    ylabel('Crossline', 'fontsize', params.fontsize,'fontweight', params.fontweight, 'fontname', 'Times New Roman');
    zlabel(params.zlabel_str, 'fontsize', params.fontsize,'fontweight', params.fontweight, 'fontname', params.fontname);
    set(gca , 'fontsize', params.fontsize,'fontweight', params.fontweight, 'fontname', params.fontname);    
    
    set(gca, 'zlim', [min(params.startTime(:)) max(params.startTime(:))+sampNum*dt]/1000);
    set(gca, 'xlim', [firstInline, firstInline+nInline-1]);
    set(gca, 'ylim', [firstCrossline, firstCrossline+nCrossline-1]);
    grid off;
%     axis([firstInline firstInline+nInline firstCrossline firstCrossline+nCrossline t(1)-dt/1000*5 t(end)+dt/1000*5]);
    view(params.view);
end

function bsChangeColorData(isurface, isShading, clim, colormap)

    if isShading
        nColor = size(colormap, 1);
        cinterval = (clim(2) - clim(1)) / nColor; 
        newLmin = clim(1) + 2*cinterval;
        d1 = isurface.CData;

        d1( find((d1 ~= inf) & (d1 < newLmin)) ) = newLmin;
        set( isurface, 'CData', d1);
    else
        set(isurface, 'AlphaDataMapping', 'direct');
        set(isurface, 'AlphaData', isurface.CData ~= 0);
        set(isurface, 'FaceColor', 'texturemap');
        set(isurface, 'FaceAlpha', 'texturemap');
    end
    

end

function [V, minTime, maxTime, newSampNum] = bsHorizontalData(data, startTime, dt)
    
    [sampNum, nInline, nCrossline] = size(data);
    
    minTime = floor((min(startTime(:)) - 2*dt)/10)*10;
    maxTime = max(startTime(:)) + 2*dt + sampNum * dt;
    
    newSampNum = round( (maxTime - minTime) / dt );
    V = zeros(newSampNum, nInline, nCrossline);
    V(:) = inf;
    
    for i = 1 : nInline
        for j = 1 : nCrossline
            
            index = round( (startTime(i, j) - minTime) / dt );  
            
            V(index+1 : index+sampNum, i, j) = data(:, i, j);
        end
    end
end