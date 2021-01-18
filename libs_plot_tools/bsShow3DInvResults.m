function bsShow3DInvResults(GInvParam, GShowProfileParam, invResult, iAtt, options, varargin)
    
    if isempty(iAtt)
        iAtt = 1;
    end
    
    [~, options] = bsGetFieldsWithDefaults(options, {'mode', 2; 'shift', 0; 'xslices', []; 'yslices', []; 'zslices', []});
    
    mode = options.mode;
    shift = options.mode;
    xslices = options.xslices;
    yslices = options.yslices;
    zslices = options.zslices;
    
    rangeInline = [min(invResult.inIds), max(invResult.inIds)];
    rangeCrossline = [min(invResult.crossIds), max(invResult.crossIds)];
    
    nCrossline = rangeCrossline(2) - rangeCrossline(1) + 1;
    nInline = rangeInline(2) - rangeInline(1) + 1;
    
    if ~iscell(invResult.type)
        type = invResult.type;
    else
        type = invResult.type{iAtt};
    end
    
    [range, scale, dataIndex, dataColorTbl, attName] = ...
            bsGetInfoByType(GShowProfileParam, GInvParam, type);
    
	if ~iscell(invResult.data)
        data = invResult.data;
    else
        data = invResult.data{iAtt};
    end
    
    data = bsReshapeDataAs3D(data, nInline, nCrossline);
    data = data / scale;
    sampNum = size(data, 1);
    
    switch mode
        case 1
            if isempty(xslices)
                xslices = rangeInline;
            end
            if isempty(yslices)
                yslices = rangeCrossline;
            end
            if isempty(zslices)
                zslices = [2202, 2200+(sampNum-3)*GInvParam.dt];
            end
            
            horiozons = [];
            startTime = ones(nInline, nCrossline) * 2200;
            
        case 2
            if isempty(xslices)
                xslices = rangeInline(1);
            end
            if isempty(yslices)
                yslices = rangeCrossline(1);
            end
            
            startTime = permute(reshape(invResult.horizon, nCrossline, nInline), [2, 1]) - GInvParam.dt*GInvParam.upNum;
            startTime = bsSmoothByGST2D(startTime, [], bsCreateGSTParam(2, 'sigma', 10));
            
            zslices = [];
            horiozons(1).horizon = startTime;
            horiozons(1).shift = (sampNum - 10) * GInvParam.dt + shift;
            
%             view_s = [121.6586   38.2342];
        case 3
            if isempty(xslices)
                xslices = round(0.5 * (rangeInline(1) + rangeInline(end)));
            end
            if isempty(yslices)
                yslices = round(0.5 * (rangeCrossline(1) + rangeCrossline(end)));
            end
            if isempty(zslices)
                zslices = 2200+sampNum/2*GInvParam.dt  + shift;
            end
            
            horiozons = [];
            startTime = ones(nInline, nCrossline) * 2200;
        case 4
            if isempty(xslices)
                xslices = round(0.5 * (rangeInline(1) + rangeInline(end)));
            end
            if isempty(yslices)
                yslices = round(0.5 * (rangeCrossline(1) + rangeCrossline(end)));
            end
            if isempty(zslices)
                zslices = [];
            end
            
            startTime = permute(reshape(invResult.horizon, nCrossline, nInline), [2, 1]) - GInvParam.dt*GInvParam.upNum;
            startTime = bsSmoothByGST2D(startTime, [], bsCreateGSTParam(2, 'sigma', 10));
            
            horiozons(1).horizon = startTime;
            horiozons(1).shift = (sampNum - 10) * GInvParam.dt + shift;
           
        case 5 
            minTime = 1;
            if isempty(xslices)
                xslices = [rangeInline(1), round((rangeInline(1)*0.6 + rangeInline(2)*0.4))];
            end
            if isempty(yslices)
                yslices = rangeCrossline(1);
            end
            if isempty(zslices)
                zslices = [minTime+shift*GInvParam.dt, minTime+(sampNum-3)*GInvParam.dt];
            end
            
            GInvParam.dt = 1000;

            horiozons = [];
            startTime = ones(nInline, nCrossline) * minTime;
            
        case 6
            minTime = 1;
            if isempty(xslices)
                xslices = [round((rangeInline(1)*0.2 + rangeInline(2)*0.8)), rangeInline(2)];
            end
            if isempty(yslices)
                yslices = rangeCrossline(1);
            end
            if isempty(zslices)
                zslices = [minTime+shift*GInvParam.dt, minTime+(sampNum-3)*GInvParam.dt];
            end
            
            GInvParam.dt = 1000;
            
            horiozons = [];
            startTime = ones(nInline, nCrossline) * minTime;
        
        case 8
            if isempty(xslices)
                xslices = rangeInline(1) + [round(0.2*nInline), round(0.8*nInline)];
            end
            if isempty(yslices)
                yslices = rangeCrossline(1) + [round(0.2*nCrossline), round(0.8*nCrossline)];
            end
            if isempty(zslices)
                zslices = [];
            end
            
            startTime = permute(reshape(invResult.horizon, nCrossline, nInline), [2, 1]) - GInvParam.dt*GInvParam.upNum;
            startTime = bsSmoothByGST2D(startTime, [], bsCreateGSTParam(2, 'sigma', 10));
            
            horiozons(1).horizon = startTime;
            horiozons(1).shift = (sampNum - 10) * GInvParam.dt + shift;
            
        otherwise
            if isempty(xslices)
                xslices = rangeInline(1);
            end
            if isempty(yslices)
                yslices = rangeCrossline(1);
            end

            startTime = permute(reshape(invResult.horizon, nCrossline, nInline), [2, 1]) - GInvParam.dt*GInvParam.upNum;
            startTime = bsSmoothByGST2D(startTime, [], bsCreateGSTParam(2, 'sigma', 10));
            horiozons(1).horizon = startTime + GInvParam.dt*GInvParam.upNum;
            horiozons(1).shift = shift;
    end
    
    
    bsShow3DVolume(data, GInvParam.dt, range/scale, xslices, yslices, zslices, horiozons, ...
        'startTime', startTime, ...
        'firstInlineId', rangeInline(1), ...
        'firstCrosslineId', rangeCrossline(1), ...
        'colormap', dataColorTbl, ...
        'fontname', GShowProfileParam.plotParam.fontname,...
        'fontsize', GShowProfileParam.plotParam.fontsize,...
        'fontweight', GShowProfileParam.plotParam.fontweight,...
        'attributeName', attName, ...
        'smoothFiltCoef', 1, varargin{:});
   
    if mode == 5
        hold on;
        % 设置透明度
        houts = bsShow3DVolume(data, GInvParam.dt, range/scale, rangeInline(2), rangeCrossline(2), minTime+GInvParam.dt*2, [], ...
            'startTime', startTime, ...
            'firstInlineId', rangeInline(1), ...
            'firstCrosslineId', rangeCrossline(1), ...
            'colormap', dataColorTbl, ...
            'fontname', GShowProfileParam.plotParam.fontname,...
            'fontsize', GShowProfileParam.plotParam.fontsize,...
            'fontweight', GShowProfileParam.plotParam.fontweight,...
            'attributeName', attName, ...
            'smoothFiltCoef', 1, varargin{:});
        
%         tc = houts{1}.
        for i = 1 : length(houts)
            CData = houts(i).CData;
            AlphaData = ones(size(CData))*1;
            if size(CData, 1) == nCrossline && size(CData, 2) == nInline
                % 切片
                AlphaData(:) = 0.4;
                AlphaData(:, xslices(2)-rangeInline(1):end) = 0;
            elseif size(CData, 1) == nInline 
                AlphaData(xslices(2)-rangeInline(1):end, 1:shift+1) = 0;
            elseif size(CData, 1) == nCrossline
                AlphaData(:, 1:shift+1) = 0;
            end
            
%             houts(i).CData = CData;
            houts(i).FaceColor = 'texturemap';
            houts(i).FaceAlpha = 'texturemap';
            houts(i).AlphaData = AlphaData;
%             set( houts(i),'FaceAlpha',  'texturemap', 'AlphaDataMapping', 'none', 'AlphaData',AlphaData);
            
        end
    elseif mode == 6
        hold on;
        % 设置透明度
        houts = bsShow3DVolume(data, GInvParam.dt, range/scale, rangeInline(1), rangeCrossline(2), minTime+GInvParam.dt, [], ...
            'startTime', startTime, ...
            'firstInlineId', rangeInline(1), ...
            'firstCrosslineId', rangeCrossline(1), ...
            'colormap', dataColorTbl, ...
            'fontname', GShowProfileParam.plotParam.fontname,...
            'fontsize', GShowProfileParam.plotParam.fontsize,...
            'fontweight', GShowProfileParam.plotParam.fontweight,...
            'attributeName', attName, ...
            'smoothFiltCoef', 1, varargin{:});
        
%         tc = houts{1}.
        for i = 1 : length(houts)
            CData = houts(i).CData;
            AlphaData = ones(size(CData))*1;
            if size(CData, 1) == nCrossline && size(CData, 2) == nInline
                % 切片
                AlphaData(:) = 1;
                AlphaData(:, 1:xslices(1)-rangeInline(1)) = 0;
            elseif size(CData, 1) == nInline 
                AlphaData(1:xslices(1)-rangeInline(1), 1:shift+1) = 0;
            elseif size(CData, 1) == nCrossline
                AlphaData(:, 1:shift+1) = 0;
            end
            
%             houts(i).CData = CData;
            houts(i).FaceColor = 'texturemap';
            houts(i).FaceAlpha = 'texturemap';
            houts(i).AlphaData = AlphaData;
%             set( houts(i),'FaceAlpha',  'texturemap', 'AlphaDataMapping', 'none', 'AlphaData',AlphaData);
            
        end
    end
    
    
end