function bsShow3DInvResults(GInvParam, GShowProfileParam, invResult, iAtt, mode, shift, varargin)
    
    if isempty(iAtt)
        iAtt = 1;
    end
    
    if isempty(shift)
        shift = 0;
    end
    
%     if isempty(rangeInline) && isempty(rangeCrossline)
%         [rangeInline, rangeCrossline] = bsGetWorkAreaRangeByParam(GInvParam);
%     end
    
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
            xslices = rangeInline;
            yslices = rangeCrossline;
            zslices = [2202, 2200+(sampNum-3)*GInvParam.dt];
            horiozons = [];
            startTime = ones(nInline, nCrossline) * 2200;
            
        case 2
            xslices = rangeInline(1);
            yslices = rangeInline(1);
            startTime = permute(reshape(invResult.horizon, nCrossline, nInline), [2, 1]) - GInvParam.dt*GInvParam.upNum;
            startTime = bsSmoothByGST2D(startTime, [], bsCreateGSTParam(2, 'sigma', 10));
            
            zslices = [];
            horiozons(1).horizon = startTime;
            horiozons(1).shift = (sampNum - 10) * GInvParam.dt + shift;
            
%             view_s = [121.6586   38.2342];
        case 3
            xslices = round(0.5 * (rangeInline(1) + rangeInline(end)));
            yslices = round(0.5 * (rangeCrossline(1) + rangeCrossline(end)));
            zslices = 2200+sampNum/2*GInvParam.dt  + shift;
            horiozons = [];
            startTime = ones(nInline, nCrossline) * 2200;
        case 4
            xslices = round(0.5 * (rangeInline(1) + rangeInline(end)));
            yslices = round(0.5 * (rangeCrossline(1) + rangeCrossline(end)));
            startTime = permute(reshape(invResult.horizon, nCrossline, nInline), [2, 1]) - GInvParam.dt*GInvParam.upNum;
            startTime = bsSmoothByGST2D(startTime, [], bsCreateGSTParam(2, 'sigma', 10));
            
            zslices = [];
            horiozons(1).horizon = startTime;
            horiozons(1).shift = (sampNum - 10) * GInvParam.dt + shift;
            
        otherwise
            xslices = rangeInline(1);
            yslices = rangeInline(1);
            startTime = permute(reshape(invResult.horizon, nCrossline, nInline), [2, 1]) - GInvParam.dt*GInvParam.upNum;
            startTime = bsSmoothByGST2D(startTime, [], bsCreateGSTParam(2, 'sigma', 10));
            
            zslices = [];
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
        
    
end