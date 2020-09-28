function bsShowSliceOf3DData(GInvParam, GShowProfileParam, invResult, iAtt, varargin)
    
    GPlotParam = GShowProfileParam.plotParam;

    if isempty(iAtt)
        iAtt = 1;
    end
    
    p = inputParser;
    
    addParameter(p, 'avg_direction', 'center');
    addParameter(p, 'avg_num', 10);
    addParameter(p, 'avg_type', 'RMS');
    addParameter(p, 'shift', 0);
    addParameter(p, 'clim', []);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    if ~iscell(invResult.type)
        type = invResult.type;
        data = invResult.data;
    else
        type = invResult.type{iAtt};
        data = invResult.data{iAtt};
    end
    
    [range, scale, dataIndex, dataColorTbl, attName] = ...
            bsGetInfoByType(GShowProfileParam, GInvParam, type);
    
    rangeInline = [min(invResult.inIds), max(invResult.inIds)];
    rangeCrossline = [min(invResult.crossIds), max(invResult.crossIds)];

    nCrossline = rangeCrossline(2) - rangeCrossline(1) + 1;
    nInline = rangeInline(2) - rangeInline(1) + 1;

    data3D = bsReshapeDataAs3D(data, nInline, nCrossline);
    
    avg_num = options.avg_num;
    pos = GInvParam.upNum + round(options.shift / GInvParam.dt);
    switch options.avg_direction
        case 'center'
            avg_num = floor(avg_num/2);
            subData = data3D(pos-avg_num:pos+avg_num, :, :);
        case 'above'
            subData = data3D(pos-avg_num+1:pos, :, :);
        case 'below'
            subData = data3D(pos:pos+avg_num-1, :, :);
    end
    
    switch lower(options.avg_type)
        case 'rms'
            avg_slice = rms(subData, 1);
        case 'mean'
            avg_slice = mean(subData, 1);
        case 'constant'
            avg_slice = data3D(pos, :, :);
        case 'max'
            avg_slice = max(subData, 1);
        case 'min'
            avg_slice = min(subData, 1);
    end
    
    slice_data = squeeze(avg_slice);
    imagesc(rangeCrossline, rangeInline, slice_data/scale);
    set(gca,'DataAspectRatio',[1 1 1]);
    
    % 获取需要旋转的角度
    angle = bsGetAngelOfCoordinates(GInvParam, rangeInline, rangeCrossline);
    view([angle, 90]);
    
    title(sprintf('Offset [%d %d]ms', options.shift-avg_num*GInvParam.dt, options.shift+avg_num*GInvParam.dt), 'fontweight', 'bold');
    if isempty(options.clim)
        clim = [prctile(slice_data(:), 5), prctile(slice_data(:), 95)];
    else
        clim = options.clim;
    end
    
    set(gca, 'colormap', dataColorTbl);
    set(gca, 'clim', clim/scale);
    
    ylabel('Inline');
    xlabel('Crossline');
    
    hc = colorbar( 'EastOutside', ...
            'fontsize', GPlotParam.fontsize,...
            'fontweight', 'bold', ...
            'fontname', GPlotParam.fontname);
    hc.Position = hc.Position .* [1.0 1 0.5 1] + [0.05 0 0 0];

    ylabel(hc, attName, ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    
    bsSetDefaultPlotSet(GPlotParam);
    
    
end


function angle = bsGetAngelOfCoordinates(GInvParam, rangeInline, rangeCrossline)
    inIds = [rangeInline(1), rangeInline(1), rangeInline(end), rangeInline(end)];
    crossIds = [rangeCrossline(1), rangeCrossline(end), rangeCrossline(1), rangeCrossline(end)];
    
    xs = zeros(1, 4);
    ys = zeros(1, 4);
    
    [~, ~, fileName, GSegyInfo] ...
        = bsGetWorkAreaRangeByParam(GInvParam);
    
    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
    sizeTrace = GSegyInfo.volHeader.sizeTrace+240;    
    
    for i = 1 : 4
        index = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inIds(i), crossIds(i));
        fseek(GSegyInfo.fid, 3600+sizeTrace*(index-1), -1);
        trHeader = bsReadTraceHeader(GSegyInfo);
            
        xs(i) = GSegyInfo.getX(trHeader.X);
        ys(i) = GSegyInfo.getY(trHeader.Y);
    end
    
    v1 = [xs(2) ys(2)] - [xs(1) ys(1)];
    v2 = [1 0];
    
    angle = acos(dot(v1, v2)/(norm(v1)*norm(v2)));
    angle = angle/pi*180;
end

