function bsShow3DSlice(GInvParam, GShowProfileParam, wellLogs, data, offset, varargin)

    [~, nInline, nCrossline] = size(data);
    
    p = inputParser;
    
    addParameter(p, 'firstInlineId', 1, @(x) (isscalar(x)&&(x>0)) );
    addParameter(p, 'firstCrosslineId', 1, @(x) (isscalar(x)&&(x>0)) );
    addParameter(p, 'colormap', bsGetColormap('velocity') );
    addParameter(p, 'avg_direction', 'center');
    addParameter(p, 'avg_num', 1);
    addParameter(p, 'avg_type', 'RMS');
    addParameter(p, 'title', '');
    addParameter(p, 'type', 'ip');
    addParameter(p, 'isNLM', 0);
    
    p.parse(varargin{:});  
    options = p.Results;
    
    [range, scale, ~, dataColorTbl, attName] = ...
            bsGetInfoByType(GShowProfileParam, GInvParam, options.type);
        
    avg_num = options.avg_num;
    pos = GInvParam.upNum + round(offset / GInvParam.dt);
    
    switch options.avg_direction
        case 'center'
            avg_num = floor(avg_num/2);
            subData = data(pos-avg_num:pos+avg_num, :, :);
        case 'above'
            subData = data(pos-avg_num+1:pos, :, :);
        case 'below'
            subData = data(pos:pos+avg_num-1, :, :);
    end
    
    switch lower(options.avg_type)
        case 'rms'
            avg_slice = rms(subData, 1);
        case 'mean'
            avg_slice = mean(subData, 1);
        case 'constant'
            avg_slice = data(pos, :, :);
    end
    
    avg_slice = reshape(avg_slice, nInline, nCrossline);
    if options.isNLM
        avg_slice = bsNLMByRef(avg_slice, avg_slice);
    end
    
    figure;
    bsSetPosition(0.47, 0.45);
    bsSubPlotFit(1, 1, 1, 0.9, 0.92, 0.01, 0.06, 0.07, 0);
    y = options.firstCrosslineId:options.firstCrosslineId+nCrossline-1;
    x = options.firstInlineId:options.firstInlineId+nInline-1;
    
    imagesc(y, x, avg_slice'/scale);
    set(gca, 'colormap', dataColorTbl);
    
    if ~isempty(range)
        set(gca, 'clim', range/scale);
    end
    
    
    GPlotParam = GShowProfileParam.plotParam;
    bsSetDefaultPlotSet(GShowProfileParam.plotParam);
    
    if ~isempty(options.title)
        title(options.title, ...
            'fontsize', GPlotParam.fontsize, ...
            'fontweight', 'bold', ...
            'fontname', GPlotParam.fontname);
    end
    
    ylabel('Inline', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    xlabel('Crossline', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    
    
    
    hc = colorbar( ...
            'fontsize', GPlotParam.fontsize,...
            'fontweight', 'bold', ...
            'fontname', GPlotParam.fontname);
    ylabel(hc, attName, ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    
    wells = cell2mat(wellLogs);
    inIds = [wells.inline];
    crossIds = [wells.crossline];
    names = {wells.name};
    
    hold on;
    plot(crossIds, inIds, 'ko', 'linewidth', 2);
    text(crossIds-0.01*nCrossline, inIds-0.03*nInline, names, 'ko', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
end