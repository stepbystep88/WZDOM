function bsShowWellLocations(GInvParam, GShowProfileParam, wellLogs)

    [rangeInline, rangeCrossline] = bsGetWorkAreaRangeByParam(GInvParam);
    nInline = rangeInline(end)-rangeInline(1)+1;
    nCrossline = rangeCrossline(end)-rangeCrossline(1)+1;
    
    figure;
    bsSetPosition(0.47, 0.45);
    rectangle('Position', [rangeInline(1) rangeCrossline(1)  nInline, nCrossline] ); hold on;
    
    GPlotParam = GShowProfileParam.plotParam;
    xlabel('Inline', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    ylabel('Crossline', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    
    title('Well locations', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    
    wells = cell2mat(wellLogs);
    inIds = [wells.inline];
    crossIds = [wells.crossline];
    names = {1, length(wellLogs)};
    for i = 1 : length(wellLogs)
        names{i} = sprintf('%d-%s', i, wellLogs{i}.name);
    end
    
    hold on;
    plot(inIds, crossIds, 'ko', 'linewidth', 2);
    text(inIds-0.03*nInline, crossIds-0.05*nCrossline, names, 'ko', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    
    grid on
    
    axis([rangeInline rangeCrossline]);
end