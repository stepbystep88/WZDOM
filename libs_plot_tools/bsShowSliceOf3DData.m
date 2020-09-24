function bsShowSliceOf3DData(GInvParam, GShowProfileParam, invResult, iAtt, clim, upOffset, downOffset)
    
    GPlotParam = GShowProfileParam.plotParam;

    if isempty(iAtt)
        iAtt = 1;
    end
    
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
    
    slice_data = squeeze((mean(data3D((upOffset:downOffset) + GInvParam.upNum, :, :), 1)));
    imagesc(rangeCrossline, rangeInline, slice_data/scale);
    set(gca,'DataAspectRatio',[1 1 1]);
    
    % 获取需要旋转的角度
    angle = bsGetAngelOfCoordinates(GInvParam, rangeInline, rangeCrossline);
    view([angle, 90]);
    
    title(sprintf('Offset [%d %d]ms', upOffset*GInvParam.dt, downOffset*GInvParam.dt), 'fontweight', 'bold');
    if isempty(clim)
        clim = [prctile(slice_data(:), 5), prctile(slice_data(:), 95)];
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
    
    GSegyInfo = GInvParam.postSeisData.segyInfo;
    fileName = GInvParam.postSeisData.fileName;
    
    GSegyInfo = bsReadVolHeader(fileName, GSegyInfo);        
    sizeTrace = GSegyInfo.volHeader.sizeTrace+240;    
    
    for i = 1 : 4
        index = bsIndexOfTraceSetOnInIdAndCrossId(GSegyInfo, inIds(i), crossIds(i));
        fseek(GSegyInfo.fid, 3600+sizeTrace*(index-1), -1);
        trHeader = bsReadTraceHeader(GSegyInfo);
            
        xs(i) = trHeader.X/1e4;
        ys(i) = trHeader.Y/1e4;
    end
    
    v1 = [xs(2) ys(2)] - [xs(1) ys(1)];
    v2 = [1 0];
    
    angle = acos(dot(v1, v2)/(norm(v1)*norm(v2)));
    angle = angle/pi*180;
end

