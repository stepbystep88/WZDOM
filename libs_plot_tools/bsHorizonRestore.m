function bsHorizonRestore(profile, options)

    GPlotParam = options.plotParam;
    
    %%
    %修正井数据
    wellData = options.wellData;
    wellPos = options.wellPos;
    [~, wellNum] = size(wellData);
    
    for i = 1 : wellNum
        sp = wellPos(i) - options.coefShowWell;
        ep = wellPos(i) + options.coefShowWell;
        
        for j = sp : ep
            profile(:, j) = wellData(:, i);
        end
    end
    
    %%
    % firstCDP, basetime, options.upNum, downNum, options.dt, timeUp, timeCenter, timeDown, wellPos, txt, cond, newWellData
    [~, traceNum] = size(profile); 
    row = options.upNum + options.downNum;
    
    sequence = linspace(0, 0 + options.dt*(row-1), row);       % 生成row行的数据
    
    time0 = options.horizon - options.upNum * options.dt;
%     traceNum = length(options.baseTime);
    
    if( size(time0, 2) == 1)
        time0 = time0';
    end
    
    time = repmat(sequence', 1, traceNum)...
        + repmat(time0, row, 1);                                % 计算出剖面的每一个点的时间
    
    t0 = min(time(1,:)) - options.dt * 2;                       % 绘图的起始时间
    interval = options.dt / 1.5;                                  % 设置插值的时间间隔

    maxoptions.dtime = max(options.horizon) + options.downNum * options.dt;
    plotRow = floor((maxoptions.dtime - t0) / interval * 1.03); % 绘图的行数
    
    horizon = zeros(plotRow, traceNum);
    horizon(:, :) = inf;
    for j = 1 : traceNum
        for i = 1 : size(horizon, 1)
            t = t0 + (i-1) * interval;      % 第i层的时间
            if t >= time(1, j) && t<= time(row, j)  
                 horizon(i,j) = bsCalVal(t, time(:,j), profile(:,j));
            end
        end
    end
    
    %层位
    min_val = min(min(horizon));
    min_val = min_val - abs(min_val)*10;
    horizon(isnan(horizon )) = min_val;
    horizon(isinf(horizon )) = min_val;
    
   
    x_horizon = (1 : traceNum)';
        
    if(~isempty(options.limits))
        % 计算倒数第二个位置的数据索引
        lmin = options.limits{2};
        lmax = options.limits{3};
        
        if isempty(options.colormap)
            load attColor.mat;
            nColor = size(attColor, 1);
        else
            %set(gcf, 'colormap', options.colormap);
            nColor = size(options.colormap, 1);
        end
    
        cinterval = (lmax - lmin) / nColor;  % 颜色间隔是64个
        lmin = lmin + cinterval + 1e-5; % + 0.1 * interval;
        horizon( (min_val < horizon) & (horizon < lmin) ) = lmin;
%         s_cplot(horizon, {'interpol','pchip'}, {'time_lines',[]}, options.limits, options.figure), hold on;%层位还原
        imagesc(horizon);
        set(gca, 'clim', [options.limits{2} options.limits{3}]);
        set(gcf, 'colormap', options.colormap);
        colorbar;
    else
%         s_cplot(horizon, {'interpol','pchip'}, {'time_lines',[]}), hold on;%层位还原
    end
    
    t0 = floor(t0);
    
    %% 设置坐标轴
    [label, data] = bsGenLabel(t0, t0+(plotRow)*interval, plotRow, options.ylabelNum);
    set(gca,'Ytick', label ) ; 
    data = floor(data / 10) / 100;
    set(gca, 'YtickLabel', data , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
    
    [label, data] = bsGenLabel(0+options.firstCDP, traceNum+options.firstCDP, traceNum, options.xlabelNum);
    set(gca,'Xtick', label ) ; 
    set(gca, 'XtickLabel', data, 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
    
   
    %% 绘制层位线
    if(~isempty(options.xlabel))
%         xlabel(options.xlabel, 'fontsize', GPlotParam.fontsize,'fontweight', 'bold' , 'fontname', GPlotParam.fontname); 
        xlabel('');
%         set(gca, 'xlable'
        title(options.xlabel, 'fontsize', GPlotParam.fontsize+1,'fontweight', 'bold' , 'fontname', GPlotParam.fontname); 
    else
        set(get(gca, 'XLabel'), 'String', []);
    end
    ylabel('Time(s)', 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);   
%     set(get(gca,'Xlabel'),'String', 'Vp/(m/s)','fontsize',options.fontsize);
%     set(get(gca,'Xlabel'),'String',{'\bf{P-Impedence}', '\bf{(g/cm^3\bullet km/s})'});
    
    if( options.isShowHorizon)
        y = 1 + (options.horizon - t0) / interval;
        plot(x_horizon, y, 'k-','LineWidth', GPlotParam.linewidth); hold on;
    end
    
%     hold on,plot(x_horizon,ones(traceNum,1)*(plotRow),'k-');
    
    %% 设置标题
    if ~isempty(options.title)
%         title(options.title, 'position',[40, 30], 'FontSize', options.fontsize);
%         title(options.title, 'FontSize', options.fontsize);
        text(12, 10, options.title, 'fontsize', GPlotParam.fontsize+3,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
    end
    
    %% 设置色标
    if isempty(options.colormap)
        load attColor;
        set(gcf, 'colormap', attColor);
    else
        set(gcf, 'colormap', options.colormap);
    end
    
    if isempty(options.wellPos) || isempty(options.wellData)
        return;
    end
    
%     hold on;
% %     plot(wellPos, 2*ones(1, wellNum), 'kv', 'linewidth', 2); hold on;
%     plot(wellPos, plotRow*ones(1, wellNum), 'k^', 'linewidth', 2); 
    
    return;
    %% 绘制井数据
    
    
    [~, wellNum] = size(wellData);
    
    newWellData = zeros(row, wellNum);
    well_plot = zeros(plotRow, wellNum);
    
    for k = 1 : wellNum
%         newWellData(:,k) = ButtLowPassFilter(newWellData(:,k), 0.4);
        
        newWellData(:,k) = wellData(end-options.upNum-options.downNum+1:end, k);
        
        minIndex = -1;
        maxIndex = -1;
        
        for i = 1 : plotRow
            t = t0 + (i-1) * interval;
            
            if(minIndex == -1 && t >= time(1, wellPos(k)))
                minIndex = i;
            end
            
            if(maxIndex == -1 && t > time(row, wellPos(k)))
                maxIndex = i - 1;
            end
            
            if ( (t >= time(1, wellPos(k))) && (t<= time(row, wellPos(k))) )
                well_plot(i, k) = bsCalVal(t,time(:,wellPos(k)), newWellData(:,k));
            end
            
        end
        
        if(minIndex > 1)
            well_plot(1 : minIndex-1, k) = well_plot(minIndex, k);
        end

        if(maxIndex < plotRow)
            well_plot(maxIndex+1 : plotRow, k) = well_plot(maxIndex, k);
        end
    end
    
%     lwell = options.coefShowWell * traceNum;
%     lwell = 0.05 * traceNum;
    
%     for i = 1 : wellNum
%         num = find(horizon(:, wellPos(i))~= minData);
%         newWellData = well_plot(num, i);
%         [x_well, ~] = plot_well(newWellData, wellPos(i), lwell);
%         y_well = num;
%         hold on; plot(x_well,y_well,'k-','LineWidth', GPlotParam.linewidth) ;%画井
%         x_well(:) = wellPos(i);
%         hold on; plot(x_well,y_well,'k--','LineWidth', GPlotParam.linewidth) ;%画井
%     end

end