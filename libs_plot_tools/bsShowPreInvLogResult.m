function bsShowPreInvLogResult(GPreInvParam, GShowProfileParam, ...
    invVals, trueLogFiltcoef, isShowSynSeisComparasion)

    if ~exist('isShowSynSeisComparasion', 'var')
        isShowSynSeisComparasion = 1;
    end
    
    GPlotParam = GShowProfileParam.plotParam;
    
    nItems = length(invVals);
    model = invVals{1}.model;
    trueLog = model.trueLog;
    
    sampNum = size(model.trueLog, 1);
    t = (model.t0 : GPreInvParam.dt : model.t0 + (sampNum-1)*GPreInvParam.dt) / 1000;
    
    
    bsShowInvLog();
    
    if isShowSynSeisComparasion
        bsShowSyntheticSeisData();
    end
    
    function bsShowInvLog()
        hf = figure;
        bsSetPreFigureSize(nItems);
        
        for iItem = 1 : nItems
        
            figure(hf);
            

            invVal = invVals{iItem};
            if trueLogFiltcoef>0 && trueLogFiltcoef<1
                trueLog = bsFiltWelllog(trueLog, trueLogFiltcoef);
            end

            bsShowPostSubInvLogResult(GPlotParam, ...
                invVal.vp/1000, trueLog(:, 2)/1000, model.initLog(:, 2)/1000, ...
                t, invVal.name, 'V_P (km/s)', GShowProfileParam.range.vp/1000, ...
                nItems, iItem, 1);
            
            bsShowPostSubInvLogResult(GPlotParam, ...
                invVal.vs/1000, trueLog(:, 3)/1000, model.initLog(:, 3)/1000, ...
                t, invVal.name, 'V_S (km/s)', GShowProfileParam.range.vs/1000, ...
                nItems, iItem, 2);
            
            bsShowPostSubInvLogResult(GPlotParam, ...
                invVal.rho, trueLog(:, 4), model.initLog(:, 4), ...
                t, invVal.name, '\rho (g/cm^3)', GShowProfileParam.range.rho, ...
                nItems, iItem, 3);

        end

        legends = {'Initial model', 'Real model', 'Inversion result'};
        bsSetLegend(GPlotParam, {'g', 'b-.', 'r'}, legends);
       
    end

    function bsShowSyntheticSeisData()
        hf = figure;
        bsSetPreFigureSize(ceil((nItems+2)/3));
        angleData = model.angleData;
        if max(angleData < 5)
            angleData = angleData / pi * 180;
        end
        angleData = round(angleData);
        
        synFromTrue = reshape(model.G * model.trueX, sampNum-1, GPreInvParam.angleTrNum);
        seisData = reshape(model.original_d, sampNum-1, GPreInvParam.angleTrNum);
        
        bsShowPreSubSynSeisData(GPlotParam, 'real', seisData, t, angleData, nItems, 1);
        bsShowPreSubSynSeisData(GPlotParam, 'synthetic from welllog', synFromTrue, t, angleData, nItems, 2);
            
        for iItem = 1 : nItems
        
            figure(hf);
            

            invVal = invVals{iItem};
            
            invLog = [bsGetDepth(invVal.vp, GPreInvParam.dt), invVal.vp, invVal.vs, invVal.rho];
            x1 = bsPreBuildModelParam(invLog, GPreInvParam.mode, model.lsdCoef);
            
            synFromInv = reshape(model.G * x1, sampNum-1, GPreInvParam.angleTrNum);
            
            
            
            bsShowPreSubSynSeisData(GPlotParam, invVal.name, synFromInv, t, angleData, nItems, iItem+2);
            
%             bsShowPreSubSynSeisData(GPlotParam, 
%                 t, invVal.name, GShowProfileParam.range.seismic, nItems, iItem)
        end

%         legends = {'Real data', 'Synthetic from true welllog', 'Synthetic from inversion result'};
%         bsSetLegend(GPlotParam, {'k', 'b', 'r'}, legends, 'Seismic');
    end
end

function bsSetLegend(GPlotParam, colors, legends, attrName)
    hL = subplot('position', [0.25    0.01    0.500    0.03]);
    poshL = get(hL, 'position');     % Getting its position

    plot(0, 0, colors{1}, 'linewidth', GPlotParam.linewidth); hold on;
    plot(0, 0, colors{2}, 'LineWidth', GPlotParam.linewidth);   hold on;
    plot(0, 0, colors{3},'LineWidth', GPlotParam.linewidth);    hold on;

    lgd = legend(legends);

    set(lgd, 'Orientation', 'horizon', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    set(lgd, 'position', poshL);      % Adjusting legend's position
    axis(hL, 'off');                 % Turning its axis off
    
    if exist('attrName', 'var')
        annotation('textbox', [0.05, 0.07, 0, 0], 'string', attrName, ...
            'edgecolor', 'none', 'fitboxtotext', 'on', ...
            'fontsize', GPlotParam.fontsize,...
            'fontweight', 'bold', ...
            'fontname', GPlotParam.fontname);
    end
end

function bsSetFigureSize(nPlot)
    switch nPlot
        case 1
            bsSetPosition(0.13, 0.45);
        case 2
            bsSetPosition(0.56, 0.45);
        case 3
            bsSetPosition(0.56, 0.45);
        otherwise
            bsSetPosition(0.34, 0.45);
    end
end

function bsSetSubPlotSize(nItems, iItem)
    switch nItems
        case {1, 2, 3}
            bsSubPlotFit(1, nItems, iItem, 0.96, 0.9, 0.08, 0.11, 0.085, 0.045);
        case 4
            bsSubPlotFit(1, nItems, iItem, 0.96, 0.92, 0.08, 0.11, 0.085, 0.045);
    end
end

function bsSetPreFigureSize(nPlot)
    switch nPlot
        case 1
            bsSetPosition(0.44, 0.27);
        case 2
            bsSetPosition(0.44, 0.48);
        case 3
            bsSetPosition(0.44, 0.68);
        otherwise
            bsSetPosition(0.44, 0.88);
    end
end

function bsSetPreSubPlotSize(nItems, iItem, k)
    
    index = 3 * (iItem - 1) + k;
    switch nItems
        case 1
            bsSubPlotFit(nItems, 3, index, 0.93, 0.75, 0.02, 0.00, 0.06, -0.05);
        case 2
            bsSubPlotFit(nItems, 3, index, 0.93, 0.88, 0.02, 0.05, 0.06, 0.01);
        case 3
            bsSubPlotFit(nItems, 3, index, 0.93, 0.91, 0.02, 0.04, 0.06, 0.01);
        case 4
            bsSubPlotFit(nItems, 3, index, 0.93, 0.92, 0.02, 0.03, 0.06, 0.00);
    end
end

%% show comparasions of welllog inversion results
function bsShowPostSubInvLogResult(GPlotParam, ...
    invVal, trueVal, initVal, ...
    t, tmethod, attName, range, nItems, iItem, k)

    bsSetPreSubPlotSize(nItems, iItem, k);
    plot(initVal, t, 'g', 'linewidth', GPlotParam.linewidth); hold on;
    plot(trueVal, t, 'b-.', 'LineWidth', GPlotParam.linewidth);   hold on;
    plot(invVal, t, 'r','LineWidth', GPlotParam.linewidth);    hold on;
        
    if k == 1
        ylabel('Time (s)');
    else
        set(gca, 'ytick', [], 'yticklabel', []);
    end

    if k == 2
        title(sprintf('(%s) %s', char( 'a' + (iItem-1) ), tmethod), ...
            'fontsize', GPlotParam.fontsize+2, 'fontweight', 'bold', 'interpreter','latex');
    end
    
    if iItem == nItems
        xlabel(attName);
    else
        set(gca, 'xtick', [], 'xticklabel', []);
    end

    set(gca, 'ydir', 'reverse');
    
    if ~isempty(range)
        set(gca, 'xlim', range);
    end
    
    set(gca, 'ylim', [t(1), t(end)]);
    
    set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
end

%% show comparasions of synthetic seismic data and real seismic data
function bsShowPostSubSynSeisData(GPlotParam, ...
    synFromTrue, synFromInv, seisData, ...
    t, tmethod, range, nItems, iItem)

    bsSetSubPlotSize(nItems, iItem);
    
    plot(synFromTrue, t(1:end-1), 'b-.', 'linewidth', GPlotParam.linewidth); hold on;
    plot(seisData, t(1:end-1), 'k','LineWidth', GPlotParam.linewidth);    hold on;
    plot(synFromInv, t(1:end-1), 'r', 'LineWidth', GPlotParam.linewidth);   hold on;
    
    ylabel('Time (s)');
    set(gca,'ydir','reverse');
    
    title(sprintf('(%s) %s', char( 'a' + (iItem-1) ), tmethod));
    if ~isempty(range)
        set(gca, 'xlim', range) ; 
    end
    set(gca, 'ylim', [t(1) t(end)]);
%     bsTextSeqIdFit(ichar - 'a' + 1, 0, 0, 12);

    set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
end

function bsShowPreSubSynSeisData(GPlotParam, tmethod, data, t, angleData, nItems, iItem)

%     bsSetPreSubPlotSize(nItems, iItem, k);
    nRow = ceil((nItems+2)/3);
%     subplot(nRow, 3, iItem);
    bsSubPlotFit(nRow, 3, iItem, 0.93, 0.92, 0.02, 0.11, 0.06, -0.02);
    
    wiggles = {'k', 'b', 'r'};
    peaks = {'k', 'b', 'r'};
    
    seismic = s_convert(data, t(1), t(2)-t(1));
    s_wplot(seismic, {'figure', 'old'}, {'xaxis', angleData});%, ...
%         {'wiggle_color', wiggles{k}}, ...
%         {'peak_fill', peaks{k}});
    
    
    if mod(iItem, 3) == 1
        ylabel('Time (s)');
    else
        set(gca, 'ytick', [], 'yticklabel', []);
    end
    
    xlabel(sprintf('(%s) %s', char( 'a' + (iItem-1) ), tmethod), ...
            'fontsize', GPlotParam.fontsize+2, 'fontweight', 'bold');

    set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
end
