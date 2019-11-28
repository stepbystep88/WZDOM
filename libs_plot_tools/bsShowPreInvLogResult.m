function bsShowPreInvLogResult(GPreInvParam, GPlotParam, GShowProfileParam, ...
    invVals, trueLogFiltcoef, isShowSynSeisComparasion)

    if ~exist('isShowSynSeisComparasion', 'var')
        isShowSynSeisComparasion = 1;
    end
    

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
        for iItem = 1 : nItems
        
            figure(hf);
            bsSetPreFigureSize(iItem);

            invVal = invVals{iItem};
            if trueLogFiltcoef>0 && trueLogFiltcoef<1
                trueLog = bsFiltWelllog(trueLog, trueLogFiltcoef);
            end

            bsShowPostSubInvLogResult(GPlotParam, ...
                invVal.vp/1000, trueLog(:, 2)/1000, model.initLog(:, 2)/1000, ...
                t, invVal.name, 'V_P (km/s)', GShowProfileParam.rangeVP/1000, ...
                nItems, iItem, 1);
            
            bsShowPostSubInvLogResult(GPlotParam, ...
                invVal.vs/1000, trueLog(:, 3)/1000, model.initLog(:, 3)/1000, ...
                t, invVal.name, 'V_S (km/s)', GShowProfileParam.rangeVS/1000, ...
                nItems, iItem, 2);
            
            bsShowPostSubInvLogResult(GPlotParam, ...
                invVal.rho, trueLog(:, 4), model.initLog(:, 4), ...
                t, invVal.name, '\rho (g/cm^3)', GShowProfileParam.rangeRho, ...
                nItems, iItem, 3);

        end

        legends = {'Initial model', 'Real model', 'Inversion result'};
        bsSetLegend(GPlotParam, {'g', 'b-.', 'r'}, legends);
       
    end

    function bsShowSyntheticSeisData()
        hf = figure;
        for iItem = 1 : nItems
        
            figure(hf);
            bsSetFigureSize(iItem);

            invVal = invVals{iItem};
            
            G = bsPostGenGMatrix(GPreInvParam.wavelet, sampNum);
            synFromInv = G * log(invVal.vp .* invVal.vs);
            synFromTrue = G * log(trueLog(:, 2) .* trueLog(:, 4));
            seisData = bsGetPostSeisData(GPreInvParam, invVal.inline, invVal.crossline,...
                model.t0, sampNum-1);
            
            bsShowPostSubSynSeisData(GPlotParam, ...
                synFromTrue, synFromInv, seisData, ...
                t, invVal.name, ...
                GShowProfileParam.rangeSeismic, nItems, iItem);

        end

        legends = {'Synthetic from true welllog', 'Real data', 'Synthetic from inversion result'};
        bsSetLegend(GPlotParam, {'b-.', 'k', 'r'}, legends, 'Seismic');
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
            set(gcf, 'position', [  336   240   247   483]);
        case 2
            set(gcf, 'position', [  336 240 1080 483]);
        case 3
            set(gcf, 'position', [  336 240 1080 483]);
        otherwise
            set(gcf, 'position', [687   134   658   543]);
    end
end

function bsSetSubPlotSize(nItems, iItem)
    switch nItems
        case {1, 2, 3}
            bsSubPlotFit(1, nItems, iItem, 0.96, 0.92, 0.08, 0.11, 0.085, 0.045);
        case 4
            bsSubPlotFit(1, nItems, iItem, 0.96, 0.92, 0.08, 0.11, 0.085, 0.045);
    end
end

function bsSetPreFigureSize(nPlot)
    switch nPlot
        case 1
            set(gcf, 'position', [  336    98   847   295]);
        case 2
            set(gcf, 'position', [  336    98   847   521]);
        case 3
            set(gcf, 'position', [   336    98   847   733]);
        otherwise
            set(gcf, 'position', [336    42   847   953]);
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
            'fontsize', GPlotParam.fontsize+2, 'fontweight', 'bold');
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
