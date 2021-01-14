function bsShowPostInvLogResult(GPostInvParam, GShowProfileParam, invVals, trueLogFiltcoef, isShowSynthetic)

    if ~exist('isShowSynthetic', 'var')
        isShowSynthetic = 0;
    end
    
    GPlotParam = GShowProfileParam.plotParam;
    
    nItems = length(invVals);
    model = invVals{1}.model;
    
    sampNum = size(model.trueLog, 1);
    t = (model.t0 : GPostInvParam.dt : model.t0 + (sampNum-1)*GPostInvParam.dt) / 1000;
    
    
    bsShowInvLog();

    if isShowSynthetic
        bsShowSyntheticSeisData();
    end
    
    function bsShowInvLog()
        hf = figure;
        bsSetFigureSize(nItems);
        
        if strcmpi(GShowProfileParam.language, 'en')
            attName = 'I_{\it{P}} (g/cc*km/s)';
            legends = {'Initial model', 'True model', 'Inversion result'};
            
        else
            attName =  '阻抗 (g/cc*km/s)';
            legends = {'初始模型', '真实模型', '反演结果'};
        end
        
        for iItem = 1 : nItems
        
            figure(hf);
            

            invVal = invVals{iItem};
            trueLog = model.trueLog;

            if trueLogFiltcoef>0 && trueLogFiltcoef<1
                trueLog = bsButtLowPassFilter(trueLog, trueLogFiltcoef);
            end

            bsShowPostSubInvLogResult(GShowProfileParam, ...
                invVal.Ip/1000, trueLog/1000, model.initLog/1000, ...
                t, invVal.name, ...
                GShowProfileParam.range.ip/1000, nItems, iItem, attName);
            
            MRRMSE = bsCalcRRSE(trueLog, model.initLog, invVal.Ip);
            fprintf('[RRMSE=%.3f]\n', MRRMSE);
        

        end
        
        
        bsSetLegend(GPlotParam, {'g', 'k-.', 'r'}, legends, attName);
        
       
    end

    function bsShowSyntheticSeisData()
        hf = figure;
        for iItem = 1 : nItems
        
            figure(hf);
            bsSetFigureSize(iItem);

            invVal = invVals{iItem};
            model = invVal.model;

            G = model.orginal_G;
            
%             G = bsPostGenGMatrix(GPostInvParam.wavelet, sampNum);
            synFromInv = G * log(invVal.Ip);
            synFromTrue = G * log(model.trueLog);
            seisData = model.dTrue;
            
            bsShowPostSubSynSeisData(GShowProfileParam, ...
                synFromTrue, synFromInv, seisData, ...
                t, invVal.name, ...
                GShowProfileParam.range.seismic, nItems, iItem);

        end
        
        if strcmpi(GShowProfileParam.language, 'en')
            legends = {'Synthetic from true welllog', 'Real data', 'Synthetic from inversion result'};
        else
            legends = {'基于测井的合成记录', '实际观测记录', '基于反演结果的合成记录'};
        end
        
        bsSetLegend(GPlotParam, {'b-.', 'k', 'r'}, legends, '');
    end
end

function bsSetLegend(GPlotParam, colors, legends, attrName)
    hL = subplot('position', [0.25    0.02    0.500    0.04]);
    poshL = get(hL, 'position');     % Getting its position

    plot(0, 0, colors{1}, 'linewidth', GPlotParam.linewidth); hold on;
    plot(0, 0, colors{2}, 'LineWidth', GPlotParam.linewidth);   hold on;
    plot(0, 0, colors{3},'LineWidth', GPlotParam.linewidth);    hold on;

    lgd = legend(legends);

    set(lgd, 'Orientation', 'horizon', ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', GPlotParam.fontweight, ...
        'fontname', GPlotParam.fontname);
    set(lgd, 'position', poshL);      % Adjusting legend's position
    axis(hL, 'off');                 % Turning its axis off
%     annotation('textbox', [0.05, 0.07, 0, 0], 'string', attrName, ...
%         'edgecolor', 'none', 'fitboxtotext', 'on', ...
%         'fontsize', GPlotParam.fontsize,...
%         'fontweight', GPlotParam.fontweight, ...
%         'fontname', GPlotParam.fontname);
end

function bsSetFigureSize(nPlot)
    switch nPlot
        case {1, 2, 3}
%             set(gcf, 'position', [  336   240   509   406]);
            set(gcf, 'position', [  336   240   818   406]);
            
        case 4
            set(gcf, 'position', [  336   240   678   406]);
        case 5
            set(gcf, 'position', [  336   240   636   406]);
        otherwise
            set(gcf, 'position', [687   134   658   543]);
    end
end

function bsSetSubPlotSize(nItems, iItem)
    switch nItems
        case {1, 2, 3}
            bsSubPlotFit(1, nItems, iItem, 0.9, 0.91, 0.05, 0.17, 0.12, 0.1);
        otherwise
            bsSubPlotFit(1, nItems, iItem, 0.96, 0.92, 0.08, 0.11, 0.085, 0.045);
    end
end

%% show comparasions of welllog inversion results
function bsShowPostSubInvLogResult(GShowProfileParam, ...
    invVal, trueVal, initVal, ...
    t, tmethod, range, nItems, iItem, attName)

    GPlotParam = GShowProfileParam.plotParam;
    bsSetSubPlotSize(nItems, iItem);
    
    plot(initVal, t, 'g', 'linewidth', GPlotParam.linewidth); hold on;
    plot(trueVal, t, 'k-.', 'LineWidth', GPlotParam.linewidth);   hold on;
    plot(invVal, t, 'r','LineWidth', GPlotParam.linewidth);    hold on;
    
    if iItem == 1
        if strcmpi(GShowProfileParam.language, 'en')
            ylabel('Time (s)');
        else
            ylabel('时间 \fontname{Times New Roman}(s)');
        end
    else
        set(gca, 'ytick', [], 'yticklabel', []);
    end
    xlabel(attName);
    
    if GShowProfileParam.showPartVert.downTime > 0 && GShowProfileParam.showPartVert.downTime > GShowProfileParam.showPartVert.upTime
        set(gca, 'ylim', [t(GShowProfileParam.showPartVert.upTime), t(GShowProfileParam.showPartVert.downTime)]);
    else
        set(gca, 'ylim', [t(1), t(end)]);
    end
    
    set(gca,'ydir','reverse');
    
    title(sprintf('(%s) %s', char( 'a' + (iItem-1) ), tmethod));
    
    if ~isempty(range)
        set(gca, 'xlim', range) ; 
    end
    
%     set(gca, 'ylim', [t(1) t(end)]);
%     bsTextSeqIdFit(ichar - 'a' + 1, 0, 0, 12);

    set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
end

%% show comparasions of synthetic seismic data and real seismic data
function bsShowPostSubSynSeisData(GShowProfileParam, ...
    synFromTrue, synFromInv, seisData, ...
    t, tmethod, range, nItems, iItem)

    bsSetSubPlotSize(nItems, iItem);
    
    GPlotParam = GShowProfileParam.plotParam;
    
    plot(synFromTrue/norm(synFromTrue), t(1:end-1), 'b-.', 'linewidth', GPlotParam.linewidth); hold on;
    plot(seisData/norm(seisData), t(1:end-1), 'k','LineWidth', GPlotParam.linewidth);    hold on;
    plot(synFromInv/norm(synFromInv), t(1:end-1), 'r', 'LineWidth', GPlotParam.linewidth);   hold on;
    
    if strcmpi(GShowProfileParam.language, 'en')
        ylabel('Time (s)');
    else
        ylabel('时间 \fontname{Times New Roman}(s)');
    end
    set(gca,'ydir','reverse');
    
    title(sprintf('(%s) %s', char( 'a' + (iItem-1) ), tmethod));
    if ~isempty(range)
        set(gca, 'xlim', range) ; 
    end
    set(gca, 'ylim', [t(1) t(end)]);
%     bsTextSeqIdFit(ichar - 'a' + 1, 0, 0, 12);

    set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
end
