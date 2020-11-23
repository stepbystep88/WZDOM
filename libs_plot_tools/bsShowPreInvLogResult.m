function bsShowPreInvLogResult(GInvParam, GShowProfileParam, ...
    invVals, trueLogFiltcoef, isShowSynSeisComparasion, timeLine)

    if ~exist('isShowSynSeisComparasion', 'var')
        isShowSynSeisComparasion = 1;
    end
    
    GPlotParam = GShowProfileParam.plotParam;
    
    nItems = length(invVals);
    model = invVals{1}.model;
    trueLog = model.trueLog;
    
    sampNum = size(model.trueLog, 1);
    t = (model.t0 : GInvParam.dt : model.t0 + (sampNum-1)*GInvParam.dt) / 1000;
    
    type = invVals{1}.type;
    ntype = length(type);
    
    % 绘制层位线
    is_plot_timeline = false;
    
    if exist('timeLine', 'var')
        is_plot_timeline = true;
        times = zeros(1, length(timeLine));
        
        for iLine = 1 : length(timeLine)
            times(iLine) = bsGetHorizonTime(timeLine{iLine}, invVals{1}.inIds, invVals{1}.crossIds, 0) / 1000;
        end
    end
        
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
            
            for k = 1 : ntype
                
                [~, ~, ~, ~, attName] = ...
                    bsGetInfoByType(GShowProfileParam, GInvParam, type{k});

                switch lower(type{k})
                    case 'vp'
                        
                        index = GInvParam.indexInWellData.vp;
                        range = GShowProfileParam.range.vp/1000;
                        
                        bsShowPostSubInvLogResult(GShowProfileParam, ...
                        invVal.data{k}/1000, trueLog(:, index)/1000, model.initLog(:, index)/1000, ...
                        t, invVal.name, attName, range, ...
                        nItems, iItem, k, ntype);
                    case 'vs'
                        index = GInvParam.indexInWellData.vs;
                        range = GShowProfileParam.range.vs/1000;
                        
                        bsShowPostSubInvLogResult(GShowProfileParam, ...
                        invVal.data{k}/1000, trueLog(:, index)/1000, model.initLog(:, index)/1000, ...
                        t, invVal.name, attName, range, ...
                        nItems, iItem, k, ntype);
            
                    case 'rho'
                        index = GInvParam.indexInWellData.rho;
                        range = GShowProfileParam.range.rho;
                        
                        bsShowPostSubInvLogResult(GShowProfileParam, ...
                            invVal.data{k}, trueLog(:, index), model.initLog(:, index), ...
                            t, invVal.name, attName, range, ...
                            nItems, iItem, k, ntype);
                    case {'vp_vs', 'vsvs_ratio', 'vp_vs_ratio'}
                        vpvs_trueLog = bsGetVp_Vs(trueLog(:, GInvParam.indexInWellData.vp), trueLog(:, GInvParam.indexInWellData.vs));
                        vpvs_initLog = bsGetVp_Vs(model.initLog(:, GInvParam.indexInWellData.vp), model.initLog(:, GInvParam.indexInWellData.vs));
                        range = GShowProfileParam.range.vpvs_ratio;
                        
                        bsShowPostSubInvLogResult(GShowProfileParam, ...
                            invVal.data{k}, vpvs_trueLog, vpvs_initLog, ...
                            t, invVal.name, attName, range, ...
                            nItems, iItem, k, ntype);
                end
                
                if is_plot_timeline
                    colors = {'c--', 'k--', 'y--'};
                    
                    for ii = 1 : length(timeLine)
                        plot(range, [times(ii), times(ii)], colors{ii}, 'linewidth', 2);
                    end
                end
        
            end
            
        end
        
        if strcmpi(GShowProfileParam.language, 'en')
            legends = {'Initial model', 'Real model', 'Inversion result'};
        else
            legends = {'初始模型', '真实模型', '反演结果'};
        end
        
        bsSetLegend(GPlotParam, {'g', 'k', 'r'}, legends);
       
        
        
    end

    function bsShowSyntheticSeisData()
        hf = figure;
        bsSetPreFigureSize(ceil((nItems+2)/2));
        angleData = model.angleData;
        if max(angleData < 5)
            angleData = angleData / pi * 180;
        end
        angleData = round(angleData);
        
        synFromTrue = reshape(model.G * model.trueX, sampNum-1, GInvParam.angleTrNum);
        seisData = reshape(model.original_d, sampNum-1, GInvParam.angleTrNum);
        
        if strcmpi(GShowProfileParam.language, 'en')
            bsShowPreSubSynSeisData(GPlotParam, 'real', seisData, t, angleData, nItems, 1);
            bsShowPreSubSynSeisData(GPlotParam, 'synthetic from welllog', synFromTrue, t, angleData, nItems, 2);
        else
            bsShowPreSubSynSeisData(GPlotParam, '实际观测道集', seisData, t, angleData, nItems, 1);
            bsShowPreSubSynSeisData(GPlotParam, '基于测井的合成道集', synFromTrue, t, angleData, nItems, 2);
        end
        
        for iItem = 1 : nItems
        
            figure(hf);
            

            invVal = invVals{iItem};
            
%             invLog = [bsGetDepth(invVal.data{1}, GInvParam.dt), invVal.data{1}, invVal.data{2}, invVal.data{3}];
%             x1 = bsPreBuildModelParam(invLog, GInvParam.mode, model.lsdCoef);
            d = invVal.model.G * invVal.xOut;
            synFromInv = reshape(d, sampNum-1, GInvParam.angleTrNum);
            
            
            
            bsShowPreSubSynSeisData(GPlotParam, sprintf('基于%s的合成道集', invVal.name), synFromInv, t, angleData, nItems, iItem+2);
            
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
%             bsSetPosition(0.44, 0.48);
            set(gcf, 'position', [308          54        1231         710]);
        case 3
            bsSetPosition(0.44, 0.68);
        otherwise
            bsSetPosition(0.44, 0.88);
    end
end

function bsSetPreSubPlotSize(nItems, iItem, k, ntype)
    if nargin < 4
        ntype = 3;
    end
    
    index = ntype * (iItem - 1) + k;
    switch nItems
        case 1
            bsSubPlotFit(nItems, ntype, index, 0.93, 0.75, 0.02, 0.00, 0.06, -0.05);
        case 2
            bsSubPlotFit(nItems, ntype, index, 0.93, 0.88, 0.02, 0.05, 0.06, 0.01);
        case 3
            bsSubPlotFit(nItems, ntype, index, 0.93, 0.88, 0.02, 0.05, 0.06, 0.01);
        case 4
            bsSubPlotFit(nItems, ntype, index, 0.93, 0.92, 0.02, 0.03, 0.06, 0.00);
    end
end

%% show comparasions of welllog inversion results
function bsShowPostSubInvLogResult(GShowProfileParam, ...
    invVal, trueVal, initVal, ...
    t, tmethod, attName, range, nItems, iItem, k, ntype)
    
    GPlotParam = GShowProfileParam.plotParam;
    
    bsSetPreSubPlotSize(nItems, iItem, k, ntype);
    plot(initVal, t, 'g', 'linewidth', GPlotParam.linewidth); hold on;
    plot(trueVal, t, 'b-.', 'LineWidth', GPlotParam.linewidth);   hold on;
    plot(invVal, t, 'r','LineWidth', GPlotParam.linewidth);    hold on;
        
    if exist('timeLine', 'var')
        
    end
    
    if k == 1
        if strcmpi(GShowProfileParam.language, 'en')
            ylabel('Time (s)');
        else
            ylabel('时间 (s)');
        end
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

function bsShowPreSubSynSeisData(GShowProfileParam, tmethod, data, t, angleData, nItems, iItem)

    GPlotParam = GShowProfileParam.plotParam;
    
%     bsSetPreSubPlotSize(nItems, iItem, k);
    nRow = ceil((nItems+2)/2);
%     subplot(nRow, 3, iItem);
    bsSubPlotFit(nRow, 2, iItem, 0.90, 0.93, 0.02, 0.11, 0.07, -0.0);
    
    wiggles = {'k', 'b', 'r'};
    peaks = {'k', 'b', 'r'};
    
    seismic = s_convert(data, t(1), t(2)-t(1));
    s_wplot(seismic, {'figure', 'old'}, {'xaxis', angleData});%, ...
%         {'wiggle_color', wiggles{k}}, ...
%         {'peak_fill', peaks{k}});
    
    
    if mod(iItem, 2) == 1
        if strcmpi(GShowProfileParam.language, 'en')
            ylabel('Time (s)');
        else
            ylabel('时间 (s)');
        end
    else
        set(gca, 'ytick', [], 'yticklabel', []);
    end
    
    xlabel(sprintf('(%s) %s', char( 'a' + (iItem-1) ), tmethod), ...
            'fontsize', GPlotParam.fontsize+2, 'fontweight', 'bold');

    set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
end
