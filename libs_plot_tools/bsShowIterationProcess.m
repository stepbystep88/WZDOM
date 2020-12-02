function bsShowIterationProcess(GInvParam, GShowProfileParam, invVals, outputs, wellLog, methods)
    
    nItems = length(invVals);
    model = invVals{1}.model;
    
    
    
    colorTbl = bsGetColormap('separate');
    names = cell(1, nItems);
    
    if strcmpi(GInvParam.flag, 'poststack')
        bsShowPost();
    else
        bsShowPre();
    end
    
    function bsShowPost()
        figure;
        for i = 1 : nItems
            midResults = outputs{i}.midResults;
            nIter = size(midResults.f, 2);
            names{i} = invVals{i}.name;
            rrsem = bsCalclPostRRSEM(midResults.x, model.trueLog);

            plot(1:nIter, rrsem, '-', 'color', colorTbl{i}, 'linewidth', 2); hold on;
        end

        xlabel('Iteration');
        ylabel('RMSE of model');
        legend(names, 'fontweight', 'bold');
        title(wellLog.name, 'fontweight', 'bold');
        set(gca, 'xlim', [1, nIter]);

        bsSetDefaultPlotSet(GShowProfileParam.plotParam);

        figure;
        for i = 1 : nItems
            midResults = outputs{i}.midResults;
            nIter = size(midResults.f, 2);
            names{i} = invVals{i}.name;
            plot(1:nIter, midResults.f, '-', 'color', colorTbl{i}, 'linewidth', 2); hold on;
        end

        xlabel('Iteration');
        ylabel('Objective function value');
        legend(names, 'fontweight', 'bold');
        title(wellLog.name, 'fontweight', 'bold');
        set(gca, 'xlim', [1, nIter]);
        bsSetDefaultPlotSet(GShowProfileParam.plotParam);
    end

    function bsShowPre()
        
%         left_color = colorTbl{4};
%         right_color = colorTbl{8};
        left_color = [0.1 0.1 1];
        right_color = [1 0.1 0.1];
        
        rrsem = cell(3, nItems);
        fs = cell(1, nItems);
        legend_strs = cell(1, nItems*2);
        nIter = 1;
        
        for i = 1 : nItems
            midResults = outputs{i}.midResults;
            nIter = max(nIter, size(midResults.f, 2));
            legend_strs{i} = sprintf('%s-相对误差', invVals{i}.name);
            legend_strs{i+nItems} = sprintf('%s-目标函数值', invVals{i}.name);
            names{i} = invVals{i}.name;
            
            [rrsem{1, i}, rrsem{2, i}, rrsem{3, i}] = bsCalclPreRRSEM(GInvParam, midResults.x, model.trueLog, model.initLog, model);
            fs{i} = sum((model.d - model.G*midResults.x).^2);
%             fs{i} = midResults.f;
        end
        
        fig = figure;
        set(gcf, 'position', [ 450   325   768   343]);
%         set(fig,'defaultAxesColorOrder',[bsHex2RGB({left_color}); bsHex2RGB({right_color})]);
%         set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        
        shapes = {':', '--', '-'};
        colors = colorTbl([4,2,1]);
        plotParam = GShowProfileParam.plotParam;
        
        loc = [0.94, 0.83, 0.07, 0.075, 0.075, 0];
        
        for i = 1 : 1
            
%             subplot(1, 2, 1);
            bsSubPlotFit(1,2,1, loc(1), loc(2), loc(3), loc(4), loc(5), loc(6));
%             yyaxis left;
            lines = [];
            for j = 1 : nItems
                line = plot(1:length(rrsem{i, j}), rrsem{i, j}, '-', 'color', colors{j}, 'linewidth', 2); hold on;
                lines = [lines, line];
                
                switch (upper(methods{j}.flag))
                    case {'SSR', 'CSR', 'DLSR'}
                        index = 1 : methods{j}.options.innerIter + 1 : length(rrsem{i, j});
                        index(1) = [];
                        x = 1 : length(rrsem{i, j});
                        y = rrsem{i, j};
                        plot(x(index), y(index), 'o', 'color', colors{j}, 'linewidth', 2); hold on;
                    otherwise
                end
                
            end
            xlabel('迭代次数');
            ylabel('模型的相对均方根误差');
            title('(a) 模型收敛情况');
            legend(lines, names, 'fontname', plotParam.fontname, 'fontsize', plotParam.fontsize);
%             title(wellLog.name);
            set(gca, 'xlim', [1, nIter-10]);
            bsSetDefaultPlotSet(GShowProfileParam.plotParam);
            
            
            %激活右侧
%             yyaxis right;
%             subplot(1, 2, 2);
            bsSubPlotFit(1,2,2, loc(1), loc(2), loc(3), loc(4), loc(5), loc(6));
            line = [];
            for j = 1 : nItems
                line = plot(1:length(fs{j}), fs{j}, '-', 'color', colors{j}, 'linewidth', 2); hold on;
                lines = [lines, line];
                
                switch (upper(methods{j}.flag))
                    case {'SSR', 'CSR', 'DLSR'}
                        index = 1 : methods{j}.options.innerIter + 1 : length(fs{i, j});
                        index(1) = [];
                        x = 1 : length(fs{i, j});
                        y = fs{i, j};
                        plot(x(index), y(index), 'o', 'color', colors{j}, 'linewidth', 2); hold on;
                    otherwise
                end
            end
            ylabel('目标函数值');
            xlabel('迭代次数');
            title('(b) 目标函数收敛情况');
            
            set(gca, 'xlim', [1, nIter-10]);
%             legend(legend_strs, 'fontname', plotParam.fontname, 'fontsize', plotParam.fontsize);
            
            bsSetDefaultPlotSet(GShowProfileParam.plotParam);
        end
    end
end




function [rrsem_vp, rrsem_vs, rrsem_rho] = bsCalclPreRRSEM(GInvParam, xs, trueLog, initLog, model)
    [vp, vs, rho] = bsPreRecoverElasticParam(xs, GInvParam.mode, model.lsdCoef);
    nIter = size(xs, 2);
    rrsem_vp = zeros(1, nIter);
    rrsem_vs = zeros(1, nIter);
    rrsem_rho = zeros(1, nIter);
    
    true_rrmse_vp = sqrt(mse(trueLog(:, 2) - initLog(:, 2)));
    true_rrmse_vs = sqrt(mse(trueLog(:, 3) - initLog(:, 3)));
    true_rrmse_rho = sqrt(mse(trueLog(:, 4) - initLog(:, 4)));
    
    for i = 1 : nIter
        rrsem_vp(i) = sqrt(mse((vp(:, i) - trueLog(:, 2)))) / true_rrmse_vp;
        rrsem_vs(i) = sqrt(mse((vs(:, i) - trueLog(:, 3)))) / true_rrmse_vs;
        rrsem_rho(i) = sqrt(mse((rho(:, i) - trueLog(:, 4)))) / true_rrmse_rho;
    end
    
end

function rrsem = bsCalclRRSEM(xs, trueLog)
    ips = exp(xs);
    nIter = size(xs, 2);
    rrsem = zeros(1, nIter);
    
    for i = 1 : nIter
        ip = ips(:, i);
        rrsem(i) = sqrt(mse((ip - trueLog)));
    end
    
end