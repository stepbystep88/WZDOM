function bsShowPreFFTResult(GInvParam, GShowProfileParam, invVals, wellInfo)

    GPlotParam = GShowProfileParam.plotParam;
    
    nItems = length(invVals);
    model = invVals{1}.model;
    dt = GInvParam.dt;


    

    trueLog = model.trueLog;
    initLog = model.initLog;
%         sampNum = floor(length(model.dTrue)/2)*2;
    attNames = {'纵波速度', '横波速度', '密度'};
    hf = figure;
	set(gcf, 'position', [18         391        1857         604]);
    for itt = 2 : 3
        
        subplot(1, 2, itt-1);
        [pTrue, f] = bsGetFrequencies(trueLog(:, itt), dt);
        pInit = bsGetFrequencies(initLog(:, itt), dt);
        d = model.d;
        seisData = reshape(d, size(initLog, 1)-1, GInvParam.angleTrNum);
        [pD, fd] = bsGetFrequencies(sum(seisData, 2), dt);
        pD = pD / norm(pD) * norm(pTrue);

        cTbl = bsGetColormap('separate');
        legends = {};
        for iItem = 1 : nItems

            invVal = invVals{iItem};
            data = invVal.data{itt - 1};

            index = find(data > 0);
            nonzero = data(min(index) : max(index));

            if length(nonzero) < length(trueLog)
    %             x = linspace(1, length(trueLog), length(nonzeroIp));
    %             newIp = interp1(x, nonzeroIp, 1:length(trueLog));
                [pInv, pf] = bsGetFrequencies(nonzero, dt);
                plot(pf, pInv, 'color', cTbl{iItem + 4}, 'linewidth', 2); hold on;
            else
                pInv = bsGetFrequencies(data, dt);
                pInv = pInv / norm(pInv) * norm(pTrue);
                plot(f, pInv, 'color', cTbl{iItem + 4}, 'linewidth', 2); hold on;
            end

            legends = [legends, invVal.name];
        end

        plot(f, pInit, '--', 'color', cTbl{1}, 'linewidth', 2); 
        plot(f, pTrue, '-.', 'color', cTbl{2}, 'linewidth', 2);
        plot(fd, pD, '--', 'color', cTbl{4}, 'linewidth', 2);

        legends = [legends, '初始模型', '实际测井', '地震资料'];

        lgd = legend(legends);
        set(lgd, ...
            'fontsize', GPlotParam.fontsize, ...
            'fontweight', 'bold', ...
            'fontname', GPlotParam.fontname);
        xlabel('Frequency');
        ylabel('Amplitude');

        if exist('wellInfo', 'var')
            title(sprintf('%s-%s', attNames{itt-1}, wellInfo.name));
        end
    %         set(gca, 'ylim', [0, 2000]);
        set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
    end
    
    

end

function [P1, f] = bsGetFrequencies(data, dt)
    
    L = length(data);
%     n = L;
    n = floor(L/2)*2;
	data = data(1:n);
    
%     n = 2^nextpow2(L);
%     data = [data; zeros(n - L, 1)];

    Y = fft(data, n);
%     P2 = abs(Y / n);
    P2 = abs(Y);
    P1 = P2(1:n/2+1);
    Fs = 1/dt*1000;
    f = Fs * (0:(n/2))/n;
    
    f = f(2:end);
    P1 = P1(2:end);
    
end

