function bsShowPostFFTResult(GPostInvParam, GShowProfileParam, invVals)

    GPlotParam = GShowProfileParam.plotParam;
    
    nItems = length(invVals);
    model = invVals{1}.model;
    dt = GPostInvParam.dt;


    hf = figure;
    set(gcf, 'position', [40         171        1387         807]);
    subplot('Position', [0.1, 0.1 0.85 0.85]);

    trueLog = model.trueLog;
    initLog = model.initLog;
%         sampNum = floor(length(model.dTrue)/2)*2;

    [pTrue, f] = bsGetFrequencies(trueLog, dt);
    pInit = bsGetFrequencies(initLog, dt);
    [pD, fd] = bsGetFrequencies(model.dTrue, dt);
    pD = pD / norm(pD) * norm(pTrue);



    cTbl = bsGetColormap('separate');
    legends = {};
    for iItem = 1 : nItems

        invVal = invVals{iItem};
        Ip = invVal.Ip;
        index = find(Ip ~= 0);
        nonzeroIp = Ip(min(index) : max(index));

        if length(nonzeroIp) < length(trueLog)
%             x = linspace(1, length(trueLog), length(nonzeroIp));
%             newIp = interp1(x, nonzeroIp, 1:length(trueLog));
            [pInv, pf] = bsGetFrequencies(nonzeroIp, dt);
            plot(pf, pInv, 'color', cTbl{iItem + 4}, 'linewidth', 2); hold on;
        else
            pInv = bsGetFrequencies(Ip, dt);
            pInv = pInv / norm(pInv) * norm(pTrue);
            plot(f, pInv, 'color', cTbl{iItem + 4}, 'linewidth', 2); hold on;
        end

        legends = [legends, invVal.name];
    end

    plot(f, pInit, '--', 'color', cTbl{1}, 'linewidth', 2); 
    plot(f, pTrue, '-.', 'color', cTbl{2}, 'linewidth', 2);
%     plot(fd, pD, '--', 'color', cTbl{4}, 'linewidth', 2);

    legends = [legends, '初始模型', '实际测井', '归一化地震数据频带'];

    lgd = legend(legends);
    set(lgd, ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'bold', ...
        'fontname', GPlotParam.fontname);
    xlabel('Frequency');
    ylabel('Amplitude');
%         set(gca, 'ylim', [0, 2000]);
    set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);

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

