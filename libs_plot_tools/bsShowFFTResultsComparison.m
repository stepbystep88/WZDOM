function bsShowFFTResultsComparison(dt, curves, titles)

    GShowProfileParam = bsCreateGShowProfileParam();
    GPlotParam = GShowProfileParam.plotParam;
    GPlotParam.fontsize = 12;
    
    nItems = size(curves, 2);


    hf = figure;
    subplot('Position', [0.1, 0.15 0.85 0.8]);

    cTbl = bsGetColormap('separate');
    shapes = {'-', '-.', '--', '-.', '--'};
    for iItem = 1 : nItems

        Ip = curves(:, iItem);
        index = find(Ip ~= 0);
        nonzeroIp = Ip(min(index) : max(index));

        if length(nonzeroIp) < size(curves, 1)
            [pInv, pf] = bsGetFrequencies(nonzeroIp, dt);
%             pInv = pInv / norm(pInv) * norm(curves(:, 1));
            plot(pf, pInv, shapes{iItem}, 'color', cTbl{iItem}, 'linewidth', 2); hold on;
        else
            [pInv, f] = bsGetFrequencies(Ip, dt);
%             pInv = pInv / norm(pInv) * norm(curves(:, 1));
            plot(f, pInv, shapes{iItem}, 'color', cTbl{iItem}, 'linewidth', 2); hold on;
        end

    end


    lgd = legend(titles);
    set(lgd, ...
        'fontsize', GPlotParam.fontsize, ...
        'fontweight', 'normal', ...
        'fontname', GPlotParam.fontname);
    
    if strcmpi(GShowProfileParam.language, 'en')
        xlabel('Frequency');
        ylabel('Amplitude');
    else
        xlabel('Ƶ�� \fontname{Times New Roman}(Hz)');
        ylabel('���');
    end
    
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

