function bsShowFFTResultsComparison(GPostInvParam, GShowProfileParam, curves, titles)

    GPlotParam = GShowProfileParam.plotParam;
    
    nItems = size(curves, 2);
    dt = GPostInvParam.dt;


    hf = figure;
    subplot('Position', [0.1, 0.1 0.85 0.85]);

    cTbl = bsGetColormap('separate');
    shapes = {'-', '-.', '--', '-.', '--'};
    for iItem = 1 : nItems

        Ip = curves(:, iItem);
        index = find(Ip ~= 0);
        nonzeroIp = Ip(min(index) : max(index));

        if length(nonzeroIp) < size(curves, 1)
            [pInv, pf] = bsGetFrequencies(nonzeroIp, dt);
            plot(pf, pInv, shapes{iItem}, 'color', cTbl{iItem}, 'linewidth', 2); hold on;
        else
            [pInv, f] = bsGetFrequencies(Ip, dt);
            pInv = pInv / norm(pInv) * norm(curves(:, 1));
            plot(f, pInv, shapes{iItem}, 'color', cTbl{iItem}, 'linewidth', 2); hold on;
        end

    end


    lgd = legend(titles);
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

