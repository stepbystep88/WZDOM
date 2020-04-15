function bsShowMultiDictionary(unionDic, Idx, colormap_K)

    figure;
    nAtom = size(unionDic, 2);
    sampNum = size(unionDic, 1) / 3;
    
%     if(nAtom >= 72)
%         nAtom = 72;
%     end
    
    D = bsGen1DDiffOperator(sampNum, 1, 1);
    diff = D * unionDic(1:sampNum, :);
    
    data = zeros(nAtom, 2);
    for i = 1 : nAtom
%         data(i, 1) = norm(diff(:, i));
        data(i, 1) = norm(diff(:, i), 1);
%         data(i, 1) = norm(unionDic(:, i), 2);
        data(i, 2) = i;
    end
    
%     l1_diff = sum(abs(diff));
%     data = [l1_diff; 1:nAtom];
    data = sortrows(data, 1);
    
    unionDic = unionDic(:, data(:, 2));
%     if(nAtom <= 25)
% %         row = round( sqrt( GSparseInvParam.nAtom / 2 ) );
%         row = 5;
%     else
%         row = 5;
%     end
    
    row = 5;
    col = 12;
%     nAtom = row * col;
    x = 0.04;
    y = 0.06;
    annotation('arrow', [x x], [1 y]);
    annotation('arrow', [x 1], [y y]);
    index = floor(linspace(1, nAtom, row * col));
    unionDic = unionDic(:, index);
%     col = ceil(nAtom / row);
    nAtom = row * col;
    colors = {'r', 'b--', 'k-.'};
    
    for i = 1 : nAtom
        bsSubPlotTightestHL(row, col, i, x+0.01, y+0.01);
        
%         tmp = abs(Dic(:, i));
        for k = 1 : 3
            ps = (k-1)*sampNum + 1;
            pe = k * sampNum;
            tmp = unionDic(ps:pe, i);
%             tmp = mapminmax(tmp', 0, 1);
%             tmp = tmp';
            plot(tmp, 1:length(tmp), colors{k}, 'linewidth', 1.5);
            hold on;
        end
        
        if nargin >= 3
            set(gca,'color', [colormap_K(Idx(i), :), 0.5]);
        end
        
%         maxVal = max(tmp);
%         minVal = min(tmp);
%         set(gca, 'xlim', [0.95*minVal, 1.05*maxVal]);
        
        
%         set(gca, 'xlim', [-0.1, 1.1]);
        maxVal = max(unionDic(:, i));
        minVal = min(unionDic(:, i));
        gap = maxVal - minVal;
        
        set(gca, 'xlim', [minVal - 0.1*gap, maxVal + 0.1*gap]);
        set(gca, 'ylim', [0, sampNum+1]);
        set(gca,'xtick',[],'xticklabel',[])
        set(gca,'ytick',[],'yticklabel',[])
        
        set(gca,'ydir','reverse');
    end
    
%     set(gcf, 'position', stpCalcBestPosition(0.3, 0.3) );
    set(gcf, 'position', [1031         121         617         303]);
%     gtext('Time increases', 'fontsize', GPlotParam.fontsize+1,'fontweight', 'bold', 'fontname', GPlotParam.fontname, 'rotation', -90);
%     gtext('Numerical value increases', 'fontsize', GPlotParam.fontsize+1,'fontweight', 'bold', 'fontname', GPlotParam.fontname);
end