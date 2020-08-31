function GShowProfileParam = bsShowMultiClassfiedDIC(GShowProfileParam, trainDICParam, MDIC, K)

    nSpecialFeat = trainDICParam.isAddLocInfo * 2 + trainDICParam.isAddTimeInfo;
    MDIC = MDIC(nSpecialFeat+1:end, :);
    
    [MDIC, Idx] = bsClusterDIC(MDIC, K);
    [colormap_K] = bsHex2RGB(bsGetColormap('separate'));
    GShowProfileParam.plotParam.colormap = colormap_K(2:K+1, :);
    GShowProfileParam.plotParam.Idx = Idx;
    
    figure; 
    imagesc(Idx); 
    title('(b) 原子类别', 'fontweight', 'bold');
    xlabel('原子序列');
    ylabel('类别');
    colorbar;
    set(gcf, 'position', [251         299        1347         130]);
    bsSetDefaultPlotSet(GShowProfileParam.plotParam);
    set(gca, 'colormap', colormap_K(2:K+1, :));
    
    nSub = 60 / K;
    
    showDIC = zeros(size(MDIC, 1), 60);
    showIdx = [];
    for iK = 1 : K
        sPos = (iK-1)*nSub + 1;
        ePos = sPos + nSub - 1;

        index = find(Idx==iK);
        showDIC(:, sPos:ePos) = MDIC(:, index(1:nSub)); 
        showIdx = [showIdx, Idx(index(1:nSub))];
    end

    bsShowMultiDictionary(showDIC, showIdx, colormap_K(2:K+1, :));
    set(gcf, 'position', [668   121   980   497]);
end