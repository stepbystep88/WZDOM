function [newData] = bsPreSparseRebuildWelllog(GInvParam, GShowProfileParam, wellInfo, GSParam)
    
    wellData = wellInfo.wellLog(:, [GInvParam.indexInWellData.vp, GInvParam.indexInWellData.vs, GInvParam.indexInWellData.rho]);
    
    [sampNum, nAtt] = size(wellData);
    GSParam = bsInitDLSRPkgs(GSParam, sampNum);
 
    sizeAtom = GSParam.sizeAtom;
    ncell = GSParam.ncell;
    patches = zeros(sizeAtom*3+GSParam.nSpecialFeat, ncell);

    newData = zeros(size(wellData));
    
    for j = 1 : ncell
        js = GSParam.index(j);

        for i = 1 : nAtt

            sPos = sizeAtom*(i-1) + 1 + GSParam.nSpecialFeat;
            ePos = sPos + sizeAtom - 1;

            iData = wellData(js : js+sizeAtom-1, i);
            patches(sPos:ePos, j) = iData;
        end
    end

    % add location or/and time information
    if GSParam.trainDICParam.isAddLocInfo && GSParam.trainDICParam.isAddTimeInfo
        patches(1:GSParam.nSpecialFeat, :) = [ones(1, ncell) * wellInfo.inline; ones(1, ncell) * wellInfo.crossline; 1 : ncell];
    elseif GSParam.trainDICParam.isAddLocInfo
        patches(1:GSParam.nSpecialFeat, :) = [ones(1, ncell) * wellInfo.inline; ones(1, ncell) * wellInfo.crossline;];
    elseif GSParam.trainDICParam.isAddTimeInfo
        patches(1:GSParam.nSpecialFeat, :) = 1 : ncell;
    end

    % normalization
    normal_patches = (patches - GSParam.min_values) ./ (GSParam.max_values - GSParam.min_values);

    if strcmp(GSParam.trainDICParam.feature_reduction, 'all')
        normal_patches = GSParam.output.B' * normal_patches;
    end

    if GSParam.isModifiedDIC
        normal_patches = GSParam.M  * normal_patches;
    end

    gammas = omp(GSParam.DIC'*normal_patches, ...
                GSParam.omp_G, ...
                GSParam.sparsity);
    new_patches = GSParam.DIC *  gammas;

    if strcmp(GSParam.trainDICParam.feature_reduction, 'all')
        new_patches = GSParam.output.B * new_patches;
    end

    new_patches = new_patches .* (GSParam.max_values - GSParam.min_values) + GSParam.min_values;



    %% reconstruct model by equations
    for i = 1 : nAtt
        sPos = sizeAtom*(i-1) + 1 +  + GSParam.nSpecialFeat;
        ePos = sPos + sizeAtom - 1;

        i_new_patches = new_patches(sPos:ePos, :);
        newData(:, i) = bsAvgPatches(i_new_patches, GSParam.index, sampNum);
    end

    showResults(GShowProfileParam.plotParam, GSParam, wellData, newData, wellInfo.wellLog(:, GInvParam.indexInWellData.time)/1000, gammas, wellInfo.name);
    
end

function showResults(GPlotParam, GSParam, oldData, newData, t, gammas, name)
    figure;
    gammas = abs(gammas);
    
    for i = 1 : 4
%         bsSetPreSubPlotSize(4, i);
        bsSubPlotFit(1, 4, i, 0.93, 0.85, 0.02, 0.07, 0.06, 0.00);
        
        if i < 4
            plot(oldData(:, i), t, 'b-.', 'LineWidth', GPlotParam.linewidth);   hold on;
            plot(newData(:, i), t, 'r','LineWidth', GPlotParam.linewidth);    hold on;
            
            switch i
                case 1
                    attName = 'Vp (m/s)';
                case 2
                    attName = 'Vs (m/s)';
                case 3
                    attName = 'Rho (g/cc)';
            end
            
            set(gca, 'xlim', [min(oldData(:, i)),max(oldData(:, i))] );
            set(gca, 'ylim', [t(1), t(end)]);
        else
%             DIC = GSParam.DIC;
            Idx = GPlotParam.Idx;
            
            labels = zeros(size(newData, 1), 1);
            for j = 1 : GSParam.ncell
                js = GSParam.index(j) + GSParam.sizeAtom / 2;
                [~, index] = max(gammas(:, j));
                labels(js) = Idx(index);
            end
            
            imagesc(labels);
            set(gca, 'colormap', [1 1 1; GPlotParam.colormap]);
            set(gca, 'clim', [0, max(Idx)]);
            set(gca, 'xtick', [], 'xticklabel', []);
        end
        
        if i == 1
            ylabel('Time (s)');
        else
            set(gca, 'ytick', [], 'yticklabel', []);
        end
        xlabel(attName);
            
        if i == 2
            title(sprintf('%s', name), ...
                'fontsize', GPlotParam.fontsize+2, 'fontweight', 'bold', 'interpreter','none');
        end
    
    end

    
%     xlabel(attName);
    set(gca, 'ydir', 'reverse');
    set(gca , 'fontsize', GPlotParam.fontsize,'fontweight', GPlotParam.fontweight, 'fontname', GPlotParam.fontname);
    
end

function GSParam = bsInitDLSRPkgs(GSParam, sampNum)
    
    validatestring(string(GSParam.reconstructType), {'equation', 'simpleAvg'});
    
    trainDICParam = GSParam.trainDICParam;
    
    sizeAtom = trainDICParam.sizeAtom;
    nAtom= trainDICParam.nAtom;
    GSParam.sizeAtom = sizeAtom;
    GSParam.nAtom = nAtom;
    GSParam.nrepeat = sizeAtom - GSParam.stride;
    
    GSParam.nSpecialFeat = trainDICParam.isAddLocInfo * 2 + trainDICParam.isAddTimeInfo;
    
    index = 1 : GSParam.stride : sampNum - sizeAtom + 1;
    if(index(end) ~= sampNum - sizeAtom + 1)
        index = [index, sampNum - sizeAtom + 1];
    end
    
    GSParam.index = index;
    GSParam.ncell = length(index);
        
    rangeCoef = GSParam.rangeCoef;
    
    GSParam.min_values = repmat(rangeCoef(:, 1), 1, GSParam.ncell);
	GSParam.max_values = repmat(rangeCoef(:, 2), 1, GSParam.ncell);
    
    GSParam.omp_G = GSParam.DIC' * GSParam.DIC;
end


