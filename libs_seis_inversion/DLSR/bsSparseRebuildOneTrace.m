function [reconstructed, gammas] = bsSparseRebuildOneTrace(GSParam, input, gamma, inline, crossline, nlm_ps)
    % sparse reconstruction
    
    [sampNum, ~] = size(input{1});
    nAtt = length(input);
    sizeAtom = GSParam.sizeAtom;

    reconstructed = zeros(sampNum, nAtt);
    
    % 首先将输入转化为小块 
    if nargin > 5 && ~isempty(nlm_ps)
        [all_patches] = bsSparseTransInput2PatchesByNLM(GSParam, input, nlm_ps, inline, crossline);
    else
        [all_patches] = bsSparseTransInput2Patches(GSParam, input, inline, crossline);
    end
    
    % 归一化
	[normal_patches, output] = bsSparseNormalization(GSParam.trainDICParam.normalizationMode, all_patches, GSParam.min_values, GSParam.max_values);
    % 稀疏表示
    gammas = omp(GSParam.ODIC'*normal_patches, GSParam.omp_G, GSParam.sparsity);
    new_patches = GSParam.ODIC *  gammas;
    % 反归一化
    denormal_patches = bsSparseDenormalization(GSParam.trainDICParam.normalizationMode,new_patches, output, GSParam.min_values, GSParam.max_values);
    
%     if length(GSParam.trainDICParam.normalizationMode) > 5
%         figure(100); 
%         subplot(1, 2, 1);
%         plot(all_patches(:,1), 'b', 'linewidth', 2);  hold on;
%         plot(denormal_patches(:,1), 'r', 'linewidth', 2); 
%         legend('原始归一化前', '重构去归一化后');
% 
%         subplot(1, 2, 2);
%         plot(normal_patches(:,1), 'k', 'linewidth', 2); hold on;
%         plot(new_patches(:,1), 'r', 'linewidth', 2); 
%         legend('原始归一化后', '重构去归一化前');
%     end
    
    % 将重构的小块整理为完整的一个信号
    for i = 1 : nAtt
        sPos = sizeAtom*(i-1) + 1 + GSParam.nSpecialFeat;
        ePos = sPos + sizeAtom - 1;

        i_new_patches = denormal_patches(sPos:ePos, :);
        reconstructed(:, i) = bsAvgPatches(i_new_patches, GSParam.index, sampNum);
        
        if gamma < 1 && gamma > 0
            reconstructed(:, i) = reconstructed(:, i) * gamma + input{i}(:, 1) * (1 - gamma);
        elseif gamma == 0
            reconstructed(:, i) = input{i}(:, 1);
        elseif gamma == 1
%             reconstructed(:, i) = reconstructed(:, i);
        end
    end
    
end