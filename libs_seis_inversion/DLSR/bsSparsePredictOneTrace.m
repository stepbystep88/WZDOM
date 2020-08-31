function [reconstructed, gammas] = bsSparsePredictOneTrace(GSParam, input, inline, crossline)
%% 根据低分辨率字典预测高分辨率
    
    [sampNum, ~] = size(input{1});
    
    % 首先将输入转化为小块
    [all_patches] = bsSparseTransInput2Patches(GSParam, input, inline, crossline);
    % 归一化
	[normal_patches, output] = bsSparseNormalization(GSParam.trainDICParam.normalizationMode, all_patches, GSParam.low_min_values, GSParam.low_max_values);
    % 稀疏表示
    gammas = omp(GSParam.lowDIC'*normal_patches, GSParam.omp_low_G, GSParam.sparsity);
    new_patches = GSParam.highDIC *  gammas;
    % 反归一化
    denormal_patches = bsSparseDenormalization(GSParam.trainDICParam.normalizationMode, new_patches, output, GSParam.high_min_values, GSParam.high_max_values);
    
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
    reconstructed = bsAvgPatches(denormal_patches, GSParam.index, sampNum);
    
end