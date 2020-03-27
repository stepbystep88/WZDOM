function [newGammas, gammas] = bsOMP(D, patches, DTD, K, nIndex)

    % 第一步，利用OMP算法获取稀疏原子的位置
    gammas = omp(D'*patches, DTD, K);
    newGammas = zeros(size(gammas));
    
    nPatch = size(patches, 2);
    for i = 1 : nPatch
        nonzero = find(gammas(:, i) ~= 0);
        pos = nonzero;
        
        for j = 1 : length(nonzero)
            neibors = nIndex(nonzero(j), :);
            pos = [pos; neibors'];
            
        end
        
        pos = unique(pos);
        subD = D(:, pos');
        
%         d = sum(subD, 2);
%         newGammas(pos, i) = d' * patches(:, i) / (d'*d);
        % |subD*a - p|_2 -> subD'*subD*al = p
        newGammas(pos, i) = inv(subD' * subD) * (subD' * patches(:, i));
%         subD = D(:, nonzero);
%         
%         newGammas(nonzero, i) = inv(subD' * subD + 0.001) * (subD' * patches(:, i));
    end
end