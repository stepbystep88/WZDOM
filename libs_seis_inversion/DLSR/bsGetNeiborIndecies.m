function indecies = bsGetNeiborIndecies(D, n)
% 查找字典D中的每一个原子的n个相似原子索引
    [~, nAtom] = size(D);
    
    similarites = zeros(nAtom, nAtom);
    indecies = zeros(nAtom, n);
    
    for i = 1 : nAtom
        for j = i+1 : nAtom
            x = D(:, i);
            y = D(:, j);
            
            tmp = sum((x-y).^2);
            similarites(i, j) = tmp;
            similarites(j, i) = similarites(i, j);
        end
    end

    
    for i = 1 : nAtom
%         [B, I] = sort(similarites(i, :), 'descend');
        [B, I] = sort(similarites(i, setdiff(1:nAtom, i)));
        indecies(i, :) = I(1:n);
    end
end