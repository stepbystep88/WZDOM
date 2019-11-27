function v = bsAvgPatches(vs, index, sampNum)
    
    v = zeros(sampNum, 1);
    count = zeros(sampNum, 1);
    [sizeAtom, ncell] = size(vs);
    
    for i = 1 : ncell
        pos = index(i);
        v(pos : pos+sizeAtom-1) = v(pos : pos+sizeAtom-1) + vs(:, i);
        count(pos : pos+sizeAtom-1) = count(pos : pos+sizeAtom-1) + 1;
    end
    v = v ./ count;
end
