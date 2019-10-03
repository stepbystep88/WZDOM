function isSame = bsIsMember(va, nests, type)
    switch type
        case 'rows'
            for i = 1 : size(nests, 1)
                if bsCompareTwoVectors(va, nests(i, :))
                    isSame = true;
                    return;
                end
            end
        case 'columns'
            for i = 1 : size(nests, 2)
                if bsCompareTwoVectors(va, nests(:, i))
                    isSame = true;
                    return;
                end
            end
    end
    
    isSame = false;
    
end