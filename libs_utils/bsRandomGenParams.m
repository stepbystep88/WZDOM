function params = bsRandomGenParams(range, isASet, n)
    if isASet
        len = length(range);
        index = randi(len, 1, n);
        params = range(index);
    else
        a = range(1);
        b = range(2);
        params = a + (b-a).*rand(1, n);
%         params = linspace(a, b, n);

    end
end