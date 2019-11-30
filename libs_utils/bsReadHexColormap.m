function colormap = bsReadHexColormap(fileName)
    data = importdata(fileName);
    
    colormap = zeros(length(data), 3);
    
    for i = 1 : length(data)
        colormap(i, :) = hex2rgb(sprintf('#%s', data{i}));
    end
    
end