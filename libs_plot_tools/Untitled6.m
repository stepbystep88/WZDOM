function [colormap_K] = bsHex2RGB(separate)
    ncolor = length(separate);
    colormap_K = zeros(ncolor, 3);
    
    for i = 1 : ncolor
        colormap_K(i, :) = hex2rgb(separate{i});
    end
    
end