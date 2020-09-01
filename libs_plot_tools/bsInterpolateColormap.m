function [res] = bsInterpolateColormap(name, nColor, colormap, position)

    if ~exist('nColor', 'var')
        nColor = 256;
    end
    
    if ~exist('colormap', 'var')
        colormap = get(gcf,'Colormap');
    end
    
    if ~exist('position', 'var')
        x1 = round(linspace(1, nColor, size(colormap, 1)));
    else
        x1 = position;
    end
    
    x = 1 : nColor;
    
    newColor = zeros(nColor, 3);
    
    for i = 1 : 3
        newColor(:, i) = interp1(x1, colormap(:, i), x, 'linear');
    end
    
    eval([name, '=', mat2str(newColor)]);
    
    res = eval(name);
    
    save(name, name);
end