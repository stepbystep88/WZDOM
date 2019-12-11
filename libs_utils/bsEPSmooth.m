function [filteredData] = bsELPSmooth(data, varargin)
%% smooth data by using the edge-preserving smoothing method
% see paper https://library.seg.org/doi/pdf/10.1190/1.2213050 for details
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% 
% Input
% options.windowSize        the size of the sliding 2D window
% options.windowShape       the type of window, see details of the paper ...
% https://library.seg.org/doi/pdf/10.1190/1.2213050
%
% Output
% filteredData              the filtered data
% -------------------------------------------------------------------------

    p = inputParser;
    
    
    addParameter(p, 'windowSize', 5);
    addParameter(p, 'windowShape', 'cubic');
    addParameter(p, 'isParallel', 1);   
    addParameter(p, 'numWorkers', bsGetMaxNumWorkers());
   
    p.parse(varargin{:});  
    options = p.Results;
    
   
    switch length(size(data))
        case 1
            
        case 2
            [panels] = bs2DPanels(options.windowSize);
            filteredData = bsELP2D(data, options.windowSize, panels);
            
        case 3
%             [filteredData] = bsELP3D(data, options);
        otherwise
            error('The input data to processed by NLM filter must be a 1D, 2D or a 3D data.');
    end

    
    
end

function [data] = bsELP2D(data, n, panels)
    [n1, n2] = size(data);
    
    for j = 1 : n1
        for i = 1 : n2
            data(j, i) = bs2DGetValue(data, i, j, n1, n2, n, panels);
        end
    end
end

function [panels] = bs2DPanels(n)
    panels = cell(1, 9);
    
    half = (n-1) / 2;
    np = (n-2)*half + 1;
    panels{1} = zeros(2, np);
    iter = 0;
    for i = -half : 0
        for j = -half+1:half-1
            if i~= 0 && j~= 0
                iter = iter + 1;
                panels{1}(:, iter) = [j; i];
            end
        end
    end
    
    np = (half+1).^2 - 2;
    panels{2} = zeros(2, np);
    iter = 0;
    for i = 0 : half
        for j = -half : 0
            if (i==0 && j==-half) || (i==half && j==0)
            else
                iter = iter + 1;
                panels{1}(:, iter) = [j; i];
            end
        end
    end
    theta = 90;
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    for i = 1 : 4
        panels{1+i*2} = R * panels{i*2-1};
        panels{2+i*2} = R * panels{i*2};
    end
    
    np = (n-2).^2;
    panels{9} = zeros(2, np);
    iter = 0;
    for i = -half : half
        for j = -half : half
            if i ~= half && j == half
                iter = iter + 1;
                panels{1}(:, iter) = [j; i];
            end
        end
    end
end

function idata = bs2DGetValue(data, i, j, n1, n2, n, panels)
    half = (n-1) / 2;
    nPanel = length(panels);
    central = inf(n*n, nPanel);
    variances = inf(n*n, nPanel);
    
    for i1 = half : half
        for j1 = half : half
            
            if ~bsCheckValidIndex(j1+j, i1+i, n1, n2)
                central(j1, i1) = data(j1+j, i1+i);
                
                for k = 1 : nPanel
                    panel = panels{k};
                    panel = [panel(1, :) + j1; panel(2, :) + i1];
                    index = bsCheckValidIndex(panel(1, :), panel(2, :), n1, n2);
                    panel(:, index) = [];

                    kdata = zeros(1, length(panel));
                    for m = 1 : size(panel, 2)
                    	kdata(m) = data(panel(1, m), panel(2, m));
                    end
                    
                    variances(j1, i1) = var(kdata);
                end
            end
        end
    end
    
    [~, index] = min(variances(:));
    idata = central(index);
end

function index = bsCheckValidIndex(j, i, n1, n2)
    index = i < 1 | j < 1 | j > n1 | i > n2;
end
