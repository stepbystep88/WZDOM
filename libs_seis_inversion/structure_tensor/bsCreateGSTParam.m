function options = bsCreateGSTParam(dim, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create options for function bsSmoothByGST2D.m and bsSmoothByGST3D.m
%
% Copyright (C) 2020. Bin She. All rights reserved.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Sep 2020
% -------------------------------------------------------------------------

% see function compute_structure_tensor.m
    
    if isempty(dim)
        dim = 2;
    end
    
    p = inputParser;
    
    % the sparsity controlling the sparse representation 
    addParameter(p, 'sub', 8); 
    addParameter(p, 'len', 4); 
    addParameter(p, 'sigma', 2); 
    addParameter(p, 'show_mid_results', false); 
    
    addParameter(p, 'm', 4); 
    addParameter(p, 'Cm', 3.31488); 
    addParameter(p, 'lambda', 1e-4); 
    addParameter(p, 'tau', 0.1); 
    addParameter(p, 'iterNum', 25);
    addParameter(p, 'ttv', 0);
    addParameter(p, 'colormap', bsGetColormap('seismic'));
    
    goptions.order = 2;
    addParameter(p, 'nabla', @(f)grad(f, goptions));
    
    
%     @(f)medfilt2(f, [2, 2])
%     @(f)filter2(fspecial('average', [2,2]), f, 'same')
    switch dim
        case 2
            addParameter(p, 'filter_fcn', @(f)imfilter(f, fspecial('gaussian', [2,2]), 'replicate'));
        case 3
            addParameter(p, 'filter_fcn', @(f)imgaussfilt3(f, 1));
    end
    
    p.parse(varargin{:});  
    options = p.Results;
    
end