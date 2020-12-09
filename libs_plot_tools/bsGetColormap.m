function colormap = bsGetColormap(name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get colormap by name
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------

    switch lower(name)
        case {'vel', 'velocity', 'pvel'}
            tmp = load('color_hrs_pvel.mat');
            colormap = tmp.color_hrs_pvel;
        case {'hrs'}
            tmp = load('color_hrs.mat');
            colormap = tmp.color_hrs;
            
        case {'jason_impedance', 'jason'}
            tmp = load('color_jason_impedance.mat');
            colormap = tmp.color_jason_impedance;
            
        case {'jason_impedance2', 'jason2'}
            tmp = load('color_jason_impedance2.mat');
            colormap = tmp.color_jason_impedance2;
            
        case 'seismic'
            tmp = load('color_seismic.mat');
            colormap = tmp.color_seismic;
            
        case {'original', [], ''}
            tmp = load('color_original.mat');
            colormap = tmp.color_original;
            
        case {'distinction', 'distinguish', 'separate'}
            tmp = load('colorTbl.mat');
            colormap = tmp.colors;
        otherwise
            tmp = load('color_original.mat');
            colormap = tmp.color_original;
    end
end