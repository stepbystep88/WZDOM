function options = bsCreateGShowProfileParam(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a struct to save the information of a segy file 
%
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% -------------------------------------------------------------------------
    p = inputParser;
    
    % how many x label ticks are shown in profile
    addParameter(p, 'xLabelNum', 10); 
    
    % how many y label ticks are shown in profile
    addParameter(p, 'yLabelNum', 8); 
    
    addParameter(p, 'emptyOffset', 2); 
    
    addParameter(p, 'language', 'zh'); 
    
    % 
%     addParameter(p, 'isLegend', 1); 
    
    % the range of different attributes, each range is a vecotr of 2
    % elements [x0, x1]
    addParameter(p, 'range', struct(...
        'seismic', [], ...
        'ip', [], ...
        'vp', [], ...
        'vs', [], ...
        'rho', [], ...
        'brittleness', [], ...
        'vpvs_ratio', [], ...
        'possion', [], ...
        'toc', [] ...
    )); 
    
    % the coloramp for showing a profile
    vel_colormap = bsGetColormap('velocity');
    addParameter(p, 'colormap', struct(...
        'seismic', bsGetColormap('seismic'), ...
        'ip', vel_colormap, ...
        'vp', vel_colormap, ...
        'vs', vel_colormap, ...
        'rho', vel_colormap, ...
        'brittleness', vel_colormap, ...
        'vpvs_ratio', vel_colormap, ...
        'possion', vel_colormap, ...
        'toc', vel_colormap, ...
        'allTheSame', [] ...
    )); 
    
    % the basic settings including fontsize, fontname, etc.
    addParameter(p, 'plotParam', bsGetDefaultPlotSet()); 
    
    % whether to display the horizon
    addParameter(p, 'isShowHorizon', 1); 
    
    % whether to display the colored well
    addParameter(p, 'isShowColorWells', 1);
    addParameter(p, 'isShowWellNames', 1);
    
    % the coeficient of low-pass fiter for showing the welllog data
    addParameter(p, 'showWellFiltCoef', 1); 
    
    % the thickness of a colored well shown in a profile
    addParameter(p, 'showWellOffset', 3); 
    
    addParameter(p, 'showWellMode', 'color');   % should be color or curve 
    
    % the scale factor of showing a profile, must be a non-negative integer
    addParameter(p, 'scaleFactor', 2);
    
    addParameter(p, 'isScaleHorizon', 1);
    
    % the scale factor of showing a profile, must be a non-negative integer
    addParameter(p, 'edgeOffsetNum', 5);
    
    % whether set the colormap as reverse mode
    addParameter(p, 'isColorReverse', 0);

    % how many traces are displayed with respect to the left-most well location
    addParameter(p, 'showLeftTrNumByWells', 10000);
    
    % how many traces are displayed with respect to the right-most well location
    addParameter(p, 'showRightTrNumByWells', 10000);

    
%     addParameter(p, 'isShowOneFigure', 0);
%     
%     addParameter(p, 'figureOptions', struct(...
%         'nRow', [],...
%         'nCol', 3, ...
%         'loc', [0.89, 0.95, 0.025, 0.05, 0.05, 0.01], ...
%         'colorbar_pos', [0.92 0.05 0.01 0.92]...
%         ));
    
    % along time, which vertical part of a profile will be display, see
    % function of bsShowVerticalPartData in bsShowInvProfiles.m for details
    addParameter(p, 'showPartVert', struct( ...
        'mode', 'off', ...
        'downTime', 0, ...  
        'upTime', 0, ...
        'horizonIds', [] ...
    ));

    p.parse(varargin{:});  
    options = p.Results;
    
end

