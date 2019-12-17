function bsShowWellLogs(wellLogs, dataIndex, types, varargin)
    
    p = inputParser;
    
    addParameter(p, 'showNum', 12);
    addParameter(p, 'isNormal', 1);
    addParameter(p, 'range', []);
    addParameter(p, 'colors', bsGetColormap('separate'));
    
    p.parse(varargin{:});  
    options = p.Results;
    
    wellNum = length(wellLogs);
    nData = length(dataIndex);
    
    figure;
    bsSetPosition(0.78, 0.56);

    assert(length(dataIndex) == length(types));
    
    nSubFigure = min(wellNum, options.showNum);
    for i = 1 : nSubFigure

        switch nSubFigure
            case {1, 2, 3, 4, 5, 6}
                subplot(1, nSubFigure, i);
            otherwise
                subplot(2, ceil(nSubFigure/2), i);
        end
        
        wellInfo = wellLogs{i};
        for j = 1 : nData
            jData = wellInfo.wellLog(:, dataIndex(j));
            if options.isNormal
                jData = normalize(jData, 'range');
            end
            
            plot(jData, 1:length(jData), 'color', options.colors{j}, 'linewidth', 2); hold on;
        end
        
        set(gca, 'ydir', 'reverse');
        if ~isempty(options.range) && ~options.isNormal
            set(gca, 'xlim', options.range());
        end
        
        if i == 1
            legend(types);
            ylabel('Sample number');
        end

        if isfield(wellInfo, 'name')
            title(wellInfo.name);
        elseif isfield(wellInfo, 'wellName')
            title(wellInfo.wellName);
        end

        bsSetDefaultPlotSet(bsGetDefaultPlotSet());
    end
end