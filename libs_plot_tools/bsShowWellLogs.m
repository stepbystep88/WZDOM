function bsShowWellLogs(wellLogs, timeIndex, dataIndex, types, varargin)
    
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
%                 subplot(1, nSubFigure, i);
                bsSubPlotFit(1, nSubFigure, i, 0.9, 0.85, 0.05, 0.05, 0.1, -0.01);
            otherwise
                subplot(2, ceil(nSubFigure/2), i);
        end
        
        wellInfo = wellLogs{i};
        time = wellInfo.wellLog(:, timeIndex)/1000;
        
        for j = 1 : nData
            jData = wellInfo.wellLog(:, dataIndex(j));
            if options.isNormal
                jData = normalize(jData, 'range');
            end
            
            
            
            plot(jData, time, 'color', options.colors{j}, 'linewidth', 2); hold on;
        end
        
        set(gca, 'ydir', 'reverse');
        if ~isempty(options.range) && ~options.isNormal
            set(gca, 'xlim', options.range);
        end
        set(gca, 'ylim', [time(1) time(end)]);
        
        if i == 1
            legend(types, 'fontweight', 'bold');
            ylabel('Time (s)');
        end

        xlabel('Amplitude');
        
%         xlabel('');
        if isfield(wellInfo, 'name')
            title(sprintf('(%s) %s', 'a'+i-1, wellInfo.name));
        elseif isfield(wellInfo, 'wellName')
            title(wellInfo.wellName);
        end

        bsSetDefaultPlotSet(bsGetDefaultPlotSet());
    end
end