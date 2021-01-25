function outResult = bsPreRebuildByCSRWithWholeProcess(GInvParam, timeLine, wellLogs, invResult, name, varargin)

    p = inputParser;
    GTrainDICParam = bsCreateGTrainDICParam(...
        'csr', ...
        'title', 'high_freq', ...
        'sizeAtom', 90, ...
        'sparsity', 5, ...
        'iterNum', 5, ...
        'nAtom', 4000, ...
        'filtCoef', 1);
    
    is3D = length(invResult.inIds) > 5000;
    
    addParameter(p, 'ratio_to_reconstruction', 1);
    addParameter(p, 'mode', 'full_freq');
    addParameter(p, 'isSaveSegy', false);
    addParameter(p, 'lowCut', 0.1);
    addParameter(p, 'highCut', 1);
    addParameter(p, 'sparsity', 5);
    addParameter(p, 'gamma', 0.5);
    addParameter(p, 'title', 'HLF');
    addParameter(p, 'trainNum', length(wellLogs));
    addParameter(p, 'exception', []);
    addParameter(p, 'mustInclude', []);
    addParameter(p, 'GTrainDICParam', GTrainDICParam);
    addParameter(p, 'is3D', is3D);
    % 相邻多个块同时稀疏表示的个数
    addParameter(p, 'nMultipleTrace', 1);
    addParameter(p, 'isInterpolation', 0);
    addParameter(p, 'invDataAtWellLocation', []);
    
    if is3D
        addParameter(p, 'gst_options', bsCreateGSTParam(3));
    else
        addParameter(p, 'gst_options', bsCreateGSTParam(2));
    end
    
    p.parse(varargin{:});  
    options = p.Results;
    GTrainDICParam = options.GTrainDICParam;
    GTrainDICParam.filtCoef = 1;
    
    switch options.mode
        case 'full_freq'
            GTrainDICParam.normailzationMode = 'off';
            DICTitle = sprintf('%s_highCut_%.2f', ...
                GTrainDICParam.title, options.highCut);
        case 'low_high'
%             GTrainDICParam.normailzationMode = 'whole_data_max_min';
            options.gamma = 1;
            DICTitle = sprintf('%s_lowCut_%.2f_highCut_%.2f', ...
                GTrainDICParam.title, options.lowCut, options.highCut);
    end
    
    [wellPos, wellIndex, wellNames] = bsFindWellLocation(wellLogs, invResult.inIds, invResult.crossIds);
    
    if isempty(options.invDataAtWellLocation)
        
        names = join(wellNames, ',');
        try
            fprintf('反演结果中检测到共%d道过井，分别为:\n\tinIds:%s...\n\tcrossIds:%s\n\twellNames=%s\n', ...
                length(wellPos), ...
                mat2str(invResult.inIds(wellPos)), ...
                mat2str(invResult.crossIds(wellPos)), ...
                names{1});

            wellLogs = wellLogs(wellIndex);
        catch
            error('There is no wells in the inverted data!!!');
        end


        try
            set_diff = setdiff(1:length(wellPos), options.exception);
            train_ids = bsRandSeq(set_diff, options.trainNum);
            train_ids = unique([train_ids, options.mustInclude]);

        catch
            options.trainNum = length(wellPos);
            train_ids = 1 : length(wellPos);
        end
    else
        % 相反，则直接运用给定的options.invDataAtWellLocation来训练字典
        options.trainNum = length(wellLogs);
        train_ids = 1 : length(wellLogs);
        
        assert(size(options.invDataAtWellLocation.data{1}, 2) == length(wellLogs), '训练井的数量应该与其对应的反演结果道数一致');
%         [outLogs] = bsGetPairOfInvAndWell(GInvParam, timeLine, wellLogs, options.invDataAtWellLocation, GInvParam.indexInWellData.ip, options);
        
    end
    
    outResult = bsSetFields(invResult, {'name', name});
    
    %% 叠前有多种属性
    for i = 1 : length(invResult.type)
        switch lower(invResult.type{i})
            case 'vp'
                dataIndex = GInvParam.indexInWellData.vp;
                GTrainDICParam.title = ['vp_', DICTitle];
                iData = 1;
            case 'vs'
                dataIndex = GInvParam.indexInWellData.vs;
                GTrainDICParam.title = ['vs_', DICTitle];
                iData = 2;
            case 'rho'
                dataIndex = GInvParam.indexInWellData.rho;
                GTrainDICParam.title = ['rho_', DICTitle];
                iData = 3;
            case {'vpvs_ratio', 'vp_vs'}
                dataIndex = GInvParam.indexInWellData.vpvs_ratio;
                GTrainDICParam.title = ['vpvs_ratio_', DICTitle];
                iData = 4;
            case 'possion'
                dataIndex = GInvParam.indexInWellData.possion;
                GTrainDICParam.title = ['possion_', DICTitle];
                iData = 5;
        end
        
        
        if isempty(options.invDataAtWellLocation)
            [outLogs] = bsGetPairOfInvAndWell(GInvParam, timeLine, wellLogs, invResult.data{iData}(:, wellPos), dataIndex, options);
        else
            [outLogs] = bsGetPairOfInvAndWell(GInvParam, timeLine, wellLogs, options.invDataAtWellLocation.data{iData}, dataIndex, options);
        end
        
        fprintf('训练联合字典中...\n');
        [DIC, train_ids, rangeCoef, output] = bsTrainDics(GTrainDICParam, outLogs, train_ids, [ 1, 2]);
        GInvWellSparse = bsCreateGSparseInvParam(DIC, GTrainDICParam, ...
            'sparsity', options.sparsity, ...
            'stride', 1);
        
        
        
        GInvWellSparse.rangeCoef = rangeCoef;
        GInvWellSparse.output = output;
        GInvWellSparse.wellPos = wellPos;
%         [testData] = bsPostReBuildByCSR(GInvParam, GInvWellSparse, wellInvResults{1}.data{iData}, wellInIds, wellCrossIds, options);
       
%         figure; 
%         plot(testData(:, 2), 'r', 'linewidth', 2); hold on; 
%         plot(outLogs{1}.wellLog(:, 2), 'k', 'linewidth', 2); 
%         plot(outLogs{1}.wellLog(:, 1), 'b', 'linewidth', 2);
%         legend('重构结果', '实际测井', '反演结果', 'fontsize', 11);
%         set(gcf, 'position', [261   558   979   420]);
%         bsShowFFTResultsComparison(GInvParam.dt, [outLogs{2}.wellLog(:, 1), outLogs{2}.wellLog(:, 2), testData(:, 2)], {'反演结果', '实际测井', '重构结果'});
        
%         startTime = invResult.horizon - GInvParam.upNum * GInvParam.dt;
%         sampNum = GInvParam.upNum + GInvParam.downNum;
%         inputData = bsStackPreSeisData(GInvParam.preSeisData.fileName, GInvParam.preSeisData.segyInfo, ...
%             invResult.inIds, invResult.crossIds, startTime, sampNum, GInvParam.dt);
    
        % 联合字典稀疏重构
        fprintf('联合字典稀疏重构: 参考数据为反演结果...\n');
        if ~options.isInterpolation
%             if options.nMultipleTrace <= 1
%                 [outputData, highData, gamma_vals, gamma_locs] = bsPostReBuildPlusInterpolationCSR(GInvParam, GInvWellSparse, invResult.data{i}, invResult.data{i},...
%                     invResult.inIds, invResult.crossIds, options);
%             else
                [outputData, highData, gamma_vals, gamma_locs] = bsPostReBuildMulTraceCSR(GInvParam, GInvWellSparse, invResult.data{i}, invResult.data{1}, invResult.inIds, invResult.crossIds, options);
%             end
        else
            gamma_locs = [];
            gamma_vals = [];
            [outputData, highData] = bsPostReBuildInterpolation(GInvParam, outLogs(train_ids), invResult.data{i}, invResult.inIds, invResult.crossIds, options);
        end
    
        

        outResult.data{i} = outputData;
        outResult.high_freq_data{i} = highData;
        outResult.gamma_vals{i} = gamma_vals;
        outResult.gamma_locs{i} = gamma_locs;
        outResult.DIC{i} = DIC;
    end
    
    if options.isSaveSegy
        try
            bsWriteInvResultsIntoSegyFiles(GInvParam, {outResult}, options.title, 0);
        catch
            fprintf('保存结果为segy文件失败...\n');
        end
    end
    
end