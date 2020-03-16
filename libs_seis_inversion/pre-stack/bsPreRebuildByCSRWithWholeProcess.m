function outResult = bsPreRebuildByCSRWithWholeProcess(GInvParam, timeLine, wellLogs, method, invResult, name, varargin)

    p = inputParser;
    GTrainDICParam = bsCreateGTrainDICParam(...
        'csr', ...
        'title', 'high_freq', ...
        'sizeAtom', 90, ...
        'sparsity', 5, ...
        'iterNum', 5, ...
        'nAtom', 4000, ...
        'filtCoef', 1);
    

    addParameter(p, 'mode', 'full_freq');
    addParameter(p, 'wellFiltCoef', 0.1);
    addParameter(p, 'lowCut', 0.1);
    addParameter(p, 'highCut', 1);
    addParameter(p, 'sparsity', 5);
    addParameter(p, 'gamma', 0.5);
    addParameter(p, 'title', 'HLF');
    addParameter(p, 'trainNum', length(wellLogs));
    addParameter(p, 'exception', []);
    addParameter(p, 'mustInclude', []);
    addParameter(p, 'GTrainDICParam', GTrainDICParam);
    
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
            GTrainDICParam.normailzationMode = 'whole_data_max_min';
            options.gamma = 1;
            DICTitle = sprintf('%s_lowCut_%.2f_highCut_%.2f', ...
                GTrainDICParam.title, options.lowCut, options.highCut);
    end
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];

    % inverse a profile
%     [~, ~, GInvParamWell] = bsBuildInitModel(GInvParam, timeLine, wellLogs, ...
%         'title', 'all_wells', ...
%         'inIds', wellInIds, ...
%         'filtCoef', options.wellFiltCoef, ...
%         'isRebuild', 1, ...
%         'crossIds', wellCrossIds ...
%     );

    fprintf('反演所有测井中...\n');
    wellInvResults = bsPreInvTrueMultiTraces(GInvParam, wellInIds, wellCrossIds, timeLine, {method});
    [wellInvResults, ~, ~] = bsPreGetOtherAttributesByInvResults(wellInvResults, GInvParam, wellLogs);
    
    set_diff = setdiff(1:length(wellInIds), options.exception);
    train_ids = bsRandSeq(set_diff, options.trainNum);
    train_ids = unique([train_ids, options.mustInclude]);
    outResult = bsSetFields(invResult, {'name', name});
    
    %% 叠前有三种属性
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
        
        [outLogs] = bsGetPairOfInvAndWell(GInvParam, timeLine, wellLogs, wellInvResults{1}.data{iData}, dataIndex, options);
        fprintf('训练联合字典中...\n');
        [DIC, train_ids, rangeCoef] = bsTrainDics(GTrainDICParam, outLogs, train_ids, ...
            [ 1, 2], GTrainDICParam.isRebuild);
        GInvWellSparse = bsCreateGSparseInvParam(DIC, GTrainDICParam, ...
            'sparsity', options.sparsity, ...
            'stride', 1);
    
        options.rangeCoef = rangeCoef;
        [testData] = bsPostReBuildByCSR(GInvParam, GInvWellSparse, wellInvResults{1}.data{iData}, options);
    
%         figure; plot(testData(:, 1), 'r', 'linewidth', 2); hold on; 
%         plot(outLogs{1}.wellLog(:, 2), 'k', 'linewidth', 2); 
%         plot(outLogs{1}.wellLog(:, 1), 'b', 'linewidth', 2);
%         legend('重构结果', '实际测井', '反演结果', 'fontsize', 11);
%         set(gcf, 'position', [261   558   979   420]);
%         bsShowFFTResultsComparison(GInvParam.dt, [outLogs{1}.wellLog, testData(:, 1)], {'反演结果', '实际测井', '重构结果'});

        % 联合字典稀疏重构
        fprintf('联合字典稀疏重构: 参考数据为反演结果...\n');
        [outputData] = bsPostReBuildByCSR(GInvParam, GInvWellSparse, invResult.data{i}, options);

        outResult.data{i} = outputData;

    end
    
    bsWriteInvResultsIntoSegyFiles(GInvParam, {outResult}, options.title, 0);

end