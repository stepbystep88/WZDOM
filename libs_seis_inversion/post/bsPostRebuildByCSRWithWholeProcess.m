function outResult = bsPostRebuildByCSRWithWholeProcess(GInvParam, timeLine, wellLogs, methods, invResult, name, varargin)

    p = inputParser;
    GTrainDICParam = bsCreateGTrainDICParam(...
        'csr', ...
        'title', '', ...
        'sizeAtom', 90, ...
        'sparsity', 5, ...
        'iterNum', 5, ...
        'nAtom', 4000, ...
        'filtCoef', 1);
    

    addParameter(p, 'mode', 'low_high');
    addParameter(p, 'lowCut', 0.1);
    addParameter(p, 'highCut', 1);
    addParameter(p, 'sparsity', 5);
    addParameter(p, 'gamma', 0.5);
    addParameter(p, 'wellFiltCoef', 0.1);
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
            GTrainDICParam.title = sprintf('%s_highCut_%.2f', ...
                GTrainDICParam.title, options.highCut);
        case 'low_high'
            GTrainDICParam.normailzationMode = 'whole_data_max_min';
            options.gamma = 1;
            GTrainDICParam.title = sprintf('%s_lowCut_%.2f_highCut_%.2f', ...
                GTrainDICParam.title, options.lowCut, options.highCut);
    end
    
    wells = cell2mat(wellLogs);
    wellInIds = [wells.inline];
    wellCrossIds = [wells.crossline];

    % inverse a profile
    [~, ~, GInvParamWell] = bsBuildInitModel(GInvParam, timeLine, wellLogs, ...
        'title', 'all_wells', ...
        'inIds', wellInIds, ...
        'isRebuild', 1, ...
        'filtCoef', options.wellFiltCoef, ...
        'crossIds', wellCrossIds ...
    );

    fprintf('反演所有测井中...\n');
    wellInvResults = bsPostInvTrueMultiTraces(GInvParamWell, wellInIds, wellCrossIds, timeLine, methods);
    
    set_diff = setdiff(1:length(wellInIds), options.exception);
    train_ids = bsRandSeq(set_diff, options.trainNum);
    train_ids = unique([train_ids, options.mustInclude]);
    [outLogs] = bsGetPairOfInvAndWell(GInvParam, timeLine, wellLogs, wellInvResults{1}.data, ...
        GInvParam.indexInWellData.ip, options);
    
%     bsShowFFTResultsComparison(GInvParam.dt, outLogs{1}.wellLog, {'反演结果', '联合字典配对'});


    %% 训练字典
    fprintf('训练联合字典中...\n');
    
    [DIC, train_ids, rangeCoef] = bsTrainDics(GTrainDICParam, outLogs, train_ids, ...
        [ 1, 2], GTrainDICParam.isRebuild);

    GInvWellSparse = bsCreateGSparseInvParam(DIC, GTrainDICParam, ...
        'sparsity', options.sparsity, ...
        'stride', 1);
    
    options.rangeCoef = rangeCoef;
    [testData] = bsPostReBuildByCSR(GInvParam, GInvWellSparse, wellInvResults{1}.data, options);
    
%     figure; plot(testData(:, 1), 'r', 'linewidth', 2); hold on; 
%     plot(outLogs{1}.wellLog(:, 2), 'k', 'linewidth', 2); 
%     plot(outLogs{1}.wellLog(:, 1), 'b', 'linewidth', 2);
%     legend('重构结果', '实际测井', '反演结果', 'fontsize', 11);
%     set(gcf, 'position', [261   558   979   420]);
%     bsShowFFTResultsComparison(GInvParam.dt, [outLogs{10}.wellLog, testData(:, 10)], {'反演结果', '实际测井', '重构结果'});
    
    %% 联合字典稀疏重构
    fprintf('联合字典稀疏重构: 参考数据为反演结果...\n');
    [outputData] = bsPostReBuildByCSR(GInvParam, GInvWellSparse, invResult.data, options);
    
    
    
    outResult = bsSetFields(invResult, {'data', outputData; 'name', name});
    
                        
    bsWriteInvResultsIntoSegyFiles(GInvParam, {outResult}, options.title);
end