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
    % ���ڶ����ͬʱϡ���ʾ�ĸ���
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
            fprintf('���ݽ���м�⵽��%d���������ֱ�Ϊ:\n\tinIds:%s...\n\tcrossIds:%s\n\twellNames=%s\n', ...
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
        % �෴����ֱ�����ø�����options.invDataAtWellLocation��ѵ���ֵ�
        options.trainNum = length(wellLogs);
        train_ids = 1 : length(wellLogs);
        
        assert(size(options.invDataAtWellLocation.data{1}, 2) == length(wellLogs), 'ѵ����������Ӧ�������Ӧ�ķ��ݽ������һ��');
%         [outLogs] = bsGetPairOfInvAndWell(GInvParam, timeLine, wellLogs, options.invDataAtWellLocation, GInvParam.indexInWellData.ip, options);
        
    end
    
    outResult = bsSetFields(invResult, {'name', name});
    
    %% ��ǰ�ж�������
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
        
        fprintf('ѵ�������ֵ���...\n');
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
%         legend('�ع����', 'ʵ�ʲ⾮', '���ݽ��', 'fontsize', 11);
%         set(gcf, 'position', [261   558   979   420]);
%         bsShowFFTResultsComparison(GInvParam.dt, [outLogs{2}.wellLog(:, 1), outLogs{2}.wellLog(:, 2), testData(:, 2)], {'���ݽ��', 'ʵ�ʲ⾮', '�ع����'});
        
%         startTime = invResult.horizon - GInvParam.upNum * GInvParam.dt;
%         sampNum = GInvParam.upNum + GInvParam.downNum;
%         inputData = bsStackPreSeisData(GInvParam.preSeisData.fileName, GInvParam.preSeisData.segyInfo, ...
%             invResult.inIds, invResult.crossIds, startTime, sampNum, GInvParam.dt);
    
        % �����ֵ�ϡ���ع�
        fprintf('�����ֵ�ϡ���ع�: �ο�����Ϊ���ݽ��...\n');
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
            fprintf('������Ϊsegy�ļ�ʧ��...\n');
        end
    end
    
end