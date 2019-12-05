function [filteredData, weightInfo] = bsNLMByRef(data, refData, varargin)
%% smooth data by using the NLM algorithm, the similarites are referenced from refData
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Dec 2019
% 
% Input
% options.p                 parameter for inverse distance weight function,
%                           its reasonable range is [0.5 3]
% options.is3D              whether the data is a 3D data
% options.nPointsUsed       the most number of points used to calculate the
%                           weight information
% options.stride            step size of sliding window
% options.searchOffset      indicating the search range
% options.windowSize        the size of the sliding 2D window

%
% Output
% filteredData              the filtered data
% -------------------------------------------------------------------------

    if nargin == 1 || isempty(refData)
        % when there is no reference data as input, we set the reference
        % data as the original data itself
        refData = data;
    end
    
    p = inputParser;
    
    addParameter(p, 'searchOffset', 3);
    addParameter(p, 'stride', [1 1 1]);
    addParameter(p, 'searchStride', [1 2 2]);
    addParameter(p, 'windowSize', 5);
    addParameter(p, 'nPointsUsed', 5);
    addParameter(p, 'weightInfo', []);
    addParameter(p, 'isParallel', 1);   
    addParameter(p, 'numWorkers', bsGetMaxNumWorkers());
    addParameter(p, 'p', 2);
   
    p.parse(varargin{:});  
    options = p.Results;
    
    assert(isequal(size(data), size(refData)), 'data and refData must have the same size.');
    
    switch length(size(data))
        case 2
            is3D = false;
        case 3
            is3D = true;
        otherwise
            error('The input data to processed by NLM filter must be a 2D profile or a 3D volume.');
    end
    
    if is3D
        [filteredData, weightInfo] = bsNLM3D(data, refData, options);
    else
        [filteredData, weightInfo] = bsNLM2D(data, refData, options);
    end
    
    
end


function [filteredData, weightInfo] = bsNLM3D(data, refData, options)
    
    ws = options.windowSize;
    
    [o1, o2, o3] = size(data);
    p1 = round(1 : options.stride(1) : o1-ws+1);
    p2 = round(1 : options.stride(2) : o2-ws+1);
    p3 = round(1 : options.stride(3) : o3-ws+1);
    
    n1 = length(p1);
    n2 = length(p2);
    n3 = length(p3);
    nPatches = n1 * n2 * n3;
    patches = cell(1, nPatches);
    refPatches = cell(1, nPatches);
    
    for i = 1 : nPatches
        if mod(i, 1000000) == 0
            fprintf('Cutting data into patches %.2f%%...\n', i/nPatches*100);
        end
        [i1, i2, i3] = bsGetId3D(i, n1, n2, n3);
        patches{i} = data(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1, p3(i3):p3(i3)+ws-1);

        if isempty(options.weightInfo)
            refPatches{i} = refData(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1, p3(i3):p3(i3)+ws-1);
        end
    end
    
    if isempty(options.weightInfo)
        [weights, index] = bsGetKNearestWeights3D(refPatches, n1, n2, n3, options);
        weightInfo.weights = weights;
        weightInfo.index = index;
    else
        weightInfo = options.weightInfo;
        weights = weightInfo.weights;
        index = weightInfo.index;
    end
    
    K = options.nPointsUsed;
    filteredData = zeros(o1, o2, o3);
    filteredNum = zeros(o1, o2, o3);
    
    % average data by calculated patches
    for i = 1 : nPatches
        
        if mod(i, 10000) == 0
            fprintf('Averaging patches as a whole volume %.2f%%...\n', i/nPatches*100);
        end
        
        [i1, i2, i3] = bsGetId3D(i, n1, n2, n3);
        
        iIndex = index(:, i);
        weight = weights(:, i);

        avgPatch = zeros(ws, ws);
        for k = 1 : K
            avgPatch = avgPatch + weight(k) * patches{iIndex(k)};
        end

        filteredData(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1, p3(i3):p3(i3)+ws-1) = ...
            filteredData(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1, p3(i3):p3(i3)+ws-1)...
            + avgPatch;
        filteredNum(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1, p3(i3):p3(i3)+ws-1) = ...
            filteredNum(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1, p3(i3):p3(i3)+ws-1) + 1;
    end
    
    filteredData = filteredData ./ filteredNum;
end

function [i1, i2, i3] = bsGetId3D(i, n1, n2, n3)

    i1 = i / (n2 * n3);
    if floor(i1) == i1
        i2 = n2;
        i3 = n3;
    else
        i1 = floor(i1) + 1;
        i = i - (i1-1) * (n2 * n3);
        
        i2 = i / n3;
        if floor(i2) == i2
            i3 = n3;
        else
            i2 = floor(i2) + 1;
            i3 = i - (i2-1) * n3;
        end
    end
    
end

function [weights, index] = bsGetKNearestWeights3D(patches, n1, n2, n3, options)
    nPatches = n1 * n2 * n3;
    
    N = options.searchOffset;
    K = options.nPointsUsed;
    
    weights = zeros(K, nPatches);
    index = zeros(K, nPatches);
    normCoef = norm(patches{1}(:));
    
    fprintf('Calculating the weight information...\n');
    p = options.p;
    searchStride = options.searchStride;
    
    if options.isParallel
        pbm = bsInitParforProgress(options.numWorkers, nPatches, 'Calculating the weight information', [], 1);
        
        parfor i = 1 : nPatches
            [weights(:, i), index(:, i)] = bsGetIWeight3D(patches, normCoef, i, n1, n2, n3, N, K, -p, searchStride);
            bsIncParforProgress(pbm, i, 2000);
        end
    else
        for i = 1 : nPatches
            if mod(i, 2000) == 0
                fprintf('Calculating the weight information %.2f%%...\n', i/nPatches*100);
            end
            [weights(:, i), index(:, i)] = bsGetIWeight3D(patches, normCoef, i, n1, n2, n3, N, K, -p, searchStride);
        end
    end
end

function [weight, index] = bsGetIWeight3D(patches, normCoef, i, n1, n2, n3, N, K, e, searchStride)
    [i1, i2, i3] = bsGetId3D(i, n1, n2, n3);
    
    D = inf((2*N+1)^2, 2);
    iter = 0;

    for c1 = i1-N : searchStride(1) : i1+N
        if c1<=0 || c1>n1
            continue;
        end

        for c2 = i2-N : searchStride(2) : i2+N
            if c2<=0 || c2>n2
                continue;
            end
            
            for c3 = i3-N : searchStride(3) : i3+N
                if c3<=0 || c3>n3 || (c2==i2 && c1==i1 && c3==i3)
                    continue;
                end
            
                iter = iter + 1;

                j = (c1-1)*(n2*n3) + (c2-1)*n3 + c3;

                difference = patches{j} - patches{i};
                D(iter, 1) = norm(difference(:)/normCoef);
                D(iter, 2) = j;
            end
        end
    end

    [B, ~] = bsMinK(D, K, 1);
    weight = B(:, 1).^e;
    weight = weight / sum(weight);
    index = B(:, 2);
end

function [filteredData, weightInfo] = bsNLM2D(data, refData, options)
    ws = options.windowSize;
    [o1, o2] = size(data);
    p1 = round(1 : options.stride(1) : o1-ws+1);
    p2 = round(1 : options.stride(2) : o2-ws+1);
    
    n1 = length(p1);
    n2 = length(p2);
    nPatches = n1 * n2;
    patches = cell(1, nPatches);
    refPatches = cell(1, nPatches);
    
    for i = 1 : nPatches
        [i1, i2] = bsGetId2D(i, n1, n2);
        patches{i} = data(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1);

        if isempty(options.weightInfo)
            refPatches{i} = refData(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1);
        end
    end
    
    if isempty(options.weightInfo)
        [weights, index] = bsGetKNearestWeights2D(refPatches, n1, n2, options);
        weightInfo.weights = weights;
        weightInfo.index = index;
    else
        weightInfo = options.weightInfo;
        weights = weightInfo.weights;
        index = weightInfo.index;
    end
    
    K = options.nPointsUsed;
    filteredData = zeros(o1, o2);
    filteredNum = zeros(o1, o2);
    
    % average data by calculated patches
    for i = 1 : nPatches
        [i1, i2] = bsGetId2D(i, n1, n2);
        
        iIndex = index(:, i);
        weight = weights(:, i);

        avgPatch = zeros(ws, ws);
        for k = 1 : K
            avgPatch = avgPatch + weight(k) * patches{iIndex(k)};
        end

        filteredData(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1) = ...
            filteredData(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1) + avgPatch;
        filteredNum(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1) = ...
            filteredNum(p1(i1):p1(i1)+ws-1, p2(i2):p2(i2)+ws-1) + 1;
    end
    
    filteredData = filteredData ./ filteredNum;
end


function [i1, i2] = bsGetId2D(i, n1, n2)
    i2 = mod((i-1), n2) + 1;
    i1 = (i - i2) / n2 + 1;
end

function [weights, index] = bsGetKNearestWeights2D(patches, n1, n2, options)
    nPatches = n1 * n2;
    
    N = options.searchOffset;
    K = options.nPointsUsed;
    
    weights = zeros(K, nPatches);
    index = zeros(K, nPatches);
    normCoef = norm(patches{1});
    
    fprintf('Calculating the weight information...\n');
    p = options.p;
    searchStride = options.searchStride;
    
    if options.isParallel
        parfor i = 1 : nPatches
            [weights(:, i), index(:, i)] = bsGetIWeight2D(patches, normCoef, i, n1, n2, N, K, -p, searchStride);
        end
    else
        for i = 1 : nPatches
            [weights(:, i), index(:, i)] = bsGetIWeight2D(patches, normCoef, i, n1, n2, N, K, -p, searchStride);
        end
    end
end

function [weight, index] = bsGetIWeight2D(patches, normCoef, i, n1, n2, N, K, e, searchStride)
    [i1, i2] = bsGetId2D(i, n1, n2);
    
    D = inf((2*N+1)^2, 2);
    iter = 0;

    for c1 = i1-N : searchStride(1) : i1+N
        if c1<=0 || c1>n1
            continue;
        end

        for c2 = i2-N : searchStride(2) : i2+N
            if c2<=0 || c2>n2 || (c2==i2 && c1==i1)
                continue;
            end

            iter = iter + 1;

            j = (c1-1)*n2 + c2;

            difference = patches{j} - patches{i};
            D(iter, 1) = norm(difference/normCoef);
            D(iter, 2) = j;
        end
    end

    [B, ~] = bsMinK(D, K, 1);
    weight = B(:, 1).^e;
    weight = weight / sum(weight);
    index = B(:, 2);
end