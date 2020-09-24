function [f1, S] = bsSmoothByGST2D(f, refData, options, S, blur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth a 2D volume data by using gradient structure tensor (GST)
% 
% Example: 
% options = bsCreateGSTParam(2);
% smoothed_data = bsSmoothByGST2D(original_data, [], options);
%
% Copyright (C) 2020. Bin She. All rights reserved.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Sep 2020
% -------------------------------------------------------------------------

%     [n1, n2] = size(f);
    
    if nargin <= 2 || isempty(options)
        options = bsCreateGSTParam(2);
    end
    
    if nargin <= 3 || isempty(S)
        if isempty(refData)
            [S, ~, ~, ~, ~,~, blur] = bsGetStructureTensor2D(f, options);
        else
            [S, ~, ~, ~, ~,~, blur] = bsGetStructureTensor2D(refData, options);
        end
    else
        if nargin <= 4 || isempty(blur)
            blur = bsCreateBlurFcn(f);
        end
    end
    
    

    %%
    % Shortcut for the multiplication \(S u\) of tensor
    % \(S\) by vector field \(u\).

    Mult = @(S,u)cat(3, S(:,:,1).*u(:,:,1) + S(:,:,3).*u(:,:,2), ...
                              S(:,:,3).*u(:,:,1) + S(:,:,2).*u(:,:,2) );

    %%
    % Step size \(\tau\).

    tau = options.tau;
    ttv = options.ttv;

    nabla = options.nabla;
    
    %%
    % First initialize the image to diffuse at time \(t=0\).

    %% Perform the full diffusion up to a large enough time.
    kdisp = round(linspace(0, options.iterNum,5)); kdisp(1) = [];
    k =1;

%     if options.show_mid_results
%         figure;
%     end
    
%     if ~isempty(options.filter_fcn)
%         if isa(options.filter_fcn, 'function_handle')
%             f1 = options.filter_fcn(f);
%         elseif  options.filter_fcn
%             f1 = blur(f, 0.01);
%         end
%         
%     end
    
    f1 = f;
    
    if isempty(tau)
        tau = 0.1 * div( Mult(S, nabla(f1) ) ) / f1;
    end
    
    maxVal = max(f1(:));
    minVal = min(f1(:));
    
    for i=1:options.iterNum
        
        dxy = nabla(f1);
        f1 = f1 + tau * div( Mult(S,  dxy) ) + ttv * div( sign(dxy) );
        if i==kdisp(k) && options.show_mid_results
            subplot(2,2,k);
            imagesc(f1);
            k = k+1;
        end
        
        f1(f1<minVal) = minVal;
        f1(f1>maxVal) = maxVal;
    end
        
%     f1 = filter2(options.h, f1, 'same');
%     if ~isempty(options.filter_fcn)
%         f1 = options.filter_fcn(f1);
%     end

    if ~isempty(options.filter_fcn)
        if isa(options.filter_fcn, 'function_handle')
            f1 = options.filter_fcn(f1);
        elseif  options.filter_fcn
            f1 = blur(f1, 0.01);
        end
        
    end
end
