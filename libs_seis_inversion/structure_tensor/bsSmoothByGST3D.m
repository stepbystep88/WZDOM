function f1 = bsSmoothByGST3D(f, refData, options, S, blur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth a 3D volume data by using gradient structure tensor (GST)
% 
% Example: 
% options = bsCreateGSTParam(3);
% smoothed_data = bsSmoothByGST3D(original_data, [], options);
%
% Copyright (C) 2020. Bin She. All rights reserved.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Sep 2020
% -------------------------------------------------------------------------

    if nargin <= 3 || isempty(S)
        if isempty(refData)
    %         bf = blur(f, 0.1);
    %         f = bf;
            [S, ~, blur] = bsGetStructureTensor3D(f, options);
        else
    %         S = T(refData / prctile(refData(:), 90), sigma);
            [S, ~, blur] = bsGetStructureTensor3D(refData, options);
        end
    else
        if nargin <= 4 || isempty(blur)
            blur = bsCreateBlurFcn(f);
        end
    end
    
    
    %%
    % Shortcut for the multiplication \(S u\) of tensor
    % \(S\) by vector field \(u\).
%     syms a1 a2 a3 a4 a5 a6
%     A = [a1 a4 a6; a4 a2 a5; a6 a5 a3];
%     syms b1 b2 b3
%     b = [b1; b2; b3];
%     A * b
    Mult = @(S,u)cat(4, S(:,:,:,1).*u(:,:,:,1) + S(:,:,:,4).*u(:,:,:, 2) + S(:,:,:,6).*u(:,:,:, 3), ...
                        S(:,:,:,2).*u(:,:,:,2) + S(:,:,:,4).*u(:,:,:, 1) + S(:,:,:,5).*u(:,:,:, 3), ...
                        S(:,:,:,3).*u(:,:,:,3) + S(:,:,:,5).*u(:,:,:, 2) + S(:,:,:,6).*u(:,:,:, 1) ...
                        );

    %%
    % Step size \(\tau\).
    tau = options.tau;
    
    %%
    % First initialize the image to diffuse at time \(t=0\).

    %% Perform the full diffusion up to a large enough time.
    maxVal = max(f(:));
    minVal = min(f(:));
    
    if ~isempty(options.filter_fcn)
        if isa(options.filter_fcn, 'function_handle')
            f1 = options.filter_fcn(f);
        elseif  options.filter_fcn
            f1 = blur(f, 0.01);
        end
        
        f1(f1<minVal) = minVal;
        f1(f1>maxVal) = maxVal;
    else
        f1 = f;
    end
    
    nabla = options.nabla;
    
    for i=1:options.iterNum
        
        f1 = f1 + tau * div( Mult(S, nabla(f1) ) );
    end
    
    
%     f1 = filter2(options.h, f1, 'same');
    if ~isempty(options.filter_fcn)
        if isa(options.filter_fcn, 'function_handle')
            f1 = options.filter_fcn(f1);
        elseif  options.filter_fcn
            f1 = blur(f1, 0.01);
        end
        
        f1(f1<minVal) = minVal;
        f1(f1>maxVal) = maxVal;
        
    end
end


