function [S, org_S, blur] = bsGetStructureTensor3D(f, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate 3D structures for a volume data
%
% Copyright (C) 2020. Bin She. All rights reserved.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Sep 2020
% -------------------------------------------------------------------------

    blur = bsCreateBlurFcn(f);

    %%
    % We use in the following a centered finite difference approximation of
    % \(\nabla f\), which is a vector field in \(\RR^{n \times n \times 2}\).
    nabla = options.nabla;

    tensorize = @(u)cat(4, u(:,:,:,1).^2, u(:,:,:,2).^2, u(:,:,:,3).^2, u(:,:,:,1).*u(:,:,:,2), u(:,:,:,2).*u(:,:,:,3), u(:,:,:,1).*u(:,:,:,3));

    %%
    % Rotate a tensor field by \(\pi/2\) (for display only).
%     A =
%     [ aa, ab, ac]
%     [ ab, bb, bc]
%     [ ac, bc, cc]
%     inv(A) = 
%     [  bc^2 - bb*cc, ab*cc - ac*bc, ac*bb - ab*bc]
%     [ ab*cc - ac*bc,  ac^2 - aa*cc, aa*bc - ab*ac]
%     [ ac*bb - ab*bc, aa*bc - ab*ac,  ab^2 - aa*bb]
    
%     aa = T(:, :, :, 1);
%     bb = T(:, :, :, 2);
%     cc = T(:, :, :, 3);
%     ab = T(:, :, :, 4);
%     bc = T(:, :, :, 5);
%     ac = T(:, :, :, 6);
    
    
    %% Structure Tensor
    % The structure tensor is a field of symetric positive matrices 
    % that encodes the local orientation and anisotropy of an image.
    T = @(f,sigma)blur( tensorize( nabla(f) ), sigma);
    

    %% Eigen-decomposition and Anisotropy
    % A symmetric tensor field \(S(x)\) can be decomposed as
    % \[ S(x) = \lambda_1(x) e_1(x) \otimes e_1(x) + \lambda_2(x) e_2(x) \otimes  e_2(x), \]
    % where \((e_1(x),e_2(x))\) are the orthogonal eigenvector fields, \(0 \leq \lambda_2(x) \leq \lambda_1(x)\)
    % are the eigenvalues.
    if options.show_mid_results
        figure;
%         tbl = bsGetColormap('seismic');
        bsShow3DVolume(f, 1, [prctile(f(:), 10), prctile(f(:), 90)], 1, 1, 1, [], 'colormap', options.colormap);
    end
    
    absbf = abs(f);
    org_S = T(f ./ prctile(absbf(:), 90), options.sigma);
        
    % 重新构建S，加强横向连续性的约束
    if options.show_mid_results
        [S, Ds] = reconstruted_tensor(org_S, options);
        bsPlotTensorField3D(f, Ds, options);
    else
        [S] = reconstruted_tensor(org_S, options);
    end
    
end

function [S, Ds] = reconstruted_tensor(T, options)
        
        
    [n1, n2, n3, ~] = size(T);

    m = options.m;
    Cm = options.Cm;
    lambda = options.lambda;

    %%
    % Define \(\phi\).
    phi = @(s,lambda)1-exp( -Cm./(s/lambda).^m );
 
    Ts = reshape(T, [], 6);
    n = n1 * n2 * n3;
    S = zeros(n, 6);
    
    if nargout > 1
        Ds = cell(1, n);
        isOutDs = true;
    else
        Ds = [];
        isOutDs = false;
    end
    
    c1 = 0.001;
    c2 = 1;
    
    pbm = bsInitParforProgress([], ...
            n, ...
            'Calculating the structure tensor...', ...
            './', ...
            0);
    
    parfor i = 1 : n
        v = Ts(i, :);

        [V, D] = eig([v(1) v(4) v(6); v(4) v(2) v(5); v(6) v(5) v(3)]);
        sd = D(1, 1) + D(2, 2) + D(3, 3);
        
        d1 = D(1,1) / sd;
        d2 = D(2,2) / sd;
        d3 = D(3,3) / sd;
        
        l1 = c1;
        l3 = c1 + (1-c1) * exp(-c2/(d3 - d1)^2);
        l2 = c1 + (1-c1) * exp(-c2/(d3 - d2)^2);

        Si = l3*(V(:, 1)*V(:, 1)') + l2*(V(:, 2)*V(:, 2)') + l1*(V(:, 3)*V(:, 3)');

        S(i, :) = [Si(1, 1), Si(2, 2), Si(3, 3), Si(1, 2), Si(2, 3), Si(1, 3)];
%         S(i, :) = v;
        
        if isOutDs
            Ds{i}= V;
        end
        
        bsIncParforProgress(pbm, i, 1000000);
    end
    
    bsDeleteParforProgress(pbm);
    
    S = reshape(S, n1, n2, n3, []);
    if isOutDs
        Ds = reshape(Ds, n1, n2, n3);
    end
end