function [S, org_S, e1, e2, lambda1,lambda2, blur] = bsGetStructureTensor2D(f, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate 2D structures for a 2D data
%
% Copyright (C) 2020. Bin She. All rights reserved.
% Programmed by: Bin She (Email: bin.stepbystep@gmail.com)
% Programming dates: Sep 2020
% -------------------------------------------------------------------------

     %%
    % We define circular convolution
    % \[ (f \star h)_i = \sum_j f_j h_{i-j}. \]
    % Note that here, \(f\) can be multi-channel, in which case each channel is
    % convolved with \(h\). This will be useful to blur tensor fields.
    [n1, n2] = size(f);
    
    blur = bsCreateBlurFcn(f);

    %%
    % We use in the following a centered finite difference approximation of
    % \(\nabla f\), which is a vector field in \(\RR^{n \times n \times 2}\).
    nabla = options.nabla;

    %%
    % We define the tensor product associated to a vector \(u = (u_1,u_2), v=(u_1,u_2)
    % \in \RR^{2}\) as the symetric matrix
    % \[ 
    %       u \otimes v = u v^* = 
    %       \begin{pmatrix} u_1 v_1 & v_1 u_2 \\ u_1 v_2 & u_2 v_2 \end{pmatrix}
    %   \in \RR^{2 \times 2}.
    % \]
    % It is extended to vector fields \( (u(x))_x \in \RR^{N \times 2} \) as
    % \[  (u \otimes v)(x) = u(x) \otimes v(x) \]

    %%
    % A tensor field \(T\) is a collection of symmetric positive definite
    % matrices \(T(x) \in \RR^{2 \times 2}\).

    %%
    % A simple way to build a tensor field is by auto-tensorization of a
    % vector field \(u(x)\), i.e. \(T = u \otimes u\). 

    %%
    % Define a shortcut for \(u \otimes u\)
    % (we make use of symmetry to only store 3 components).

    tensorize = @(u)cat(3, u(:,:,1).^2, u(:,:,2).^2, u(:,:,1).*u(:,:,2));

    %%
    % Rotate a tensor field by \(\pi/2\) (for display only).

    rotate = @(T)cat(3, T(:,:,2), T(:,:,1), -T(:,:,3));

    %% Structure Tensor
    % The structure tensor is a field of symetric positive matrices 
    % that encodes the local orientation and anisotropy of an image.

    %%
    % It was initially introduced for corner detection <#biblio [HarSteph88]> <#biblio [Forstner86]>
    % and oriented texture analysis <#biblio [KassWit85]>.

    %%
    % Given an image \(f\), its structure tensor with scale \( \sigma>0 \) is
    % defined as
    % \[ T_\si = h_\si \star T_0 \qwhereq T_0 = \nabla f \otimes \nabla f. \]
    % For each location \(x\), \(T_\si(x)\) is thus a positive definite matrix.

    T = @(f,sigma)blur( tensorize( nabla( (f) ) ), sigma);

    %%
    % The matrix \(T_\si(x)\) can be understood as the local covariance matrix
    % of the set of gradient vector around \(x\).

    %%
    % Another way to get some insight about this tensor field is to consider a
    % localized version \(f_x\) of the image around point \(x\), defined by 
    % \(f_x(y) = h_\si(x-y)^{1/2} f(y)\), which is close to zero when \(y\) is
    % far away from \(x\).
    % One has the following Taylor expansion of the \(L^2\) norm between two
    % close enough localizations:
    % \[ \norm{f_x - f_{x+\de}}^2 = \de^* T_\si(x) \de + O(\norm{\de}^3). \]

    %%
    % To better understand the behavior of \(T_\si\) as a function of \(\si\),
    % one can computes its Taylor expansion for small \(\si\)
    % \[ T_\si(x) = T_0(x) + \si^2 Hf(x)^2 + O(\si^3), \]
    % where \(Hf(x) \in \RR^{2 \times 2}\) is the Hessian matrix of \(f\) at point \(x\).
    % This shows that when \(\si\) increases, the intial rank-1 tensor \(T_0(x)\) 
    % becomes full rank because it integrates energy from \(Hf(x)^2\).

    %%
    % A convenient way to display a tensor field 
    % such as \(T_\si\) is to draw 
    % an ellispe \(\Ee_x\) at each pixel \(x\) as the (scaled and translated) unit ball of the tensor
    % \[ \Ee_x = \enscond{\de \in \RR^2}{ \de^* T_\si(x) \de \leq 1 }. \]
    % This allows one to visualize the anisotropy and orientation encoded in
    % the tensor field.

    %% 
    % Display \(T_\si\) for \(\si=0.1\) (the tensors are almost rank-1):

%     [T, rotate] = bsCreateTensorFcn(f);

%     options.sub = 16;
    if options.show_mid_results
        figure;
        subplot('position', [0 0 1 1]);
        plot_tensor_field(rotate(T(f, options.sigma)), f, options);
%         title(['\sigma=' num2str(options.sigma)]);
    end

    %% Eigen-decomposition and Anisotropy
    % A symmetric tensor field \(S(x)\) can be decomposed as
    % \[ S(x) = \lambda_1(x) e_1(x) \otimes e_1(x) + \lambda_2(x) e_2(x) \otimes  e_2(x), \]
    % where \((e_1(x),e_2(x))\) are the orthogonal eigenvector fields, \(0 \leq \lambda_2(x) \leq \lambda_1(x)\)
    % are the eigenvalues.

    %%
    % Compute the eigenvalues of \(S \in \RR^{2 \times 2}\) as
    % \[ \la_i =  \frac{1}{2} \pa{ S_{1,1}+S_{2,2} \pm \sqrt{\Delta(S)} }
    %       \qwhereq \Delta(S) = (S_{1,1}-S_{2,2})^2 + 4 S_{1,2}^2, \]
    % where one should use the \(+\) sign for \(i=1\).

    delta = @(S)(S(:,:,1)-S(:,:,2)).^2 + 4*S(:,:,3).^2;
    eigenval = @(S)deal( ...
        (S(:,:,1)+S(:,:,2)+sqrt(delta(S)))/2,  ...
        (S(:,:,1)+S(:,:,2)-sqrt(delta(S)))/2 );

    %%
    % Compute (at each location \(x\)) the leading eigenvector as
    % \[ e_1 = \frac{1}{Z} \begin{pmatrix}
    %       2 S_{1,2} \\
    %       S_{2,2}-S_{1,1} + \sqrt{\Delta(S)}
    %   \end{pmatrix} \]
    % where \(Z\) is a normalization factor ensuring \(\norm{e_1}=1\).

    normalize = @(u)u./repmat(sqrt(sum(u.^2,3))+1e-8, [1 1 2]);
    eig1 = @(S)normalize( cat(3,2*S(:,:,3), S(:,:,2)-S(:,:,1)+sqrt(delta(S)) ) );

    %%
    % Vector \(e_2\) is obtained by applying a \(\pi/2\) rotation to \(e_1\), 
    % which defines the eigenbasis.

    ortho = @(u)deal(u, cat(3,-u(:,:,2), u(:,:,1)));
    eigbasis = @(S)ortho(eig1(S));

    %%
    % Compute the eigendecomposition of \(T_\si\).

    sigma = options.sigma;
    org_S = T(f, sigma);
    
%     if isempty(refData)
%         
%     else
%         S = T(refData, sigma);
%     end
    
    [lambda1,lambda2] = eigenval(org_S);
    [e1,e2] = eigbasis(org_S);
    
    %%
    % Compute the energy and anisotropy
    % \[ E(x) = \sqrt{\lambda_1(x)+\lambda_2(x)} \qandq
    %   A(x) =  \frac{\lambda_1(x)-\la_2(x)}{\la_1(x) + \lambda_2(x)} \in [0,1]. \]
    %%
    % Display it.
    if options.show_mid_results
        figure;
        E = sqrt(lambda1+lambda2);
        A = (lambda1-lambda2)./(lambda1+lambda2);
        imageplot({E A}, {'E', 'A'});
    end
    
    %%
    % Implement the reconstruction formula
    % \[ S = \la_1 (e_1 \otimes e_1) + \la_2 (e_2 \otimes e_2). \]
    recompose = @(lambda1,lambda2,e1,e2)repmat(lambda1,[1 1 3]).*tensorize(e1) + repmat(lambda2,[1 1 3]).*tensorize(e2);
    
     m = options.m;
    Cm = options.Cm;
    
    %%
    % Define \(\phi\).
    phi = @(s,lambda)1-exp( -Cm./(s/lambda).^m );
    lambda = options.lambda;
        %%
    % Display \(\phi(s)\) and \(s\phi(s)\) for \(\la=1\).
%     if options.show_mid_results
%         s = linspace(0,0.05)';
%         figure;
%         plot(s, [phi(s,0.01) s.*phi(s,0.01)], 'LineWidth', 2); 
%         legend('\phi(s)', 's \phi(s)');
%     end
    lambda1 = lambda1 / (norm(lambda1) + 1e-8);
    S = recompose(phi(lambda1,lambda),ones(n1, n2),e1,e2);
end