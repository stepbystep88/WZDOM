%% Tensor-driven Diffusion Flows
% This numerical tour explores the structure tensor to represent the
% geometry of images and textures. It applies it to perform anisotropic
% image diffusion.
% A good reference for diffusion flows in image processing is <#biblio [Weickert98]>.

% perform_toolbox_installation('signal', 'general');

%% Helpers Functions
% We define here a few features (convolution, gradient, etc.) that will be
% used in the sequel.

%%
% Size of the image of \(N=n \times n\) pixels.

n = 256;
close all;
%%
% Load an image \(f\).

name = 'hibiscus';
f0 = load_image(name,n);
f0 = rescale( sum(f0,3) );
% f = profileData(:, :, 2);
% f = invResults{3}.data;

% f0 = reshape(profileData(:,:,2), size(profileData, 1), size(profileData, 2));
f = f0;

[n1,n2] = size(f);

% [X,Y] = meshgrid(1:n2, 1:n1);
% [Xq,Yq] = meshgrid(1:1:n2, 1:0.25:n1);
% 
% Vq = interp2(X,Y,f,Xq,Yq);
% % 
% f = Vq;

[n1,n2] = size(f);
% 
% f = f(1:700, 1:700);
%%
% Display it.

figure;
imagesc(f);

%%
% We define circular convolution
% \[ (f \star h)_i = \sum_j f_j h_{i-j}. \]
% Note that here, \(f\) can be multi-channel, in which case each channel is
% convolved with \(h\). This will be useful to blur tensor fields.

cconv = @(f,h)real(ifft2(fft2(f).*repmat(fft2(h),[1 1 size(f,3)])));

%%
% Define a Gaussian blurring kernel of width \(\si\):
% \[ h_\si(x) = \frac{1}{Z} e^{ -\frac{x_1^2+x_2^2}{2\si^2} }\]
% where \(Z\) ensures that \(\hat h_\si(0)=1\).

t = [0:n/2 -n/2+1:-1];
[X2,X1] = meshgrid(t,t);
normalize = @(h)h/sum(h(:));
h = @(sigma)normalize( exp( -(X1.^2+X2.^2)/(2*sigma^2) ) );

%%
% Define the convolution with \(h_\si\).

blur = @(f,sigma)cconv(f,h(sigma));

%%
% We use in the following a centered finite difference approximation of
% \(\nabla f\), which is a vector field in \(\RR^{n \times n \times 2}\).

options.order = 2;
nabla = @(f)grad(f,options);

%%
% We define the tensor product associated to a vector \(u = (u_1,u_2), v=(u_1,u_2)
% \in \RR^{2}\) as the symetric matrix
% \[ 
%       u \otimes v = u v^* = 
%       \begin{pmatrix} u_1 v_1 & v_1 u_2 \\爑_1 v_2 & u_2 v_2 \end{pmatrix}
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

T = @(f,sigma)blur( tensorize( nabla(f) ), sigma);

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

[T, rotate] = bsCreateTensorFcn(f);

options.sub = 8;
figure; sigma = 0.1;
plot_tensor_field(rotate(T(f,sigma)), f, options);
title(['\sigma=' num2str(sigma)]);
set(gcf, 'position', [680   223   800   755]);

%%
% For \(\si=4\):

figure; sigma = 2;
plot_tensor_field(rotate(T(f,sigma)), f, options);
title(['\sigma=' num2str(sigma)]);
set(gcf, 'position', [680   223   800   755]);

%%
% For \(\si=10\):

figure; sigma = 10;
plot_tensor_field(rotate(T(f,sigma)), f, options);
title(['\sigma=' num2str(sigma)]);
set(gcf, 'position', [680   223   800   755]);

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

normalize = @(u)u./repmat(sqrt(sum(u.^2,3)), [1 1 2]);
eig1 = @(S)normalize( cat(3,2*S(:,:,3), S(:,:,2)-S(:,:,1)+sqrt(delta(S)) ) );

%%
% Vector \(e_2\) is obtained by applying a \(\pi/2\) rotation to \(e_1\), 
% which defines the eigenbasis.

ortho = @(u)deal(u, cat(3,-u(:,:,2), u(:,:,1)));
eigbasis = @(S)ortho(eig1(S));

%%
% Compute the eigendecomposition of \(T_\si\).

sigma = 2;
S = T(f,sigma);
[lambda1,lambda2] = eigenval(S);
[e1,e2] = eigbasis(S);

sub = 4;
figure;

for k = 1 :2
    subplot(1, 2, k);
    imagesc(f0); hold on;
    
    for i = sub : sub : n1- sub
        for j = sub : sub :n2 - sub

            if k == 1
                v = e1(i, j, :);
            else
                v = e2(i, j, :);
            end
            
            vv = [v(1) v(2)] * 4;

            px = [i, j] - 0.5 * vv;
            py = [i, j] + 0.5 * vv;

            plot([px(1) py(1)], [px(2), py(2)], 'r-', 'linewidth', 2);

        end
    end
end

%%
% Implement the reconstruction formula
% \[ S = \la_1 (e_1 \otimes e_1) + \la_2 (e_2 \otimes e_2). \]

recompose = @(lambda1,lambda2,e1,e2)repmat(lambda1,[1 1 3]).*tensorize(e1) + repmat(lambda2,[1 1 3]).*tensorize(e2);

%%
% Check that the recomposition is exact.

mynorm = @(x)norm(x(:));
S1 = recompose(lambda1,lambda2,e1,e2);
fprintf('Should be 0: %.3f\n', mynorm(S-S1));

%%
% The eigenvalues of \(T_\si\) can be used to detect interest point in the
% image:

%%
% * A flat region is composed of pixels \(x\) with \(\la_1(x) \approx \la_2(x) \approx 0\).
% * A straight edge is composed of pixels \(x\) with \(0 \approx \la_2(x) \ll \la_1(x)\).
% * A corner is composed of pixels \(x\) with \(0 \ll \la_2(x) \approx \la_1(x)\).

%%
% This idea is at the heart of the Forstner/Harris corner detector <#biblio [HarSteph88]> <#biblio [Forstner86]>.

%%
% Compute the energy and anisotropy
% \[ E(x) = \sqrt{\lambda_1(x)+\lambda_2(x)} \qandq
%   A(x) =  \frac{\lambda_1(x)-\la_2(x)}{\la_1(x) + \lambda_2(x)} \in [0,1]. \]

E = sqrt(lambda1+lambda2);
A = (lambda1-lambda2)./(lambda1+lambda2);

%%
% Display it.

figure;
imageplot({E A}, {'E', 'A'});
set(gcf, 'position', [680   223   800   755]);
%% Tensor Driven Anisotropic Diffusion
% A tensor field \(S\) can be used as anisotropic metric to drive a diffusion PDE flow.
% The good reference for such a flow is <#biblio [Weickert98]>.

%%
% This defines an anisotropic diffusion flow \(t \mapsto f_t\)
% \[ \pd{f_t}{t}(x) = \text{div}\pa{ S(x) \nabla f_t(x) } \]
% where \(f_0\) is a given data at time \(t=0\).

%%
% Note that this is actually a linear PDE, since \(S\) does not evolve in
% time. But in practice, \(S\) is usually computed from \(f_0\), so that
% the mapping \(f_0 \mapsto f_t\) is actually non-linear.

%%
% This PDE is discretized in time using a explicit time stepping
% \[ f^{(\ell+1)}(x) = f^{(\ell)}(x) + \tau \text{div}\pa{ S(x) \nabla f^{(\ell)}(x) } \]

%%
% The time step \(\tau\) should be small enough for the diffusion to be stable.

%%
% To produce edge-enhancing diffusion, we define \(S\) from the structure
% tensor field \(T_\si\) by re-normalizing the eigenvalues. 
% \[ S(x) = \phi(\lambda_1(x)) e_1(x)e_1(x)^* + e_2(x)e_2(x)^*, \]
% where \(\phi : \RR^+ \rightarrow \RR^+\) is defined, following
% <#biblio [Weickert98]>, as
% \[ \phi(s) = 1 - \text{exp}\pa{
%      -\frac{C_m}{(s/\la)^m} }. \]
% Here \(m\) is a given exponent, and 
% the constant \(C_m\) ensures that \(s \phi(s)\) is increasing for 
% \(s < \la\) and decreasing for \(s > \la\), which produces the
% edge-enhancing effect.

%%
% Set the values of \(m\) and \(C_m\).

m = 4;
Cm = 3.31488;

%%
% Define \(\phi\).

phi = @(s,lambda)1-exp( -Cm./(s/lambda).^m );

%%
% Display \(\phi(s)\) and \(s\phi(s)\) for \(\la=1\).

s = linspace(0,0.05)';
figure;
plot(s, [phi(s,0.01) s.*phi(s,0.01)], 'LineWidth', 2); 
legend('\phi(s)', 's \phi(s)');

%%
% Select \(\lambda\).

lambda = 1e-4;

%%
% Select \(\si\).

sigma = 2;

%%
% Compute the eigen-decomposition of \(T_\si\).

S = T(f,sigma);
[lambda1,lambda2] = eigenval(S);
[e1,e2] = eigbasis(S);


%%
% Compute \(S\).
S = recompose(phi(lambda1,lambda),ones(n1, n2),e1,e2);
% S = recompose(ones(n1, n2)*0,ones(n1, n2)*1,e1,e2);

%%
% Note that this remapping of the eigenvalues of \(T\) to the eigenvalues
% of \(S\) exchanges the roles of the eigenaxes. This causes the diffusion
% to be stronger along the edges, and to be small perpenticular to it.

%%
% This flow can thus be seen as an anisotropic version of the famous
% Perona-Malick flow <#biblio [PerMal90]>. Note that the Perona-Malick flow
% is often refered to as an _anisotropic diffusion_, but it is actually
% incorrect, because the diffusion tensor associated to is is actually
% isotropic, since it corresponds to using a time-dependent tensor field
% \[ S(x) = \phi(\norm{\nabla f_t(x)}) \text{Id}_2 . \] 

%%
% Shortcut for the multiplication \(S u\) of tensor
% \(S\) by vector field \(u\).

Mult = @(S,u)cat(3, S(:,:,1).*u(:,:,1) + S(:,:,3).*u(:,:,2), ...
                          S(:,:,3).*u(:,:,1) + S(:,:,2).*u(:,:,2) );
                        
%%
% Step size \(\tau\).
                        
tau = .1;

%%
% First initialize the image to diffuse at time \(t=0\).

% f1 = f;
% f1 = invResults{2}.data(:, 1:n);

%%
% Perform one step of the diffusion.

% f1 = f + tau * div( Mult(S, nabla(f) ) );

%EXO
%% Perform the full diffusion up to a large enough time.
f = bsAddNoise(f0, 1, 20, [], [], [], []);

f1 = filter2(fspecial('average', [5,5]), f, 'same');
f2 = medfilt2(f, [5 5]);
% f3 = filter2(fspecial('gaussian', [5,5]), f, 'same');

f3 = bsNLMByRef(f, [], 'windowSize', [5, 5], 'searchOffset', 5, 'nPointsUsed', 5);

GST_options = bsCreateGSTParam();
GST_options.sigma = 2;
GST_options.iterNum = 50;
GST_options.ttv = 0.001;
GST_options.tau = 0.;
GST_options.lambda = 1e-4;
GST_options.show_mid_results = false;
GST_options.sub = 8;
GST_options.filter_fcn = @(f)medfilt2(f, [2, 2]);

f4 = bsSmoothByGST2D(f, [],  GST_options);

GST_options.ttv = 0.0001;
GST_options.tau = 0.05;

f5 = bsSmoothByGST2D(f, [],  GST_options);
f6 = bsSmoothByGST2D(f, f0,  GST_options);

tbl = bsGetColormap('original');
figure;
set(gcf, 'position', [269         117        1602         789]);

range = [prctile(f0(:), 5), prctile(f0(:), 95)];
subplot(2, 4, 1); imagesc(f0); set(gca, 'clim', range); title('原始图像');  colormap(tbl);
subplot(2, 4, 2); imagesc(f); set(gca, 'clim', range); title('含噪图像');   colormap(tbl);
subplot(2, 4, 3); imagesc(f1); set(gca, 'clim', range); title('均值滤波');  colormap(tbl);
subplot(2, 4, 4); imagesc(f2); set(gca, 'clim', range); title('中值滤波');  colormap(tbl);
subplot(2, 4, 5); imagesc(f3); set(gca, 'clim', range); title('非局部均值滤波');          colormap(tbl);
subplot(2, 4, 6); imagesc(f4); set(gca, 'clim', range); title('全变差光滑化');            colormap(tbl);
subplot(2, 4, 7); imagesc(f5); set(gca, 'clim', range); title('各向异性结构张量滤波 (噪声张量)');                colormap(tbl);
subplot(2, 4, 8); imagesc(f6); set(gca, 'clim', range); title('各向异性结构张量滤波 (干净张量)');    colormap(tbl);

figure;
set(gcf, 'position', [269         117        1602         789]);

res = cell(2, 4);
iterNum = [10 25 50 100];
for i = 1 : 4
    GST_options.iterNum = iterNum(i);
    GST_options.ttv = 0.001;
    GST_options.tau = 0.;   
    res{1, i} = bsSmoothByGST2D(f, f0,  GST_options);
    
    GST_options.ttv = 0.0001;
    GST_options.tau = 0.05;
    res{2, i} = bsSmoothByGST2D(f, f0,  GST_options);
    
    subplot(2, 4, i); imagesc(res{1, i}); set(gca, 'clim', range); title(sprintf('全变差 (iter:%d)', iterNum(i)));  colormap(tbl);
    subplot(2, 4, 4 + i); imagesc(res{2, i}); set(gca, 'clim', range); title(sprintf('各向异性结构张量滤波 (iter:%d)', iterNum(i)));  colormap(tbl);
end
