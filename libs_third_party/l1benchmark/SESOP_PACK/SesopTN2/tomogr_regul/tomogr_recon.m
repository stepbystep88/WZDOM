% profile on; tomogr_recon; profile viewer;
clear all


par.max_iter=30; % Max iterations of Conjugate Gradients
res = 200;              % Image resolution in pixels (square image)

par.x00=zeros(res);
%par.x00(ceil(res/4),ceil(res/4))=1;
par.x00= phantom('Modified Shepp-Logan', res);




Nangles=2*res;
angles = linspace(-45,135,Nangles);   % Our fast Radon works only for angles from -45 to 135 degrees
%angles = linspace(-45,45,Nangles);;

par.k_small_bins =4; % number of small bins per interval between pixels (used in fastradon.m)
par.flag_fastradon=1;

par.angles=angles;
par.res=res;

y=myradon(par.x00,par); 
%y=fastradon(par.x00,par);
figure; subplot(121);imagesc(par.x00);subplot(122);imagesc(y);colorbar

%tic; y1=radon(par.x00,angles); toc
%figure; subplot(121);imagesc(par.x00);subplot(122);imagesc(y1(3:end-2,:)-y);colorbar
%return


[Nbins,Nang]=size(y);

par.Nbins=Nbins;
par.Nang=Nang;

Rty=myadjradon(y,par);
%Rty=fastadjradon(y,par);



x0=zeros(res);
par.report=@report_recon;     % function for plotting CG progress results

figure;
tic
x = CGmz(@RtR,x0,Rty,par);  % Gonjugate Gradient solution of R'Rx=R'y;  RtR.m computes R'Rx
toc

% figure
% tic
% x = CGmz(@RtRfast,x0,Rty,par);  % Gonjugate Gradient solution of R'Rx=R'y;  RtR.m computes R'Rx
% toc

return




















B = [];
angles = linspace(0,180,Nangles);;
for angle = angles,
    B = [B;fastfp(A,angle,4)];
end

%load dataset1.mat


x0 = zeros(res);
c = zeros(res);

% for i = 1:length(angles)
%     %c = c + fastbpc1(B(i,:),angles(i),4);
%     c = c + fastbp(B(i,:),angles(i),4);
% end

c=myadjradon(B',angles);

% c1=myadjradon(B',angles);
% figure;imagesc(c);colorbar
% figure;imagesc(c1);colorbar
% return

tic
%x = CG(@func1c,x0,-c,angles);
x = CGmz1(@forw_back_proj,x0, c,angles);
%x = CGmz1(@func1c,x0, c,angles);
toc