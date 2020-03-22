% In order to approximately invert Radon transfotm, we compute
% reconstruction functional "w", which has minimal L2-norm given bluring 
% window "wind". Penalty "lam" compromises between ||w|| and the
% influence of pixels located out of the window.

%---------------------- By Michael Zibulevsky, 28.03.2006; version 14.04.2006 ------------------------------%

clear all

Nangles=400;
%angles = linspace(0,180,Nangles);;
%angles = linspace(0,90,Nangles);;
angles = linspace(-45,135,Nangles);   % Our fast Radon works only for angles from -45 to 135 degrees


par.max_iter=100; % Max iterations of Conjugate Gradients
res = 200;               % Resolution of square image in pixels

lam=1e0;  % penalty par. for norm w 
k1=15;         % horizontal window size in image space;
k2=k1;       %vertical window size in image space;;

ky=150;       % window size in projection space
iiy=[1:ky]';

mu=20;   %penalty par. for the constraint 1'A_1*w=1

par.k_small_bins =4; % numper of small bins per interval between pixels (used in fastradon.m)
par.flag_fastradon=1;


i1=ceil(res/4);    %window corner
i2=ceil(res/4);
ii=[i1:i1+k1-1; i2:i2+k2-1]; % window indeces

i1c=i1+ceil(k1/2);  %center of the window
i2c=i2+ceil(k2/2);

xc=zeros(res); xc(i1c,i2c)=1; % One at the central pixel of the window 



wind=zeros(res);
wind(ii(1,:),ii(2,:))=1;


par.wind=wind;
par.lam=lam;
par.mu=mu;
par.angles=angles;
par.res=res;

b=myradon(wind,par);
[Nbins,Nang]=size(b);

par.b=b;
par.Nbins=Nbins;
par.Nang=Nang;

yc=myradon(xc,par);


[ttt,ind_ymax]=max(yc);
y1=zeros(Nbins,Nang);
iiy=zeros(ky,Nang);
for i=1:Nang, 
	iiw(:,i) = ind_ymax(i)-ceil(ky/2)+[1:ky]'; 
	y1( iiw(:,i)  ,i) =1;
end

par.iiw=iiw; % indeces of nonzero elements in w

for i=1:par.Nang, b_local(:,i)= b(par.iiw(:,i),i);end

%y1(iiy)=100;


figure; subplot(221);imagesc(xc);colorbar;
subplot(222);imagesc(yc);colorbar;
subplot(223);imagesc(y1);colorbar;

%return

%w0=zeros(size(b));
%w = CGmz1(@fw,w0,mu*b,par);

figure;
w0_local=zeros(size(iiw));
par.report=@report_local;  % function for plot  CG progress results
w_local = CGmz(@fw_local,w0_local,mu*b_local,par); % Gonjugate Gradient optimization
w=zeros(par.Nbins,par.Nang);
for i=1:Nang, w(par.iiw(:,i),i)=w_local(:,i);end

Atw=myadjradon(w,par);

w=w/(wind(:)'*Atw(:));  % normalization 
Atw=Atw/(wind(:)'*Atw(:));

%figure;imagesc((Atw)); colorbar;
%title(sprintf('NormW = %g  normA2tw=%g lam=%g',norm(w(:)), norm(Atw(:).*(1-wind(:))),  lam));
%figure;imagesc(indicator_wind);colorbar
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















