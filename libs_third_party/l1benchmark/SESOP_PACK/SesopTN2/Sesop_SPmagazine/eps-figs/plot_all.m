%plot_all.m
%
%called by show_all_results_FINAL1.m

i_run=numel(reports);
%ind=1:i_run;              % which curves to show


ColorPlots=1;

if ColorPlots
	%colorvec = [{'k'},{'b'},{'r'},{'g'},{'y'},{'c'},{'m'},{'.'},{'--'},{'r.'}];
	%colorvec = [{'k'},{'k:'},{'b'},{'b:'},{'r'},{'r:'},{'g'},{'g:'},{'y'},{'c'},{'m'},{'.'},{'--'},{'r.'}];
	%colorvec = [{'k'},{'k--'},{'b'},{'b--'},{'r'},{'r--'},{'g'},{'g--'},{'c'},{'c--'},{'m'},{'m--'},{'.'},{'--'},{'r.'}];
   colorvec = [{'k'},{'k--'},{'b'},{'b--'},{'r'},{'r--'},{'g'},{'g--'},{'c'},{'c--'},{'m'},{'m--'},{'k-.'},{'b-.'},{'r-.'},{'g-.'},{'c-.'},{'m-.'}];

	
else
	%colorvec = [{'k:'},{'k-.'},{'k--'},{'k.'},{'k'},{'k+'},{'k^'},{'ks'},{'k'}];
	colorvec = [{'k'},{'k:'},{'k-.'},{'k--'},{'k.'},{'k+'},{'k^'},{'ko'},{'ks'} {'k*'} {'kv-'}];
	%colorvec = [{'k'},{'k:'},{'k-.'},{'k.'},{'k+'},{'k^'},{'ko'},{'ks'} {'k*'} {'kd'},{'k--'}];
	%colorvec = [{'k'},{'k:'},{'k.'},{'k+'},{'k-.'},{'k^'},{'ko'},{'ks'},{'k--'}];
	%colorvec = [{'k'},{'k.'},{'k:'},{'k+'},{'k-.'},{'k^'},{'ko'},{'ks'},{'k--'}];

	%colorvec = [{'k'},{'k:'},{'k-.'},{'ko'},{'k+'},{'k^'},{'ko'},{'ks'},{'k--'}];
	%colorvec = [{'k'},{'k.'},{'ko'},{'k*'},{'k^'},{'ks'},{'kd'},{'k>'},{'kv'}];
end


% if isfield(reports(1),'FlagNewReport'), FlagNewReport=1; %new report format
% else                                    FlagNewReport=0;
% end
% 

axislimits_ft  =axislimits_f;   axislimits_ft(2)  =tmax;
axislimits_snrt=axislimits_snr; axislimits_snrt(2)=tmax;

fbest=1e100; for i=1:i_run, fbest=min(fbest,min(reports(i).func_values));end
f0=reports(1).func_values(1)-fbest;

f_figure_handle=figure('Position',[00 40 400 600],'name',['Objective function' report_name] );
mylegends=[];
k=0;
for i=ind,
	if strcmp(reports(i).method,'L1-LS-IntPoint'),
		h=semilogy(reports(i).nniter_fg, (reports(i).func_values-fbest)/f0,colorvec{i});
	else
		h=semilogy(reports(i).nniter, (reports(i).func_values(reports(i).nniter+1)-fbest)/f0,colorvec{i});
   end
   set(h,'LineWidth',2); 
   k=k+1; mylegends{k}=reports(i).method;
   h=xlabel('Iteration');set(h,'FontSize',18); 
   h=ylabel('f - fbest');set(h,'FontSize',18); 
   hold on;
end
set(gca,'FontSize',18); 
legend(mylegends);
axis(axislimits_f);
print( '-depsc2', [name '_ObjFunc']);


% Timing f

f_figure_handle=figure('Position',[400 40 400 600],'name',['Objective function with time' report_name] );
mylegends=[];
k=0;
for i=ind,
   ttt=reports(i).times;ttt=ttt-ttt(1);
	if strcmp(reports(i).method,'L1-LS-IntPoint'),
		h=semilogy(ttt, (reports(i).func_values-fbest)/f0,colorvec{i});
	else
		h=semilogy(ttt(reports(i).nniter+1), (reports(i).func_values(reports(i).nniter+1)-fbest)/f0,colorvec{i});
   end
   set(h,'LineWidth',2); 
   k=k+1; mylegends{k}=reports(i).method;
   h=xlabel('CPU Time, Sec');set(h,'FontSize',18); 
   h=ylabel('f - fbest');set(h,'FontSize',18); 
   hold on;
end
set(gca,'FontSize',18); 
legend(mylegends);
axis(axislimits_ft);
print( '-depsc2', [name '_ObjFuncTime']);





% figure('Position',[400 40 500 800],'name',['Objective with time' report_name] );
% for i=ind,
%    ttt=reports(i).times;ttt=ttt-ttt(1);
%    semilogy(ttt, reports(i).func_values-fbest,colorvec{i});
%    xlabel('CPU time, Sec');
%    ylabel('f - fbest');
%    hold on;
% end
% legend(mylegends);
% % 

figure('Position',[800 40 400 600],'name',['SNR ' report_name] );
for i=ind,
   h=plot(reports(i).nniter, reports(i).Xsnr,colorvec{i});
   set(h,'LineWidth',2); hold on;
   h=xlabel('Iteration');set(h,'FontSize',18); 
   h=ylabel('SNR, dB');set(h,'FontSize',18); 
end
set(gca,'FontSize',18); 
legend(mylegends,'Location','SouthEast'); %grid
axis(axislimits_snr);
print( '-depsc2', [name '_snr']);


figure('Position',[1200 40 400 600],'name',['SNR withtime ' report_name] );
for i=ind,
   %if isfield(reports(i),'SNRtime') , ttt=reports(i).SNRtime;ttt=ttt-ttt(1);
   if strcmp(reports(i).method,'L1-LS-IntPoint'),
      ttt=reports(i).SNRtime; ttt=ttt-ttt(1);
   else
      ttt=reports(i).times(reports(i).nniter+1); ttt=ttt-ttt(1);
   end
   h=plot(ttt, reports(i).Xsnr,colorvec{i});
   set(h,'LineWidth',2); hold on;
   h=xlabel('CPU time, Sec');set(h,'FontSize',18); 
   h=ylabel('SNR, dB');set(h,'FontSize',18); 
end
set(gca,'FontSize',18); 
legend(mylegends,'Location','SouthEast'); %grid
axis(axislimits_snrt);
print( '-depsc2', [name '_snrTime']);

 
% figure('Position',[1200 40 500 800],'name',['SNR with time ' report_name] );
% for i=1:i_run,
%    if isfield(reports(i),'SNRtime') , ttt=reports(i).SNRtime;ttt=ttt-ttt(1);
%    else ttt=reports(i).times(reports(i).nniter+1);ttt=ttt-ttt(1);
%    end
%    plot(ttt, reports(i).Xsnr,colorvec{i});
%    xlabel('CPU time, Sec');ylabel('SNR, Db')
%    hold on;
% end
% legend(mylegends,'Location','SouthEast');grid

% for i=1:i_run, reports(i).options=[];reports(i).options1=[];end
% save reports_ComprSens_128aaaa.mat reports

