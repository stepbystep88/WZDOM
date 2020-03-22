function nnreport(vcwb,report,par,Xtrain,ytrain,Xtest,ytest),
% Plot Neural Net error during learning SESOP_TN iterations
%
% Call: nnreport(vcwb,report,par,Xtrain,ytrain,Xtest,ytest)
%
%global GlobalErrTrain GlobalErrTest GlobalNNiter

% Michael Zibulevsky, 05.08.2008        
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 

global GlobalErrTrain GlobalErrTest GlobalNNiter

GlobalNNiter=[GlobalNNiter report.Niter];

ynntrain=nnet(vcwb,par,Xtrain);
GlobalErrTrain = [GlobalErrTrain sumsqr(ytrain-ynntrain)];

ynntest=nnet(vcwb,par,Xtest);
GlobalErrTest = [GlobalErrTest sumsqr(ytest-ynntest)];

figure(par.nnreport_figure_handle);
subplot(311);
plot(GlobalNNiter,[GlobalErrTrain' GlobalErrTest']);
%semilogy(GlobalNNiter,[GlobalErrTrain' GlobalErrTest']);
legend('Train Error','Test Error'); grid
subplot(312);plot([ytrain(:) ynntrain(:)]);grid;legend('ytrain', 'ynntrain')
subplot(313);plot([ytest(:) ynntest(:)]);grid;legend('ytest', 'ynntest')
drawnow
