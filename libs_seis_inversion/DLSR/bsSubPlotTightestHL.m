function bsSubPlotTightestHL( M, N, index, sx, sy)
%% 本程序用于非常紧凑的显示字典学习的结果
%
% 输入
% M             行数
% N             列数
% index         第几个
%
% 输出           无

    dw = 0.93/N;  dh = 0.9/M;
    w = dw; h = dh ;
    
    y = floor( (index-1) / N);
    x = mod(index-1, N);
    subplot('Position', [x*dw+sx (M-y-1)*dh+sy w h]);
end
