function bsSaveFigure(path, name)

    if isempty(name) || strcmp(name, 'tmp') || name(1) == '_'
        return;
    end

    warning('off');
    mkdir(sprintf('%s/figure', path));
    mkdir(sprintf('%s/jpg', path));
    mkdir(sprintf('%s/emf', path));
    mkdir(sprintf('%s/eps', path));
    mkdir(sprintf('%s/pdf', path));
    mkdir(sprintf('%s/tif', path));
    mkdir(sprintf('%s/png', path));
    mkdir(sprintf('%s/mat', path));
    
    
    set(gcf, 'paperPositionMode', 'auto');
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    
    if(fig_pos(1) < 1)
        set(gcf,'Units','normalized');
    end
    
    fig.PaperSize = [fig_pos(3) fig_pos(4)];

    fileName = sprintf('%s/mat/%s.mat', path, name);
    save(fileName);
    warning('on');
    
    fileName = sprintf('%s/figure/%s.fig', path, name);
    savefig(gcf, fileName);
    
%     fileName = sprintf('%s/jpg/%s.jpg', path, name);
%     print('-djpeg', '-r600', fileName);
    
%     fileName = sprintf('%s/tif/%s.tif', path, name);
%     print('-dtiff', '-r600', fileName);
%     
%     fileName = sprintf('%s/emf/%s.emf', path, name);
%     print('-dmeta', fileName);
%     
    fileName = sprintf('%s/png/%s.png', path, name);
    print('-dpng', '-r600', fileName);
    
    fileName = sprintf('%s/jpg/%s.jpg', path, name);
    print('-djpeg', '-r600', fileName);
%     
    fileName = sprintf('%s/eps/%s.eps', path, name);
    print('-depsc', fileName);
% %     
    fileName = sprintf('%s/pdf/%s.pdf', path, name);
    print('-dpdf', '-bestfit', fileName);
    saveas(gcf, fileName);
    
end