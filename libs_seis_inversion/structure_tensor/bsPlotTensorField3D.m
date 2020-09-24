function bsPlotTensorField3D(M, Ds, options)

    if ~isempty(M)
        figure; 
        range = [prctile(M(:), 10), prctile(M(:), 90)];
        [n1, n2, n3] = size(Ds);

        xx = 1 : n1;
        yy = 1 : n2;
        zz = 1 : n3;
        
        for i = 1 : 3
            subplot(1, 3, i);
            
            bsShow3DVolume(M, 1000, range, 1, 1, (n1-2)*1000, [], 'colormap', colormap(gray), 'view', [133.2431   39.8213], 'attributeName', []);
            
            plot_direction(Ds, i, xx, yy, zz, options);
            
            zlabel('Sample number');
        end
    end
end



function plot_direction(Ds, iatt, xx, yy, zz, options)
    
    sub = options.sub;
    len = options.len;
    hold on;
    
    axeschild = get(gca,'children');
    for k = 1 : length(axeschild)
        x = get(axeschild(k),'XData');
        y = get(axeschild(k),'YData');
        z = get(axeschild(k),'ZData');
        
        [n1, n2] = size(x);
        for i = sub : sub : n1- sub
            for j = sub : sub :n2 - sub
                [~, pi] = min(abs(x(i, j) - yy));
                [~, pj] = min(abs(y(i, j) - zz));
                [~, pk] = min(abs(z(i, j) - xx));
                
                V = Ds{pk, pi, pj};
                vv = V(:, iatt)' * len;
                
                % z x y
                % x y z
                vv = [vv(2)  vv(3) vv(1)];
                px = [pi, pj, pk] - 0.5 * vv;
                py = [pi, pj, pk] + 0.5 * vv;
                
                plot3([px(1) py(1)], [px(2), py(2)], [px(3), py(3)], 'r-', 'linewidth', 2);
                
            end
        end
    end
    
end

