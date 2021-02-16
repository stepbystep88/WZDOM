function bsPlotTensorField3D(M, Ds, lens, options)

    if ~isempty(M)
%         figure; 
        range = [prctile(M(:), 10), prctile(M(:), 90)];
        [n1, n2, n3] = size(Ds);

        xx = 1 : n1;
        yy = 1 : n2;
        zz = 1 : n3;
        
        strs = {'first', 'second', 'third'};
        plotParam = bsGetDefaultPlotSet();
        plotParam.fontsize = 18;
        
        for i = 1 : 4
%             subplot(1, 3, i);
            figure; 
            set(gcf, 'position', [ 46         189        750         500]);
            subplot('position', [0.116666666666667 0.153999572753906 0.874 0.836000427246094]);
            bsShow3DVolume(M, 1000, range, 1, 1, (n1-18)*1000, [], 'colormap', bsGetColormap('seismic'), 'view', [ 124.7150   60.9999], 'attributeName', []);
            
            if i < 4
                plot_direction(Ds, lens, 3 - i + 1, xx, yy, zz, options);
                title(sprintf('The %s eigenvector field', strs{i}), 'fontweight', 'bold', 'fontsize', 12);
            end
            
%             zlabel('Sample number');
            zlabel('采样点数');
            bsSetDefaultPlotSet(plotParam);
        end
    end
end



function plot_direction(Ds, lens, iatt, xx, yy, zz, options)
    
    sub = options.sub;
    len_coef = options.len;
    hold on;
    
    axeschild = get(gca,'children');
    for k = 1 : 1
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
                len = lens(pk, pi, pj) * len_coef;
                vv = V(:, iatt)' * len_coef;
                
                % z x y
                % x y z
                vv = [vv(2)  vv(3) vv(1)];
                px = [pi, pj, pk] - 0.5 * vv;
                py = [pi, pj, pk] + 0.5 * vv;
                
                plot3([px(1) py(1)], [px(2), py(2)], [px(3), py(3)], 'k-', 'linewidth', 2);
                
            end
        end
    end
    
end

