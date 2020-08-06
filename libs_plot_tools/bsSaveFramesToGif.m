function bsSaveFramesToGif(frames, fileName)
    nFrame = size(frames, 2);
    
%     set(gcf, 'color', 'white');
%     set(gca, 'color', 'white');
 
    for i = 1 : nFrame
        
        im = frame2im(frames(:, i));
        [I, map] = rgb2ind(im, 256);
        
        
        if i == 1
            imwrite(I, map, fileName, 'gif', 'DelayTime',0.1, 'loopcount', inf);
        else
            if(mod(i, 2) == 0)
                imwrite(I, map, fileName, 'gif', 'DelayTime',0.05, 'writemode','append');
            end
        end
    end
end