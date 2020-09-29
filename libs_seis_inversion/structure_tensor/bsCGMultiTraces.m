function xk = bsCGMultiTraces(xs, b, bsAx, maxIter, S, gst_options, is3D, nInline, nCrossline)
    xk = xs;
    rk = b - bsAx(xk);
    pk = rk;
    
    if is3D
        Mult = @(S,u)cat(4, S(:,:,:,1).*u(:,:,:,1) + S(:,:,:,4).*u(:,:,:, 2) + S(:,:,:,6).*u(:,:,:, 3), ...
                        S(:,:,:,2).*u(:,:,:,2) + S(:,:,:,4).*u(:,:,:, 1) + S(:,:,:,5).*u(:,:,:, 3), ...
                        S(:,:,:,3).*u(:,:,:,3) + S(:,:,:,5).*u(:,:,:, 2) + S(:,:,:,6).*u(:,:,:, 1) ...
                        );
    else
        
        Mult = @(S,u)cat(3, S(:,:,1).*u(:,:,1) + S(:,:,3).*u(:,:,2), ...
                      S(:,:,3).*u(:,:,1) + S(:,:,2).*u(:,:,2) );
    end
    
    nabla =  gst_options.nabla;
    tau = gst_options.tau;
    
    for k = 1 : maxIter
%         alphak = 0.01;
        Apk = bsAx(pk);
        try
            t1 = rk.^2;
            t2 = pk .* Apk;
            
            alphak = sum(t1(:)) / sum(t2(:));
            if alphak < -5 || alphak > 5
                alphak = 0.01;
            end
            
        catch ME
%             warning(ME.msg);
            warning('Getting alphak encouters an error.');
            alphak = 0.001;
        end
        
%         if alphak < 0 || alphak > 1
%             disp(alphak);
%             alphak = 0.001;
%         end
%         Ip = exp(xk);
        Ip = xk;
        
        if is3D
            fprintf('Excuting the %d-th/%d iteration...\n', k, maxIter);
            xk3D = bsReshapeDataAs3D(Ip, nInline, nCrossline);
            delta = bsReshapeDataAs2D(div(Mult(S, nabla(xk3D))));
        else
            delta = div(Mult(S, nabla(Ip)));
        end
        
        xk = xk + alphak * pk + tau * delta;
%         xk = log(exp(xk + alphak * pk) + tau * delta);
        rk_new = rk - alphak * Apk;
        
        t1 = rk_new.^2;
        t2 = rk.^2;
        
        betak = sum(t1(:)) / sum(t2(:));
        pk = rk_new + betak * pk;
        rk = rk_new;
        
    end
    
    
end