function [R, Rp, Rs, Rd] = bsCreateRMatrix(index, sizeAtom, sampNum)

    ncell = length(index);
    
    R = cell(1, ncell);
    Rp = cell(1, ncell);
    Rs = cell(1, ncell);
    Rd = cell(1, ncell);
    
    for i = 1 : ncell
        dr = sparse(sizeAtom, sampNum);
        dp = sparse(sizeAtom, sampNum*3);
        ds = sparse(sizeAtom, sampNum*3);
        dd = sparse(sizeAtom, sampNum*3);
        
        s = index(i);
        for j = 1 : sizeAtom
            dr(j, s + j - 1) = 1;
            dp(j, s + j - 1) = 1;
            ds(j, sampNum + s + j - 1) = 1;
            dd(j, sampNum*2 + s + j - 1) = 1;
        end
        
        R{i} = dr;
        Rp{i} = dp;
        Rs{i} = ds;
        Rd{i} = dd;
    end
   
end