function data = bsJointData(cdata)

    data = [];
    
    for i = 1 : length(cdata)
        icdata = cdata(i);
        data = [data, icdata{1}];
    end
end