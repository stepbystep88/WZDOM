function data = bsSpmdModule(fcn, inputData, numWorkers)
    [groupData, ~] = bsSeparateData(inputData, numWorkers);
    
    spmd
        cdata = fcn(groupData);
    end
    
    data = bsJointData(cdata);
end