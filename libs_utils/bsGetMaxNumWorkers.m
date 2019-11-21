function numWorkers = bsGetMaxNumWorkers()
    myCluster = parcluster('local');
    numWorkers = myCluster.NumWorkers;
    
    fprintf('The number of workers of local cluster is %d\n', numWorkers);
end