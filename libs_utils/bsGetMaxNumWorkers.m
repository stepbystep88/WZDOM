function numWorkers = bsGetMaxNumWorkers()
    myCluster = parcluster('local');
    numWorkers = myCluster.NumWorkers;
end