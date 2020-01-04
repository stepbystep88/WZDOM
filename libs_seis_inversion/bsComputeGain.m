function [waveletGain] = bsComputeGain(post, syn)
    waveletGain = post'*syn / (syn'*syn);
%     return;
    [row, column] = size(syn);
    if column < row
        column = row;
    end
    synDiff = zeros(column-1, 1);
    postDiff = zeros(column-1, 1);
    for i= 2 : 1 : column
        synDiff(i) = syn(i) - syn(i-1);
        postDiff(i) = post(i)-post(i-1);
    end

    synPeaTroValue = zeros(column-1, 1);
    postPeaTroValue = zeros(column-1, 1);
    k = 1;
    m = 1;
    for i = 1: 1 : column-1
        if (synDiff(i) > 0) && (synDiff(i+1) < 0)
            synPeaTroValue(k) = abs(syn(i+1));
            k = k+1;
        end
        if (synDiff(i) < 0) && (synDiff(i+1) > 0)
            synPeaTroValue(k) = abs(syn(i+1));
            k = k+1;
        end
        if (postDiff(i) > 0) &&( postDiff(i+1) < 0)
            postPeaTroValue(k) = abs(post(i+1));
            m = m+1;
        end
        if (postDiff(i) < 0) && (postDiff(i+1) > 0)
            postPeaTroValue(k) = abs(post(i+1));
            m = m+1;
        end
    end

    synPeaTroValue = sort(synPeaTroValue(synPeaTroValue >0), 'descend');
    postPeaTroValue = sort(postPeaTroValue(postPeaTroValue >0), 'descend');

    waveletGain  = mean(postPeaTroValue(1:5)) / mean(synPeaTroValue(1:5));
end
