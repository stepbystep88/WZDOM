function result = bsButtBandPassFilter(data, lowf, highf, dt)

    n = 5;
    
    if nargin >= 4
        ft = 1000 / dt / 2;
        
    %     [b, a] = butter(n, [lowf/ft, highf/ft], 'bandpass');
    %     result = filtfilt(b, a, data);

        [b, a]= ellip(n, 3, 30, [lowf/ft, highf/ft]);         
        result = filtfilt(b, a, data);
    else
        
        [b, a]= ellip(n, 3, 30, [lowf, highf]);         
        result = filtfilt(b, a, data);
        
    end
    

end