function str = bsGetCurentDateAsString()
    formatOut = 'yyyy-mm-dd';
    str = datestr(datetime('now'), formatOut);
end