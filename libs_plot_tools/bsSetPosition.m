function bsSetPosition(pw, ph)
    screenSize = get(0, 'ScreenSize');
    ws = screenSize(3);
    hs = screenSize(4);
    
    set(gcf, 'position', [ws*.05, hs*.05, ws*pw, hs*ph]);
end