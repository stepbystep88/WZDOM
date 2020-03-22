clear X map;
imglist = {'flujet', ... Fluid Jet
           'spine', ... Bone
           'gatlin', ... Gatlinburg
           'durer', ... Durer
           'detail', ... Durer Detal
           'cape', ... Cape Cod
           'clown', ... Clown
           'earth', ... Earth
           'mandrill', ... Mandrill
           'spiral'};

colorlabels = {'default', 'hsv','hot','pink',...
               'cool','bone','prism','flag',...
               'gray','rand'};
;

for i=1:10,
   load(imglist{i},'X','map');
   figure;imagesc(X); colorbar;colormap('gray')
end
%colormap(map);
%colormap(colorlabels{1});
%axis off;