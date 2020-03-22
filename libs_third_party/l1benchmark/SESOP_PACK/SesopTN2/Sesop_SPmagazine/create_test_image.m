% Create/read test image X00 and resize it to [m0,n0]
%
% create_test_image.m  

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %     Read/Create test image
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

   if FlagSparseImage,
      X00=10*(sprandn(m0,n0,image_sparsity));%save X00 X00;
%       if m0==64 && n0==128, load X00_sparse64_128; 
%       else                  load X00_sparse64_64; 
%       end
      X00=full(X00);
   else 
      if strcmp(par.ImageName,'Phantom'), X00=phantom(n0);end 
      if strcmp(par.ImageName,'Lena'),
         X00=imread('..\..\..\images\lena.png');
      end
      if strcmp(par.ImageName,'Peppers')
         X00=imread('..\..\..\images\peppers.png');
      end
      X00=double(X00);X00=X00/max(X00(:)); %figure;imagesc(X00);colorbar
   end
   
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   %           Resize the image to [m0,n0]
   %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   [m,n]=size(X00);  resize_factor=1/ceil(m/m0);
   if resize_factor~=1, X00=imresize(X00,resize_factor);  [m,n]=size(X00);end
   
