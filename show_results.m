%% Recover images 

% Image recovery 
Recover_input                       = X.no_noise;
Recover_input_noise                 = X.noise;

% recover output by doing idct2 transform
Recover_output_FLIPS                = recover_vectorised_image_from_dct_transform(F, dictionary, patch_dimensions) ;


% final recovered images
im_input            = patch2image(Recover_input,im.noise,store_inp);
im_noisy            = patch2image(Recover_input_noise,im.noise,store_inp);
im_output           = patch2image(Recover_output_FLIPS,im.noise,store_inp);


%% Show images and export pdfs

% original image
fig_original_image = tiledlayout(1,1,'Padding','tight');
fig_original_image.Units = 'inches';
fig_original_image.OuterPosition = [0.25 0.25 3 3];
nexttile;

imagesc(reordered_im_input) ;
colormap(gray) ;

im = gcf;
exportgraphics(im,'lena_original.pdf','ContentType','vector') ;


% noisy image
fig_original_image = tiledlayout(1,1,'Padding','tight');
fig_original_image.Units = 'inches';
fig_original_image.OuterPosition = [0.25 0.25 3 3];
nexttile;

imagesc(reordered_im_input_noise) ;
colormap(gray) ;

im = gcf;
exportgraphics(im,'lena_noisy.pdf','ContentType','vector') ;


% recovered image
fig_original_image = tiledlayout(1,1,'Padding','tight');
fig_original_image.Units = 'inches';
fig_original_image.OuterPosition = [0.25 0.25 3 3];
nexttile;

imagesc(reordered_im_output_FLIPSe) ;
colormap(gray) ;

im = gcf;
exportgraphics(im,'lena_recovered.pdf','ContentType','vector') ;
% Set paper size to figure size
fig.PaperSize = [fig_width, fig_height];



%% Different functions needed


function x_output = patch2image(x_input,im,store_inp)
% Reconstructing image from the extracted/reconstructed (sliding) image patches

[imdim1,imdim2]=size(im);
indices_store_patches_inp= reshape(1:imdim1*imdim2,[imdim1 imdim2]);

rebuild_inp_rec     =  zeros(imdim1,imdim2);
for i=1:imdim1
    for j=1:imdim2
        idx_inp = store_inp==indices_store_patches_inp(i,j);
        rebuild_inp_rec(i,j)=mean(x_input(idx_inp),'all');
    end
end
x_output = rebuild_inp_rec;

end


% computing idct2 transform of F, i.e., getting images from their dct transforms

% phiF = phi(F), where phi = idct2 matrix

function phiF    = recover_vectorised_image_from_dct_transform(F, dictionary, patch_dimensions)

[n, N] = size(F) ;

phiF  = zeros(n,N) ;

for i = 1:N

    patch = reshape(F(:,i), patch_dimensions) ; % for fast computation of dct/basis coefficients

    if strcmp(dictionary,'dct')
        idct2_im     = idct2(patch) ;
        phiF(:,i) = reshape(idct2_im, [n 1]) ;
    end

end

end

