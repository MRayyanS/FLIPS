%% Recovering output images

% recover output F_star

Recover_output_F_star   = recover_vectorised_image_from_dct_transform(F_star, dictionary, patch_dimensions) ;
imout_F_star            = patch2image(Recover_output_F_star,im.noise,store_inp);

% recover output FLIPS

Recover_output_FLIPS   = recover_vectorised_image_from_dct_transform(F_FLIPS, dictionary, patch_dimensions) ;
imout_FLIPS            = patch2image(Recover_output_FLIPS,im.noise,store_inp);

% recover output CP

Recover_output_CP   = recover_vectorised_image_from_dct_transform(F_CP, dictionary, patch_dimensions) ;
imout_CP            = patch2image(Recover_output_CP,im.noise,store_inp);

% recover output csalsa

Recover_output_csalsa   = recover_vectorised_image_from_dct_transform(F_csalsa, dictionary, patch_dimensions) ;
imout_csalsa            = patch2image(Recover_output_csalsa,im.noise,store_inp);

% show all the images

imshow([imout_F_star imout_FLIPS imout_CP imout_csalsa]) ;







%% function: recover_vectorised_image_from_dct_transform

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



%% function: ptch2image function

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
