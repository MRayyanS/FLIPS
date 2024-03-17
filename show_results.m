%% Recover images 

% Image recovery 
Recover_input             = X_no_noise ;
Recover_input_noise       = X_noise ;

l2_recovered_output       = D*F(:,1) ;
l2_output_measurement     = C*l2_recovered_output ;
l2_recovered_output       = ((l2_output_measurement'*X_noise)/norm(l2_output_measurement, 2)^2)*l2_recovered_output ;
l2_recovered_image        = reshape(l2_recovered_output, input_im_size) ;

FLIPS_recovered_output    = D*F(:,end) ;
FLIPS_output_measurement  = C*FLIPS_recovered_output ;
FLIPS_recovered_output    = ( (FLIPS_output_measurement'*X_noise)/norm(FLIPS_output_measurement,2)^2 ).*FLIPS_recovered_output ;

Recovered_image       = reshape(FLIPS_recovered_output, input_im_size) ;




%% Show images and export pdf

% original image
fig_original_image = tiledlayout(1,1,'Padding','tight');
fig_original_image.Units = 'inches';
fig_original_image.OuterPosition = [0.25 0.25 3 3];
nexttile;

imagesc(im_original) ;
colormap(gray) ;

im = gcf;
exportgraphics(im,'lena_original.pdf','ContentType','vector') ;


% recovered image image
fig_recovered_image = tiledlayout(1,1,'Padding','tight');
fig_recovered_image.Units = 'inches';
fig_recovered_image.OuterPosition = [0.25 0.25 3 3];
nexttile;

imagesc(Recovered_image) ;
colormap(gray) ;

im = gcf;
exportgraphics(im,'lena_FLIPS_recovered.pdf','ContentType','vector') ;


%% Plotting relevant things

% Plotting sub-optimality of FLIPS
sub_opt                   = eta - eta(end).*ones(size(eta)) ;
fig_sub_opt               = tiledlayout(1,1,'Padding','tight');
fig_sub_opt.Units         = 'inches';
fig_sub_opt.OuterPosition = [0.25 0.25 3 3];
nexttile;

t = 1:maxiter ; % iteration index
semilogy(t,sub_opt,'-o','MarkerSize',4);
grid on;
axis padded;
xlabel('Iterations, $k$', 'FontSize',10,'Interpreter','latex');
legend('$ \eta(h_t) - \eta(h^*) $ ', 'FontSize',10,'Interpreter','latex','Location','southwest');
% ylabel('$ \eta(h_t) - \eta(h^*) $', 'FontSize',10,'Interpreter','latex');

im = gcf;
exportgraphics(im,'lena_sub_optimality.pdf','ContentType','vector') ;


% Plotting the distance to the true solution
distance = zeros(1, maxiter) ;
for t = 1:maxiter
    distance(t) = norm(F(:,t) - F_true, 2) ;
end

fig_distance = tiledlayout(1,1,'Padding','tight');
fig_distance.Units = 'inches';
fig_distance.OuterPosition = [0.25 0.25 3 3];
nexttile;

t = 1:maxiter ; % iteration index
semilogy(t,distance,'-x','MarkerSize',4);
grid on;
axis padded;
xlabel('Iterations, $k$', 'FontSize',10,'Interpreter','latex');
legend('$ || f_t - f_{tr} ||_2 $ ', 'FontSize',10,'Interpreter','latex','Location','northeast');
% ylabel('$ || f_t - f_{tr} ||_2 $', 'FontSize',10,'Interpreter','latex');

im = gcf;
exportgraphics(im,'lena_distance.pdf','ContentType','vector') ;


% Plotting the step-size
fig_stepsize = tiledlayout(1,1,'Padding','tight');
fig_stepsize.Units = 'inches';
fig_stepsize.OuterPosition = [0.25 0.25 3 3];
nexttile;

t = 1:maxiter ; % iteration index
plot(t,gamma,'-*','MarkerSize',4);
grid on;
axis padded;
xlabel('Iterations, $k$', 'FontSize',10,'Interpreter','latex');
legend('$ \gamma_t $ ', 'FontSize',10,'Interpreter','latex','Location','southeast');
% ylabel('$ \gamma(h_t) $', 'FontSize',20,'Interpreter','latex');

im = gcf;
exportgraphics(im,'lena_stepsize.pdf','ContentType','vector') ;