%% Plotting relevant things

% Plotting recovered solutions

fig_solutions = tiledlayout(1,1,'Padding','tight');
fig_solutions.Units = 'inches';
fig_solutions.OuterPosition = [0.25 0.25 3 3];
nexttile;
i = 1:n ;
plot(i,F_true,'-o','MarkerSize',5);
grid on;
hold on;
plot(i,F_rec,'-x','MarkerSize',5) ;
axis padded;
xlabel('$i = 1,2,\ldots, n = 500$', 'FontSize',10,'Interpreter','latex');
legend('$f_{tr}$', '$ f^* $ ', 'FontSize',15,'Interpreter','latex','Location','northeast');
im = gcf;
exportgraphics(im,'n500_solution.pdf','ContentType','vector') ;


% Plotting sub-optimality of FLIPS

sub_opt = eta - eta(end).*ones(size(eta)) ;

fig_sub_opt = tiledlayout(1,1,'Padding','tight');
fig_sub_opt.Units = 'inches';
fig_sub_opt.OuterPosition = [0.25 0.25 3 3];
nexttile;

t = 1:maxiter ; % iteration index
semilogy(t,sub_opt,'-o','MarkerSize',4);
grid on;
axis padded;
xlabel('Iterations, $k$', 'FontSize',10,'Interpreter','latex');
legend('$ \eta(h_t) - \eta(h^*) $ ', 'FontSize',10,'Interpreter','latex','Location','southwest');
% ylabel('$ \eta(h_t) - \eta(h^*) $', 'FontSize',20,'Interpreter','latex');

im = gcf;
exportgraphics(im,'n500_sub_optimality.pdf','ContentType','vector') ;


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
exportgraphics(im,'n500_distance.pdf','ContentType','vector') ;




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
legend('$ \gamma_t $ ', 'FontSize',10,'Interpreter','latex','Location','northwest');
% ylabel('$ \gamma(h_t) $', 'FontSize',20,'Interpreter','latex');

im = gcf;
exportgraphics(im,'n500_stepsize.pdf','ContentType','vector') ;



% Plotting first-order sub-optimality

fig_first_order_optimality = tiledlayout(1,1,'Padding','tight');
fig_first_order_optimality.Units = 'inches';
fig_first_order_optimality.OuterPosition = [0.25 0.25 3 3];
nexttile;

t = 1:maxiter ; % iteration index
semilogy(t,optimality_check,'-o','MarkerSize',4);
grid on;
axis padded;
xlabel('Iterations, $k$', 'FontSize',10,'Interpreter','latex');
legend('$ \frac{1}{2} + \frac{ \langle \nabla \eta (h) , h \rangle }{ 2 || \nabla \eta (h) ||_1 } $ ', 'FontSize',10,'Interpreter','latex','Location','southwest');
% ylabel('$ 1 + \frac{ \langle \nabla \eta (h) , h \rangle }{ || \nabla \eta (h) ||_1 } $', 'FontSize',20,'Interpreter','latex');

im = gcf;
exportgraphics(im,'n500_first_order_opt_check.pdf','ContentType','vector') ;