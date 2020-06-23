
%% Active subspace plots 
% Change to this directory to before running

clear variables
close all
clc

addpath(genpath('../../'))

num_local_linear_grad_points = 30;
num_bootstrap = 10000;

%% 
%load('Large_FS_cost_scenarios.mat')
file_id = '6D'
p = 6;
load_name = ['CEP_PCE_', file_id, '_p',num2str(p), '.mat']
load(load_name)

load_name = ['CEP_test_run_', file_id, '_output.mat']
load(load_name)
%load('Small_FS_cost_scenarios.mat')

Xi = Xi_WSPGL1_L;

c1_pre = Xi(:,1);
c2_pre = Xi(:,2);
c_max_wind = Xi(:,3);
c_avg_wind = Xi(:,4);
c_max_gas = Xi(:,5);

%% Ishigami 

[M, m] = size(Xi_CEP_test_run);
%f = FS_cost; title_string =  'Expansion Cost \$'; c = c1_pre; c_name = 'c1_pre'; print_name = ['FS_AS_', file_id, '_p',num2str(p) ]; n = 1;
f = SS_cost; title_string =  'Operations Cost \$'; c = c2_pre; c_name = 'c2_pre';  print_name = ['SS_AS_', file_id, '_p',num2str(p) ]; n = 2;
%f = max_wind; title_string =  'Max Installed Wind MW'; c = c_max_wind; c_name = 'c_max_wind';  print_name = ['max_wind_AS_', file_id, '_p',num2str(p) ]; n = 1;
%f = avg_wind; title_string =  'Avg Wind Avail. MW'; c = c_avg_wind; c_name = 'c_avg_wind';  print_name = ['avg_wind_AS_', file_id, '_p',num2str(p) ]
%f = max_gas; title_string =  'Max Installed Gas MW'; c = c_max_gas; c_name = 'c_max_gas'; print_name = ['max_gas_AS_', file_id, '_p',num2str(p) ]; n = 2;
X = Xi_CEP_test_run;
df = local_linear_gradients(X,f,num_local_linear_grad_points);
% M = 10000;
% X = 2*rand(M, m) - 1;
% X(:,1) = rand(M,1)/4+0.75;
% %X(:,1) = zeros;
% %X(:,1) = rand(M,1)/4-1;
% %X(:,2) = rand(M,1)/2-0.5;
% %X(:,3) = rand(M,1)/2-0.5;
% %X(:,4) = rand(M,1)/4+0.75;
% f = zeros(M,1);
% df = finite_difference_gradients(X,@FS_PCE);
% f = FS_PCE(X);
%df = local_linear_gradients(X,f,6)

sub = compute(df,num_bootstrap);
opts = plot_opts([]);
opts.fontsize = 22;
opts.markersize = 13;
%%
Y = X*sub.eigenvectors(:, 1:m);
%load('ishigami_sensitivity.mat')
%save('ishigami_as')

%%
%close all
% figure 
% eigenvalues(sub.eigenvalues, sub.e_br,opts)
% figure
% eigenvectors(sub.eigenvectors(:,1:2), opts)
% figure
% subspace_errors(sub.sub_br, opts)
%%
%close all
n=1
FF = figure('Position',[1006           1        2355        1367])
subplot(2,3,4)
% Make 1D sufficient summary plot.
scatter(Y(:,1), f,'filled', ...
     'markeredgecolor', 'none', ...'markerfacecolor', opts.color, ...
     'cdata', f, ...
     'marker', opts.marker, ...
     'sizedata', 10*opts.markersize, ...
     'linewidth', opts.linewidth)

% Format plot
% if isempty(opts.title)
%     title('Output', 'fontsize', opts.fontsize)
% else
%     title(opts.title, 'fontsize', opts.fontsize)
% end
%xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
xlabel('Active Variable 1 $(y_1 = \mathbf{w}_1^T\mathbf{\xi})$', 'fontsize', opts.fontsize,'interpreter','latex')

set(gca,  'TickLabelInterpreter','latex',...
    'fontsize', opts.fontsize)
grid('on')
box on
colormap('parula')
T = title(title_string)
set(T,'interpreter','latex')
xlim([-2,2])

subplot(2,3,5)

% Make 2D sufficient summary plot.
scatter(Y(:,1), Y(:,2), 'filled', ...
        'cdata', f, ....
        'marker', opts.marker, ...
        'sizedata', 10*opts.markersize, ...
        'linewidth', opts.linewidth)

% Format plot
% if isempty(opts.title)
%     title('Output', 'fontsize', opts.fontsize)
% else
%     title(opts.title, 'fontsize', opts.fontsize)
% end
%xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
%ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
xlabel('Active Variable 1  $(y_1 = \mathbf{w}_1^T\mathbf{\xi})$', 'fontsize', opts.fontsize,'interpreter','latex')
ylabel('Active Variable 2 $(y_2 = \mathbf{w}_2^T\mathbf{\xi})$', 'fontsize', opts.fontsize,'interpreter','latex')


T = title(title_string)
set(T,'interpreter','latex')
set(gca, 'TickLabelInterpreter','latex',...
    'fontsize', opts.fontsize)
grid('on')
cb = colorbar('TickLabelInterpreter','latex')
colormap('parula')
%T= title(cb,'f(x)');
%set(T,'interpreter','latex')
a =get(cb);
%set(cb,'Position',[a.Position(1) a.Position(2) a.Position(3) a.Position(4)])
%set(cb,'Position',[1.06*a.Position(1) a.Position(2) a.Position(3) 1*a.Position(4)])
xlim([-2,2])
ylim([-2,2])
box on 

subplot(2,3,3)
F1 = eigenvalues(sub.eigenvalues, sub.e_br,opts);
disp(F1)

subplot(2,3,2)
F1 = eigenvectors(sub.eigenvectors(:,1:2), opts);
disp(F1)

subplot(2,3,6)
subspace_errors(sub.sub_br, opts)
disp(F1)
% e = sub.eigenvalues;
% hold on
% 
% 
% % Plot eigenvalues.
% semilogy(e, ...
%          'markeredgecolor', 'k', ...
%          'markerfacecolor', opts.color, ...
%          'color', opts.color, ...
%          'marker', opts.marker, ...
%          'markersize', 13, ...
%          'linewidth', opts.linewidth)
% 
% % Format plot.
% title(opts.title, 'fontsize', opts.fontsize)
% 
% % if isempty(opts.xlabel)
% %     xlabel('Index', 'fontsize', opts.fontsize)
% % else
% %     xlabel(opts.xlabel, 'fontsize', opts.fontsize)
% % end
% 
% ylabel('Eigenvalues $\lambda_i$', 'fontsize', opts.fontsize,'interpreter','latex')
% xlabel('$i$', 'fontsize', opts.fontsize,'interpreter','latex')
% 
% 
% set(gca, ...
%     'XLim', [1, m], ...
%     'XTick', 1:m, ...
%     'XScale', 'Linear', ...
%     'YScale', 'Log', ...
%     'fontsize', opts.fontsize)
% grid on
% box on
subplot(2,3,1)

%load('ishigami_sensitivity.mat')
%load('ishigami_sensitivity.mat')
%as = as./norm(as,2);

%
W = sub.eigenvectors;
e = sub.eigenvalues;
alpha = zeros(m,1);


% n = m;
% for i = 1:m
%     nu(i) =  W(i,1:n).^2*e(1:n);
% end
% nu = alpha/norm(alpha)

nu = diag(W*diag(e)*W');
nu = nu/norm(nu)

plot(1:m,nu, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'b', ...
         'marker', opts.marker, ...
         'markersize', opts.markersize, ...
         'linewidth', opts.linewidth)
 hold on
 alpha = [];

for i = 1:m
    alpha(i) =  W(i,1:n).^2*e(1:n);
end



alpha = alpha/norm(alpha)

plot(1:m,alpha, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'r', ...
         'marker', 's', ...
         'markersize', 0.8*opts.markersize, ...
         'linewidth', opts.linewidth)
 hold on
 

 %load('c1_pre_SP.mat')
 %load('c2_pre_SP.mat')
 %load('c1_pre_Small_FS.mat')
 %load('c2_pre_Small_FS.mat')

%[Stot,S1] = get_sobol_indices(c1_pre,index_pc);
%[Stot,S1] = get_sobol_indices(c2_pre,index_pc); 

%load(ref_load_name,c_name,'index_pc')
[Stot,S1] = get_sobol_indices(c,index_pc);

 plot(1:m,Stot, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'g', ...
         'marker', opts.marker, ...
         'markersize',opts.markersize, ...
         'linewidth', opts.linewidth)
hold on;
plot(1:m,S1, ...
 'markeredgecolor', 'k', ...
 'markerfacecolor', 'm', ...
 'marker', 's', ...
 'markersize', 0.8*opts.markersize, ...
 'linewidth', opts.linewidth)

L = legend('$\nu$', ['$\alpha (n=$',num2str(n),')'] ,'$\tau$','$\tau_1$','location','best')
set(L,'interpreter','latex','fontsize',20)
set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'YScale', 'linear', ...
    'xticklabels',{'$R^{cap}$', '$c^{loss}$' , '$c^{exce}$', '$c^{ng}$','$c^{wind}$','$t^{cap}$'},...
    'TickLabelInterpreter','latex',...
    'fontsize', opts.fontsize)
grid on
ylim([0 1])      
%'xticklabels',{'$S_w$','$W_{fw}$','$A$','$\Lambda$','$q$','$\lambda$','$t_c$','$N_Z$','$W_{dg}$','$W_p$'},...
%'xticklabels',{'$\beta_1$','$\beta_2$','$\beta_3$','$\rho_1$','$\gamma_1$','$\gamma_2$','$\omega$','$\psi$'},...
%xticklabels({'S_w','W_{fw}','A','\Sigma','q','\lambda','t_c','N_Z','W_{dg}','W_p'})

%close all
W1 = W(:,1)

%%
%pause(1)
%print(print_name,'-dpng','-r300')
%saveas(FF,print_name,'epsc')

%%
% 
% F = figure('Position',[1987         347         834         748])
% %close all
% % Make 3D sufficient summary plot.
% scatter3(Y(:,1), Y(:,2), Y(:,3), 'filled', ...
%         'cdata', f, ....
%         'marker', opts.marker, ...
%         'sizedata', 50*opts.markersize, ...
%         'linewidth', opts.linewidth)
%     
%     
%     
% xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
% ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
% zlabel('Active Variable 3 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
% T = title(title_string)
% set(T,'interpreter','latex')
% set(gca, 'fontsize', opts.fontsize)
% grid('on')
% cb = colorbar();
% a = cb.Position;
% set(cb,'Position',[0.99*a(1) a(2) a(3) a(4)])
% colormap('parula')
% xlim([-2,2])
% ylim([-2,2])

%%
%saveas(F,[print_name,'3DSS'],'epsc')


%%


