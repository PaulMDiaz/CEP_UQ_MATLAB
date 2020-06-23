%% Active subspace plots 
%clear variables
close all
%clc

addpath(genpath('../../'))

%% 
%load('Large_FS_cost_scenarios.mat')

file_id = '6D'

p_ref = 6;
ref_load_name = ['CEP_PCE_', file_id, '_p',num2str(p_ref), '.mat']
p  = 6;
load_name = ['CEP_test_run_', file_id, '_output.mat']
load(load_name)
load_name = ['CEP_PCE_', file_id, '_p',num2str(p), '.mat']
load(load_name)
%load('Small_FS_cost_scenarios.mat')`   
%%

Xi = Xi_WSPGL1_L;
c1_pre = Xi(:,1);
c2_pre = Xi(:,2);
c_max_wind = Xi(:,3);
c_avg_wind = Xi(:,4);
c_max_gas = Xi(:,5);
%%
%save(load_name,'c1_pre','c2_pre','c_max_wind','c_avg_wind','c_max_gas','index_pc','Xi_WSPGL1_L', 'Xi_WSP_L', 'Xi_RR_L', 'Xi_IT_L', 'Sys_info')


%% Ishigami 
m = 6; n = 6;
%p = 4
%M = size(Xi_CEP_test_run,1);
% fs = FS_cost;
% X = Xi_CEP_test_run;
%  

%X = zeros(M,m)

%%
M = 10000; X = 2*rand(M,m)-1; 

%X = randn(M,m)
%X = Xi_CEP_test_run; M = size(Xi_CEP_test_run,1);
clc
tic


%X(:,1) = zeros;
%X(:,1) = rand(M,1)/4-1;
%X(:,1) = zeros;
%X(:,3) = 2*rand(M,1)-1;
%X(:,2) = rand(M,1)*0.05-1;
%X(:,2) = 2*rand(M,1)-1;
%X(:,3) = rand(M,1)*0.05+0.95;
%X(:,4) = 0.1*rand(M,1)-1;
%X(:,5) = 0.1*rand(M,1)+0.9;
%f = zeros(M,1);
%df = finite_difference_gradients(X,@FS_PCE);
fs = FS_PCE(X,c1_pre,index_pc);
ss = SS_PCE(X,c2_pre,index_pc);
mg = MG_PCE(X,c_max_gas,index_pc);
mw = MW_PCE(X,c_max_wind,index_pc);
toc


%% %%%%%%%%%%%%%%%%%%%%%%%%%% Expansion Costs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df = local_linear_gradients(X,fs,20)
sub = compute(df);
opts = plot_opts([]);
opts.fontsize = 36;
Y = X*sub.eigenvectors(:, 1:m);
figure()
set(gcf, 'Position', [100, 100, 1024, 1024])
subplot(2,2,3)
% Make 1D sufficient summary plot.
scatter(Y(:,1), fs,'filled', ...
     'markeredgecolor', 'none', ...'markerfacecolor', opts.color, ...
     'cdata', fs, ...
     'marker', opts.marker, ...
     'sizedata', 10*opts.markersize, ...
     'linewidth', opts.linewidth)

xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
box on
colormap('parula')
T = title('Expansion Cost \$')

set(T,'interpreter','latex')
xlim([-2,2])

subplot(2,2,4)

% Make 2D sufficient summary plot.
scatter(Y(:,1), Y(:,2), 'filled', ...
        'cdata', fs, ....
        'marker', opts.marker, ...
        'sizedata', 10*opts.markersize, ...
        'linewidth', opts.linewidth)

xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
T = title('Expansion Cost \$')
set(T,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
cb = colorbar();
colormap('parula')
a =get(cb);
xlim([-2,2])
ylim([-2,2])
box on 

subplot(2,2,2)
e = sub.eigenvalues;
hold on
% Plot eigenvalues.
semilogy(e, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', opts.color, ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)

% Format plot.
title(opts.title, 'fontsize', opts.fontsize)

ylabel('Eigenvalues $\lambda_i$', 'fontsize', opts.fontsize,'interpreter','latex')
xlabel('$i$', 'fontsize', opts.fontsize,'interpreter','latex')


set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'XScale', 'Linear', ...
    'YScale', 'Log', ...
    'fontsize', opts.fontsize)
grid on
box on
subplot(2,2,1)

subplot(2,3,1)

%load('ishigami_sensitivity.mat')
%load('ishigami_sensitivity.mat')
%as = as./norm(as,2);

%
W = sub.eigenvectors;
e = sub.eigenvalues;
alpha = zeros(m,1);

mm = m;

for i = 1:mm
    alpha(i) = W(i,1:n).^2*e(1:n);
end
alpha = alpha/norm(alpha)

plot(1:m,alpha, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'b', ...
         'marker', opts.marker, ...
         'markersize', opts.markersize, ...
         'linewidth', opts.linewidth)
 hold on
 
 n = 1
for i = 1:mm
    alpha(i) = W(i,1:n).^2*e(1:n);
end
alpha = alpha/norm(alpha)

plot(1:m,alpha, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'r', ...
         'marker', opts.marker, ...
         'markersize', opts.markersize, ...
         'linewidth', opts.linewidth)
 hold on
 

 %load('c1_pre_SP.mat')
 %load('c2_pre_SP.mat')
 %load('c1_pre_Small_FS.mat')
 %load('c2_pre_Small_FS.mat')

%[Stot,S1] = get_sobol_indices(c1_pre,index_pc);
%[Stot,S1] = get_sobol_indices(c2_pre,index_pc); 

load(ref_load_name,'c1_pre','index_pc')
[Stot,S1] = get_sobol_indices(c1_pre,index_pc);

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
 'marker', opts.marker, ...
 'markersize', opts.markersize, ...
 'linewidth', opts.linewidth)

L = legend('$\nu$', ['$\alpha (n=$',num2str(n),')'] ,['$\tau(p=$',num2str(p_ref),')'],['$\tau_1(p=$',num2str(p_ref),')'],'location','best')
set(L,'interpreter','latex','fontsize',20)
set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'xticklabels',{'$R^{cap}$', '$c^{loss}$' , '$c^{oload}$', '$c^{ng}$','$c^{wind}$','$t^{avail}$'},...
    'TickLabelInterpreter','latex',...
    'fontsize', opts.fontsize)
grid on
axis square
ylim([0 1])      
% Use for 5D
% set(gca, ...
%     'XLim', [1, m], ...
%     'XTick', 1:m, ...
%     'xticklabels',{'$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
%     'TickLabelInterpreter','latex',...
%     'fontsize', opts.fontsize)
% grid on

%%
% print_name = [file_id, 'PCE_FS_resamp_M',num2str(M)]
% print(print_name,'-dpng','-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Operation Costs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%df = finite_difference_gradients(X,@SS_PCE);

%ss = SS_cost;
df = local_linear_gradients(X,ss,20);

sub = compute(df);
opts = plot_opts([]);
opts.fontsize = 20;
Y = X*sub.eigenvectors(:, 1:m);
figure()
set(gcf, 'Position', [100, 100, 1024, 1024])
subplot(2,2,3)
% Make 1D sufficient summary plot.
scatter(Y(:,1), ss,'filled', ...
     'markeredgecolor', 'none', ...'markerfacecolor', opts.color, ...
     'cdata', ss, ...
     'marker', opts.marker, ...
     'sizedata', 10*opts.markersize, ...
     'linewidth', opts.linewidth)

xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
box on
colormap('parula')
T = title('Operation Cost \$')
set(T,'interpreter','latex')
xlim([-2,2])

subplot(2,2,4)

% Make 2D sufficient summary plot.
scatter(Y(:,1), Y(:,2), 'filled', ...
        'cdata', ss, ....
        'marker', opts.marker, ...
        'sizedata', 10*opts.markersize, ...
        'linewidth', opts.linewidth)

xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
T = title('Operation Cost \$')
set(T,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
cb = colorbar();
colormap('parula')
a =get(cb);
xlim([-2,2])
ylim([-2,2])
box on 

subplot(2,2,2)
e = sub.eigenvalues;
hold on


% Plot eigenvalues.
semilogy(e, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', opts.color, ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)

% Format plot.
title(opts.title, 'fontsize', opts.fontsize)

ylabel('Eigenvalues $\lambda_i$', 'fontsize', opts.fontsize,'interpreter','latex')
xlabel('$i$', 'fontsize', opts.fontsize,'interpreter','latex')


set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'XScale', 'Linear', ...
    'YScale', 'Log', ...
    'fontsize', opts.fontsize)
grid on
box on
subplot(2,2,1)
W = sub.eigenvectors;
e = sub.eigenvalues;
alpha = zeros(m,1);



for i = 1:m
    alpha(i) = W(i,1:m).^2*e(1:m);
end
alpha = alpha/norm(alpha)

plot(1:n,alpha, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'b', ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)
 hold on

load(ref_load_name,'c2_pre','index_pc')
[Stot,S1] = get_sobol_indices(c2_pre,index_pc); 

pause(0.1)

 plot(1:n,Stot, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'g', ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)
ylim([0 1])     
     
L = legend('$\nu$',['$\tau(p=$',num2str(p_ref),')'],'location','NorthWest')
set(L,'interpreter','latex','fontsize',20)

set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'xticklabels',{'$R^{cap}$', '$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
    'TickLabelInterpreter','latex',...
    'fontsize', opts.fontsize)
grid on

% Use for 6D

% Use for 5D
% set(gca, ...
%     'XLim', [1, m], ...
%     'XTick', 1:m, ...
%     'xticklabels',{'$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
%     'TickLabelInterpreter','latex',...
%     'fontsize', opts.fontsize)
% grid on

 %%
% print_name = [file_id, 'PCE_SS_resamp_M',num2str(M)]
% print(print_name,'-dpng','-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%% MAX GAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%df = finite_difference_gradients(X,@MG_PCE);
df = local_linear_gradients(X,mg,30);
sub = compute(df);
opts = plot_opts([]);
opts.fontsize = 20;
%
Y = X*sub.eigenvectors(:, 1:m);
figure()
set(gcf, 'Position', [100, 100, 1024, 1024])
subplot(2,2,3)
% Make 1D sufficient summary plot.
scatter(Y(:,1), mg,'filled', ...
     'markeredgecolor', 'none', ...'markerfacecolor', opts.color, ...
     'cdata', fs, ...
     'marker', opts.marker, ...
     'sizedata', 10*opts.markersize, ...
     'linewidth', opts.linewidth)

xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
box on
colormap('parula')
T = title('Installed Gas MW')
set(T,'interpreter','latex')
xlim([-2,2])

subplot(2,2,4)

% Make 2D sufficient summary plot.
scatter(Y(:,1), Y(:,2), 'filled', ...
        'cdata', mg, ....
        'marker', opts.marker, ...
        'sizedata', 10*opts.markersize, ...
        'linewidth', opts.linewidth)
    
xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
T = title('Installed Gas MW')
set(T,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
cb = colorbar();
colormap('parula')
a =get(cb);
xlim([-2,2])
ylim([-2,2])
box on 

subplot(2,2,2)
e = sub.eigenvalues;
hold on

% Plot eigenvalues.
semilogy(e, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', opts.color, ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)

% Format plot.
title(opts.title, 'fontsize', opts.fontsize)

ylabel('Eigenvalues $\lambda_i$', 'fontsize', opts.fontsize,'interpreter','latex')
xlabel('$i$', 'fontsize', opts.fontsize,'interpreter','latex')


set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'XScale', 'Linear', ...
    'YScale', 'Log', ...
    'fontsize', opts.fontsize)
grid on
box on
subplot(2,2,1)
W = sub.eigenvectors;
e = sub.eigenvalues;
alpha = zeros(m,1);
n = m

for i = 1:m
    alpha(i) = W(i,1:n).^2*e(1:n);
end
alpha = alpha/norm(alpha)

plot(1:m,alpha, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'b', ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)
 hold on
 
load(ref_load_name,'c_max_gas','index_pc') 
[Stot,S1] = get_sobol_indices(c_max_gas,index_pc);
 plot(1:m,Stot, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'g', ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)
ylim([0 1])     
     
L = legend('$\nu$',['$\tau(p=$',num2str(p_ref),')'],'location','NorthWest')
set(L,'interpreter','latex','fontsize',20)
% Use for 6D
set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'xticklabels',{'$R^{cap}$', '$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
    'TickLabelInterpreter','latex',...
    'fontsize', opts.fontsize)
grid on
% Use for 5D
% set(gca, ...
%     'XLim', [1, m], ...
%     'XTick', 1:m, ...
%     'xticklabels',{'$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
%     'TickLabelInterpreter','latex',...
%     'fontsize', opts.fontsize)
% grid on

%%
% print_name = [file_id, 'PCE_max_gas_resamp_M',num2str(M)]
% print(print_name,'-dpng','-r300')

%% %%%%%%%%%%%%%%%%%%%%%%%%%% MAX WIND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%df = finite_difference_gradients(X,@MW_PCE);
df = local_linear_gradients(X,mw,20);
sub = compute(df);
opts = plot_opts([]);
opts.fontsize = 20;
Y = -1*X*sub.eigenvectors(:, 1:m);
figure()
set(gcf, 'Position', [100, 100, 1024, 1024])
subplot(2,2,3)
% Make 1D sufficient summary plot.
scatter(Y(:,1), mw,'filled', ...
     'markeredgecolor', 'none', ...'markerfacecolor', opts.color, ...
     'cdata', mw, ...
     'marker', opts.marker, ...
     'sizedata', 10*opts.markersize, ...
     'linewidth', opts.linewidth)
xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
box on
colormap('parula')
T = title('Installed Wind MW')
set(T,'interpreter','latex')
xlim([-2,2])

subplot(2,2,4)

% Make 2D sufficient summary plot.
scatter(Y(:,1), Y(:,2), 'filled', ...
        'cdata', mw, ....
        'marker', opts.marker, ...
        'sizedata', 10*opts.markersize, ...
        'linewidth', opts.linewidth)

xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
T = title('Installed Wind MW')
set(T,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
cb = colorbar();
colormap('parula')
a =get(cb);
xlim([-2,2])
ylim([-2,2])
box on 

subplot(2,2,2)
e = sub.eigenvalues;
hold on
% Plot eigenvalues.
semilogy(e, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', opts.color, ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)

% Format plot.
title(opts.title, 'fontsize', opts.fontsize)

ylabel('Eigenvalues $\lambda_i$', 'fontsize', opts.fontsize,'interpreter','latex')
xlabel('$i$', 'fontsize', opts.fontsize,'interpreter','latex')


set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'XScale', 'Linear', ...
    'YScale', 'Log', ...
    'fontsize', opts.fontsize)
grid on
box on
subplot(2,2,1)

W = sub.eigenvectors;
e = sub.eigenvalues;
alpha = zeros(m,1);
n = m

for i = 1:m
    alpha(i) = W(i,1:n).^2*e(1:n);
end
alpha = alpha/norm(alpha)

plot(1:m,alpha, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'b', ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)
 hold on
load(ref_load_name,'c_max_wind','index_pc')
[Stot,S1] = get_sobol_indices(c_max_wind,index_pc);
 plot(1:m,Stot, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'g', ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)
ylim([0 1])     
     
L = legend('$\nu$',['$\tau(p=$',num2str(p_ref),')'],'location','NorthWest')
set(L,'interpreter','latex','fontsize',20)
% Use for 6D
set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'xticklabels',{'$R^{cap}$', '$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
    'TickLabelInterpreter','latex',...
    'fontsize', opts.fontsize)
grid on
% Use for 5D
% set(gca, ...
%     'XLim', [1, m], ...
%     'XTick', 1:m, ...
%     'xticklabels',{'$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
%     'TickLabelInterpreter','latex',...
%     'fontsize', opts.fontsize)
% grid on

%%
%print_name = [file_id, 'PCE_max_wind_resamp_M',num2str(M)]
%print(print_name,'-dpng','-r300')

%%
%close all
opts.fontsize = 30;
FF = figure('Position',[1         721        3360         647])
subplot(1,4,1)
%load('Small_FS_cost_scenarios.mat')

%histogram(FS_cost,'NumBins',10,'FaceColor',[250,70,22]/255)
nbins = 50;

load('6D_p4_hist_limits')
% fs_limits = [min(min(fs),min(FS_cost)),max(max(fs),max(FS_cost))]
% ss_limits = [min(min(ss),min(SS_cost)),max(max(ss),max(SS_cost))]
% g_limits = [min(min(mg),min(max_gas)),max(max(mg),max(max_gas))]
% w_limits = [min(min(mw),min(max_wind)),max(max(mw),max(max_wind))]
%  save('6D_p4_hist_limits','fs_limits','ss_limits','g_limits','w_limits')
%w_limits = [6000,9000]





[counts,edges] = histcounts(fs,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts, 1,'FaceColor',[250,70,22]/255);
c_vec = c1_pre;
str = {['$\mu = $',num2str(c_vec(1),'%0.3g')],[ '$ \sigma = $', num2str(sqrt(sum(c_vec(2:end).^2)),'%0.3g')]}
dim = [0.135 0.83 0.0416 0.0836];
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',opts.fontsize);
xlim(fs_limits)
ylim([0 600])
%xticks([0 5 10])

box on
xlabel('Expansion Cost \$', 'fontsize',16,'interpreter','latex')
set(gca,'fontsize',opts.fontsize)


%
subplot(1,4,2)
nbins = 50;

[counts,edges] = histcounts(ss,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts,1,'FaceColor',[0,20,137]/255);
c_vec = c2_pre;
str = {['$\mu = $',num2str(c_vec(1),'%0.3g')],[ '$ \sigma = $', num2str(sqrt(sum(c_vec(2:end).^2)),'%0.3g')]}
dim = [0.42 0.83 0.0416 0.0836];
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',opts.fontsize);
xlim(ss_limits)
ylim([0 500])
%xticks([0 5 10])

box on
xlabel('Operation Cost \$', 'fontsize',16,'interpreter','latex')
set(gca,'fontsize',opts.fontsize)

%
subplot(1,4,3)
nbins = 50;

[counts,edges] = histcounts(mg,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts,1,'FaceColor','k');
c_vec = c_max_gas;
str = {['$\mu = $',num2str(c_vec(1),'%0.3g')],[ '$ \sigma = $', num2str(sqrt(sum(c_vec(2:end).^2)),'%0.3g')]}
dim = [0.65 0.83 0 0.0836];
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',opts.fontsize);
xlim(g_limits)
ylim([0 800])
%xticks([1000 2000 3000 4000 5000])

box on
xlabel('Installed Gas MW', 'fontsize',16,'interpreter','latex')
set(gca,'fontsize',opts.fontsize)

%
subplot(1,4,4)
nbins = 50;

[counts,edges] = histcounts(mw,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts,1,'FaceColor',[.8118    0.7216    0.4863]);
c_vec = c_max_wind;
str = {['$\mu = $',num2str(c_vec(1),'%0.3g')],[ '$ \sigma = $', num2str(sqrt(sum(c_vec(2:end).^2)),'%0.3g')]}
dim = [0.84 0.83 0 0.0836];
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',opts.fontsize);
xlim(w_limits)
ylim([0 700])
%xticks([0 5 10])

box on
xlabel('Installed Wind MW', 'fontsize',16,'interpreter','latex')
set(gca,'fontsize',opts.fontsize)

%%
% 
print_name = [file_id, '_p',num2str(p),'_M',num2str(M), '_PCE_QoIs_resamp']
print(print_name,'-dpng','-r300')
saveas(FF,print_name,'epsc')

%%
%close all

nbins = 50;
F = figure('Position',[1         721        3360         647])
subplot(1,4,1)
%load('Small_FS_cost_scenarios.mat')

%histogram(FS_cost,'NumBins',10,'FaceColor',[250,70,22]/255)
[counts,edges] = histcounts(FS_cost,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts, 1,'FaceColor',[250,70,22]/255);
c_vec = c1_pre;
str = {['$\mu = $',num2str(c_vec(1),'%0.3g')],[ '$ \sigma = $', num2str(sqrt(sum(c_vec(2:end).^2)),'%0.3g')]}
dim = [0.135 0.83 0.0416 0.0836];
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',opts.fontsize);
xlim(fs_limits)
ylim([0 400])
%xticks([0 5 10])

box on
xlabel('Expansion Cost \$', 'fontsize',16,'interpreter','latex')
set(gca,'fontsize',opts.fontsize)


%
subplot(1,4,2)
nbins = 50;

[counts,edges] = histcounts(SS_cost,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts,1,'FaceColor',[0,20,137]/255);
c_vec = c2_pre;
str = {['$\mu = $',num2str(c_vec(1),'%0.3g')],[ '$ \sigma = $', num2str(sqrt(sum(c_vec(2:end).^2)),'%0.3g')]}
dim = [0.42 0.83 0.0416 0.0836];
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',opts.fontsize);
xlim(ss_limits)
%xticks([0 5 10])

box on
xlabel('Operation Cost \$', 'fontsize',16,'interpreter','latex')
set(gca,'fontsize',opts.fontsize)

%
subplot(1,4,3)
nbins = 50;

[counts,edges] = histcounts(max_gas,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts,1,'FaceColor','k');
c_vec = c_max_gas;
str = {['$\mu = $',num2str(c_vec(1),'%0.3g')],[ '$ \sigma = $', num2str(sqrt(sum(c_vec(2:end).^2)),'%0.3g')]}
dim = [0.65 0.83 0 0.0836];
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',opts.fontsize);
xlim(g_limits)
%xticks([1000 2000 3000 4000 5000])

box on
xlabel('Installed Gas MW', 'fontsize',16,'interpreter','latex')
set(gca,'fontsize',opts.fontsize)


%
subplot(1,4,4)
nbins = 50;

[counts,edges] = histcounts(max_wind,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts,1,'FaceColor',[.8118    0.7216    0.4863]);
c_vec = c_max_wind;
str = {['$\mu = $',num2str(c_vec(1),'%0.3g')],[ '$ \sigma = $', num2str(sqrt(sum(c_vec(2:end).^2)),'%0.3g')]}
dim = [0.84 0.83 0 0.0836];
ann = annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',opts.fontsize);
ylim([0 700])
xlim(w_limits)


box on
xlabel('Installed Wind MW', 'fontsize',16,'interpreter','latex')
set(gca,'fontsize',opts.fontsize)

return

%%
% print_name = [file_id, '_SPEED_QoIs']
% print(print_name,'-dpng','-r300')
% saveas(F,print_name,'epsc')




%% %%%%%%%%%%%%%%%%%%%%%%%%% MAX GAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%df = finite_difference_gradients(X,@MG_PCE);
df = local_linear_gradients(X,mg,20);
sub = compute(df);
opts = plot_opts([]);
opts.fontsize = 20;
%
Y = X*sub.eigenvectors(:, 1:m);
figure()
set(gcf, 'Position', [100, 100, 1024, 1024])
subplot(2,2,3)
% Make 1D sufficient summary plot.
scatter(Y(:,1), mg,'filled', ...
     'markeredgecolor', 'none', ...'markerfacecolor', opts.color, ...
     'cdata', fs, ...
     'marker', opts.marker, ...
     'sizedata', 10*opts.markersize, ...
     'linewidth', opts.linewidth)

xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
box on
colormap('parula')
T = title('Installed Gas MW')
set(T,'interpreter','latex')
xlim([-2,2])

subplot(2,2,4)

% Make 2D sufficient summary plot.
scatter(Y(:,1), Y(:,2), 'filled', ...
        'cdata', mg, ....
        'marker', opts.marker, ...
        'sizedata', 10*opts.markersize, ...
        'linewidth', opts.linewidth)
    
xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
T = title('Installed Gas MW')
set(T,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
cb = colorbar();
colormap('parula')
a =get(cb);
xlim([-2,2])
ylim([-2,2])
box on 

subplot(2,2,2)
e = sub.eigenvalues;
hold on

% Plot eigenvalues.
semilogy(e, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', opts.color, ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)

% Format plot.
title(opts.title, 'fontsize', opts.fontsize)

ylabel('Eigenvalues $\lambda_i$', 'fontsize', opts.fontsize,'interpreter','latex')
xlabel('$i$', 'fontsize', opts.fontsize,'interpreter','latex')


set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'XScale', 'Linear', ...
    'YScale', 'Log', ...
    'fontsize', opts.fontsize)
grid on
box on
subplot(2,2,1)
W = sub.eigenvectors;
e = sub.eigenvalues;
alpha = zeros(m,1);
n = m

for i = 1:m
    alpha(i) = W(i,1:n).^2*e(1:n);
end
alpha = alpha/norm(alpha)

plot(1:m,alpha, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'b', ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)
 hold on
 
load(ref_load_name,'c_max_gas','index_pc') 
[Stot,S1] = get_sobol_indices(c_max_gas,index_pc);
 plot(1:m,Stot, ...
         'markeredgecolor', 'k', ...
         'markerfacecolor', 'g', ...
         'color', opts.color, ...
         'marker', opts.marker, ...
         'markersize', 13, ...
         'linewidth', opts.linewidth)
ylim([0 1])     
     
L = legend('$\nu$',['$\tau(p=$',num2str(p_ref),')'],'location','NorthWest')
set(L,'interpreter','latex','fontsize',20)
% Use for 6D
set(gca, ...
    'XLim', [1, m], ...
    'XTick', 1:m, ...
    'xticklabels',{'$R^{cap}$', '$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
    'TickLabelInterpreter','latex',...
    'fontsize', opts.fontsize)
grid on
% Use for 5D
% set(gca, ...
%     'XLim', [1, m], ...
%     'XTick', 1:m, ...
%     'xticklabels',{'$c^{loss}$' , '$c^{oload}$', '$c^{ng}$', '$c^{wind}$' ,'$t^{avail}$'},...
%     'TickLabelInterpreter','latex',...
%     'fontsize', opts.fontsize)
% grid on

%%

figure
%f = fs;  title_string =  'Expansion Cost \$'; c = c1_pre; c_name = 'c1_pre'; print_name = ['FS_3D_', file_id, '_p',num2str(p) ]
%f = ss;  title_string =  'Operations Cost \$'; c = c2_pre; c_name = 'c2_pre';  print_name = ['SS_3D_', file_id, '_p',num2str(p) ]
f = mg;  title_string =  'Max Installed Gas MW'; c = c_max_gas; c_name = 'c_max_gas'; print_name = ['max_gas_3D_', file_id, '_p',num2str(p) ]
%f = mw;  title_string =  'Max Installed Wind MW'; c = c_max_wind; c_name = 'c_max_wind';  print_name = ['max_wind_3D_', file_id, '_p',num2str(p) ]


df = local_linear_gradients(X,f,20)
sub = compute(df);
opts = plot_opts([]);
opts.fontsize = 20;
Y = X*sub.eigenvectors(:, 1:m);


%close all
% Make 3D sufficient summary plot.
scatter3(Y(:,1), Y(:,2), Y(:,3), 'filled', ...
        'cdata', f, ....
        'marker', opts.marker, ...
        'sizedata', 10*opts.markersize, ...
        'linewidth', opts.linewidth)
    
    
    
xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
T = title(title_string)
set(T,'interpreter','latex')
set(gca, 'fontsize', opts.fontsize)
grid('on')
cb = colorbar();
colormap('parula')
a =get(cb);
xlim([-2,2])
ylim([-2,2])
view(-37.5,30);
%%
% close all
% %% Structured domain
% 
% Nshp = 100;
% Nlvl = 20;
% [Y1,Y2,Y3] = meshgrid(linspace(-1,1,Nshp),...
%                       linspace(-1,1,Nshp),...
%                       linspace(-1,1,Nshp));
% y = [reshape(Y1,Nshp^3,1), reshape(Y2,Nshp^3,1), reshape(Y3,Nshp^3,1)];
% % evaluate along ridge (or 3d surrogate)
%  XX = y*W(:,1:3)';  
%  
%  %f = fs; F = FS_PCE(XX,c1_pre,index_pc); title_string =  'Expansion Cost \$'; c = c1_pre; c_name = 'c1_pre'; print_name = ['FS_3D_', file_id, '_p',num2str(p) ]
% f = ss; F = SS_PCE(XX,c2_pre,index_pc); title_string =  'Operations Cost \$'; c = c2_pre; c_name = 'c2_pre';  print_name = ['SS_3D_', file_id, '_p',num2str(p) ]
% %f= mg; F = MG_PCE(XX,c_max_gas,index_pc); title_string =  'Max Installed Gas MW'; c = c_max_gas; c_name = 'c_max_gas'; print_name = ['max_gas_3D_', file_id, '_p',num2str(p) ]
% %f = mw; F = MW_PCE(XX,c_max_wind,index_pc); title_string =  'Max Installed Wind MW'; c = c_max_wind; c_name = 'c_max_wind';  print_name = ['max_wind_3D_', file_id, '_p',num2str(p) ]
% 
% 
%  
% % start by partitioning the range of F
% maxF = max(F);
% minF = min(F);
% % assign appropriate grid manually
% YY(:,:,:,1) = Y1;
% YY(:,:,:,2) = Y2;
% YY(:,:,:,3) = Y3;
% % build uniform level-set values at 90% of min and max
% lvlsets = linspace(0.9*minF,0.9*maxF,Nlvl);
% 
% FV = cell(Nlvl,1); figure; 
% for i=1:Nlvl
%     fprintf('Building isosurface %i of %i...',i,Nlvl);
%     FV{i} = isosurface(YY(:,:,:,1), YY(:,:,:,2), YY(:,:,:,3),reshape(F(1:Nshp^3),Nshp,Nshp,Nshp),lvlsets(i));
%     Ptch = patch(FV{i},'FaceVertexCData',lvlsets(i)*ones(size(FV{i}.vertices,1),1),'FaceColor','flat'); 
%     colorbar; grid on;
%     Ptch.FaceAlpha = 0.5; Ptch.EdgeColor = 'none';
%     clc;
% end
% 
% xlabel('Active Variable 1 ($\mathbf{w}_1^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
% ylabel('Active Variable 2 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
% zlabel('Active Variable 3 ($\mathbf{w}_2^T\mathbf{\xi}$)', 'fontsize', opts.fontsize,'interpreter','latex')
% T = title(title_string)
% set(T,'interpreter','latex')
% set(gca, 'fontsize', opts.fontsize)
% grid('on'); view(-37.5,30);
%%

% F = mg; Y = X*W;
% % plot level sets
% figure; IsoShadow(Y,F,Nlvl); view(-37.5,30);
% title 'Random data'