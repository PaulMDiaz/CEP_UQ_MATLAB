clear
clc
close all
addpath(genpath('./'))
file_id = '6D_UC'

load_name = ['CEP_test_run_', file_id, '_output.mat']


load(load_name)
%load('Small_FS_cost_scenarios.mat')
%%
[N,d] = size(Xi_CEP_test_run)
p  = 2
save_name = ['CEP_PCE_', file_id, '_p',num2str(p), '.mat']

index_pc = nD_polynomial_array(d,p);
P = size(index_pc,1)

Psi = zeros(N,P);
for i = 1:size(Xi_CEP_test_run)
   Psi(i,:) = piset(Xi_CEP_test_run(i,:),index_pc);
    
end
Lambda_cands = 0.5*10.^(linspace(-2,0,3)) %0.02:0.02:1 %10.^(-linspace(6,0,30))
plotting = false
%%

%c1_pre = Psi\FS_cost
%c2_pre = Psi\SS_cost
% 
% c1 = SINDy(Psi,FS_cost,'SPGL1');
% c2 = SINDy(Psi,SS_cost,'SPGL1');
% norm(Psi*c1-FS_cost)/norm(FS_cost)
% norm(Psi*c2-SS_cost)/norm(SS_cost)

%
%load(save_name)
%%

[Xi_WSPGL1_L, Xi_WSP_L, Xi_RR_L, Xi_IT_L, Sys_info] = L_Adapt_SINDy(Psi,[FS_cost,SS_cost,max_wind,avg_wind,max_gas],Lambda_cands,plotting)


Xi = Xi_WSPGL1_L; 

c1_pre = Xi(:,1); 
c2_pre = Xi(:,2); %Xi(:,2);
c_max_wind = Xi(:,3);       
c_avg_wind = Xi(:,4);
c_max_gas = Xi(:,5);
%
%save(save_name,'c1_pre','c2_pre','c_max_wind','c_avg_wind','c_max_gas','index_pc','Xi_WSPGL1_L', 'Xi_WSP_L', 'Xi_RR_L', 'Xi_IT_L', 'Sys_info')

%save(save_name,'c1_pre','c2_pre','c_max_wind','c_avg_wind','c_max_gas','index_pc','Xi','Sys_info')


%%
FS_rel_err = norm(Psi*c1_pre-FS_cost)/norm(FS_cost)
FS_l0_norm = length(find(c1_pre ~= 0))
FS_l1_norm = sum(abs(c1_pre))

SS_rel_err = norm(Psi*c2_pre-SS_cost)/norm(SS_cost)
SS_l0_norm = length(find(c2_pre ~= 0))
SS_l1_norm = sum(abs(c2_pre))

MW_rel_err = norm(Psi*c_max_wind-max_wind)/norm(max_wind)
MW_l0_norm = length(find(c_max_wind ~= 0))
MW_l1_norm = sum(abs(c_max_wind))
%norm(Psi*c_avg_wind-avg_wind)/norm(avg_wind)
MG_rel_err = norm(Psi*c_max_gas-max_gas)/norm(max_gas)
MG_l1_norm = sum(abs(c_max_gas))
MG_l0_norm = length(find(c_max_gas ~= 0))







%%
%close all
%load('CEP_PCE.mat')

[Tau_FS,S_1_FS] = get_sobol_indices(c1_pre,index_pc)

[Tau_SS, S_1_SS] = get_sobol_indices(c2_pre,index_pc)


[Tau_max_wind, S_1_max_wind] = get_sobol_indices(c_max_wind,index_pc)
[Tau_avg_wind, S_1_avg_wind] = get_sobol_indices(c_avg_wind,index_pc)
[Tau_max_gas, S_1_max_gas] = get_sobol_indices(c_max_gas,index_pc)



      
%%

%weights = get_matrix_weights(Psi);
%Psiweights = Psi*weights;
%sigma =  cross_val_sigma(Psiweights,FS_cost)
%opts = spgSetParms('iterations',50000*size(Psi,2),'verbosity',0); %20*size(D,2)
%c1_pre = weights*spg_bpdn(Psiweights,FS_cost,sigma*norm(FS_cost),opts);

%sigma =  cross_val_sigma(Psiweights,SS_cost)
%c2_pre = weights*spg_bpdn(Psiweights,SS_cost,sigma*norm(SS_cost),opts);



%c1 = Psi\FS_cost
%c2 = Psi\SS_cost

%%
%close all
figure('Position', [283,138,1019,800])
subplot(2,2,1)
%labels = {'$R^{cap}$', '$c^{loss}$' , '$c^{oload}$', '$c^{ng}$','$c^{wind}$','$line_{ex}$'};
title_string = { 'Expansion Cost', ['Sobol Indices  M = ' , num2str(N) ]}
%

labels = {'$R^{cap}$', '$c^{loss}$' , '$c^{oload}$', '$c^{ng}$','$c^{wind}$','$t^{avail}$'};
%labels = {'$c^{loss}$' , '$c^{oload}$', '$c^{ng}$','$c^{wind}$','$t^{avail}$'};
gen_sobol_plot(Tau_FS,S_1_FS, labels, title_string)

%
%
%gen_pie_chart(Tau_FS/norm(Tau_FS),labels,title_string)
%
%print_name = [file_id, '_FS']
%print(print_name,'-dpng','-r300')

% %%
% close all
% title_string = { 'Expansion Cost', ['1$^{st}$ Order Sobol Indices    M = ' , num2str(N) ]}
% gen_pie_chart(S_1_FS,labels,title_string)
% %%
% print_name = [file_id, '_S1_FS']
% print(print_name,'-dpng','-r300')

%
%close all
subplot(2,2,2)
title_string = { 'Operation Cost', ['Sobol Indices    M = ' , num2str(N) ]};
gen_sobol_plot(Tau_SS,S_1_SS, labels, title_string)
%
%print_name = [file_id, '_SS']
%print(print_name,'-dpng','-r300')

%
%close all
subplot(2,2,3)
title_string = { 'Max Wind', [' Sobol Indices    M = ' , num2str(N) ]};
gen_sobol_plot(Tau_max_wind,S_1_max_wind, labels, title_string)

%
%print_name = [file_id, '_max_wind']
%print(print_name,'-dpng','-r300')



%
%close all
%title_string = { 'Avg Wind', ['Sobol Indices    M = ' , num2str(N) ]}
%gen_sobol_plot(Tau_avg_wind,S_1_avg_wind, labels, title_string)
%
%print_name = [file_id, '_avg_wind']
%print(print_name,'-dpng','-r300')
%
%close all
subplot(2,2,4)
title_string = { 'Max Installed Gas', ['Sobol Indices    M = ' , num2str(N) ]}
gen_sobol_plot(Tau_max_gas,S_1_max_gas, labels, title_string)

%%
% print_name = [file_id, '_p',num2str(p) '_Sobol_Inds']
% print(print_name,'-dpng','-r300')




%%
% hold on
% close all
% %load('Small_FS_cost_scenarios.mat')
% 
% %histogram(FS_cost,'NumBins',10,'FaceColor',[250,70,22]/255)
% nbins = 50;
% 
% [counts,edges] = histcounts(SS_cost,nbins);
% center = 0.5*(edges(1:end-1)+edges(2:end));
% %bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
% bar(center, counts,1,'FaceColor',[0,20,137]/255);
% %xlim([1.3, 2.1].*10^14)
% %xticks([0 5 10])
% 
% set(gca,'fontsize',20)
% box on
% xlabel('Operation Cost $u(\mathbf{\xi})$ in dollars', 'fontsize',20,'interpreter','latex')
% %%
% print('SS_cost_hist','-dpng','-r300')
% 
% %%
% 
% close all
% %load('Small_FS_cost_scenarios.mat')
% 
% %histogram(FS_cost,'NumBins',10,'FaceColor',[250,70,22]/255)
% nbins = 50;
% 
% [counts,edges] = histcounts(FS_cost,nbins);
% center = 0.5*(edges(1:end-1)+edges(2:end));
% %bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
% bar(center, counts, 1,'FaceColor',[250,70,22]/255);
% %xlim([1.3, 2.1].*10^14)
% %xticks([0 5 10])
% 
% set(gca,'fontsize',20)
% box on
% xlabel('Expansion Cost $u(\mathbf{\xi})$ in dollars', 'fontsize',20,'interpreter','latex')
% 
% %%
% %print('FS_cost_hist','-dpng','-r300')
% 



%%
% %close all
% figure
% subplot(1,4,1)
% %load('Small_FS_cost_scenarios.mat')
% 
% %histogram(FS_cost,'NumBins',10,'FaceColor',[250,70,22]/255)
% nbins = 50;
% 
% [counts,edges] = histcounts(FS_cost,nbins);
% center = 0.5*(edges(1:end-1)+edges(2:end));
% %bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
% bar(center, counts, 1,'FaceColor',[250,70,22]/255);
% %xlim([4, 7].*10^9)
% %xticks([0 5 10])
% 
% set(gca,'fontsize',20)
% box on
% xlabel('Expansion Cost \$', 'fontsize',16,'interpreter','latex')
% 
% 
% %
% subplot(1,4,2)
% nbins = 50;
% 
% [counts,edges] = histcounts(SS_cost,nbins);
% center = 0.5*(edges(1:end-1)+edges(2:end));
% %bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
% bar(center, counts,1,'FaceColor',[0,20,137]/255);
% %xlim([1, 1.6].*10^8)
% %xticks([0 5 10])
% 
% set(gca,'fontsize',20)
% box on
% xlabel('Operation Cost \$', 'fontsize',16,'interpreter','latex')
% 
% %
% subplot(1,4,3)
% nbins = 50;
% 
% [counts,edges] = histcounts(max_gas,nbins);
% center = 0.5*(edges(1:end-1)+edges(2:end));
% %bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
% bar(center, counts,1,'FaceColor','b');
% %xlim([1800, 4500])
% %xticks([0 5 10])
% 
% set(gca,'fontsize',20)
% box on
% xlabel('Installed Gas MW', 'fontsize',16,'interpreter','latex')
% 
% 
% %
% subplot(1,4,4)
% nbins = 50;
% 
% [counts,edges] = histcounts(max_wind,nbins);
% center = 0.5*(edges(1:end-1)+edges(2:end));
% %bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
% bar(center, counts,1,'FaceColor','r');
% %xlim([0,4000])
% %xticks([0 5 10])
% 
% set(gca,'fontsize',20)
% box on
% xlabel('Installed Wind MW', 'fontsize',16,'interpreter','latex')

%%
% print_name = [file_id, 'SPEED_QoIs']
% print(print_name,'-dpng','-r300')
