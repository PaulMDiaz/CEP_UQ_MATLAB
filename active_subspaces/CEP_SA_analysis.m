clear
clc
addpath(genpath('./'))
load('CEP_test_run_output.mat')
%load('Small_FS_cost_scenarios.mat')

[N,d] = size(Xi_CEP_test_run)
p  = 2


index_pc = nD_polynomial_array(d,p)
P = size(index_pc,1)

Psi = zeros(N,P);
for i = 1:size(Xi_CEP_test_run)
   Psi(i,:) = piset(Xi_CEP_test_run(i,:),index_pc);
    
end
Lambda_cands = linspace(0.1,.35,5)

%%

%c1_pre = Psi\FS_cost
%c2_pre = Psi\SS_cost
% 
% c1 = SINDy(Psi,FS_cost,'SPGL1');
% c2 = SINDy(Psi,SS_cost,'SPGL1');
% norm(Psi*c1-FS_cost)/norm(FS_cost)
% norm(Psi*c2-SS_cost)/norm(SS_cost)

%

%load('CEP_PCE_p6.mat')
%%
[Xi_WSPGL1_L, ~, ~, ~, Sys_info] = adapt_SINDy(Psi,[FS_cost,SS_cost,max_wind,avg_wind,max_gas],Lambda_cands)
                
 
c1_pre = Xi_WSPGL1_L(:,1),
c2_pre = Xi_WSPGL1_L(:,2),
c_max_wind = Xi_WSPGL1_L(:,3),
c_avg_wind = Xi_WSPGL1_L(:,4),
c_max_gas = Xi_WSPGL1_L(:,5),

norm(Psi*c1_pre-FS_cost)/norm(FS_cost)
norm(Psi*c2_pre-SS_cost)/norm(SS_cost)
norm(Psi*c_max_wind-max_wind)/norm(max_wind)
norm(Psi*c_avg_wind-avg_wind)/norm(avg_wind)
norm(Psi*c_max_gas-max_gas)/norm(max_gas)





%%



%c1 = SP(Psi,0,FS_cost,true,'fold');
%c2 = SP(Psi,0,SS_cost,true,'fold');
%%
%save('CEP_PCE.mat','c1_pre','c2_pre','c_max_wind','c_avg_wind','c_max_gas','index_pc')

      
%%

%weights = get_matrix_weights(Psi);
%Psiweights = Psi*weights;
%sigma =  cross_val_sigma(Psiweights,FS_cost)
%opts = spgSetParms('iterations',50000*size(Psi,2),'verbosity',0); %20*size(D,2)
%c1_pre = weights*spg_bpdn(Psiweights,FS_cost,sigma*norm(FS_cost),opts);

%sigma =  cross_val_sigma(Psiweights,SS_cost)
%c2_pre = weights*spg_bpdn(Psiweights,SS_cost,sigma*norm(SS_cost),opts);

%%
close all
%load('CEP_PCE.mat')

[Tau_FS,S_1_FS] = get_sobol_indices(c1_pre,index_pc)

[Tau_SS, S_1_SS] = get_sobol_indices(c2_pre,index_pc)


[Tau_max_wind, S_1_max_wind] = get_sobol_indices(c_max_wind,index_pc)
[Tau_avg_wind, S_1_avg_wind] = get_sobol_indices(c_avg_wind,index_pc)
[Tau_max_gas, S_1_max_gas] = get_sobol_indices(c_max_gas,index_pc)

%c1 = Psi\FS_cost
%c2 = Psi\SS_cost

%%

close all
labels = {'$R^{cap}$', '$c^{loss}$' , '$c^{oload}$', '$c^{ng}$','$c^{wind}$'};
title_string = { 'First Stage Cost', ['Total Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(Tau_FS/norm(Tau_FS),labels,title_string)
%%
print('Tau_FS','-dpng','-r300')

%%
close all
title_string = { 'First Stage Cost', ['1$^{st}$ Order Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(S_1_FS,labels,title_string)
%%
print('S1_FS','-dpng','-r300')


%%
close all
title_string = { 'Second Stage Cost', ['Total Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(Tau_SS,labels,title_string)
%%
print('Tau_SS','-dpng','-r300')

%%
close all
title_string = { 'Second Stage Cost', ['1$^{st}$ Order Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(S_1_SS,labels,title_string)
%%
print('S1_SS','-dpng','-r300')

%%
close all
title_string = { 'Max Installed Wind', ['Total Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(Tau_max_wind,labels,title_string)
%%
print('Tau_max_wind','-dpng','-r300')

%%
close all
title_string = { 'Max Wind', ['1$^{st}$ Order Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(S_1_max_wind,labels,title_string)
%%
print('S1_max_wind','-dpng','-r300')

%%
close all
title_string = { 'Avg Wind', ['Total Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(Tau_avg_wind,labels,title_string) 

%%
print('Tau_avg_wind','-dpng','-r300')

%%
close all
title_string = { 'Avg Wind', ['1$^{st}$ Order Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(S_1_avg_wind,labels,title_string)
%%
print('S1_avg_wind','-dpng','-r300')
%%
close all
title_string = { 'Max Installed Gas', ['Total Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(Tau_max_gas,labels,title_string)
%%
print('Tau_max_gas','-dpng','-r300')

%%
close all
title_string = { 'Max Installed Gas', ['1$^{st}$ Order Sobol Indices    M = ' , num2str(N) ]}
gen_pie_chart(S_1_max_gas,labels,title_string)
%%
print('S1_max_gas','-dpng','-r300')



%%
hold on
close all
%load('Small_FS_cost_scenarios.mat')

%histogram(FS_cost,'NumBins',10,'FaceColor',[250,70,22]/255)
nbins = 50;

[counts,edges] = histcounts(SS_cost,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts,1,'FaceColor',[0,20,137]/255);
%xlim([1.3, 2.1].*10^14)
%xticks([0 5 10])

set(gca,'fontsize',20)
box on
xlabel('Operation Cost $u(\mathbf{\xi})$ in dollars', 'fontsize',20,'interpreter','latex')
%%
print('SS_cost_hist','-dpng','-r300')

%%

close all
%load('Small_FS_cost_scenarios.mat')

%histogram(FS_cost,'NumBins',10,'FaceColor',[250,70,22]/255)
nbins = 50;

[counts,edges] = histcounts(FS_cost,nbins);
center = 0.5*(edges(1:end-1)+edges(2:end));
%bar(center, counts, 0.5*10000000000000,'FaceColor',[250,70,22]/255);
bar(center, counts, 1,'FaceColor',[250,70,22]/255);
%xlim([1.3, 2.1].*10^14)
%xticks([0 5 10])

set(gca,'fontsize',20)
box on
xlabel('Expansion Cost $u(\mathbf{\xi})$ in dollars', 'fontsize',20,'interpreter','latex')

%%
%print('FS_cost_hist','-dpng','-r300')



