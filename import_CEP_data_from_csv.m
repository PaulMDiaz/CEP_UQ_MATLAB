clear 
clc

filename = '6D_UC_output.csv'

d = 6;

 f = fopen (filename);
 C = textscan(f, '%n %n %n %n %n %n %n %n %n %n %n %n %s ', 'Delimiter', ',', 'HeaderLines', 1);
 
 
 for i = 1:d
    Xi_CEP_test_run(:,i) = C{i+1};
 end
 
 FS_cost = C{8};
 SS_cost = C{9};
 max_wind = C{10};
 avg_wind = C{11};
 max_gas = C{12};
 fclose(f);
 
 good_inds = find(~isnan(Xi_CEP_test_run(:,1)));
 Xi_CEP_test_run = Xi_CEP_test_run(good_inds,:);
 FS_cost = FS_cost(good_inds);
SS_cost = SS_cost(good_inds);
max_wind = max_wind(good_inds);
max_gas = max_gas(good_inds);
avg_wind = avg_wind(good_inds);



 
 %%
clear C d f i filename ans
save('CEP_test_run_6D_UC_output.mat')