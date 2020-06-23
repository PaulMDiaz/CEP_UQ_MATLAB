function [Tau, S_1] = get_sobol_indices(c, index_pc)

%c = Psi\U;
var_f = 0;

for j = 2:size(index_pc,1)
    var_f = var_f+  c(j)^2; %*mean(Psi(:,j).^2); %var_f + (U(j)- mu_f)^2;
end
%var_f = var_f/size(index_pc,2);

Tau = zeros(size(index_pc,2),1);
S_1 = Tau;
for d = 1:size(index_pc,2)
   coeff_inds = find(index_pc(:,d) ~= 0);
   for j = 1:size(coeff_inds)
        Tau(d) = Tau(d) + c(coeff_inds(j))^2;
   end
   
end

D = size(index_pc,2);
for d = 1:D
    %coeff_inds = []
    for n = 1:size(index_pc,1)
       if index_pc(n,d) ~= 0 && sum(index_pc(n,setdiff(1:D,d))) == 0
            S_1(d) = S_1(d) + c(n)^2;        
            %coeff_inds = [coeff_inds n];
       end
    end
    
    %coeff_inds
end

Tau = Tau./var_f;
S_1 = S_1./var_f;

end

