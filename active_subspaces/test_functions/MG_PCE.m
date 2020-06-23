function f = MG_PCE(Xi, c_max_gas,index_pc)
psi = zeros(size(Xi,1),size(index_pc,1));
for i = 1:size(Xi,1)
    psi(i,:) = piset(Xi(i,:),index_pc);
end

f = psi*c_max_gas;
inds = find(f < 0);
f(inds) = zeros;


end

