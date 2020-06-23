function f = SS_PCE(Xi,c2_pre,index_pc)

psi = zeros(size(Xi,1),size(index_pc,1));
for i = 1:size(Xi,1)
    psi(i,:) = piset(Xi(i,:),index_pc);
end

f = psi*c2_pre;


end

