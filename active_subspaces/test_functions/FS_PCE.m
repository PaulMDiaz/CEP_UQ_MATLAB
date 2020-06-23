function f = FS_PCE(Xi,c1_pre, index_pc)


psi = zeros(size(Xi,1),size(index_pc,1));
for i = 1:size(Xi,1)
    psi(i,:) = piset(Xi(i,:),index_pc);
end

f = psi*c1_pre;


end

