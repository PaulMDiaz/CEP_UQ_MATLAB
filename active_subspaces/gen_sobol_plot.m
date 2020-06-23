function gen_sobol_plot(Stot,S1, labels, title_string)

%figure
m = length(Stot);
fontsize = 20;

w1 = 0.5; w2 = 0.25;
bar(1:m,Stot,w1,'FaceColor',[250,70,22]/255);
hold on
bar(1:m,S1,w2,'FaceColor',[0,20,137]/255);


L = legend('Total Index $\tau$','$1^{st}$ Order Index $S_1$','location','best');
set(L,'interpreter','latex','fontsize',fontsize)
set(gca, ...
    'XTick', 1:m, ...
    'xticklabels',labels,...
    'TickLabelInterpreter','latex',...
    'fontsize', fontsize,...
    'ylim',[0,1]);
grid on
title(title_string,'interpreter','latex')
end

