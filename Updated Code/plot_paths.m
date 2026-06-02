function plot_paths(figno,Lx,Ly,npaths,posnode)
%
% Plot particle paths
figure(figno)
plot([0 Lx Lx 0 0],[0 0 Ly Ly 0],'k-')
hold on
axis equal
axis tight
axis([0 Lx 0 Ly])
%
Np = length(npaths);
alpha_min = 0.01;
alpha_value = alpha_min + (1-alpha_min)*1/Np;
maxlinewidth = 3;
for ip = 1:Np
    plot(posnode(npaths{ip},1),posnode(npaths{ip},2),...
            'color',[1 0 0 alpha_value],'LineWidth',maxlinewidth)
end

return
