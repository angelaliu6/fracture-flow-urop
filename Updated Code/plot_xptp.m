function plot_xptp(figno,npaths,tpaths,posnode)
%
% Plot particle space-time plot
figure(figno)
%
Np = length(npaths);
for ip = 1:Np
    np = npaths{ip};
    xp = posnode(np,1);
    tp = tpaths{ip};
    stairs(xp,tp,'LineWidth',0.25)
    semilogy(xp,tp,'o','MarkerSize',3,...
        'MarkerEdgeColor','k','MarkerFaceColor','w')
    hold on
end
xlabel('x (m)')
ylabel('t (s)')

return
