function plot_flow(figno,Lx,Ly,linkflow,linknods,backbonenods,posnode,pvec,tol)
%
% Plot pressure at nodes and flowrate along links

figure(figno)
plot([0 Lx Lx 0 0],[0 0 Ly Ly 0],'k-')
hold on
axis equal
axis tight
axis([0 Lx 0 Ly])
%
maxlinewidth = 4;
linkwidth = abs(linkflow)/max(linkflow)*maxlinewidth;
for ilink = 1:length(linknods)
    if abs(linkflow(ilink)) > tol
        plot(posnode(linknods(ilink,:),1),posnode(linknods(ilink,:),2),...
            'k-','LineWidth',linkwidth(ilink))
    end
end
%
scatter(posnode(backbonenods,1),posnode(backbonenods,2),40,pvec,'filled')
colorbar

return
