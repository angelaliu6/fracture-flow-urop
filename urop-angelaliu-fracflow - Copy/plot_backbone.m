function plot_backbone(figno,backbonenods,posnode)
%
% plot backbone nodes

figure(figno)
plot(posnode(backbonenods,1),posnode(backbonenods,2),'g*','MarkerSize',5)

return
