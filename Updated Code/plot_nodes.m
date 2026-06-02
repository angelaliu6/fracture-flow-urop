function plot_nodes(figno,posnode)
%
% plot fracture intersections (nodes)

figure(figno)
plot(posnode(:,1),posnode(:,2),'ko','MarkerSize',4)

return
