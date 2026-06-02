function plot_fractures(figno,Lx,Ly,xfini,xfend,yfini,yfend)
%
% plot fractures

figure(figno)
plot([0 Lx Lx 0 0],[0 0 Ly Ly 0],'k-')
hold on
plot([xfini(:,1),xfend(:,1)]',[xfini(:,2),xfend(:,2)]','m-','LineWidth',1)
plot([yfini(:,1),yfend(:,1)]',[yfini(:,2),yfend(:,2)]','c-','LineWidth',1)
%
axis equal
axis tight
axis([0 Lx 0 Ly])

return