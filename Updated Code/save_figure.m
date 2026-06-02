function save_figure(figno,figname)

% save figure
figure(figno)
set(gcf,'PaperPosition',[1 1 6 3])
print(gcf,'-dpng','-r600',[figname,'.png'])
% print(gcf,'-depsc','-loose',[figname,'.eps'])
