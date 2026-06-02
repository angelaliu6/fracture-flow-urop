function save_figure(figno, foldername, figname)
    save figure
    figure(figno)
    set(figno,'PaperPosition',[1 1 6 3])
    % print(figno,'-dpdf','-loose',[figname,'.pdf'])
    print(figno,'-dpdf','-loose',fullfile(foldername, [figname,'.pdf'])) 