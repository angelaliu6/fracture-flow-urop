function save_figure(figno, foldername, figname)
    fig = figure(figno);
    set(fig,'PaperPosition',[1 1 6 3])
    print(fig,'-dpdf','-loose',fullfile(foldername, [char(figname),'.pdf']))