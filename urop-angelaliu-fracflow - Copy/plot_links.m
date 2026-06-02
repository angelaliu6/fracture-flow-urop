function plot_links(figno,nlink,linknods,posnode)
%
% plot links of fracture network

figure(figno);
for ilink = 1:nlink
    plot(posnode(linknods(ilink,:),1),posnode(linknods(ilink,:),2),...
        'k-','LineWidth',1.5)
end
return
