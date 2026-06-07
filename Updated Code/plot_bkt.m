function plot_bkt(figno,tpaths)
%
% Plot breakthrough times
figure(figno)
%
Np = length(tpaths);
t_bkt = [];
for ip = 1:Np
    if ~isempty(tpaths{ip})
        t_bkt(end+1) = tpaths{ip}(end);
    end
end
histogram(log10(t_bkt))
ax = gca;
ax.XTick = unique(round(ax.XTick));
ax.XTickLabel = 10.^(ax.XTick);
xlabel('Breakthrough time (s)')
