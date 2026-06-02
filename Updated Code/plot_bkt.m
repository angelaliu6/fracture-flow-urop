function plot_bkt(figno,tpaths)
%
% Plot breakthrough times
figure(figno)
%
Np = length(tpaths);
t_bkt = zeros(Np,1);
for ip = 1:Np
    t_bkt(ip) = tpaths{ip}(end);
end
histogram(log10(t_bkt))
ax = gca;
ax.XTick = unique(round(ax.XTick));
ax.XTickLabel = 10.^(ax.XTick);
xlabel('Breakthrough time (s)')
