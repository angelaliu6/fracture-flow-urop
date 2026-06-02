function [h] = create_plot_and_histogram(fulldata, subdata_length, color_c, color_h, translucencyConstant, numBins)
    % [f,xi] = ksdensity(fulldata(1:subdata_length));
    % c = plot(xi, f, color_c);
    [~, bin_edges] = histcounts(fulldata(1:subdata_length), numBins);
    bin_edges = sort(unique([bin_edges, 0, 0]));
    h = histogram(fulldata(1:subdata_length), bin_edges, 'Normalization','pdf', 'FaceColor', color_h, 'FaceAlpha', translucencyConstant);
    hold on
