function [Adjmat,backbonenods,isolatednods,inletnods,outletnods] = ...
    backbone(linknods,linktrans,posnode,nnode,Lx)

% Determine backbone of the network

% build weighted adjacency matrix as a sparse matrix
% (transmissibility matrix without the diagonal)
Adjmat = sparse([linknods(:,1);linknods(:,2)],...
           [linknods(:,2);linknods(:,1)],...
           [linktrans;linktrans]);

% build weighted graph
G = graph(linknods(:,1), linknods(:,2), linktrans, nnode);

% compute shortestpath between each node and the boundary nodes
% if shortestpath = Inf, node is isolated and does not belong to backbone
%
% inlet
inletnods = find(posnode(:,1) == 0);
inletdist = Inf*ones(nnode,length(inletnods));
for iinlet = 1:length(inletnods)
%     inletdist(:,iinlet) = graphshortestpath(Adjmat,inletnods(iinlet));
    inletdist(:,iinlet) = distances(G,inletnods(iinlet))';
end
inletdist = min(inletdist,[],2);
%
% outlet
outletnods = find(posnode(:,1) == Lx);
outletdist = Inf*ones(nnode,length(outletnods));
for ioutlet = 1:length(outletnods)
%     outletdist(:,ioutlet) = graphshortestpath(Adjmat,outletnods(ioutlet));
    outletdist(:,ioutlet) = distances(G,outletnods(ioutlet))';
end
outletdist = min(outletdist,[],2);
%
bounddist = inletdist + outletdist;
%
% backbone nodes: those connected to inlet and outlet
backbonenods = find(bounddist < Inf);
isolatednods = find(bounddist == Inf);

return

