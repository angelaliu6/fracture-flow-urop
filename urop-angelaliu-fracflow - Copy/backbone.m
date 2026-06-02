function [Adjmat,backbonenods,isolatednods,inletnods,outletnods] = ...
    backbone(linknods,linktrans,posnode,nnode,Lx)

% Determine backbone of the network

% build weighted adjacency matrix as a sparse matrix
% (transmissibility matrix without the diagonal)
try
    Adjmat = sparse([linknods(:,1);linknods(:,2)],...
           [linknods(:,2);linknods(:,1)],...
           [linktrans;linktrans]);

catch
a=0
end

% build weighted graph
G = graph(linknods(:,1), linknods(:,2), linktrans);
% highest node number in graph
mnode = max(unique(linknods));

% compute shortestpath between each node and the boundary nodes
% if shortestpath = Inf, node is isolated and does not belong to backbone

% inlet
inletnods = find(posnode(:,1) == 0);
inletdist = Inf*ones(nnode,length(inletnods)); % why is this nnode instead of mnode?
% disp(size(distances(G)))

% I believe the distances(G,inletnods(iinlet)) needs to be transposed, and
% that inletdist's number-of-rows should be set to the number of vertices in G instead of being left as :.
for iinlet = 1:length(inletnods)
    inletdist(1:mnode,iinlet) = distances(G,inletnods(iinlet)); % the matrix has d(i,j) = distance between nodes inletnods(iinlet), j
end
inletdist = min(inletdist,[],2);
%
% outlet
outletnods = find(posnode(:,1) == Lx);
outletdist = Inf*ones(nnode,length(outletnods));
for ioutlet = 1:length(outletnods)
    outletdist(1:mnode,ioutlet) = distances(G,outletnods(ioutlet));
end
outletdist = min(outletdist,[],2);
%
bounddist = inletdist + outletdist;
%
% backbone nodes: those connected to inlet and outlet
backbonenods = find(bounddist < Inf);
isolatednods = find(bounddist == Inf);

return

