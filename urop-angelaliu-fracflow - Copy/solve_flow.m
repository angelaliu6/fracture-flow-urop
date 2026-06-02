 % conninletnodes ("connected inlet nodes") is the list of indices of backbone inlet nodes, WRT the full list of nodes
 % pvec ("pressure vector")[n] is the pressure at the nth backbone node
 % precondition: backbonenods is nonempty
function [conninletnods,connoutletnods,pvec] = ...
    solve_flow(Adjmat,backbonenods,isolatednods,inletnods,outletnods,pinlet,poutlet)

% Solve flow (node pressures)
%
% transmissibility matrix
Tmat = Adjmat(backbonenods,backbonenods);
Tdiag = full(sum(Tmat,2));
Tmat = spdiags(-Tdiag,0,Tmat);
Tmat = -Tmat;
%
rhsvec = zeros(length(Tmat),1);
%
% relabel inlet/outlet nodes and impose boundary conditions
conninletnods = setdiff(inletnods,isolatednods);
for iinlet = 1:length(conninletnods)
    bbinletnods(iinlet,1) = find(backbonenods == conninletnods(iinlet)); % "backbone inlet nodes" is the list of indices of backbone inlet nodes, WRT the backbone list
end
connoutletnods = setdiff(outletnods,isolatednods);
for ioutlet = 1:length(connoutletnods)
    bboutletnods(ioutlet,1) = find(backbonenods == connoutletnods(ioutlet));
end
%
% modify transmissibility matrix and RHS vector
%
% inlet nodes
for iinlet = 1:length(bbinletnods)
    inode = bbinletnods(iinlet);
    jnode = setdiff(find(Tmat(inode,:)),inode); % jnode is the set of indices of other nodes with nonzero adjacency to the node with index inode (indices are WRT the full node list)
    rhsvec(jnode) = rhsvec(jnode) - Tmat(jnode,inode)*pinlet;
    rhsvec(inode) = Tmat(inode,inode)*pinlet;
    Tmat(jnode,inode) = 0;
    Tmat(inode,jnode) = 0;
end

%
% outlet nodes
for ioutlet = 1:length(bboutletnods)
    inode = bboutletnods(ioutlet);
    jnode = setdiff(find(Tmat(inode,:)),inode);
    rhsvec(jnode) = rhsvec(jnode) - Tmat(jnode,inode)*poutlet;
    rhsvec(inode) = Tmat(inode,inode)*poutlet;
    Tmat(jnode,inode) = 0;
    Tmat(inode,jnode) = 0;
end

%
% solve system
pvec = Tmat\rhsvec;

return
