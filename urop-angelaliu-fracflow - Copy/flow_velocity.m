% linkflow[l] stores the flow Q through the lth link (the volume of liquid flow per unit depth per unit time)
% linkvel[l] stores the velocity q through the lth link (the volume of liquid flow per unit area per unit time)
function [linkflow,linkvel,flow_in] = ...
    flow_velocity(linktrans,linkwidth,linknods,backbonenods,isolatednods,...
    conninletnods,connoutletnods,pvec,tol)
%
% Compute flowrate and velocity through each link

linkflow = zeros(size(linktrans));
linkvel = linkflow;
for ilink = 1:length(linktrans)
    %
    if (isempty(intersect(isolatednods,linknods(ilink,:))))
        inode = find(backbonenods == linknods(ilink,1));
        jnode = find(backbonenods == linknods(ilink,2));
        linkflow(ilink) = linktrans(ilink)*(pvec(inode) - pvec(jnode));
        linkvel(ilink) = linkflow(ilink)/linkwidth(ilink);
    end
end
%
% find the indices of the links at the boundaries
conninletlinks = zeros(size(conninletnods));
for iinlet = 1:length(conninletnods)
    conninletlinks(iinlet) = find(linknods(:,1) == conninletnods(iinlet));
end
connoutletlinks = zeros(size(connoutletnods));
for ioutlet = 1:length(connoutletnods)
    connoutletlinks(ioutlet) = find(linknods(:,2) == connoutletnods(ioutlet));
end
% check that the total computed flow in from left and out from right are the same
flow_in = sum(linkflow(conninletlinks));
flow_out = sum(linkflow(connoutletlinks));
if any(linkflow(conninletlinks) < 0)
    disp(' Warning: Flowrate in some inlet nodes < 0\n')
end
if any(linkflow(connoutletlinks) < 0)
    disp(' Warning: Flowrate in some outlet nodes < 0\n')
end
if abs(flow_in - flow_out) > tol
    disp(' Warning: Flowrate in not equal to flowrate out\n')
end

return
