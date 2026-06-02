function [linknods,linkapert,linktrans,nlink] = ...
    build_network(nnode,posnode,nfracx,nfracy,Nx,Ny,bfx,bfy,visc)

% Build fracture network
% Determine links between nodes and their transmissibilities

%
% x-links (x-fracture number already sorted in ascending order)
nlink = 0;
linknods = zeros(1,2);
linktrans = zeros(2,1);
linkapert = zeros(2,1);
for inode = 1:nnode-1
    %
    ifracx = nfracx(inode);
    jnode = inode + 1;
    jfracx = nfracx(jnode);
    if (jfracx == ifracx)
        nlink = nlink + 1;
        linknods(nlink,1) = inode;
        linknods(nlink,2) = jnode;
        linkapert(nlink) = bfx(ifracx);
        linktrans(nlink) = bfx(ifracx)^3/(12*visc)/...
            (posnode(jnode,1)-posnode(inode,1));
    end
end

% y-links (need to sort for y-fracture number in ascending order)
[Y,ids] = sort(nfracy);
for inode = 1:nnode-1
    %
    ifracy = nfracy(ids(inode));
    jnode = inode + 1;
    jfracy = nfracy(ids(jnode));
    if (jfracy == ifracy) & (ifracy >= 1) & (ifracy <= Ny)
        nlink = nlink + 1;
        linknods(nlink,1) = ids(inode);
        linknods(nlink,2) = ids(jnode);
        linkapert(nlink) = bfy(ifracy);
        linktrans(nlink) = bfy(ifracy)^3/(12*visc)/...
            (posnode(ids(jnode),2)-posnode(ids(inode),2));
    end
end

return

