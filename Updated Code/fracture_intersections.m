function [posnode,nfracx,nfracy,nnode] = ...
    fracture_intersections(Lx,Ly,Nx,Ny,xfmid,xfini,xfend,yfmid,yfini,yfend)

% Compute fracture intersections (nodes)
nnode = 0;
posnode = zeros(1,2);
nfracx = zeros(2,1);
nfracy = zeros(2,1);
%
for ix = 1:Nx
    xini = xfini(ix,1);
    xend = xfend(ix,1);
    ynode = xfmid(ix,2);
    %
    % left boundary (must be connected)
    if (xini < 0)
        for iy = 1:Ny
            xmid = yfmid(iy,1);
            yini = yfini(iy,2);
            yend = yfend(iy,2);
            %
            if (xmid < xend) & (ynode > yini) & (ynode < yend)
                nnode = nnode + 1;
                posnode(nnode,:) = [0, ynode];
                nfracx(nnode) = ix;
                nfracy(nnode) = 0;
                break
            end
        end
    end
    %
    % fracture-fracture intersections
    for iy = 1:Ny
        xnode = yfmid(iy,1);
        yini = yfini(iy,2);
        yend = yfend(iy,2);
        %
        if (xnode > xini) & (xnode < xend) & (ynode > yini) & (ynode < yend)
            nnode = nnode + 1;
            posnode(nnode,:) = [xnode, ynode];
            nfracx(nnode) = ix;
            nfracy(nnode) = iy;
        end
    end
    %
    % right boundary (must be connected)
    if (xend > Lx)
        for iy = 1:Ny
            xmid = yfmid(iy,1);
            yini = yfini(iy,2);
            yend = yfend(iy,2);
            %
            if (xmid > xini) & (ynode > yini) & (ynode < yend)
                nnode = nnode + 1;
                posnode(nnode,:) = [Lx, ynode];
                nfracx(nnode) = ix;
                nfracy(nnode) = Ny+1;
                break
            end
        end
    end
end

return
