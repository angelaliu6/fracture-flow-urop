function [is_connected, flow_in] = run_trial(Lx, Ly, Nx, Ny, lx, ly, bmin, bmax, visc, pinlet, poutlet, tol, seed)
% using the given parameters, generates one fracture network and computes
% the total flow rate for that network
    is_connected = false;
    flow_in = 0;

    % below, xf__ is a list corresponding to the Nx x-fractures; ini/end/mid
    % are (x,y) coordinates, l is the length, and b is the aperture
    [xfmid,xfini,xfend,yfmid,yfini,yfend,lfx,lfy,bfx,bfy,seed] = generate_fractures(Lx,Ly,Nx,Ny,lx,ly,bmin,bmax,seed);
    % plot_fractures(1,Lx,Ly,xfini,xfend,yfini,yfend)
    % save_figure(1,'Trial Plots',Task 1a')

    [posnode,nfracx,nfracy,nnode] = ...
        fracture_intersections(Lx,Ly,Nx,Ny,xfmid,xfini,xfend,yfmid,yfini,yfend);
    % plot_nodes(1,posnode)
    % save_figure(1,'Trial Plots','Task 2a')

    if nnode == 0
        warning('No fracture intersections\n ')
        return
    end
    [linknods,linkapert,linktrans,nlink] = ...
        build_network(nnode,posnode,nfracx,nfracy,Nx,Ny,bfx,bfy,visc);
    % plot_links(1,nlink,linknods,posnode)
    % save_figure(1,'Trial Plots','Task 2b')
    
    if nlink <= 1
        warning('No links between nodes\n ')
        return
    end
    [Adjmat,backbonenods,isolatednods,inletnods,outletnods] = ...
        backbone(linknods,linktrans,posnode,nnode,Lx);
    % plot_backbone(1,backbonenods,posnode)
    % save_figure(1,'Trial Plots','Task 3')
    
    if isempty(backbonenods)
        warning('Fracture network does not connect inlet with outlet\n ')
        return
    end
    

    [conninletnods,connoutletnods,pvec] = ...
        solve_flow(Adjmat,backbonenods,isolatednods,inletnods,outletnods,pinlet,poutlet);





    % scatter(posnode(backbonenods,1),posnode(backbonenods,2),40,pvec,'filled')
    % colorbar
    % save_figure(1,'Trial Plots','Task 4a')

    is_connected = true;
    [linkflow,linkvel,flow_in] = ...
        flow_velocity(linktrans,linkapert,linknods,backbonenods,isolatednods,...
        conninletnods,connoutletnods,pvec,tol);

    % plot_flow(2,Lx,Ly,linkflow,linknods,backbonenods,posnode,pvec,tol)
    % save_figure(2,'Trial Plots','Task 4b')

    % Task 5a
    
    % total_flow_rate = flow_in;
    % what is effective aperture?
    % pressure_and_flowrate_stats = [max(pvec),min(pvec),mean(pvec),std(pvec),max(linkflow),min(linkflow),mean(linkflow),std(linkflow)] % across different nodes, then across different links
    % figure(3)
    % plot(1:length(pvec),sort(pvec),'k')
    % title('Pressures by Node Index')
    % 
    % figure(4)
    % plot(1:length(linkflow),sort(linkflow))
    % title('Flow Rates by Link Index')
    % close(1);
    % close(2);
    % close(3);
    % close(4);
end