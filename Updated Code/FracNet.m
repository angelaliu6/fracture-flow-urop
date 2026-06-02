% Flow and Transport through a Fracture Network
clear all
close all

%% Parameters
%
% domain dimensions
Lx = 200; % [m]
Ly = 100; % [m]
%
% number of fractures for each family
Nx = 100; % [m]
Ny = 200; % [m]
%
% maximum fracture length
lx = 60; % [m]
ly = 30; % [m]
%
% min and max fracture gap width
bmin = 1e-3; % [m]
bmax = 5e-3; % [m]
%
% fluid viscosity
visc = 1e-3; % [Pa.s]
%
% boundary pressures
pinlet = 1e5; % [Pa]
poutlet = 0; % [Pa]
%
% number of particles
Np = 1000;
%
% small tolerance
tol = 1e-12;
%
% plot or not / save figure or not
ioplot = 1;
iosave = 0;


%% Generate fractures
%
seed = 1616917526;
% seed = NaN;
%
[xfmid,xfini,xfend,yfmid,yfini,yfend,lfx,lfy,bfx,bfy,seed] = ...
    generate_fractures(Lx,Ly,Nx,Ny,lx,ly,bmin,bmax,seed);

% plot fractures
if ioplot == 1
    figno = 1;
    plot_fractures(figno,Lx,Ly,xfini,xfend,yfini,yfend)
end
%
if iosave == 1
    figname = 'fig_fractures';
    save_figure(figno,figname)
end

%% Determine fracture intersections
[posnode,nfracx,nfracy,nnode] = ...
    fracture_intersections(Lx,Ly,Nx,Ny,xfmid,xfini,xfend,yfmid,yfini,yfend);
%
% plot fracture intersections
if ioplot == 1
    figno = 1;
    plot_nodes(figno,posnode)
end
%
if iosave == 1
    figname = 'fig_intersections';
    save_figure(figno,figname)
end

% if there are no fracture intersections, return
if nnode == 0
    warning('No fracture intersections\n ')
    flow_in = 0;
    return
end

%% Build fracture network
% Determine links between nodes and their transmissibilities
[linknods,linkapert,linktrans,nlink] = ...
    build_network(nnode,posnode,nfracx,nfracy,Nx,Ny,bfx,bfy,visc);
%
if ioplot == 1
    figno = 1;
    plot_links(figno,nlink,linknods,posnode)
end
%
if iosave == 1
    figname = 'fig_network';
    save_figure(figno,figname)
end

% if there are no links between nodes, return
if nlink <= 1
    warning('No links between nodes\n ')
    flow_in = 0;
    return
end

%% Determine backbone of the fracture network
[Adjmat,backbonenods,isolatednods,inletnods,outletnods] = ...
    backbone(linknods,linktrans,posnode,nnode,Lx);
%
% plot nodes in backbone
if ioplot == 1
    figno = 1;
    plot_backbone(figno,backbonenods,posnode)
end
%
if iosave == 1
    figname = 'fig_backbone';
    save_figure(figno,figname)
end

% if backbone is empty, return
if isempty(backbonenods)
    warning('Fracture network does not connect inlet with outlet\n ')
    flow_in = 0;
    return
end

%% Solve flow equation for node pressures
[conninletnods,connoutletnods,pvec] = ...
    solve_flow(Adjmat,backbonenods,isolatednods,inletnods,outletnods,pinlet,poutlet);


%% Compute flowrate and velocity through each link
[linkflow,linkvel,conninletlinks,flow_in] = ...
    flow_velocity(linktrans,linkapert,linknods,backbonenods,isolatednods,...
    conninletnods,connoutletnods,pvec,tol);
%
% plot flowrate along links and pressure at nodes
if ioplot == 1
    figno = 2;
    plot_flow(figno,Lx,Ly,linkflow,linknods,backbonenods,posnode,pvec,tol)
end
%
if iosave == 1
    figname = 'fig_flow';
    save_figure(figno,figname)
end

%% Simulate transport
[npaths,tpaths] = ...
    simul_transport(linkflow,linkvel,linknods,posnode,...
    conninletnods,conninletlinks,connoutletnods,tol,Np);
%
% plot particle paths
if ioplot == 1
    figno = 3;
    plot_paths(figno,Lx,Ly,npaths,posnode)
end
%
if iosave == 1
    figname = 'fig_paths';
    save_figure(figno,figname)
end
%
% plot particle space-time plot
if ioplot == 1
    figno = 4;
    plot_xptp(figno,npaths,tpaths,posnode)
end
%
if iosave == 1
    figname = 'fig_xptp';
    save_figure(figno,figname)
end
%
% plot breakthrough times
if ioplot == 1
    figno = 5;
    plot_bkt(figno,tpaths)
end
%
if iosave == 1
    figname = 'fig_bkt';
    save_figure(figno,figname)
end

return
