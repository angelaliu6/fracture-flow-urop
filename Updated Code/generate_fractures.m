function [xfmid,xfini,xfend,yfmid,yfini,yfend,lfx,lfy,bfx,bfy,seed] = ...
    generate_fractures(Lx,Ly,Nx,Ny,lx,ly,bmin,bmax,seed)
%
% Generate fractures -- two families

%
% control random number generation
if isnan(seed)
    rng('shuffle');
    state = rng();
    seed = state.Seed;
else
    rng(seed);
end

%
% generate positions of fracture midpoints
xfmid = [Lx*rand(Nx,1),sort(Ly*rand(Nx,1))];
yfmid = [sort(Lx*rand(Ny,1)),Ly*rand(Ny,1)];
%
% generate fracture lengths (uniform distrib)
lfx = lx*rand(Nx,1);
lfy = ly*rand(Ny,1);
%
% generate fracture widths (uniform distrib)
bfx = bmin + (bmax - bmin)*rand(Nx,1);
bfy = bmin + (bmax - bmin)*rand(Ny,1);
%
% initial and end points of fractures
xfini = [xfmid(:,1)-0.5*lfx,xfmid(:,2)];
xfend = [xfmid(:,1)+0.5*lfx,xfmid(:,2)];
yfini = [yfmid(:,1),yfmid(:,2)-0.5*lfy];
yfend = [yfmid(:,1),yfmid(:,2)+0.5*lfy];

return

