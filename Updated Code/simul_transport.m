function [npaths,tpaths] = ...
    simul_transport(linkflow,linkvel,linknods,posnode,...
    conninletnods,conninletlinks,connoutletnods,tol,Np)

% initialize cell arrays with node and time paths
npaths = cell(Np,1);
tpaths = cell(Np,1);

% loop over particles
for ip = 1:Np
    %
    % simulate one path
    [npath,tpath] = ...
        transport_one(linkflow,linkvel,linknods,posnode,...
        conninletnods,conninletlinks,connoutletnods,tol);
    %
    % store path
    npaths{ip} = npath;
    tpaths{ip} = tpath;

end

return
