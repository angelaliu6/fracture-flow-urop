function [npaths,tpaths] = ...
    simul_transport(linkflow,linkvel,linknods,posnode,...
    conninletnods,conninletlinks,connoutletnods,tol,Np, timeout)

% initialize cell arrays with node and time paths
npaths = cell(Np,1);
tpaths = cell(Np,1);

if isempty(conninletlinks) || isempty(conninletnods)
    warning('simul_transport: network not connected, returning Inf paths');
    for ip = 1:Np
        tpaths{ip} = [Inf];
    end
    return
end

% loop over particles
for ip = 1:Np
    %
    % simulate one path
    [npath,tpath] = ...
        transport_one(linkflow,linkvel,linknods,posnode,...
        conninletnods,conninletlinks,connoutletnods,tol, timeout);
    %
    % store path
    npaths{ip} = npath;
    tpaths{ip} = tpath;

end

return
