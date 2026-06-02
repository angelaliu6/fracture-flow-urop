function [npath,tpath] = ...
    transport_one(linkflow,linkvel,linknods,posnode,...
    conninletnods,conninletlinks,connoutletnods,tol)

% initialize path: nodes and times
npath = [];
tpath = [];

% select a connected inlet node to place the particle
% probability proportional to flowrate of inlet link
qq_in = linkflow(conninletlinks);
QQ_in = cumsum(qq_in)/sum(qq_in);
P = rand();
idx = find(QQ_in > P);
n_inlet = conninletnods(idx(1));

% place the particle at inlet node
it = 1;
npath(it) = n_inlet;
tpath(it) = 0;

keepgoing = 1;
while keepgoing ~= 0

    % find neighboring links and their (outward) flow rate
    nt = npath(it);
    neighlinks = [];
    neighflows = [];
    % links starting at the node
    ll = find(linknods(:,1) == nt);
    neighlinks = [neighlinks; ll];
    neighflows = [neighflows; linkflow(ll)];
    % links ending at the node
    ll = find(linknods(:,2) == nt);
    neighlinks = [neighlinks; ll];
    neighflows = [neighflows; -linkflow(ll)];
    % check mass balance at the node (not inlet node)
    Q_net = sum(neighflows);
    if (abs(Q_net) > tol && it > 1)
        warning('Nodal net flowrate not zero\n')
    end
    
    % choose link according to outward flow rate
    qq_out = neighflows.*(neighflows>0);
    QQ_out = cumsum(qq_out)/sum(qq_out);
    P = rand();
    idx = find(QQ_out > P);
    link_out = neighlinks(idx(1));
    
    % determine new node and transition time
    if linknods(link_out,1) == nt
        ntp1 = linknods(link_out,2);
    else
        ntp1 = linknods(link_out,1);
    end
    linklen = abs(posnode(ntp1,1)-posnode(nt,1)) ...
              + abs(posnode(ntp1,2)-posnode(nt,2));
    dt = linklen/abs(linkvel(link_out));
    
    % move particle
    npath(it+1) = ntp1;
    tpath(it+1) = tpath(it) + dt;
    
    % stop if new node is at outlet, otherwise continue
    if ~isempty(find(connoutletnods == ntp1))
        keepgoing = 0;
    else
        it = it+1;
    end

end

return
