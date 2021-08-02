function [jlist,md,dlist]   = SDL_Rand(node,tbl,N)

% create N sets of new nodes
% Input:
% ---      node,      p nodes from original data, 1xp
% ---      tbl,       table containing the XYZ coordinates of all of the nodes, 148x4
% ---      N,         number of sets of nodes should be made
% Output:
% ---      jlist,     randomly generated N sets of nodes, Nxp
% ---      md,        mean distance between paris of roginal nodes
% ---      dlist,     mean distance between pairs of randomly chosen nodes


nn = size(tbl,1);   % number of nodes in the whole brain, i.e. 148
pn = length(node);  % number of nodes within each set
nl = sum(node<=74); % number of nodes in the left hemisphere, left = 1:74, right = 75:148
md = pdist(tbl{node,2:4}); md = mean(md)/size(md,1); % mean distance between pairs of nodes
jlist = []; dlist = [];

for i = 1:N
    flag = 1; % 1=need to search for new nodes, 0=found
    while flag
        node1 = randperm(nn);
        node1 = node1(1:pn);
        if sum(node1<=74) == nl % if the number of nodes in either hemisphere is consistent with the original nodes
            flag = 0;
            jlist(i,:) = node1;
            
            md1 = pdist(tbl{node1,2:4}); md1 = mean(md1)/size(md1,1); % mean distance between pairs of nodes
            dlist(i) = md1; % mean distance between pairs of randomly chosen nodes
        end
    end
end


%% End
end