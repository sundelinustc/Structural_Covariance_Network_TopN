function SDL_SCA_TopN4(SDL)

% SCA within the top-N areas with CT or SA reduction


%% Loading information
% Residuals
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149}; % residules
if strcmp(SDL.data_type{1}(4:end),'Gender')
    idx11 = strcmp(T.Group,'PTSD') & strcmp(T.Gender,'M');  % group 11, PTSD & Male
    idx12 = strcmp(T.Group,'PTSD') & strcmp(T.Gender,'F');  % group 12, PTSD & Female
    idx21 = strcmp(T.Group,'CONT') & strcmp(T.Gender,'M');  % group 21, CONT & Male
    idx22 = strcmp(T.Group,'CONT') & strcmp(T.Gender,'F');  % group 22, CONT & Female
elseif strcmp(SDL.data_type{1}(4:end),'Dep')
    %     try
    %         idx11 = strcmp(T.Group,'PTSD') & strcmp(T.Dep,'1');  % group 11, PTSD & depression
    %         idx12 = strcmp(T.Group,'PTSD') & strcmp(T.Dep,'0');  % group 12, PTSD & no-depression
    %         idx21 = strcmp(T.Group,'CONT') & strcmp(T.Dep,'1');  % group 21, CONT & depression
    %         idx22 = strcmp(T.Group,'CONT') & strcmp(T.Dep,'0');  % group 22, CONT & no-depression
    %     catch
    idx11 = strcmp(T.Group,'PTSD') & (T.Dep==1);  % group 11, PTSD & depression
    idx12 = strcmp(T.Group,'PTSD') & (T.Dep==0);  % group 12, PTSD & no-depression
    idx21 = strcmp(T.Group,'CONT') & (T.Dep==1);  % group 21, CONT & depression
    idx22 = strcmp(T.Group,'CONT') & (T.Dep==0);  % group 22, CONT & no-depression
    %     end
else
end

[sum(idx11),sum(idx12),sum(idx21),sum(idx22)]

% Effect size
fdir = fullfile(SDL.out,SDL.data_type{1}(1:2)); 
fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{1}(1:2),'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn,'kk','d'); fprintf('Loaded: Effect Size <- %s\n\n\n',fn);
% T    --- table containing raw CT values and sbj info
% R    --- residules of CT after lme regression
% idx1 --- row index of PTSD in T & R
% idx2 --- row index of CONT in T & R
% d    --- Cohen's d for PTSD vs CONT
% kk   --- the list of ranked area No. for PTSD vs CONT based on d, i.e. kk(1) showes the area with the largest CT reduction  

% Region's label & XYZ coordinates
fdir = fullfile(SDL.path,'Original');
load(fullfile(fdir,'RegionLabel.mat')); % load struct D which contains the 148 areas name & XYZ locations
tbl = cell2table(D,'VariableNames',{'Area','X','Y','Z'});


%% Structural Covariance Matrix -- Original
if strcmp(SDL.ana_type{1},'corr')
    func = @corr;
elseif strcmp(SDL.ana_type{1},'partialcorr')
    func = @partialcorr;
else
end
S0 = SDL_SCA4(R,idx11,idx12,idx21,idx22,func,0);     % SC of original PTSD & CONT data

% %% Generate multiple sets of random nodes with matching hemispherical distribution & distance
% SC0 = []; SC1 = []; % mean SC values of top-N, original(0) & permutation(1) data
% NR = 5000; % number of sets of N random nodes
% for i = 2:length(kk) % per top-N
%     % original nodes
%     SC0(i).node   = kk(1:i); % No. of the top-N areas
%     
%     % random nodes
%     [jlist,md,dlist] = SDL_Rand(SC0(i).node,tbl,NR); % the No. of the top-N areas
%     
%     % matching Euclidian distance
%     % it is weired that can't use ttest(dlist-md);
%     k = 0;
%     while (k < 3000) && (signrank(dlist-md) < 0.05) % loop<3000 times if the mean distance of random nodes is different from original nodes distance
%         k = k + 1;
%         if mean(dlist) < md % if mean distance of random nodes is smaller than md
%             d1 = mean(dlist);
%             while d1 < md
%                 [j1,~,d1] = SDL_Rand(SC0(i).node,tbl,1);
%             end
%             idx = find(dlist==min(dlist)); % replace the set with minimal distance with a new set with distance > md
%             idx = idx(1);
%             jlist(idx,:) = j1;
%             dlist(idx) = d1;
%         elseif mean(dlist) > md % if mean distance of random nodes is larger than md
%             d1 = mean(dlist);
%             while d1 > md
%                 [j1,~,d1] = SDL_Rand(SC0(i).node,tbl,1);
%             end
%             idx = find(dlist==max(dlist)); % replace the set with maximal distance with a new set with distance < md
%             idx = idx(1);
%             jlist(idx,:) = j1;
%             dlist(idx) = d1;
%         else
%         end
%     end
%     
%     SC1(i).node   = jlist; % No. of the top-N areas
%     fprintf('Generated: randomly chosen %d sets of %3d nodes with matching hemisphere and distance\n',NR,i);
%     
% end
% %save
% fdir = fullfile(SDL.out,SDL.data_type{1}); 
% fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% save(fn,'S0','SC0','SC1','-v7.3');
% fprintf('Saved: Top-N analyses midway results saved in \n-->%s\n',fn);

%% SC, CI & p values of top-N areas with CT/SA reduction
% loading
fdir = fullfile(SDL.out,SDL.data_type{1}(1:2)); % load SC0 and SC1 from existed results, i.e. CleanData, to avoid wasting time
fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1}(1:2),'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');
fprintf('Loaded: Top-N analyses midway results loaded from \n<--%s\n',fn);

% mean SC, CI and p values
CIType = 0.95; % 95% confidene interval (CI) 
for i = 2:size(SC0,2) % per top-N
    % original data
    node0 = SC0(i).node;
    SC0(i).Z11 = triu(S0.Y11(node0,node0),1); SC0(i).Z12 = triu(S0.Y12(node0,node0),1); SC0(i).Z21 = triu(S0.Y21(node0,node0),1); SC0(i).Z22 = triu(S0.Y22(node0,node0),1);% original data, values above diagonal
    SC0(i).Z11c = nonzeros(SC0(i).Z11); SC0(i).Z12c = nonzeros(SC0(i).Z12); SC0(i).Z21c = nonzeros(SC0(i).Z21); SC0(i).Z22c = nonzeros(SC0(i).Z22);% non-zero values, Num of connections X 1
    
    SC0(i).mean11raw = mean(SC0(i).Z11c); % mean of all connections, 1X1, original data, group 11
    SC0(i).mean12raw = mean(SC0(i).Z12c); % mean of all connections, 1X1, original data, group 12
    SC0(i).mean21raw = mean(SC0(i).Z21c); % mean of all connections, 1X1, original data, group 21
    SC0(i).mean22raw = mean(SC0(i).Z22c); % mean of all connections, 1X1, original data, group 22
    
    SC0(i).mean11pos = mean(SC0(i).Z11c(SC0(i).Z11c>0));
    SC0(i).mean12pos = mean(SC0(i).Z12c(SC0(i).Z12c>0));
    SC0(i).mean21pos = mean(SC0(i).Z21c(SC0(i).Z21c>0));
    SC0(i).mean22pos = mean(SC0(i).Z22c(SC0(i).Z22c>0));
    
    SC0(i).mean11neg = mean(SC0(i).Z11c(SC0(i).Z11c<0));
    SC0(i).mean12neg = mean(SC0(i).Z12c(SC0(i).Z12c<0));
    SC0(i).mean21neg = mean(SC0(i).Z21c(SC0(i).Z21c<0));
    SC0(i).mean22neg = mean(SC0(i).Z22c(SC0(i).Z22c<0));
    
    SC0(i).mean11abs = mean(abs(SC0(i).Z11c));
    SC0(i).mean12abs = mean(abs(SC0(i).Z12c));
    SC0(i).mean21abs = mean(abs(SC0(i).Z21c));
    SC0(i).mean22abs = mean(abs(SC0(i).Z22c));
    
    SC0(i).Z11 = []; SC0(i).Z12 = []; SC0(i).Z21 = []; SC0(i).Z22 = [];
    SC0(i).Z11c = []; SC0(i).Z12c = [];SC0(i).Z21c = []; SC0(i).Z22c = [];
    
    % permuted data
    for j = 1:size(SC1(i).node,1) % per permutation (totally 5000 times)
        node1 = SC1(i).node(j,:);
        SC1(i).Z11(:,:,j) = triu(S0.Y11(node1,node1),1); SC1(i).Z12(:,:,j) = triu(S0.Y12(node1,node1),1); SC1(i).Z21(:,:,j) = triu(S0.Y21(node1,node1),1); SC1(i).Z22(:,:,j) = triu(S0.Y22(node1,node1),1);% original data, values above diagonal
        SC1(i).Z11c(:,j)  = nonzeros(SC1(i).Z11(:,:,j)); SC1(i).Z12c(:,j) = nonzeros(SC1(i).Z12(:,:,j)); SC1(i).Z21c(:,j)  = nonzeros(SC1(i).Z21(:,:,j)); SC1(i).Z22c(:,j) = nonzeros(SC1(i).Z22(:,:,j));% non-zero values, Num of connections X Num of permutations
   
        SC1(i).mean11pos(j) = mean(nonzeros(SC1(i).Z11c(:,j).*(SC1(i).Z11c(:,j)>0)));
        SC1(i).mean12pos(j) = mean(nonzeros(SC1(i).Z12c(:,j).*(SC1(i).Z12c(:,j)>0)));
        SC1(i).mean21pos(j) = mean(nonzeros(SC1(i).Z21c(:,j).*(SC1(i).Z21c(:,j)>0)));
        SC1(i).mean22pos(j) = mean(nonzeros(SC1(i).Z22c(:,j).*(SC1(i).Z22c(:,j)>0)));
        
        SC1(i).mean11neg(j) = mean(nonzeros(SC1(i).Z11c(:,j).*(SC1(i).Z11c(:,j)<0)));
        SC1(i).mean12neg(j) = mean(nonzeros(SC1(i).Z12c(:,j).*(SC1(i).Z12c(:,j)<0)));
        SC1(i).mean21neg(j) = mean(nonzeros(SC1(i).Z21c(:,j).*(SC1(i).Z21c(:,j)<0)));
        SC1(i).mean22neg(j) = mean(nonzeros(SC1(i).Z22c(:,j).*(SC1(i).Z22c(:,j)<0)));
    end
    SC1(i).mean11raw = mean(SC1(i).Z11c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 11
    SC1(i).mean12raw = mean(SC1(i).Z12c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 12
    SC1(i).mean21raw = mean(SC1(i).Z21c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 21
    SC1(i).mean22raw = mean(SC1(i).Z22c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 22
    
    SC1(i).mean11abs = mean(abs(SC1(i).Z11c),1);
    SC1(i).mean12abs = mean(abs(SC1(i).Z12c),1);
    SC1(i).mean21abs = mean(abs(SC1(i).Z21c),1);
    SC1(i).mean22abs = mean(abs(SC1(i).Z22c),1);
    
    SC1(i).Z11 = []; SC1(i).Z12 = [];SC1(i).Z21 = []; SC1(i).Z22 = [];
    SC1(i).Z11c = []; SC1(i).Z12c = [];SC1(i).Z21c = []; SC1(i).Z22c = [];
    
%     % confidence interval (CI)
%     SC1(i).CI11(:)    = SDL_CI(SC1(i).mean11',CIType); % 95% CI
%     SC1(i).CI12(:)    = SDL_CI(SC1(i).mean12',CIType); % 95% CI
%     SC1(i).CI21(:)    = SDL_CI(SC1(i).mean21',CIType); % 95% CI
%     SC1(i).CI22(:)    = SDL_CI(SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_11v21(:) = SDL_CI(SC1(i).mean11' - SC1(i).mean21',CIType); % 95% CI
%     SC1(i).diffCI_12v22(:) = SDL_CI(SC1(i).mean12' - SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_11v12(:) = SDL_CI(SC1(i).mean11' - SC1(i).mean12',CIType); % 95% CI
%     SC1(i).diffCI_21v22(:) = SDL_CI(SC1(i).mean21' - SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_main1(:) = SDL_CI(SC1(i).mean11' + SC1(i).mean12' - SC1(i).mean21' - SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_main2(:) = SDL_CI(SC1(i).mean11' + SC1(i).mean21' - SC1(i).mean12' - SC1(i).mean22',CIType); % 95% CI
%     SC1(i).diffCI_inter(:) = SDL_CI(SC1(i).mean11' - SC1(i).mean12' - SC1(i).mean21' + SC1(i).mean22',CIType); % 95% CI
%     
%     % p values
%     SC1(i).p11     = SDL_p_permutation(SC0(i).mean11,SC1(i).mean11); % group1 vs permutation
%     SC1(i).p12     = SDL_p_permutation(SC0(i).mean12,SC1(i).mean12); % group2 vs permutation
%     SC1(i).p21     = SDL_p_permutation(SC0(i).mean21,SC1(i).mean21); % group1 vs permutation
%     SC1(i).p22     = SDL_p_permutation(SC0(i).mean22,SC1(i).mean22); % group2 vs permutation
%     SC1(i).p_diff_11v21 = SDL_p_permutation(SC0(i).mean11-SC0(i).mean21,SC1(i).mean11-SC1(i).mean21); % group 11-21 vs permutation
%     SC1(i).p_diff_12v22 = SDL_p_permutation(SC0(i).mean12-SC0(i).mean22,SC1(i).mean12-SC1(i).mean22); % group 12-22 vs permutation
%     SC1(i).p_diff_11v12 = SDL_p_permutation(SC0(i).mean11-SC0(i).mean12,SC1(i).mean11-SC1(i).mean12); % group 11-12 vs permutation
%     SC1(i).p_diff_21v22 = SDL_p_permutation(SC0(i).mean21-SC0(i).mean22,SC1(i).mean21-SC1(i).mean22); % group 21-22 vs permutation
%     SC1(i).p_diff_main1 = SDL_p_permutation(SC0(i).mean11+SC0(i).mean12-SC0(i).mean21-SC0(i).mean22,...
%         SC1(i).mean11+SC1(i).mean12-SC1(i).mean21-SC1(i).mean22); % group 11+12-21-22 vs permutation
%     SC1(i).p_diff_main2 = SDL_p_permutation(SC0(i).mean11-SC0(i).mean12+SC0(i).mean21-SC0(i).mean22,...
%         SC1(i).mean11-SC1(i).mean12+SC1(i).mean21-SC1(i).mean22); % group 11-12+21-22 vs permutation
%     SC1(i).p_diff_inter = SDL_p_permutation(SC0(i).mean11-SC0(i).mean12-SC0(i).mean21+SC0(i).mean22,...
%         SC1(i).mean11-SC1(i).mean12-SC1(i).mean21+SC1(i).mean22); % group 11-12-21+22 vs permutation
    
    fprintf('Calculated: SC, CI & p-val for Top-%3d Areas\n',i);
end

%save
fdir = fullfile(SDL.out,SDL.data_type{1}); mkdir(fdir);
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
save(fn,'SC0','SC1','-v7.3');
fprintf('Saved: Top-N analyses results saved in \n-->%s\n',fn);

% saved into a .mat fil, to be read by R
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn);
data = zeros(5001,147,16); % 5001(1 actual + 5000 random) x 147(network size: 2~148) x 16(PTSD-male,PTSD-female,non-PTSD-male,non-PTSD-female groups x raw/pos/neg/abs)
for j = 1:147 % per network size (2~148)
    data(1,j,1) = SC0(j+1).mean11raw; data(2:5001,j,1) = SC1(j+1).mean11raw; % mean11raw
    data(1,j,2) = SC0(j+1).mean12raw; data(2:5001,j,2) = SC1(j+1).mean12raw; % mean12raw
    data(1,j,3) = SC0(j+1).mean21raw; data(2:5001,j,3) = SC1(j+1).mean21raw; % mean21raw
    data(1,j,4) = SC0(j+1).mean22raw; data(2:5001,j,4) = SC1(j+1).mean22raw; % mean22raw  
    
    data(1,j,5) = SC0(j+1).mean11pos; data(2:5001,j,5) = SC1(j+1).mean11pos; % mean11pos
    data(1,j,6) = SC0(j+1).mean12pos; data(2:5001,j,6) = SC1(j+1).mean12pos; % mean12pos
    data(1,j,7) = SC0(j+1).mean21pos; data(2:5001,j,7) = SC1(j+1).mean21pos; % mean21pos
    data(1,j,8) = SC0(j+1).mean22pos; data(2:5001,j,8) = SC1(j+1).mean22pos; % mean22pos    
    
    data(1,j,9) = SC0(j+1).mean11neg; data(2:5001,j,9) = SC1(j+1).mean11neg; % mean11neg
    data(1,j,10)= SC0(j+1).mean12neg; data(2:5001,j,10)= SC1(j+1).mean12neg; % mean12neg
    data(1,j,11)= SC0(j+1).mean21neg; data(2:5001,j,11)= SC1(j+1).mean21neg; % mean21neg
    data(1,j,12)= SC0(j+1).mean22neg; data(2:5001,j,12)= SC1(j+1).mean22neg; % mean22neg    
    
    data(1,j,13)= SC0(j+1).mean11abs; data(2:5001,j,13)= SC1(j+1).mean11abs; % mean11abs
    data(1,j,14)= SC0(j+1).mean12abs; data(2:5001,j,14)= SC1(j+1).mean12abs; % mean12abs
    data(1,j,15)= SC0(j+1).mean21abs; data(2:5001,j,15)= SC1(j+1).mean21abs; % mean21abs
    data(1,j,16)= SC0(j+1).mean22abs; data(2:5001,j,16)= SC1(j+1).mean22abs; % mean22abs
end
fn = fullfile(fdir,'Results_for_R.mat');
save(fn,'data');


% %% Plot SC values of Top-N areas with reduction in PTSD, CONT and PTSD-CONT 
% fdir = fullfile(SDL.out,SDL.data_type{1}); 
% fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% % load(fn,'SC0','SC1');
% 
% X = 2:size(SC0,1);
% for i = 2:size(SC0,2)
%     X11(i-1) = SC0(i).mean11; X12(i-1) = SC0(i).mean12; X21(i-1) = SC0(i).mean21; X22(i-1) = SC0(i).mean22;
%     
%     p11(i-1) = SC1(i).p11; p12(i-1) = SC1(i).p12; p21(i-1) = SC1(i).p21; p22(i-1) = SC1(i).p22;
%     
%     X11R(i-1) = mean(SC1(i).mean11); X12R(i-1) = mean(SC1(i).mean12); X21R(i-1) = mean(SC1(i).mean21); X22R(i-1) = mean(SC1(i).mean22);
%     X11RA(i-1,:) = SC1(i).mean11; X12RA(i-1,:) = SC1(i).mean12; X21RA(i-1,:) = SC1(i).mean21; X22RA(i-1,:) = SC1(i).mean22;
%     
%     CI11L(i-1) = SC1(i).CI11(1); CI11R(i-1) = SC1(i).CI11(2); CI12L(i-1) = SC1(i).CI12(1); CI12R(i-1) = SC1(i).CI12(2);
%     CI21L(i-1) = SC1(i).CI21(1); CI21R(i-1) = SC1(i).CI21(2); CI22L(i-1) = SC1(i).CI22(1); CI22R(i-1) = SC1(i).CI22(2);
% 
%     diff_11v21L(i-1) = SC1(i).diffCI_11v21(1); diff_11v21R(i-1) = SC1(i).diffCI_11v21(2);
%     diff_12v22L(i-1) = SC1(i).diffCI_12v22(1); diff_12v22R(i-1) = SC1(i).diffCI_12v22(2);
%     diff_11v12L(i-1) = SC1(i).diffCI_11v12(1); diff_11v12R(i-1) = SC1(i).diffCI_11v12(2);
%     diff_21v22L(i-1) = SC1(i).diffCI_21v22(1); diff_21v22R(i-1) = SC1(i).diffCI_21v22(2);
%     diff_main1L(i-1)  = SC1(i).diffCI_main1(1); diff_main1R(i-1) = SC1(i).diffCI_main1(2);
%     diff_main2L(i-1)  = SC1(i).diffCI_main2(1); diff_main2R(i-1) = SC1(i).diffCI_main2(2);
%     diff_interL(i-1)  = SC1(i).diffCI_inter(1); diff_interR(i-1) = SC1(i).diffCI_inter(2);
%     
%     
%     p_diff_11v21(i-1) = SC1(i).p_diff_11v21; 
%     p_diff_12v22(i-1) = SC1(i).p_diff_12v22;
%     p_diff_11v12(i-1) = SC1(i).p_diff_11v12;
%     p_diff_21v22(i-1) = SC1(i).p_diff_21v22;
%     p_diff_main1(i-1) = SC1(i).p_diff_main1;
%     p_diff_main2(i-1) = SC1(i).p_diff_main2;
%     p_diff_inter(i-1) = SC1(i).p_diff_inter;
% end
% 
% 
% %% FDR corrected p values per top-area
% fprintf('\n\nFDR correction: group 11 vs random\n')
% [~, ~, ~, adj_p11] = fdr_bh(p11,0.05,'pdep','yes'); 
% find(adj_p11<=0.05)
% fprintf('\n\nFDR correction: group 12 vs random\n')
% [~, ~, ~, adj_p12] = fdr_bh(p12,0.05,'pdep','yes');
% find(adj_p12<=0.05)
% fprintf('\n\nFDR correction: group 21 vs random\n')
% [~, ~, ~, adj_p21] = fdr_bh(p21,0.05,'pdep','yes'); 
% find(adj_p21<=0.05)
% fprintf('\n\nFDR correction: group 22 vs random\n')
% [~, ~, ~, adj_p22] = fdr_bh(p22,0.05,'pdep','yes');
% find(adj_p22<=0.05)
% 
% fprintf('\n\nFDR correction: group 11v21 vs random\n')
% [~, ~, ~, adj_p_diff_11v21] = fdr_bh(p_diff_11v21,0.05,'pdep','yes'); 
% find(adj_p_diff_11v21<=0.05)
% fprintf('\n\nFDR correction: group 12v22 vs random\n')
% [~, ~, ~, adj_p_diff_12v22] = fdr_bh(p_diff_12v22,0.05,'pdep','yes'); 
% find(adj_p_diff_12v22<=0.05)
% fprintf('\n\nFDR correction: group 11v12 vs random\n')
% [~, ~, ~, adj_p_diff_11v12] = fdr_bh(p_diff_11v12,0.05,'pdep','yes'); 
% find(adj_p_diff_11v12<=0.05)
% fprintf('\n\nFDR correction: group 21v22 vs random\n')
% [~, ~, ~, adj_p_diff_21v22] = fdr_bh(p_diff_21v22,0.05,'pdep','yes'); 
% find(adj_p_diff_21v22<=0.05)
% % fprintf('\n\nFDR correction: main1 (PTSD vs CONT) vs random\n')
% % [~, ~, ~, adj_p_diff_main1] = fdr_bh(p_diff_main1,0.05,'pdep','yes'); 
% % find(adj_p_diff_main1<=0.05)
% % fprintf('\n\nFDR correction: main2 vs random\n')
% % [~, ~, ~, adj_p_diff_main2] = fdr_bh(p_diff_main2,0.05,'pdep','yes'); 
% % find(adj_p_diff_main2<=0.05)
% fprintf('\n\nFDR correction: interaction vs random\n')
% [~, ~, ~, adj_p_diff_inter] = fdr_bh(p_diff_inter,0.05,'pdep','yes'); 
% find(adj_p_diff_inter<=0.05)
% 
% 
% %% p value for area under curve
% p11a = SDL_p_permutation(mean(X11),mean(X11RA,1))
% p12a = SDL_p_permutation(mean(X12),mean(X12RA,1))
% p21a = SDL_p_permutation(mean(X21),mean(X21RA,1))
% p22a = SDL_p_permutation(mean(X22),mean(X22RA,1))
% 
% p_diff_11v21a = SDL_p_permutation(mean(X11)-mean(X21),mean(X11RA,1)-mean(X21RA,1))
% p_diff_12v22a = SDL_p_permutation(mean(X12)-mean(X22),mean(X12RA,1)-mean(X22RA,1))
% p_diff_11v12a = SDL_p_permutation(mean(X11)-mean(X12),mean(X11RA,1)-mean(X12RA,1))
% p_diff_21v22a = SDL_p_permutation(mean(X21)-mean(X22),mean(X21RA,1)-mean(X22RA,1))
% p_diff_11v22a = SDL_p_permutation(mean(X11)-mean(X22),mean(X11RA,1)-mean(X22RA,1))
% p_diff_12v21a = SDL_p_permutation(mean(X12)-mean(X21),mean(X12RA,1)-mean(X21RA,1))
% % p_diff_main1a = SDL_p_permutation(mean(X11)+mean(X12)-mean(X21)-mean(X22),...
% %     mean(X11RA,1)+mean(X12RA,1)-mean(X21RA,1)-mean(X22RA,1))
% % p_diff_main2a = SDL_p_permutation(mean(X11)-mean(X12)+mean(X21)-mean(X22),...
% %     mean(X11RA,1)-mean(X12RA,1)+mean(X21RA,1)-mean(X22RA,1))
% p_diff_intera = SDL_p_permutation(mean(X11)-mean(X12)-mean(X21)+mean(X22),...
%     mean(X11RA,1)-mean(X12RA,1)-mean(X21RA,1)+mean(X22RA,1))
% 
% plist = [p_diff_11v21a,p_diff_12v22a,p_diff_11v12a,p_diff_21v22a,p_diff_11v22a,p_diff_12v21a];
% [~, ~, ~, adj_p_mul] = fdr_bh(plist,0.05,'pdep','yes')
% 
% 
% 
% 
% %% Plot curves & heatmap of SC values
% if strcmp(SDL.data_type{1}(4:end),'Gender')
%     txt1 = 'Male'; txt2 = 'Female';
% elseif strcmp(SDL.data_type{1}(4:end),'Dep')
%     txt1 = 'Depress'; txt2 = 'No-Depress';
% else
% end
% 
% get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
% figure; % original vs permutation
% subplot(3,4,1);
% X = 2:size(S0.Y11,1); x2 = [X, fliplr(X)]; inBetween = [CI11L, fliplr(CI11R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,X11,'r-',X,X11R,'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['PTSD & ',txt1]);
% 
% subplot(3,4,2);
% X = 2:size(S0.Y12,1); x2 = [X, fliplr(X)]; inBetween = [CI12L, fliplr(CI12R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,X12,'r-',X,X12R,'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['PTSD & ',txt2]);
% 
% subplot(3,4,3);
% X = 2:size(S0.Y21,1); x2 = [X, fliplr(X)]; inBetween = [CI21L, fliplr(CI21R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,X21,'r-',X,X21R,'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['CONT & ',txt1]);
% 
% subplot(3,4,4);
% X = 2:size(S0.Y22,1); x2 = [X, fliplr(X)]; inBetween = [CI22L, fliplr(CI22R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,X22,'r-',X,X22R,'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['CONT & ',txt2]);
% 
% subplot(3,4,5);
% X = 2:size(S0.Y11,1); x2 = [X, fliplr(X)]; inBetween = [diff_11v21L, fliplr(diff_11v21R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X11-X21),'r-',X,(X11R-X21R),'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['PTSD-CONT, ',txt1]);
% 
% subplot(3,4,6);
% X = 2:size(S0.Y12,1); x2 = [X, fliplr(X)]; inBetween = [diff_12v22L, fliplr(diff_12v22R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X12-X22),'r-',X,(X12R-X12R),'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['PTSD-CONT, ',txt2]);
% 
% subplot(3,4,7);
% X = 2:size(S0.Y11,1); x2 = [X, fliplr(X)]; inBetween = [diff_11v12L, fliplr(diff_11v12R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X11-X12),'r-',X,(X11R-X12R),'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['PTSD, ',txt1,'-',txt2]);
% 
% subplot(3,4,8);
% X = 2:size(S0.Y21,1); x2 = [X, fliplr(X)]; inBetween = [diff_21v22L, fliplr(diff_21v22R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X21-X22),'r-',X,(X21R-X22R),'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['CONT, ',txt1,'-',txt2]);
% 
% subplot(3,4,9);
% X = 2:size(S0.Y11,1); x2 = [X, fliplr(X)]; inBetween = [diff_main1L, fliplr(diff_main1R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X11+X12-X21-X22),'r-',X,(X11R+X12R-X21R-X22R),'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title('PTSD-CONT');
% 
% subplot(3,4,10);
% X = 2:size(S0.Y11,1); x2 = [X, fliplr(X)]; inBetween = [diff_main2L, fliplr(diff_main2R)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X11-X12+X21-X22),'r-',X,(X11R-X12R+X21R-X22R),'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title([txt1,'-',txt2]);
% 
% subplot(3,4,11);
% X = 2:size(S0.Y11,1); x2 = [X, fliplr(X)]; inBetween = [diff_interL, fliplr(diff_interR)];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
% plot(X,(X11-X12-X21+X22),'r-',X,(X11R-X12R-X21R+X22R),'b-','LineWidth',2);
% xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title(['PTSD(',txt1,'-',txt2,')-CONT(',txt1,'-',txt2,')']);
% 
% savefig(fullfile(SDL.out,SDL.data_type{1},['Results_Curves_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
% 
% % MN1 = []; MN2 = []; MNd = [];
% % % re-organize matrix based on kk
% % for i = 1:length(kk)
% %     for j = 1:length(kk)
% %         MN1(i,j) = S0.Y1(kk(i),kk(j));
% %         MN2(i,j) = S0.Y2(kk(i),kk(j));
% %     end
% % end
% % MNd = MN1 - MN2;
% % % figure; 
% % subplot(3,2,2);imagesc(MN1); colormap('hot'); colorbar; title('PTSD'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([0,1.2]);
% % subplot(3,2,4);imagesc(MN2); colormap('hot'); colorbar; title('CONT'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([0,1.2]);
% % subplot(3,2,6);imagesc(MNd); colormap('hot'); colorbar; title('PTSD-CONT'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([-0.1,0.1]);
% % savefig(fullfile(SDL.out,SDL.data_type,['Results_2Dmap_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
% 
% 
% % %% Plot correlation
% % % Both pos and neg efect size
% % figure('DefaultAxesFontSize',18)
% % ax1 = subplot(1,2,1); % PTSD
% % a = triu(S0.Y1,1); a = a + a'; % structural covariance after removing Inf in the diagnose
% % scatter(ax1,d,mean(a),45,'ko');h1 = lsline(ax1); h1.Color = 'r'; h1.LineWidth = 2; axis([-0.13,0.13,-0.025,0.025]);
% % title('PTSD'); xlabel(['Effect Size of ',SDL.data_type{1},' Difference']); ylabel('Structural Covariance');
% % [r,p] = corr(d',mean(a)');
% % txt = sprintf('R=%1.3f, p=%1.3f',r,p);
% % text(0,0.015,txt,'FontSize',18);
% % 
% % ax2 = subplot(1,2,2); % CONT
% % a = triu(S0.Y2,1); a = a + a'; % structural covariance after removing Inf in the diagnose
% % scatter(ax2,d,mean(a),45,'ko');h2 = lsline(ax2); h2.Color = 'r'; h2.LineWidth = 2;  axis([-0.13,0.13,-0.025,0.025]);
% % title('CONT'); xlabel(['Effect Size of ',SDL.data_type{1},' Difference']); ylabel('Structural Covariance');
% % [r,p] = corr(d',mean(a)');
% % txt = sprintf('R=%1.3f, p=%1.3f',r,p);
% % text(0,0.015,txt,'FontSize',18);
% % 
% % savefig(fullfile(SDL.out,SDL.data_type,['Results_Corr_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
% % 
% % % Only neg efect size
% % figure('DefaultAxesFontSize',18)
% % ax1 = subplot(1,2,1); % PTSD
% % a = triu(S0.Y1,1); a = a + a'; % structural covariance after removing Inf in the diagnose
% % a = [d',mean(a)']; a = a(a(:,1)<0,:);
% % scatter(ax1,a(:,1),a(:,2),45,'ko');h1 = lsline(ax1); h1.Color = 'r'; h1.LineWidth = 2; axis([-0.11,0.01,-0.025,0.025]);
% % title('PTSD'); xlabel(['Effect Size of ',SDL.data_type{1},' Difference']); ylabel('Structural Covariance');
% % [r,p] = corr(a(:,1),a(:,2));
% % txt = sprintf('R=%1.3f, p=%1.3f',r,p);
% % text(-0.1,0.022,txt,'FontSize',18);
% % 
% % ax2 = subplot(1,2,2); % CONT
% % a = triu(S0.Y2,1); a = a + a'; % structural covariance after removing Inf in the diagnose
% % a = [d',mean(a)']; a = a(a(:,1)<0,:);
% % scatter(ax2,a(:,1),a(:,2),45,'ko');h2 = lsline(ax2); h2.Color = 'r'; h2.LineWidth = 2;  axis([-0.11,0.01,-0.025,0.025]);
% % title('CONT'); xlabel(['Effect Size of ',SDL.data_type{1},' Difference']); ylabel('Structural Covariance');
% % [r,p] = corr(a(:,1),a(:,2));
% % txt = sprintf('R=%1.3f, p=%1.3f',r,p);
% % text(-0.1,0.022,txt,'FontSize',18);
% % 
% % savefig(fullfile(SDL.out,SDL.data_type,['Results_Corr_NegEffectSize_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));

%% End
end