function SDL_SCA_TopN_Age10(SDL)

% SCA within the top-N areas with CT or SA reduction


%% Loading information
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Residules <- %s\n\n\n',fn);
R = T{:,2:149}; % residules

a = prctile(T.Age,[20 40 60 80]); % divide into 5 groups with compariable sample size
[sum(T.Age<a(1)),sum(T.Age>=a(1)&T.Age<a(2)),sum(T.Age>=a(2)&T.Age<a(3)),sum(T.Age>=a(3)&T.Age<a(4)),sum(T.Age>=a(4))]


idx11 = strcmp(T.Group,'PTSD') & (T.Age<a(1));               % group 11, PTSD & 0<=Age<23
idx12 = strcmp(T.Group,'PTSD') & (T.Age>=a(1) & T.Age<a(2)); % group 12, PTSD & 23<=Age<30
idx13 = strcmp(T.Group,'PTSD') & (T.Age>=a(2) & T.Age<a(3)); % group 13, PTSD & 30<=Age<38
idx14 = strcmp(T.Group,'PTSD') & (T.Age>=a(3) & T.Age<a(4)); % group 14, PTSD & 38<=Age<49
idx15 = strcmp(T.Group,'PTSD') & (T.Age>=a(4));              % group 15, PTSD & 49<=Age<199

idx21 = strcmp(T.Group,'CONT') & (T.Age<a(1));               % group 21, CONT & 0<=Age<23
idx22 = strcmp(T.Group,'CONT') & (T.Age>=a(1) & T.Age<a(2)); % group 22, CONT & 23<=Age<30
idx23 = strcmp(T.Group,'CONT') & (T.Age>=a(2) & T.Age<a(3)); % group 23, CONT & 30<=Age<38
idx24 = strcmp(T.Group,'CONT') & (T.Age>=a(3) & T.Age<a(4)); % group 24, CONT & 38<=Age<49
idx25 = strcmp(T.Group,'CONT') & (T.Age>=a(4));              % group 25, CONT & 49<=Age<199

fn = fullfile(fdir,['Results_EffectSize_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: Effect Size <- %s\n\n\n',fn);
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
S0 = SDL_SCA10(R,idx11,idx12,idx13,idx14,idx15,idx16,idx21,idx22,idx23,idx24,idx25,idx26,func,0);     % SC of original PTSD & CONT data

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
% fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
% fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% save(fn,'S0','SC0','SC1','-v7.3');
% fprintf('Saved: Top-N analyses midway results saved in \n-->%s\n',fn);

%% SC, CI & p values of top-N areas with CT/SA reduction
% loading
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}(1:2)); % load SC0 and SC1 from existed results, i.e. CleanData, to avoid wasting time
fn = fullfile(fdir,['Results_TopN_PTSD_vs_CONT_',SDL.data_type{1}(1:2),'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');
fprintf('Loaded: Top-N analyses midway results loaded from \n<--%s\n',fn);

% mean SC, CI and p values
CIType = 0.95; % 95% confidene interval (CI) 
for i = 2:size(SC0,2) % per top-N
    % original data, values above diagonal
    node0 = SC0(i).node;
    SC0(i).Z11 = triu(S0.Y11(node0,node0),1); 
    SC0(i).Z12 = triu(S0.Y12(node0,node0),1); 
    SC0(i).Z13 = triu(S0.Y13(node0,node0),1);
    SC0(i).Z14 = triu(S0.Y14(node0,node0),1);
    SC0(i).Z15 = triu(S0.Y15(node0,node0),1);
    SC0(i).Z16 = triu(S0.Y16(node0,node0),1);
    
    SC0(i).Z21 = triu(S0.Y21(node0,node0),1); 
    SC0(i).Z22 = triu(S0.Y22(node0,node0),1);
    SC0(i).Z23 = triu(S0.Y23(node0,node0),1);
    SC0(i).Z24 = triu(S0.Y24(node0,node0),1);
    SC0(i).Z25 = triu(S0.Y25(node0,node0),1);
    SC0(i).Z26 = triu(S0.Y26(node0,node0),1);
    
    % non-zero values, Num of connections X 1
    SC0(i).Z11c = nonzeros(SC0(i).Z11); 
    SC0(i).Z12c = nonzeros(SC0(i).Z12); 
    SC0(i).Z13c = nonzeros(SC0(i).Z13);
    SC0(i).Z14c = nonzeros(SC0(i).Z14);
    SC0(i).Z15c = nonzeros(SC0(i).Z15);
    SC0(i).Z16c = nonzeros(SC0(i).Z16);
    
    SC0(i).Z21c = nonzeros(SC0(i).Z21); 
    SC0(i).Z22c = nonzeros(SC0(i).Z22);% 
    SC0(i).Z23c = nonzeros(SC0(i).Z23);
    SC0(i).Z24c = nonzeros(SC0(i).Z24);
    SC0(i).Z25c = nonzeros(SC0(i).Z25);
    SC0(i).Z26c = nonzeros(SC0(i).Z26);
    
    % mean of all connections, 1X1, original data
    SC0(i).mean11 = mean(SC0(i).Z11c); % mean of all connections, 1X1, original data, group 11
    SC0(i).mean12 = mean(SC0(i).Z12c); % mean of all connections, 1X1, original data, group 12
    SC0(i).mean13 = mean(SC0(i).Z13c);
    SC0(i).mean14 = mean(SC0(i).Z14c);
    SC0(i).mean15 = mean(SC0(i).Z15c);
    SC0(i).mean16 = mean(SC0(i).Z16c);
    
    SC0(i).mean21 = mean(SC0(i).Z21c); % mean of all connections, 1X1, original data, group 21
    SC0(i).mean22 = mean(SC0(i).Z22c); % mean of all connections, 1X1, original data, group 22
    SC0(i).mean23 = mean(SC0(i).Z23c);
    SC0(i).mean24 = mean(SC0(i).Z24c);
    SC0(i).mean25 = mean(SC0(i).Z25c);
    SC0(i).mean26 = mean(SC0(i).Z26c);
    
    
    SC0(i).Z11 = []; 
    SC0(i).Z12 = []; 
    SC0(i).Z13 = [];
    SC0(i).Z14 = [];
    SC0(i).Z15 = [];
    SC0(i).Z16 = [];
    
    SC0(i).Z21 = []; 
    SC0(i).Z22 = [];
    SC0(i).Z23 = []; 
    SC0(i).Z24 = []; 
    SC0(i).Z25 = []; 
    SC0(i).Z26 = []; 
    
    SC0(i).Z11c = []; 
    SC0(i).Z12c = [];
    SC0(i).Z13c = [];
    SC0(i).Z14c = [];
    SC0(i).Z15c = [];
    SC0(i).Z16c = [];
    
    SC0(i).Z21c = []; 
    SC0(i).Z22c = [];
    SC0(i).Z23c = [];
    SC0(i).Z24c = [];
    SC0(i).Z25c = [];
    SC0(i).Z26c = [];
    
    % permuted data
    for j = 1:size(SC1(i).node,1) % per permutation (totally 5000 times)
        node1 = SC1(i).node(j,:);
        SC1(i).Z11(:,:,j) = triu(S0.Y11(node1,node1),1); 
        SC1(i).Z12(:,:,j) = triu(S0.Y12(node1,node1),1); 
        SC1(i).Z13(:,:,j) = triu(S0.Y13(node1,node1),1);
        SC1(i).Z14(:,:,j) = triu(S0.Y14(node1,node1),1);
        SC1(i).Z15(:,:,j) = triu(S0.Y15(node1,node1),1);
        SC1(i).Z16(:,:,j) = triu(S0.Y16(node1,node1),1);
        
        SC1(i).Z21(:,:,j) = triu(S0.Y21(node1,node1),1); 
        SC1(i).Z22(:,:,j) = triu(S0.Y22(node1,node1),1);
        SC1(i).Z23(:,:,j) = triu(S0.Y23(node1,node1),1);
        SC1(i).Z24(:,:,j) = triu(S0.Y24(node1,node1),1);
        SC1(i).Z25(:,:,j) = triu(S0.Y25(node1,node1),1);
        SC1(i).Z26(:,:,j) = triu(S0.Y26(node1,node1),1);
        
        SC1(i).Z11c(:,j) = nonzeros(SC1(i).Z11(:,:,j)); 
        SC1(i).Z12c(:,j) = nonzeros(SC1(i).Z12(:,:,j)); 
        SC1(i).Z13c(:,j) = nonzeros(SC1(i).Z13(:,:,j));
        SC1(i).Z14c(:,j) = nonzeros(SC1(i).Z14(:,:,j));
        SC1(i).Z15c(:,j) = nonzeros(SC1(i).Z15(:,:,j));
        SC1(i).Z16c(:,j) = nonzeros(SC1(i).Z16(:,:,j));
        
        SC1(i).Z21c(:,j) = nonzeros(SC1(i).Z21(:,:,j)); 
        SC1(i).Z22c(:,j) = nonzeros(SC1(i).Z22(:,:,j));% non-zero values, Num of connections X Num of permutations
        SC1(i).Z23c(:,j) = nonzeros(SC1(i).Z23(:,:,j));
        SC1(i).Z24c(:,j) = nonzeros(SC1(i).Z24(:,:,j));
        SC1(i).Z25c(:,j) = nonzeros(SC1(i).Z25(:,:,j));
        SC1(i).Z26c(:,j) = nonzeros(SC1(i).Z26(:,:,j));
    end
    SC1(i).mean11 = mean(SC1(i).Z11c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 11
    SC1(i).mean12 = mean(SC1(i).Z12c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 12
    SC1(i).mean13 = mean(SC1(i).Z13c,1);
    SC1(i).mean14 = mean(SC1(i).Z14c,1);
    SC1(i).mean15 = mean(SC1(i).Z15c,1);
    SC1(i).mean16 = mean(SC1(i).Z16c,1);
    
    SC1(i).mean21 = mean(SC1(i).Z21c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 21
    SC1(i).mean22 = mean(SC1(i).Z22c,1); % mean of all connections per permutation, 1 X Num of permutations, random data, group 22
    SC1(i).mean23 = mean(SC1(i).Z23c,1);
    SC1(i).mean24 = mean(SC1(i).Z24c,1);
    SC1(i).mean25 = mean(SC1(i).Z25c,1);
    SC1(i).mean26 = mean(SC1(i).Z26c,1);
    
    SC1(i).Z11 = []; 
    SC1(i).Z12 = [];
    SC1(i).Z13 = [];
    SC1(i).Z14 = [];
    SC1(i).Z15 = [];
    SC1(i).Z16 = [];
    
    SC1(i).Z21 = []; 
    SC1(i).Z22 = [];
    SC1(i).Z23 = [];
    SC1(i).Z24 = [];
    SC1(i).Z25 = [];
    SC1(i).Z26 = [];
    
    SC1(i).Z11c = []; 
    SC1(i).Z12c = [];
    SC1(i).Z13c = [];
    SC1(i).Z14c = [];
    SC1(i).Z15c = [];
    SC1(i).Z16c = [];
    
    SC1(i).Z21c = []; 
    SC1(i).Z22c = [];
    SC1(i).Z23c = [];
    SC1(i).Z24c = [];
    SC1(i).Z25c = [];
    SC1(i).Z26c = [];
    
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
    % p values
%     SC1(i).p11     = SDL_p_permutation(SC0(i).mean11,SC1(i).mean11); % group1 vs permutation
%     SC1(i).p12     = SDL_p_permutation(SC0(i).mean12,SC1(i).mean12); % group2 vs permutation
%     SC1(i).p21     = SDL_p_permutation(SC0(i).mean21,SC1(i).mean21); % group1 vs permutation
%     SC1(i).p22     = SDL_p_permutation(SC0(i).mean22,SC1(i).mean22); % group2 vs permutation
SC1(i).p_diff_11v21 = SDL_p_permutation(SC0(i).mean11-SC0(i).mean21,SC1(i).mean11-SC1(i).mean21); % group 11-21 vs permutation
SC1(i).p_diff_12v22 = SDL_p_permutation(SC0(i).mean12-SC0(i).mean22,SC1(i).mean12-SC1(i).mean22); % group 12-22 vs permutation
SC1(i).p_diff_13v23 = SDL_p_permutation(SC0(i).mean13-SC0(i).mean23,SC1(i).mean13-SC1(i).mean23); % group 13-23 vs permutation
SC1(i).p_diff_14v24 = SDL_p_permutation(SC0(i).mean14-SC0(i).mean24,SC1(i).mean14-SC1(i).mean24); % group 14-24 vs permutation
SC1(i).p_diff_15v25 = SDL_p_permutation(SC0(i).mean15-SC0(i).mean25,SC1(i).mean15-SC1(i).mean25); % group 15-25 vs permutation
SC1(i).p_diff_16v26 = SDL_p_permutation(SC0(i).mean16-SC0(i).mean26,SC1(i).mean16-SC1(i).mean26); % group 16-26 vs permutation

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
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
save(fn,'SC0','SC1','-v7.3');
fprintf('Saved: Top-N analyses results saved in \n-->%s\n',fn);


%% Plot SC values of Top-N areas with reduction in PTSD, CONT and PTSD-CONT 
fdir = fullfile(SDL.path,'Outputs',SDL.data_type{1}); 
fn = fullfile(fdir,['Results_TopN_SC_CI_p_PTSD_vs_CONT_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
load(fn,'SC0','SC1');

X = 2:size(SC0,1);
for i = 2:size(SC0,2)
    
    X11(i-1) = SC0(i).mean11; 
    X12(i-1) = SC0(i).mean12; 
    X13(i-1) = SC0(i).mean13; 
    X14(i-1) = SC0(i).mean14; 
    X15(i-1) = SC0(i).mean15; 
    X16(i-1) = SC0(i).mean16; 
    
    X21(i-1) = SC0(i).mean21; 
    X22(i-1) = SC0(i).mean22;
    X23(i-1) = SC0(i).mean23;
    X24(i-1) = SC0(i).mean24;
    X25(i-1) = SC0(i).mean25;
    X26(i-1) = SC0(i).mean26;
    
%     p11(i-1) = SC1(i).p11; p12(i-1) = SC1(i).p12; p21(i-1) = SC1(i).p21; p22(i-1) = SC1(i).p22;
    
    X11R(i-1) = mean(SC1(i).mean11); 
    X12R(i-1) = mean(SC1(i).mean12); 
    X13R(i-1) = mean(SC1(i).mean13);
    X14R(i-1) = mean(SC1(i).mean14);
    X15R(i-1) = mean(SC1(i).mean15);
    X16R(i-1) = mean(SC1(i).mean16);
    
    X21R(i-1) = mean(SC1(i).mean21); 
    X22R(i-1) = mean(SC1(i).mean22);
    X23R(i-1) = mean(SC1(i).mean23);
    X24R(i-1) = mean(SC1(i).mean24);
    X25R(i-1) = mean(SC1(i).mean25);
    X26R(i-1) = mean(SC1(i).mean26);
    
    
    
    X11RA(i-1,:) = SC1(i).mean11; 
    X12RA(i-1,:) = SC1(i).mean12; 
    X13RA(i-1,:) = SC1(i).mean13; 
    X14RA(i-1,:) = SC1(i).mean14; 
    X15RA(i-1,:) = SC1(i).mean15; 
    X16RA(i-1,:) = SC1(i).mean16; 
    
    X21RA(i-1,:) = SC1(i).mean21; 
    X22RA(i-1,:) = SC1(i).mean22;
    X23RA(i-1,:) = SC1(i).mean23;
    X24RA(i-1,:) = SC1(i).mean24;
    X25RA(i-1,:) = SC1(i).mean25;
    X26RA(i-1,:) = SC1(i).mean26;
    
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
    p_diff_11v21(i-1) = SC1(i).p_diff_11v21; 
    p_diff_12v22(i-1) = SC1(i).p_diff_12v22;
    p_diff_13v23(i-1) = SC1(i).p_diff_13v23;
    p_diff_14v24(i-1) = SC1(i).p_diff_14v24;
    p_diff_15v25(i-1) = SC1(i).p_diff_15v25;
    p_diff_16v26(i-1) = SC1(i).p_diff_16v26;
%     p_diff_11v12(i-1) = SC1(i).p_diff_11v12;
%     p_diff_21v22(i-1) = SC1(i).p_diff_21v22;
%     p_diff_main1(i-1) = SC1(i).p_diff_main1;
%     p_diff_main2(i-1) = SC1(i).p_diff_main2;
%     p_diff_inter(i-1) = SC1(i).p_diff_inter;
end


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
fprintf('\n\nFDR correction: group 11v21 vs random\n')
[~, ~, ~, adj_p_diff_11v21] = fdr_bh(p_diff_11v21,0.05,'pdep','yes'); 
find(adj_p_diff_11v21<=0.05)
fprintf('\n\nFDR correction: group 12v22 vs random\n')
[~, ~, ~, adj_p_diff_12v22] = fdr_bh(p_diff_12v22,0.05,'pdep','yes'); 
find(adj_p_diff_12v22<=0.05)
fprintf('\n\nFDR correction: group 13v23 vs random\n')
[~, ~, ~, adj_p_diff_13v23] = fdr_bh(p_diff_13v23,0.05,'pdep','yes'); 
find(adj_p_diff_13v23<=0.05)
fprintf('\n\nFDR correction: group 14v24 vs random\n')
[~, ~, ~, adj_p_diff_14v24] = fdr_bh(p_diff_14v24,0.05,'pdep','yes'); 
find(adj_p_diff_14v24<=0.05)
fprintf('\n\nFDR correction: group 15v25 vs random\n')
[~, ~, ~, adj_p_diff_15v25] = fdr_bh(p_diff_15v25,0.05,'pdep','yes'); 
find(adj_p_diff_15v25<=0.05)
fprintf('\n\nFDR correction: group 16v26 vs random\n')
[~, ~, ~, adj_p_diff_16v26] = fdr_bh(p_diff_16v26,0.05,'pdep','yes'); 
find(adj_p_diff_16v26<=0.05)
% fprintf('\n\nFDR correction: group 11v12 vs random\n')
% [~, ~, ~, adj_p_diff_11v12] = fdr_bh(p_diff_11v12,0.05,'pdep','yes'); 
% find(adj_p_diff_11v12<=0.05)
% fprintf('\n\nFDR correction: group 21v22 vs random\n')
% [~, ~, ~, adj_p_diff_21v22] = fdr_bh(p_diff_21v22,0.05,'pdep','yes'); 
% find(adj_p_diff_21v22<=0.05)
% fprintf('\n\nFDR correction: main1 (PTSD vs CONT) vs random\n')
% [~, ~, ~, adj_p_diff_main1] = fdr_bh(p_diff_main1,0.05,'pdep','yes'); 
% find(adj_p_diff_main1<=0.05)
% fprintf('\n\nFDR correction: main2 vs random\n')
% [~, ~, ~, adj_p_diff_main2] = fdr_bh(p_diff_main2,0.05,'pdep','yes'); 
% find(adj_p_diff_main2<=0.05)
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
p_diff_11v21a = SDL_p_permutation(mean(X11)-mean(X21),mean(X11RA,1)-mean(X21RA,1))
p_diff_12v22a = SDL_p_permutation(mean(X12)-mean(X22),mean(X12RA,1)-mean(X22RA,1))
p_diff_13v23a = SDL_p_permutation(mean(X13)-mean(X23),mean(X13RA,1)-mean(X23RA,1))
p_diff_14v24a = SDL_p_permutation(mean(X14)-mean(X24),mean(X14RA,1)-mean(X24RA,1))
p_diff_15v25a = SDL_p_permutation(mean(X15)-mean(X25),mean(X15RA,1)-mean(X25RA,1))
p_diff_16v26a = SDL_p_permutation(mean(X16)-mean(X26),mean(X16RA,1)-mean(X26RA,1))
% p_diff_11v12a = SDL_p_permutation(mean(X11)-mean(X12),mean(X11RA,1)-mean(X12RA,1))
% p_diff_21v22a = SDL_p_permutation(mean(X21)-mean(X22),mean(X21RA,1)-mean(X22RA,1))
% p_diff_main1a = SDL_p_permutation(mean(X11)+mean(X12)-mean(X21)-mean(X22),...
%     mean(X11RA,1)+mean(X12RA,1)-mean(X21RA,1)-mean(X22RA,1))
% p_diff_main2a = SDL_p_permutation(mean(X11)-mean(X12)+mean(X21)-mean(X22),...
%     mean(X11RA,1)-mean(X12RA,1)+mean(X21RA,1)-mean(X22RA,1))
% p_diff_intera = SDL_p_permutation(mean(X11)-mean(X12)-mean(X21)+mean(X22),...
%     mean(X11RA,1)-mean(X12RA,1)-mean(X21RA,1)+mean(X22RA,1))


% %% Plot curves & heatmap of SC values
% if strcmp(SDL.data_type{1}(4:end),'Age')
%     txt1 = 'Young'; txt2 = 'Old';
% elseif strcmp(SDL.data_type{1}(4:end),'Gender')
%     txt1 = 'Male'; txt2 = 'Female';
% elseif strcmp(SDL.data_type{1}(4:end),'Dep')
%     txt1 = 'Depress'; txt2 = 'No-Depress';
% else
% end

get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure; % original vs permutation
L = 1:length(X11);
plot(L,(X11-X11R)-(X21-X21R),...
    L,(X12-X12R)-(X22-X22R),...
    L,(X13-X13R)-(X23-X23R),...
    L,(X14-X14R)-(X24-X24R),...
    L,(X15-X15R)-(X25-X25R),...
    L,(X16-X16R)-(X26-X26R),'LineWidth', 2);
legend('<20','<30&>=20','<40&>=30','<50&>=40','<60&>=50','>=60');
xlabel('Top-N Regions'); ylabel('Structural Covariance');title('(PTSD-CONT)Actual - (PTSD-CONT)Random')


% b = mean([(X1-X1R)',...
%     (X2-X2R)',...
%     (X3-X3R)',...
%     (X4-X4R)',...
%     (X5-X5R)'])';
% bar(0,b);
% xlabel('PTSD Severity','FontName', 'Arial','FontSize',16); 
% ylabel('AUC','FontName', 'Arial','FontSize',16);
% xticklabels({' '});
% box off




% [x,y]=meshgrid(6,147);
% z=[ (X11-X11R)-(X21-X21R);...
%     (X12-X12R)-(X22-X22R);...
%     (X13-X13R)-(X23-X23R);...
%     (X14-X14R)-(X24-X24R);...
%     (X15-X15R)-(X26-X26R);...
%     (X16-X16R)-(X26-X26R)];
% mesh(1:x,1:y,z)
% colormap(jet)
% title('Default colors for mesh BEFORE 2014b')
% 
% 
% imagesc(z); colormap('hot'); colorbar;
% xlabel('Top-N Regions'); ylabel('Age');title('Structural Covariance of (PTSD-CONT)Actual - (PTSD-CONT)Random')


savefig(fullfile(SDL.path,'Outputs',SDL.data_type{1},['Results_Curves_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));


%% End
end