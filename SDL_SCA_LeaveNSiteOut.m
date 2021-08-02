function SDL_SCA_LeaveNSiteOut(SDL)
% Leave-N-Site-Out analyses.
% CT and SA data: PTSD-vs-random, CONT-vs-random, PTSD-vs-random; not for
% severity or interactions with other variables.
% Input
% - SDL, a structure conatining the basic paths
% Parameters (in this script)
% - SDL.data_type, data type, i.e. CT or SA
% - nleft, the number of sites to be left out
% Output
% - M, a matrix of SC, 3 (network types: atrophic, hyperatrophic & stable) x (2 groups [PTSD & CONT] * top-1~148 areas, top-1=0) x nloop

%% Parameters
SDL.data_type = {
    'CT';
    'SA'
    };

nleft = 3; % number of sites to be left out

addpath(fullfile(SDL.path,'Scripts','ComBat')); % add path of combat scripts

% %% Calculation
% for j = 1:size(SDL.data_type,1) % per data type
%     dtype = SDL.data_type{j}; % data type
%
%     % Load data
%     fname = SDL.raw;
%     T0 = readtable(fname,'sheet',dtype);
%     T0 = T0(T0.Included==1,:); % remove the red rows
%     T0 = T0(T0.Gender==1 | T0.Gender==2,:); % keep only clear males & females
%     fprintf('Completed: Load<- %s\n', SDL.raw);
%
%     nsite = unique(T0.Site); % index of site (1-29 in this study)
%
%     % Leave Site(s) Out Analyses
%     if nleft == 1 % leave-one-site-out
%         nloop = length(nsite)+1; % total number of loops for calculation, N of leave-one-site-out + 1 total sites = 29 + 1
%     else % leave-sites-out
%         nloop = 5000 + 1; % % N of leave-sites-out + 1 total sites
%     end
%     M = zeros(3,148*2,nloop); % matrix to store SC, 3 (network types) x (2 groups * top-1~148 areas, top-1=0) x nloop
%
%     for i = 1:nloop % per loop
%         T = T0; % T0 for raw data (can't change), T for analyses (can change)
%         if i == nloop % including all sites
%             fprintf('\nIncluding all sites\n');
%         else % leave site(s) out
%             if nleft == 1 % leave-one-site-out
%                 SiteToLeave = i; % for leave-one-site-out analyses
%             else
%                 SiteToLeave = randperm(length(nsite));
%                 SiteToLeave = SiteToLeave(1:3); % remove rnadomly selected 3 sites
%             end
%             idx = logical(sum((T.Site==SiteToLeave),2)); % index of site(s) to be removed
%             T(idx,:) = []; % remove particular sites
%             fprintf('\nRemoving Site=');fprintf(' %02d',SiteToLeave);fprintf('\n');
%         end
%
%         T = SDL_in_Clean(T); % clean & re-organize data
%         T = SDL_in_ComBat(T);% ComBat harmonization
%         T = SDL_in_residuals(T); % Residuals after linear regression
%         d = SDL_in_RankOrder(T); % Cohen's d, PTSD vs. CONT
%         M(:,:,i) = SDL_in_SCA_TopN(T,d); % mean SC of top-n regions; col 1,2&3 = atrophic, hypertrophic & stable networs, respectively
%         fprintf('Matrix of SC: M calculated at time %s\n', datestr(now,'HH:MM:SS.FFF'))
%     end
%
%     fot = ['Results_Leave_',num2str(nleft),'_Sites_Out_',dtype,'.mat'];
%     save(fot,'M','-v7.3'); fprintf('Saved: %s\n\n\n',fot);
%
% end
%
% fprintf('\n=============Completed! Leave %02d Site(s) Out Analyses=============\n',num2str(nleft));


%% Plots
% SDL_Plot(SDL); % Plot the distributions of AUC of 5,000 leave-3-sites-out
SDL_individual_CI(SDL); % whether any value of all-sites calculation beyond CI of leave-3-sites-out calculation



%% End
end

%% CLean and re-organize data
function T = SDL_in_Clean(T)
% Rename Groups
T.Grp(T.Group==1) = {'PTSD'};
T.Grp(T.Group==2) = {'CONT'};
T.Group = T.Grp; T.Grp = [];
fprintf('Completed: Rename PTSD & CONT\n');

%% Re-organize data table into ['CONT';'PTSD']
T1 = T(strcmp('PTSD',T.Group),:);
T2 = T(strcmp('CONT',T.Group),:);
T = [T1;T2];
fprintf('Completed: Re-organize data table into PTSD and CONT\n');

% Rename Gender
T.Gen(T.Gender==1) = {'M'};
T.Gen(T.Gender==2) = {'F'};
T.Gender = T.Gen; T.Gen = [];
T = T(strcmp('M',T.Gender) | strcmp('F',T.Gender),:); % exclude the subjects missing gender values
fprintf('Completed: Re-organize values of Gender\n');

% Replace cortical thickness 0 with NaN
col_CT     = 2:149; % the columns containing cortical thickness data
T1         = T{:,col_CT};
T1(T1==0)  = NaN; % the cortical thickness data are in columns 2-149
T(:,col_CT)= array2table(T1);
fprintf('Completed: Replace all cortical values of 0s with NaN\n');

% Replace NaN with Column Mean per group
T1 = T(strcmp('PTSD',T.Group),:);
T2 = T(strcmp('CONT',T.Group),:);
for j = [col_CT,150] % per column of CT/SA (2~149), Age (150)
    T1(isnan(T1{:,j}),j) = array2table(nanmean(T1{:,j}));
    T2(isnan(T2{:,j}),j) = array2table(nanmean(T2{:,j}));
end
T = [T1;T2];
fprintf('Completed: Replace cortical thickness of NaN with Column Mean per group\n');
end


%% ComBat
function T = SDL_in_ComBat(T)%% ComBat protocal
tic
bat = T.Site; % site
nsite = unique(bat); % unique values of site
for k = 1:size(bat,1) % per site value
    batch(k) = find(nsite==bat(k)); % batch value is in fact the unique index of a particular site
end
dat = T{:,2:149}'; % values, Num_voxel x Num-sbj
age = T.Age;
T.Sex(strcmp(T.Gender,'M')) = 1; T.Sex(strcmp(T.Gender,'F')) = 2;
sex = dummyvar(T.Sex); %sex, 1=male, 2=female, dummyvar only accepts positive integer
T.Grp(strcmp(T.Group,'PTSD')) = 2; T.Grp(strcmp(T.Group,'CONT')) = 1;
disease = dummyvar(T.Grp); % diagnosis, 1=CONT, 2=PTSD, dummyvar only accepts positive integer
mod = [age sex(:,2) disease(:,2)];
% disease 1st column is control, 2nd is patients
% sex 1st column is male, 2nd is female
% dep 2nd column is subjects with depression

% data_harmonized = combat(dat, batch, mod, 0); % non-parametric adjustments
data_harmonized = combat(dat, batch, mod, 1); % parametric adjustments

T{:,2:149} = data_harmonized'; % replace the raw values in the table
fprintf('Completed: ComBat harmonization, ');toc
end



%% Residuals after linear regression
function T = SDL_in_residuals(T)
R = []; % matrix containing residuals after regression
T1 = T;
T.Gender = nominal(T.Gender);
T.mean = mean(T{:,2:149},2); % mean across areas

for j = 1:148 % per column of data values
    tbl = T(:,{'Age','Gender','mean'});
    tbl.data = T{:,j+1}; % data are in column 2-149
    lme = fitlme(tbl,'data ~ Age + Age^2 + Gender + mean');
    R(:,j) = residuals(lme);
end
T = T1;
T{:,2:149} = R; % replace raw values with residules
fprintf('Completed: Residuals after linear regression\n');
end

%% Cohen's d
function d = SDL_in_RankOrder(T)
% Cohen's d
idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT');
for j = 1:148 % per area
    d(j) = computeCohen_d(T{idx1,j+1},T{idx2,j+1}); % Cohen's d, PTSD vs. CONT
end
fprintf('Completed: Cohen''s d\n');
end

%% Mean SC of top-n regions
function M = SDL_in_SCA_TopN(T,d)
% Input
% -- T, table containing the variables of interest, col 2~149 are regional CT or SA values
% -- d, Cohen's d of PTSD vs. CONT

% ranks according to Cohen's d
kk = zeros(3,size(d,2)); % 3x148 matrix, row 1,2&3 correspond to atrophic, hypertrophic and stable networks, respectively
[~,kk(1,:)] = sort(d,'ascend');      % rank order for atrophic network, i.e. from PTSD<CONT to PTSD>CONT
[~,kk(2,:)] = sort(d,'descend');     % rank order for hypertrophic network, i.e. from PTSD>CONT to PTSD<CONT
[~,kk(3,:)] = sort(abs(d),'ascend'); % rank order for stable network, i.e. from minimal difference to larger difference
fprintf('Completed: Rank orders according to Cohen''s d\n');

R = T{:,2:149}; % residules
idx1 = strcmp(T.Group,'PTSD'); idx2 = strcmp(T.Group,'CONT'); % index of two groups

% SC per pair of regions in two groups
S0 = SDL_SCA(R,idx1,idx2,@corr,0);
% S0.Y1 --- R-to-Z transformation of corr matrix of PTSD
% S0.Y2 --- R-to-Z transformation of corr matrix of CONT
fprintf('Completed: SC per pair of regions in both groups\n');

% Mean SC of top-n regions
M = zeros(3,2*length(d)); % 3x(2x148) matrix for mean SC per top-n network; col 1&149 = 0; col 2-148, PTSD; col 150-296, CONT
for i = 2:length(d) % per top-N
    for k = 1:3 % per network type
        tic
        % original data
        node0 = kk(k,1:i); % index of nodes of top-n regions
        
        % SC per pair of regions, values above diagonal
        Z1 = triu(S0.Y1(node0,node0),1); % PTSD
        Z2 = triu(S0.Y2(node0,node0),1); % CONT
        
        % non-zero values (Num of connections x 1 vector)
        Z1c = nonzeros(Z1);
        Z2c = nonzeros(Z2);
        
        % mean SC (positive connections only !!! 1x1 value)
        M(k,i)           = mean(Z1c(Z1c>0)); % PTSD
        M(k,i+length(d)) = mean(Z2c(Z2c>0)); % CONT
    end
    fprintf('Completed: Mean SC of top-%03d regions, \n',i);toc
end


end


%% distribution of AUC of 5,000 leave-3-sites-out
function SDL_Plot(SDL)
figure;
set(gcf,'color','w'); % white background
% CT-based network
load('Results_Leave_3_Sites_Out_CT.mat'); % load M matrix
for i = 1:3 % per tyep of network, i.e. atrophic, hypertrophic, stable
    % PTSD
    a = M(i,2:148,:); % ith network, top-2~148 regions in PTSD, 5000 leave-N-sites-out + 1 all-sites
    AUC = nansum(a,2); % AUC, summed up from top-2 to top-148 networks, NaN removed
    nloop = size(AUC,3); % number of loop, 5000+1
    v = reshape(AUC,nloop,1); % extratc value per round of loop
    v1 = v;
    subplot(3,6,6*i-5); % left most column 3 panels
    histogram(v(1:nloop-1)); xlabel('AUC'); ylabel('Count');
    xline(v(nloop),'-r','LineWidth',1.5); % value of all sites
    CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
    xline(CI(1),'--b','LineWidth',1); xline(CI(2),'--b','LineWidth',1); % lines for 95% CI
    box off; % remove border around the axes
    fprintf('CT_network_type=%d, group=PTSD, [actual,average_random,CI(1),CI(2)]:\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n',...
        i,v(nloop),nanmean(v(1:nloop-1)),CI(1),CI(2));
    
    % CONT
    a = M(i,150:296,:); % ith network, top-2~148 regions in CONT, 5000 leave-2-sites-out + 1 all-sites
    AUC = nansum(a,2); % AUC, summed up from top-2 to top-148 networks, NaN removed
    nloop = size(AUC,3); % number of loop, 5000+1
    v = reshape(AUC,nloop,1); % extratc value per round of loop
    v2 = v;
    subplot(3,6,6*i-4); % left 2nd column 3 panels
    histogram(v(1:nloop-1)); xlabel('AUC'); ylabel('Count');
    xline(v(nloop),'-r','LineWidth',1.5); % value of all sites
    CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
    xline(CI(1),'--b','LineWidth',1); xline(CI(2),'--b','LineWidth',1); % lines for 95% CI
    box off; % remove border around the axes
    fprintf('CT_network_type=%d, group=CONT, [actual,average_random,CI(1),CI(2)]:\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n',...
        i,v(nloop),nanmean(v(1:nloop-1)),CI(1),CI(2));
    
    % PTSD-CONT
    v = v1-v2; % PTSD-CONT
    subplot(3,6,6*i-3); % left 2nd column 3 panels
    histogram(v(1:nloop-1)); xlabel('AUC'); ylabel('Count');
    xline(v(nloop),'-r','LineWidth',1.5); % value of all sites
    CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
    xline(CI(1),'--b','LineWidth',1); xline(CI(2),'--b','LineWidth',1); % lines for 95% CI
    box off; % remove border around the axes
    fprintf('CT_network_type=%d, group=PTSD-CONT, [actual,average_random,CI(1),CI(2)]:\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n',...
        i,v(nloop),nanmean(v(1:nloop-1)),CI(1),CI(2));
end

% SA-based network
load('Results_Leave_3_Sites_Out_SA.mat'); % load M matrix
for i = 1:3 % per tyep of network, i.e. atrophic, hypertrophic, stable
    % PTSD
    a = M(i,2:148,:); % ith network, top-2~148 regions in PTSD, 5000 leave-2-sites-out + 1 all-sites
    AUC = nansum(a,2); % AUC, summed up from top-2 to top-148 networks, NaN removed
    nloop = size(AUC,3); % number of loop, 5000+1
    v = reshape(AUC,nloop,1); % extratc value per round of loop
    v1 = v;
    subplot(3,6,6*i-2); % right 2nd column 3 panels
    histogram(v(1:nloop-1)); xlabel('AUC'); ylabel('Count');
    xline(v(nloop),'-r','LineWidth',1.5); % value of all sites
    CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
    xline(CI(1),'--b','LineWidth',1); xline(CI(2),'--b','LineWidth',1); % lines for 95% CI
    box off; % remove border around the axes
    fprintf('SA_network_type=%d, group=PTSD, [actual,average_random,CI(1),CI(2)]:\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n',...
        i,v(nloop),nanmean(v(1:nloop-1)),CI(1),CI(2));
    
    % CONT
    a = M(i,150:296,:); % ith network, top-2~148 regions in CONT, 5000 leave-2-sites-out + 1 all-sites
    AUC = nansum(a,2); % AUC, summed up from top-2 to top-148 networks, NaN removed
    nloop = size(AUC,3); % number of loop, 5000+1
    v = reshape(AUC,nloop,1); % extratc value per round of loop
    v2 = v;
    subplot(3,6,6*i-1); % right most column 3 panels
    histogram(v(1:nloop-1)); xlabel('AUC'); ylabel('Count');
    xline(v(nloop),'-r','LineWidth',1.5); % value of all sites
    CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
    xline(CI(1),'--b','LineWidth',1); xline(CI(2),'--b','LineWidth',1); % lines for 95% CI
    box off; % remove border around the axes
    fprintf('SA_network_type=%d, group=CONT, [actual,average_random,CI(1),CI(2)]:\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n',...
        i,v(nloop),nanmean(v(1:nloop-1)),CI(1),CI(2));
    
    % PTSD-CONT
    v = v1-v2; % PTSD-CONT
    subplot(3,6,6*i); % right most column 3 panels
    histogram(v(1:nloop-1)); xlabel('AUC'); ylabel('Count');
    xline(v(nloop),'-r','LineWidth',1.5); % value of all sites
    CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
    xline(CI(1),'--b','LineWidth',1); xline(CI(2),'--b','LineWidth',1); % lines for 95% CI
    box off; % remove border around the axes
    fprintf('SA_network_type=%d, group=PTSD-CONT, [actual,average_random,CI(1),CI(2)]:\t%1.4f\t%1.4f\t%1.4f\t%1.4f\n',...
        i,v(nloop),nanmean(v(1:nloop-1)),CI(1),CI(2));
end

end



%% whether any value of all-sites calculation beyond the CI of 5,000 leave-3-sites-out calculations
function SDL_individual_CI(SDL)

% CT-based network
load('Results_Leave_3_Sites_Out_CT.mat'); % load M matrix
for i = 1:3 % per tyep of network, i.e. atrophic, hypertrophic, stable
    % PTSD
    a = M(i,2:148,:); % ith network, top-2~148 regions in PTSD, 5000 leave-2-sites-out + 1 all-sites
    a1 = a;
    nloop = size(a,3); % number of loop, i.e. 5000+1
    nsize = size(a,2); % number of network size, i.e. 147 correspond to top-2~148 networks
    for j = 1:nsize % per network size
        v = reshape(a(1,j,:),nloop,1); % extratc value per round of loop, 5001x1 vector
        CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
        v0 = v(nloop); % value of all-sites calculation
        if (v0 < CI(1)) | (v0 > CI(2))
            switch i
                case 1
                    txt = 'atrophic';
                case 2
                    txt = 'hypertrophic';
                otherwise
                    txt = 'stable';
            end
            fprintf('Beyond CI !: data_type=CT,\tgroup=PTSD,\tnetwork=%s,\tsize=top-%03d\n',txt,j);
        end
    end
    
    % CONT
    a = M(i,150:296,:); % ith network, top-2~148 regions in CONT, 5000 leave-2-sites-out + 1 all-sites
    a2 = a;
    nloop = size(a,3); % number of loop, i.e. 5000+1
    nsize = size(a,2); % number of network size, i.e. 147 correspond to top-2~148 networks
    for j = 1:nsize % per network size
        v = reshape(a(1,j,:),nloop,1); % extratc value per round of loop, 5001x1 vector
        CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
        v0 = v(nloop); % value of all-sites calculation
        if (v0 < CI(1)) | (v0 > CI(2))
            switch i
                case 1
                    txt = 'atrophic';
                case 2
                    txt = 'hypertrophic';
                otherwise
                    txt = 'stable';
            end
            fprintf('Beyond CI !: data_type=CT,\tgroup=CONT,\tnetwork=%s,\tsize=top-%03d\n',txt,j);
        end
    end
    
    % PTSD-CONT
    a = a1-a2;
    nloop = size(a,3); % number of loop, i.e. 5000+1
    nsize = size(a,2); % number of network size, i.e. 147 correspond to top-2~148 networks
    for j = 1:nsize % per network size
        v = reshape(a(1,j,:),nloop,1); % extratc value per round of loop, 5001x1 vector
        CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
        v0 = v(nloop); % value of all-sites calculation
        if (v0 < CI(1)) | (v0 > CI(2))
            switch i
                case 1
                    txt = 'atrophic';
                case 2
                    txt = 'hypertrophic';
                otherwise
                    txt = 'stable';
            end
            fprintf('Beyond CI !: data_type=CT,\tgroup=PTSD-CONT,\tnetwork=%s,\tsize=top-%03d\n',txt,j);
        end
    end
    
end


% SA-based network
load('Results_Leave_3_Sites_Out_SA.mat'); % load M matrix
for i = 1:3 % per tyep of network, i.e. atrophic, hypertrophic, stable
    % PTSD
    a = M(i,2:148,:); % ith network, top-2~148 regions in PTSD, 5000 leave-2-sites-out + 1 all-sites
    a1 = a;
    nloop = size(a,3); % number of loop, i.e. 5000+1
    nsize = size(a,2); % number of network size, i.e. 147 correspond to top-2~148 networks
    for j = 1:nsize % per network size
        v = reshape(a(1,j,:),nloop,1); % extratc value per round of loop, 5001x1 vector
        CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
        v0 = v(nloop); % value of all-sites calculation
        if (v0 < CI(1)) | (v0 > CI(2))
            switch i
                case 1
                    txt = 'atrophic';
                case 2
                    txt = 'hypertrophic';
                otherwise
                    txt = 'stable';
            end
            fprintf('Beyond CI !: data_type=SA,\tgroup=PTSD,\tnetwork=%s,\tsize=top-%03d\n',txt,j);
        end
    end
    
    % CONT
    a = M(i,150:296,:); % ith network, top-2~148 regions in CONT, 5000 leave-2-sites-out + 1 all-sites
    a2 = a;
    nloop = size(a,3); % number of loop, i.e. 5000+1
    nsize = size(a,2); % number of network size, i.e. 147 correspond to top-2~148 networks
    for j = 1:nsize % per network size
        v = reshape(a(1,j,:),nloop,1); % extratc value per round of loop, 5001x1 vector
        CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
        v0 = v(nloop); % value of all-sites calculation
        if (v0 < CI(1)) | (v0 > CI(2))
            switch i
                case 1
                    txt = 'atrophic';
                case 2
                    txt = 'hypertrophic';
                otherwise
                    txt = 'stable';
            end
            fprintf('Beyond CI !: data_type=SA,\tgroup=CONT,\tnetwork=%s,\tsize=top-%03d\n',txt,j);
        end
    end
    
    % PTSD-CONT
    a = a1-a2;
    nloop = size(a,3); % number of loop, i.e. 5000+1
    nsize = size(a,2); % number of network size, i.e. 147 correspond to top-2~148 networks
    for j = 1:nsize % per network size
        v = reshape(a(1,j,:),nloop,1); % extratc value per round of loop, 5001x1 vector
        CI = prctile(v(1:nloop-1),[2.5 97.5]); % 95% CI (in fact 95% percentile)
        v0 = v(nloop); % value of all-sites calculation
        if (v0 < CI(1)) | (v0 > CI(2))
            switch i
                case 1
                    txt = 'atrophic';
                case 2
                    txt = 'hypertrophic';
                otherwise
                    txt = 'stable';
            end
            fprintf('Beyond CI !: data_type=SA,\tgroup=PTSD-CONT,\tnetwork=%s,\tsize=top-%03d\n',txt,j);
        end
    end
    
end

end