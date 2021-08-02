
% Mean SC values of top-N areas with CT/SA reduction
% Input data
% S0 --- SC original data
% S1 --- SC permutation data
% Y1 --- group 1
% Y2 --- group2


%% Load data
fn = fullfile(SDL.path,'Outputs','CleanData',['Data_SC_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn);
fprintf('Completed: Structural covariance matrix loaded from <- %s\n\n\n',fn);

%% Mean SC values of top-N areas with CT/SA reduction
SC0 = []; SC1 = []; % mean SC values of top-N, original(0) & permutation(1) data
% SC0.Y1 = zeros(1,size(S0.Y1,1));       SC0.Y2 = zeros(1,size(S0.Y2,1));      % original data for group 1&2
% SC1.Y1 = zeros(1,size(S1.Y1,1),SDL.N); SC1.Y2 = zeros(1,size(S1.Y2,1),SDL.N);% permutations for group 1&2
% ST0 = [];
for i = 2:size(S0.Y1,1) % per top-N
    j  = kk(1:i); % top-N areas
    jp = randperm(size(S0.Y1,1)); jp = jp(1:i); % randperm areas
    
    SC0.Y1 = S0.Y1(j,j); SC0.Y2 = S0.Y2(j,j); % original data, SC values of the top-N areas,group 1 & 2
    SC1.Y1 = S1.Y1(jp,jp,:); SC1.Y2 = S1.Y2(jp,jp,:); % permutation data, SC values of the top-N areas,group 1 & 2
    
    % original data
    SC0.Z1 = triu(SC0.Y1,1); SC0.Z2 = triu(SC0.Y2,1); % original data, values above diagonal
    SC0.Z1c = nonzeros(SC0.Z1(:)); SC0.Z2c = nonzeros(SC0.Z2(:)); % non-zero values
    SC0.mean1(i) = mean(SC0.Z1c); % mean of nonzeros, original data, group 1
    SC0.mean2(i) = mean(SC0.Z2c); % mean of nonzeros, original data, group 2
    
    % permutaton data
    SC1.Z1 = []; SC1.Z2 = []; SC1.Z1c = []; SC1.Z2c = [];
    for k = 1:SDL.N % per permutation
        % permutation data, values above diagonal
        SC1.Z1(:,:,k) = triu(SC1.Y1(:,:,k),1);
        SC1.Z1c(:,k)  = nonzeros(SC1.Z1(:,:,k)); % non-zero values
        SC1.Z2(:,:,k) = triu(SC1.Y2(:,:,k),1);
        SC1.Z2c(:,k)  = nonzeros(SC1.Z2(:,:,k)); % non-zero values
    end
    SC1.mean1(i) = mean(mean(SC1.Z1c)); % mean of nonzeros, permutation data, group 1
    SC1.mean2(i) = mean(mean(SC1.Z2c)); % mean of nonzeros, permutation data, group 2
    SC1.diff(i)  = mean(mean(SC1.Z1c)) - mean(mean(SC1.Z2c)); % difference between group 1&2
    
    CIType = 0.95; % 95% confidene interval (CI) 
    SC1.CI1(:,i) = SDL_CI(mean(SC1.Z1c,1)',CIType); % 95% CI
    SC1.CI2(:,i) = SDL_CI(mean(SC1.Z2c,1)',CIType); % 95% CI
    SC1.diffCI(:,i) = SDL_CI(mean(SC1.Z1c,1)' - mean(SC1.Z2c,1)',0.95); % 95% CI
    
    SC0.p1(i)   = SDL_p_permutation(SC0.mean1(i),mean(SC1.Z1c)); % group1 vs permutation
    SC0.p2(i)   = SDL_p_permutation(SC0.mean2(i),mean(SC1.Z2c)); % group2 vs permutation
    SC0.p_diff(i) = SDL_p_permutation(SC0.mean1(i)-SC0.mean2(i),mean(SC1.Z1c)-mean(SC1.Z2c)); % group1-2 vs permutation
    
    fprintf('Completed: SC for Top-%d Areas\n',i);
end


%% Plot SC values of Top-N areas with reduction in PTSD, CONT and PTSD-CONT 
X = 2:size(S0.Y1,1);

figure; % PTSD vs CONT
plot(X,SC0.mean1(X),'r-',X,SC0.mean2(X),'b-','LineWidth',2);
xlabel('Nmber of Areas'); ylabel('Structural Covariance');legend('PTSD','CONT');
savefig(fullfile(SDL.path,'Outputs','CleanData',['Results_Line_PTSD_vs_CONT',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));

figure; % original vs permutation
subplot(3,1,1);
X = 3:size(S0.Y1,1);x2 = [X, fliplr(X)]; inBetween = [SC1.CI1(1,X)+SC1.mean1(X), fliplr(SC1.CI1(2,X)+SC1.mean1(X))];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
X = 2:size(S0.Y1,1);plot(X,SC0.mean1(X),'r-',X,SC1.mean1(X),'b-','LineWidth',2);
xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title('PTSD');
subplot(3,1,2);
X = 3:size(S0.Y1,1);x2 = [X, fliplr(X)]; inBetween = [SC1.CI2(1,X)+SC1.mean2(X), fliplr(SC1.CI2(2,X)+SC1.mean2(X))];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
X = 2:size(S0.Y1,1);plot(X,SC0.mean2(X),'r-',X,SC1.mean2(X),'b-','LineWidth',2);
xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title('CONT');
subplot(3,1,3);
X = 3:size(S0.Y1,1);x2 = [X, fliplr(X)]; inBetween = [SC1.diffCI(1,X)+SC1.diff(X), fliplr(SC1.diffCI(2,X)+SC1.diff(X))];fill(x2, inBetween, [91, 207, 244] / 255); hold on;
X = 2:size(S0.Y1,1);plot(X,(SC0.mean1(X)-SC0.mean2(X)),'r-',X,SC1.diff(X),'b-','LineWidth',2);
xlabel('Number of Areas'); ylabel('Structural Covariance');legend('95% CI','Areas with Reduction','Random Areas'); title('PTSD-CONT');
savefig(fullfile(SDL.path,'Outputs','CleanData',['Results_Line_Original_vs_Random',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));

%% Plot heatmap of SC values
MN1 = []; MN2 = []; MNd = [];
% re-organize matrix based on kk
for i = 1:length(kk)
    for j = 1:length(kk)
        MN1(i,j) = S0.Y1(kk(i),kk(j));
        MN2(i,j) = S0.Y2(kk(i),kk(j));
    end
end
MNd = MN1 - MN2;
figure; 
subplot(1,3,1);pcolor(MN1); colormap('hot'); colorbar; title('PTSD'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([0,1.2]);
subplot(1,3,2);pcolor(MN2); colormap('hot'); colorbar; title('CONT'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([0,1.2]);
subplot(1,3,3);pcolor(MNd); colormap('hot'); colorbar; title('PTSD-CONT'); xlabel('Top-N Areas'); ylabel('Top-N Areas');caxis([-0.1,0.1]);
savefig(fullfile(SDL.path,'Outputs','CleanData',['Results_2Dmap_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));

% % calculate CI (CI = mean(x)+- t * (s / square(n)))
% function CI = SDL_CI(y,CIType)
% % Input
% %    y      --- data
% %    CIType --- 0.95 or 0.99
% % Output
% %    CI     --- confidence interval
% N = size(y,1); % Number of ?Experiments? In Data Set
% yMean = mean(y); % Mean Of All Experiments At Each Value Of ?x?
% ySEM = std(y)/sqrt(N); % Compute ?Standard Error Of The Mean? Of All Experiments At Each Value Of ?x? . !!!??? May be Wrong!!!
% 
% a = (1-CIType)/2;
% ts = tinv([a 1-a], N-1);          % Calculate Probability Intervals Of t-Distribution
% CI = bsxfun(@times, ySEM, ts(:)); % Calculate Confidence Intervals Of All Experiments At Each Value Of ?x?
% end
% 
% 
% % calculate p values of permutation
% function p = SDL_p_permutation(v0,v1)
% % Input
% %    v0 --- value of original data, 1x1
% %    v1 --- value of permutations, 1xN
% % Output
% %    p  --- statistical significance
% 
% p = sum((v1-v0)>0)/length(v1);
% end