function S = SDL_SCA(R,idx1,idx2,func,N)
% To calculate the mean SC and R-to-Z transformation
% Input
% R    --- Residules matrix after lme to regress out nusiance regressors
% idx1 --- row index of group1, i.e. PTSD
% idx2 --- row index of group2, i.e. CONT
% func --- matlab functions for calculation, i.e. @corr or @partialcorr
% N    --- number of permutations, no permutation if 0
% Output
% S    --- a structure containing
% S.M1 --- corr/partialcorr matrix of group1
% S.M2 --- corr/partialcorr matrix of group2
% S.Y1 --- R-to-Z transformation of corr/partialcorr matrix of group1
% S.Y2 --- R-to-Z transformation of corr/partialcorr matrix of group2

S = [];
if N==0 % no permutation (for original data only)
    %% correlation or partial correlation
    S.M1 = func(R(idx1,:));
    S.M2 = func(R(idx2,:));
    %% R-to-Z transformation, i.e. 0.5*log((1+x)/(1-x))
    S.Y1 = atanh(S.M1);
    S.Y2 = atanh(S.M2);
    
%     get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
%     figure;
%     subplot(2,2,1);histogram(S.M1);  xlabel('Coeff.'); ylabel('Count');title('Group = PTSD');
%     subplot(2,2,2);histogram(S.M2);  xlabel('Coeff.'); ylabel('Count');title('Group = CONT');
%     subplot(2,2,3);histogram(S.Y1);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
%     subplot(2,2,4);histogram(S.Y2);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
    fprintf('Completed: SC matrix and R-to-Z transformation, N=0\n\n');
else % permutation
    idx0 = [idx1;idx2]; % shuffling group labels 
    for i = 1:N
        %% Shuffle labels
        idx = idx0(randperm(length(idx0))); % shuffle the labels of two groups
        %% correlation or partial correlation
        S.M1(:,:,i) = func(R(idx(idx1),:));
        S.M2(:,:,i) = func(R(idx(idx2),:));
        %% R-to-Z transformation, i.e. 0.5*log((1+x)/(1-x))
        S.Y1 = atanh(S.M1);
        S.Y2 = atanh(S.M2);
        fprintf('Completed: SC matrix and R-to-Z transformation, N=%4d\n',i);
    end
end
end