function S = SDL_SCA4(R,idx11,idx12,idx21,idx22,func,N)
% To calculate the mean SC and R-to-Z transformation
% Input
% R    --- Residules matrix after lme to regress out nusiance regressors
% idx11 --- row index of group11, i.e. PTSD & e.g. Age<33
% idx12 --- row index of group12, i.e. PTSD & e.g. Age>=33
% idx21 --- row index of group21, i.e. CONT & e.g. Age<33
% idx22 --- row index of group22, i.e. CONT & e.g. Age>=33
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
    S.M11 = func(R(idx11,:));
    S.M12 = func(R(idx12,:));
    S.M21 = func(R(idx21,:));
    S.M22 = func(R(idx22,:));
    %% R-to-Z transformation, i.e. 0.5*log((1+x)/(1-x))
    S.Y11 = atanh(S.M11);
    S.Y12 = atanh(S.M12);
    S.Y21 = atanh(S.M21);
    S.Y22 = atanh(S.M22);
    
%     get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
%     figure;
%     subplot(2,4,1);histogram(S.M11);  xlabel('Coeff.'); ylabel('Count');title('Group = PTSD &1');
%     subplot(2,4,2);histogram(S.M12);  xlabel('Coeff.'); ylabel('Count');title('Group = PTSD &2');
%     subplot(2,4,3);histogram(S.M21);  xlabel('Coeff.'); ylabel('Count');title('Group = CONT &1');
%     subplot(2,4,4);histogram(S.M22);  xlabel('Coeff.'); ylabel('Count');title('Group = CONT &2');
%     subplot(2,4,5);histogram(S.Y11);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
%     subplot(2,4,6);histogram(S.Y12);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
%     subplot(2,4,7);histogram(S.Y21);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
%     subplot(2,4,8);histogram(S.Y22);  xlabel('R-to-Z Trans. of Coeff.'); ylabel('Count');
    fprintf('Completed: SC matrix and R-to-Z transformation, N=0\n\n');
else % permutation
    idx0 = [idx11;idx12;idx21;idx22]; % shuffling group labels 
    for i = 1:N
        %% Shuffle labels
        idx = idx0(randperm(length(idx0))); % shuffle the labels of two groups
        %% correlation or partial correlation
        S.M11(:,:,i) = func(R(idx(idx11),:));
        S.M12(:,:,i) = func(R(idx(idx12),:));
        S.M21(:,:,i) = func(R(idx(idx21),:));
        S.M22(:,:,i) = func(R(idx(idx22),:));
        %% R-to-Z transformation, i.e. 0.5*log((1+x)/(1-x))
        S.Y11 = atanh(S.M11);
        S.Y12 = atanh(S.M12);
        S.Y21 = atanh(S.M21);
        S.Y22 = atanh(S.M22);
        fprintf('Completed: SC matrix and R-to-Z transformation, N=%4d\n',i);
    end
end
end