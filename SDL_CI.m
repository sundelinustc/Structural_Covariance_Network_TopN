function CI = SDL_CI(y,CIType)
% calculate CI (CI = mean(x)+- t * (s / square(n)))
% Input
%    y      --- data
%    CIType --- 0.95 or 0.99
% Output
%    CI     --- confidence interval


% %% Based on Normal Distribution (not good for my data based on histogram display)
% % The output should + mean to gain the actual CI 
% N = size(y,1); % Number of ?Experiments? In Data Set
% yMean = mean(y); % Mean Of All Experiments At Each Value Of ?x?
% ySEM = std(y)/sqrt(N); % Compute ?Standard Error Of The Mean? Of All Experiments At Each Value Of ?x? . !!!??? May be Wrong!!!
% 
% a = (1-CIType)/2;
% ts = tinv([a 1-a], N-1);          % Calculate Probability Intervals Of t-Distribution
% CI = bsxfun(@times, ySEM, ts(:)); % Calculate Confidence Intervals Of All Experiments At Each Value Of ?x?
% 
% 
% 
% %% Another method based on Normal Distribution (not good for my data based on histogram display)
% https://www.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
% x = y;
% CIFcn = @(x,p)std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * tinv(abs([0,1]-(1-p/100)/2),sum(~isnan(x(:)))-1) + mean(x(:),'omitnan'); 
% CI = CIFcn(x,95); 


%% Based on not normal distribution
% https://www.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
x = y;
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));
CI = CIFcn(x,95); 


end