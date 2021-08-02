function p = SDL_p_permutation(v0,v1)
% calculate p values of permutation
% Input
%    v0 --- value of original data, 1x1
%    v1 --- value of permutations, 1xN or Nx1
% Output
%    p  --- statistical significance, i.e. the probability that random data exceedes the original data

if v0 > mean(v1)
    p = sum((v1-v0)>0)/length(v1);
else
    p = sum((v1-v0)<0)/length(v1);
end

    
end