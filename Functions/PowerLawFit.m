function beta = PowerLawFit(y,x);
% function beta = PowerLawFit(y,x);
% fits the model y = beta(1) * x .^ beta(2)
% using a reweighted linearized form

logy = y.*log(y);
sol = pinv(cat(1,y.*ones(size(x)),y.*log(x))')*logy';
beta = [exp(sol(1)), sol(2)];
