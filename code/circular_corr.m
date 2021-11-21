
function [r, pval] = circular_corr(a,b)

n = length(a);

a2 = sum(exp(1i*a),1);
b2 = sum(exp(1i*b),1);


% obtain mean by
a3 = angle(a2);
b3 = angle(b2);


% compute correlation coeffcient from p. 176
num = sum(sin(a - a3) .* sin(b - b3));
den = sqrt(sum(sin(a - a3).^2) .* sum(sin(b - b3).^2));
r = num / den;	

% compute pvalue
l20 = mean(sin(a - a3).^2);
l02 = mean(sin(b - b3).^2);
l22 = mean((sin(a - a3).^2) .* (sin(b - b3).^2));

ts = sqrt((n * l20 * l02)/l22) * r;
pval = 2 * (1 - normcdf(abs(ts)));