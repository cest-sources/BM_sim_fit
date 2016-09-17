function string = errorbar_str(val, errorbar)
%ERRORBAR_STR format a value and error bar as a reasonable-looking string.
%
% string = errorbar_str(val, errorbar)
% string = errorbar_str(vals)
%
% If a vector of values is supplied then first compute mean and standard error
% and then format them. Note that no sanity checking is done to verify that the
% mean and variance give sensible estimators.
%
% I always give error bars to 2sf (even if not justified). I give val the same
% level of precision and format it using the same exponent. I avoid scientific
% notation for moderate numbers and if it avoids implying more precision than
% specified above. If the error bar is bigger than the value then I just
% represent each to 2sf in some reasonable way.

% SOOO much longer than I intended. I am *sure* there's a better way of doing
% this!

% Iain Murray, November 2007, NaN/Inf handling October 2014
% Tweaks for more corner cases October/November 2015.

if nargin == 1
    assert(isvector(val));
    errorbar = std(val)/sqrt(length(val));
    val = mean(val);
end

assert(isscalar(val) && isscalar(errorbar));
if isinf(val) || isnan(val)
    string = sprintf('%g', val);
    return
end

if errorbar == 0
    string = [two_sf(val) ' ± 0'];
    return
end

if isnan(errorbar)
    string = [two_sf(val) ' ± 0'];
    return
end

oom = @(x) floor(log10(abs(x)));
oom_errorbar = oom(errorbar);

% Do 2sf rounding now, so I don't get confused about what precision the most
% significant digit has:
errorbar = round(errorbar/(10^(oom_errorbar-1))) * (10^(oom_errorbar-1));
oom_errorbar = oom(errorbar);

oom_val = oom(val);
exponent = max(oom(val), oom(errorbar));
precision = oom(errorbar) - 1; % oom of last digit (errorbars always 2sf)

if oom_errorbar > oom_val
    % Errorbar is bigger than value, which makes comparisons silly, just format
    % each side independently, ensuring not implying more than 2sf.
    % TODO could improve this part in some cases.
    string = [two_sf(val), ' ± ', two_sf(errorbar)];
elseif (abs(exponent) > 4) || (precision > 0)
    % Use scientific notation for large exponents or if I would need trailing
    % zeros before the decimal point to avoid giving distracting levels of
    % precision.
    decimals = exponent - precision;
    eb_mantissa = errorbar/(10^exponent);
    val_mantissa = val/(10^exponent);
    string = [my_format(val_mantissa, exponent, decimals), ' ± ', ...
        my_format(eb_mantissa, exponent, decimals)];
else
    % Easy case: can use normal decimal notation with no exponents
    decimals = abs(precision);
    format = sprintf(['%%0.%df ± %%0.%df'], decimals, decimals);
    string = sprintf(format, val, errorbar);
end


function string = my_format(mantissa, exponent, digits)
% Outputs number as mantissa with digits decimal digits and e+exponent bunged
% on the end. Unusually the mantissa is allowed to be >=10 or <1, as it's
% more important for the exponents to match for ease of comparison of value
% and errorbar.

sign_strings = ['-+'];
signstr = @(x) sign_strings((x>=0)+1);
format = sprintf(['%%0.%dfe',signstr(exponent),'%02d'], digits, abs(exponent));
string = sprintf(format, mantissa);


function string = two_sf(num)
% Nice formatting for a number that should be displayed to two sig figs.

if num == 0
    string = '0';
    return
end
oom = @(x) floor(log10(abs(x)));
oom_num = oom(num);
if (oom_num == 1)
    string = sprintf('%d', round(num));
elseif (oom_num < -4) || (oom_num > 0)
    string = sprintf('%0.1e', num);
else
    string = sprintf('%0.2g', num);
end