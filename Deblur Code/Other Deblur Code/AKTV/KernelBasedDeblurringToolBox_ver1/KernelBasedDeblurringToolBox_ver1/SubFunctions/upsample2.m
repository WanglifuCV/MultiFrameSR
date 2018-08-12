function y = upsample2(x, u)
% UPSAMPLE2 [y = upsample2(x, u)]
% 2 dimensional upsampling function
% x : input to be upsampled
% u : upsampling factor
%
% coded by hiro on Nov 21, 2004

if u == 1
    y = x;
    return;
end

y = upsample(upsample(x.', u).', u);