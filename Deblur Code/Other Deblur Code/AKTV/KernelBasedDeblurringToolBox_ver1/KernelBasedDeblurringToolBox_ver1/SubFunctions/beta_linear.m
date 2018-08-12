function [z, zx, zy] = beta_linear

global s11;
global s12;
global s22;
global f1;
global f2;

[N, M] = size(s11);

s22inv = inv22(s22);
s12s22inv = zeros(N, M, 2);
s12s22inv(:,:,1) = s12(:,:,1) .* s22inv(:,:,1) + s12(:,:,2) .* s22inv(:,:,2);
s12s22inv(:,:,2) = s12(:,:,1) .* s22inv(:,:,3) + s12(:,:,2) .* s22inv(:,:,4);
denom = s11 - (s12s22inv(:,:,1) .* s12(:,:,1) + s12s22inv(:,:,2) .* s12(:,:,2));
z = (f1 - (s12s22inv(:,:,1) .* f2(:,:,1) + s12s22inv(:,:,2) .* f2(:,:,2))) ./ (denom + (denom == 0));

s21s11inv = zeros(N, M, 2);
ES = zeros(N, M, 4);
s21s11inv(:,:,1) = s12(:,:,1) ./ s11;
S21s11inv(:,:,2) = s12(:,:,2) ./ s11;
ES(:,:,1) = s22(:,:,1) - s21s11inv(:,:,1) .* s12(:,:,1);
ES(:,:,2) = s22(:,:,2) - s21s11inv(:,:,1) .* s12(:,:,2);
ES(:,:,3) = s22(:,:,3) - s21s11inv(:,:,2) .* s12(:,:,1);
ES(:,:,4) = s22(:,:,4) - s21s11inv(:,:,2) .* s12(:,:,2);
ESinv = inv22(ES);
tmp1 = - s21s11inv(:,:,1) .* f1 + f2(:,:,1);
tmp2 = - s21s11inv(:,:,2) .* f1 + f2(:,:,2);
zx = ESinv(:,:,1) .* tmp1 + ESinv(:,:,3) .* tmp2;
zy = ESinv(:,:,2) .* tmp1 + ESinv(:,:,4) .* tmp2;