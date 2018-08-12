% function y = ishift2(img, sx, sy)
% % ISHIFT2 [y = ishift2(img, sx, sy)]
% % 2 dimensional integer shift function
% % img : image to be shifted
% % sx  : x direction shift distance
% % sy  : y direction shift distance
% %
% % coded by hiro on Nov 22, 2004
% 
% sx = round(sx);
% sy = round(sy);
% [N M] = size(img);
% y = zeros(N, M);
% if sx >= 0
%     if sy >= 0
%         y(sx+1:N, sy+1:M) = img(1:N-sx, 1:M-sy);
%     else
%         y(sx+1:N, 1:M+sy) = img(1:N-sx, -sy+1:M);
%     end
% else
%     if sy >= 0
%         y(1:N+sx, sy+1:M) = img(-sx+1:N, 1:M-sy);
%     else
%         y(1:N+sx, 1:M+sy) = img(-sx+1:N, -sy+1:M);
%     end
% end
function y = ishift2(y, sx, sy)
% ISHIFT2 [y = ishift2(img, sx, sy)]
% 2 dimensional integer shift function
% img : image to be shifted
% sx  : x direction shift distance
% sy  : y direction shift distance
%
% coded by hiro on Nov 22, 2004

sx = round(sx);
sy = round(sy);
[N M] = size(y);
if sx >= 0
    if sy >= 0
        y(sx+1:N, sy+1:M) = y(1:N-sx, 1:M-sy);
        y(1:sx, :) = 0;
        y(:, 1:sy) = 0;
    else
        y(sx+1:N, 1:M+sy) = y(1:N-sx, -sy+1:M);
        y(1:sx, :) = 0;
        y(:, M+sy+1:end) = 0;
    end
else
    if sy >= 0
        y(1:N+sx, sy+1:M) = y(-sx+1:N, 1:M-sy);
        y(N+sx+1:end, :) = 0;
        y(:, 1:sy) = 0;
    else
        y(1:N+sx, 1:M+sy) = y(-sx+1:N, -sy+1:M);
        y(N+sx+1:end, :) = 0;
        y(:, M+sy+1:end) = 0;
    end
end