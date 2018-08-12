function [mvi, mvd] = convertmv(mvs, r, ref)

% motion
mv = -[0 mvs(1, :); 0 mvs(2, :)] * r;
mv(1,:) = mv(1,:) - mv(1,ref);
mv(2,:) = mv(2,:) - mv(2,ref);
% integer motion
mvi = round(mv);
% decimal motion
mvd = mv - mvi;