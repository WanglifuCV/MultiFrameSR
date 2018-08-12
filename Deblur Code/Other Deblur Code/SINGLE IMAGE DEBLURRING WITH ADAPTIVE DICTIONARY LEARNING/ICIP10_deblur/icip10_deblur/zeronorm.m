function n = zeronorm(A)
% ZERONORM function is to find how many zeros in a vector.

[m,l] = size(A);
n=0;

for i=1:m
    for j=1:l
        if A(i,j)~=0
            n = n+1;
        end
    end
end