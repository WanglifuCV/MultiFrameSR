function  DCT  =  Get_DCT( b )

K     =  b*b;

n     =  ceil(sqrt(K));
DCT   =  zeros(b,n);
for k =  0 : 1 : n-1
    a   =  cos([0:1:b-1]'*k*pi/n);
    if k>0 
        a  =  a-mean(a); 
    end;
    DCT(:,k+1)  =  a/norm(a);
end;

DCT  =  kron(DCT,DCT);
