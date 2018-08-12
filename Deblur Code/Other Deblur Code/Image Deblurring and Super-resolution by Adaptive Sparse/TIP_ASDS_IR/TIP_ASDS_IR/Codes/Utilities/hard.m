% define the soft threshold function, which is used above.
function y = hard(x,tau)

y = x.*double((x>=tau)|(-x>=tau));