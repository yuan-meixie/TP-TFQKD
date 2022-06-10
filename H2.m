function [h]=H2(x)

if x==0|imag(x)~=0
    h=0;
else
h=-x*log2(x)-(1-x)*log2(1-x);
end
end