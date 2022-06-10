function [xbl,xbu] =  xtoxbzmg(x,ebs)


deta1=(-log(ebs)+sqrt((log(ebs))^2-8*log(ebs)*x))/2;
deta2=-log(ebs)+sqrt((log(ebs))^2-2*x*log(ebs));
xbl=x-deta1;
paicu=find(imag(xbl)==0);
xbl(paicu)=max(zeros(size(xbl(paicu))),xbl(paicu));
xbu=x+deta2;

end

