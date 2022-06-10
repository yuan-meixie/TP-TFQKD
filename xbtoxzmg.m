function [xl,xu] = xbtoxzmg(xb,ebs)
    deta1=[];
    deta2=[];
    paicu=[];
    deta1=(-log(ebs)+sqrt((log(ebs))^2-8*log(ebs)*xb))./(2*xb);
    deta2=sqrt(-2*log(ebs)./xb);
    xl=xb.*(1-deta2);
    paicu=find(imag(xl)==0);
    xl(paicu)=max(zeros(size(xl(paicu))),xl(paicu));
    xu=xb.*(1+deta1);
end