function [gamw] = gazmg(ebs,lamda,n,k)


g=2*(n+k)./(n.*k).*log(sqrt(n+k).*cc(n,k,lamda)./(sqrt(2*pi*n.*k.*lamda.*(1-lamda))*ebs));

pacu1=find(k<n);
pacu2=find(k>=n);
gam(pacu1)=((1-2*lamda(pacu1)).*n(pacu1).*g(pacu1)./(n(pacu1)+k(pacu1))+sqrt(n(pacu1).^2.*g(pacu1).^2./(n(pacu1)+k(pacu1)).^2+4*lamda(pacu1).*g(pacu1)-4*lamda(pacu1).^2.*g(pacu1)))./(2+2*n(pacu1).^2.*g(pacu1)./(n(pacu1)+k(pacu1)).^2);
gam(pacu2)=((1-2*lamda(pacu2)).*k(pacu2).*g(pacu2)./(n(pacu2)+k(pacu2))+sqrt(k(pacu2).^2.*g(pacu2).^2./(n(pacu2)+k(pacu2)).^2+4*lamda(pacu2).*g(pacu2)-4*lamda(pacu2).^2.*g(pacu2)))./(2+2*k(pacu2).^2.*g(pacu2)./(n(pacu2)+k(pacu2)).^2);
gamw=gam';
function [ccm] = cc(n,k,lamda)
ccm=[];
ccm=exp(1./(8*(n+k))+1./(12*k)-1./(12*k.*lamda+1)-1./(12*k.*(1-lamda)+1));
end
end

