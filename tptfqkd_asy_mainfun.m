function [objv] = tptfqkd_asy_mainfun(X,l,sigma,N)
%The main function for calculating the key rate of TP-TFQKD under asymmetric channel

la=(l-100)/2;
lb=(l+100)/2;
pd=10^-8;
fec=1.1;
edz=0;
alfa=0.165;
etad=0.7;
taba=10^(-alfa*la/10);
tabb=10^(-alfa*lb/10);
etaa=taba*etad;
etab=tabb*etad;
ebssec=36/24*10^-10;
ebs=36/24*10^-10;
ebscor=36/24*10^-10;


mua=X(:,1)/1000;
nua=X(:,2)/1000;
pmua=X(:,3)/1000;
pnua=X(:,4)/1000;
poad=X(:,5)/1000;
mub=X(:,6)/1000;
nub=X(:,7)/1000;
pmub=X(:,8)/1000;
pnub=X(:,9)/1000;
pobd=X(:,10)/1000;
delta=pi/X(:,11);

poap=1-pmua-pnua-poad;
pobp=1-pmub-pnub-pobd;
poa=poad+poap;
pob=pobd+pobp;

kap=@(ka,kb)etaa*ka+etab*kb;
yka=@(ka,kb)exp(-kap(ka,kb)/2)*(1-pd);
xka=@(ka,kb)sqrt(etaa*ka*etab*kb);
Qkakb=@(ka,kb)2*yka(ka,kb)*(besseli(0,xka(ka,kb))-yka(ka,kb));
Qnuaob=2*exp(-etaa*nua/2)*(1-pd)*(1-exp(-etaa*nua/2)*(1-pd));
Qmuaob=2*exp(-etaa*mua/2)*(1-pd)*(1-exp(-etaa*mua/2)*(1-pd));
Qoanub=2*exp(-etab*nub/2)*(1-pd)*(1-exp(-etab*nub/2)*(1-pd));
Qoamub=2*exp(-etab*mub/2)*(1-pd)*(1-exp(-etab*mub/2)*(1-pd));
Qoaob=2*pd*(1-pd);

xoapobp=N*poap*pobp*Qoaob;
xoapmub=N*poap*pmub*Qoamub;
xmuaobp=N*pmua*pobp*Qmuaob;
xmuamub=N*pmua*pmub*Qkakb(mua,mub);
xmaxA=max(xoapobp+xoapmub,xmuaobp+xmuamub);

xoaobd=N*(poa*pob-poap*pobp)*Qoaob;
xmuaobd=N*pmua*pobd*Qmuaob;
xoadmub=N*poad*pmub*Qoamub;
xnuaob=N*pnua*pob*Qnuaob;
xoanub=N*poa*pnub*Qoanub;


Ez=((1-edz)*xoapobp*xmuamub+edz*xmuaobp*xoapmub)/(xoapobp*xmuamub+xmuaobp*xoapmub);
nz=(xoapobp*xmuamub+xmuaobp*xoapmub)/xmaxA;

[Exoaobdl,Exoaobdu] = xtoxbzmg(xoaobd,ebs);
[Exnuaobl,~] = xtoxbzmg(xnuaob,ebs);
[~,Exmuaobdu] = xtoxbzmg(xmuaobd,ebs);
[Exoanubl,~] = xtoxbzmg(xoanub,ebs);
[Exoadmubl,Exoadmubu] = xtoxbzmg(xoadmub,ebs);

Ey01l=mub/(mub*nub-nub^2)/N*max(0,exp(nub)*Exoanubl/(poa*pnub)-nub^2/mub^2*exp(mub)*Exoadmubu/(poad*pmub)-(mub^2-nub^2)/mub^2*Exoaobdu/(poa*pob-poap*pobp));
Ey10l=mua/(mua*nua-nua^2)/N*max(0,exp(nua)*Exnuaobl/(pnua*pob)-nua^2/mua^2*exp(mua)*Exmuaobdu/(pmua*pobd)-(mua^2-nua^2)/mua^2*Exoaobdu/(poa*pob-poap*pobp));


Ez01l=N*poap*pmub*mub*exp(-mub)*Ey01l;
Ez10l=N*pmua*pobp*mua*exp(-mua)*Ey10l;
Ez00l=pmua*pobp*exp(-mua)*Exoaobdl/(poa*pob-poap*pobp);
Ez0mubl=pmua*pmub*exp(-mua)*Exoadmubl/(poad*pmub);
Exoapmubl=poap*Exoadmubl/poad;
Exoapobpl=poap*pobp*Exoaobdl/(poa*pob-poap*pobp);


Es11zl=Ez01l*Ez10l/xmaxA;
[s11zl,~]=xbtoxzmg(Es11zl,ebs);

Es0mubzl=(Exoapmubl*Ez00l+Exoapobpl*Ez0mubl)/xmaxA;
[s0mubzl,~]=xbtoxzmg(Es0mubzl,ebs);


funmyLR= @(x)(1-pd)^2*(exp(-(etaa*nua+etab*nub)/2 - sqrt(etaa*nua*etab*nub)*cos(x)) - (1 - pd)*exp(-(etaa*nua+etab*nub)))...
    .*(exp(-(etaa*nua+etab*nub)/2 + sqrt(etaa*nua*etab*nub)*cos(x)) - (1 - pd)*exp(-(etaa*nua+etab*nub)))...
    ./((1-pd)*(exp(-(etaa*nua+etab*nub)/2 + sqrt(etaa*nua*etab*nub)*cos(x)) - (1 - pd)*exp(-(etaa*nua+etab*nub)))+...
    (1-pd)*(exp(-(etaa*nua+etab*nub)/2 - sqrt(etaa*nua*etab*nub)*cos(x)) - (1 - pd)*exp(-(etaa*nua+etab*nub))));
funmyp=@(x)1./((1-pd)*(exp(-(etaa*nua+etab*nub)/2 + sqrt(etaa*nua*etab*nub)*cos(x)) - (1 - pd)*exp(-(etaa*nua+etab*nub)))+...
    (1-pd)*(exp(-(etaa*nua+etab*nub)/2 - sqrt(etaa*nua*etab*nub)*cos(x)) - (1 - pd)*exp(-(etaa*nua+etab*nub))));
    
barqnuanubLR=integral(funmyLR, sigma, delta+sigma);
qnuanubp=integral(funmyp, sigma, delta+sigma);

E11xl=2*N*pnua*pnub*nua*nub*exp(-2*(nua+nub))*Ey01l*Ey10l*qnuanubp/pi;
[s11xl,~]=xbtoxzmg(E11xl,ebs);

EQ00u=Exoaobdu/(N*(poa*pob-poap*pobp));
EQ00l=Exoaobdl/(N*(poa*pob-poap*pobp));
Em2nu0l=1/2*N*pnua*pnub*delta*exp(-nua-nub)*EQ00l/pi;
Em02nul=1/2*N*pnua*pnub*delta*exp(-nua-nub)*EQ00l/pi;
[m02nul,~]=xbtoxzmg(Em02nul+Em2nu0l,ebs);

Em00u=1/2*N*pnua*pnub*exp(-2*(nua+nub))*EQ00u^2*qnuanubp/pi;
[~,m00u]=xbtoxzmg(Em00u,ebs);

mx=2*N*pnua*pnub*barqnuanubLR/pi;

m11xu=mx-m02nul+m00u;
if m11xu>mx
    m11xu=mx;
end
    
e11x=m11xu/s11xl;

e1zzu=e11x+gazmg(ebs,e11x,s11zl,s11xl);
if isempty(e1zzu)
objv=0;
else 


ll=(s0mubzl+s11zl.*(1-H2(e1zzu))-nz*fec*H2(Ez)-log2(2/ebscor)-2*log2(2/(ebssec*ebssec))-2*log2(1/(2*ebssec)))/N;

objv=-ll;
objv(e1zzu<0|Ez>1/2|e1zzu>1/2|e11x>1/2)=0;
end

end

