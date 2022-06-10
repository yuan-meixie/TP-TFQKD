clear all
clc
tic

resolution=0.001; 
upperbound=1/resolution;
A=[-1,1,0,0,0,0,0,0,0,0,0;              %mua>nua
   0,0,0,0,0,-1,1,0,0,0,0;              %mub>nub
    0,0,1,1,1,0,0,0,0,0,0;              %pmua+pnua+poda<1
     0,0,0,0,0,0,0,1,1,1,0; ];          %pmub+pnub+podb<1
b=[-1;-1;upperbound-1;upperbound-1];


lb=[90,1,20,30,1,500,70,100,30,1,2];
ub=[150,5,100,500,100,1000,150,500,450,70,80];  
gilen=length(lb);

IntCon=[1,2,3,4,5,6,7,8,9,10,11];  
nonlcon=[]; 
Ldata=530;
% Ldata=[0:10:700];
pare_result=zeros(length(Ldata),gilen);
para_resultbest=zeros(length(Ldata),gilen);
Rtptf_asy=zeros(length(Ldata),1);

sigma= pi/36; %set the value of the angle of misalignment.
kk=13;
N= 10^kk; %set the value of date size

for j=1:length(Ldata)
    
    d=Ldata(j);
      tot=5;
fun=@(x)tptfqkd_asy_mainfun(x,d,sigma,N);

for i=1:tot

[para_result(j,:)]=ga(fun,gilen,A,b,[],[],lb,ub,nonlcon,IntCon);
mmp=-fun(para_result(j,:));
if mmp>Rtptf_asy(j)
    Rtptf_asy(j)=mmp;
    para_resultbest(j,:)=para_result(j,:);
end
end

end


Etad=0.7;
Rband=-log2(1-Etad*10.^(-0.165*Ldata/10));
semilogy(Ldata,Rtptf_asy,Ldata, Rband,'r')    

Rbest=strcat('Rtptf_asy_',num2str(kk));
eval([Rbest,'=Rtptf_asy;']);
save(strcat('Rtptf_asy_',num2str(kk),'.mat'),strcat('Rtptf_asy_',num2str(kk)));


toc