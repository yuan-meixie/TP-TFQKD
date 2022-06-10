clear all
clc
tic

resolution=0.001; 
upperbound=1/resolution;
A=[-1,1,0,0,0,0;                  %mu>nu
     0,0,1,1,1,0];                %pmu+pnu+pod<1
b=[-1;upperbound-1];

lb=[100,0,0,0,0,2];
ub=[1000,100,500,400,100,50];      
IntCon=[1,2,3,4,5,6];
nonlcon=[]; 
Ldata=530;
% Ldata=[0:10:700];
pare_result=zeros(length(Ldata),6);
para_resultbest=zeros(length(Ldata),6);
Rtptf_sy=zeros(length(Ldata),1);

sigma= pi/36; %set the value of the angle of misalignment.
kk=13;
N= 10^kk; %set the value of date size

for j=1:length(Ldata)
    
    d=Ldata(j);
      tot=5;
fun=@(x)tptfqkd_mainfun(x,d,sigma,N);

for i=1:tot

[para_result(j,:)]=ga(fun,6,A,b,[],[],lb,ub,nonlcon,IntCon);
mmp=-fun(para_result(j,:));
if mmp>Rtptf_sy(j)
    Rtptf_sy(j)=mmp;
    para_resultbest(j,:)=para_result(j,:);
end
end

end


Etad=0.7;
Rband=-log2(1-Etad*10.^(-0.165*Ldata/10));
semilogy(Ldata,Rtptf_sy,Ldata, Rband,'r')    

Rbest=strcat('Rtptf_sy_',num2str(kk));
eval([Rbest,'=Rtptf_sy;']);
save(strcat('Rtptf_sy_',num2str(kk),'.mat'),strcat('Rtptf_sy_',num2str(kk)));


toc