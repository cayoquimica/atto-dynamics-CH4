clear all; close all

mat=load('e1-q1.txt');
e1=mat(:,2);
e2=mat(:,3);
e3=mat(:,4);
r=mat(:,1);
step=0.04;
r1=-5.0:step:5.0;r1=r1';
e1int=interp1(r,e1,r1,'spline');

%e1int=e1int';
e0=-39.8385;
for i=1:length(e1int)
if e1int(i) < e0
e0=e1int(i);i
end
end
e1int=e1int-e0;
save('e1int.txt','e1int','-ascii')

for i=1:length(r)-32
    r0(i)=r(i+8);
    e10(i)=e1(i+8);
end
step2=0.04;
r11=-5.0-step2/2:step2:5.0+step2/2;r11=r11';
e1int0=interp1(r0,e10,r11,'spline');

e2int=interp1(r,e2,r11,'spline');
e3int=interp1(r,e3,r11,'spline');

%e1int=e1int';
e0=-39.8385;
for i=1:length(e1int0)
if e1int0(i) < e0
e0=e1int0(i);i
end
end
e1int0=e1int0-e0;
e2int=e2int-e0;
e3int=e3int-e0;
save('e1int0.txt','e1int0','-ascii')
save('e2int.txt','e2int','-ascii')
save('e3int.txt','e3int','-ascii')

a0=5.2917721092d-11;
v0=2.1876912633d6;
c=2.997924580d8;
h=1.05457168d-34;
lam=800d-9;
energy=2*pi*137/(lam/a0)
freq=c/lam*2*pi/(v0/a0)

% sum=0.0d0;
% for i=1:length(x)
% sum=sum+x(i)^2*e1int(i)*step;
% end
% sum

pdm1=load('fort.1041.outpdm1x.dat');
pdm10=interp1(r,pdm1,r11,'spline');
save('pdm1x.txt','pdm10','-ascii')

pdm2=load('fort.1041.outpdm2x.dat');
pdm20=interp1(r,pdm2,r11,'spline');
save('pdm2x.txt','pdm20','-ascii')

pdm3=load('fort.1041.outpdm3x.dat');
pdm30=interp1(r,pdm3,r11,'spline');
save('pdm3x.txt','pdm30','-ascii')
%figure(1),plot(pdm1,'b','LineWidth',2),hold on,plot(pdm2,'k','LineWidth',2),hold on,plot(pdm3,'r','LineWidth',2)
%figure(2),plot(pdm10,'b','LineWidth',2),hold on,plot(pdm20,'k','LineWidth',2),hold on,plot(pdm30,'r','LineWidth',2)

E0=0.01;
t0=400;
phi=0;
w=0.056937;
sig=50;
i=0.d0;
t1=350;
phi=w*(t1-t0)+pi/2;
sig2=30;
%for i=1:100
for j=1:80000+1
    %for i=1:length(r1)
    t=(j-1)/100;
        %yt(j)= E0/w * (-(t-t0)/sig^2*sin(w*(t-t0)+phi*i/5) + w*cos(w*(t-t0)+phi*i/5) ) * exp(-(t-t0)^2/(2*sig^2)) ;
        yt(j)= E0/w * (cos(w*(t-t0)+pi/2) ) * exp(-(t-t0)^2/(2*sig^2)) ;
        yt2(j)= E0/w * (cos(w*(t-t1)+phi) ) * exp(-(t-t1)^2/(2*sig2^2)) ;
    %end
end
yt=yt';
yt2=yt2';
plot(yt)
pause
hold on,plot(yt2,'r')
pause
%end

save('field.txt','yt','-ascii')

x21=load('fort.1041.outtdm21x.dat');
for i=1:length(r)-0
    r01(i)=r(i+0);
    x211(i)=x21(i+0);
end
tdm21x=interp1(r01,x211,r11,'pchip');
save('tdm21x.txt','tdm21x','-ascii')




x31=load('fort.1041.outtdm31x.dat');
for i=1:length(r)-4
    r02(i)=r(i+4);
    x311(i)=x31(i+4);
end
tdm31x=interp1(r02,x311,r11,'pchip');
save('tdm31x.txt','tdm31x','-ascii')

x32=load('fort.1041.outtdm32x.dat');
for i=1:length(r)
    r03(i)=r(i);
    x321(i)=x32(i);
end
tdm32x=interp1(r03,x321,r11,'pchip');
save('tdm32x.txt','tdm32x','-ascii')

%figure(3),plot(x21,'b','LineWidth',2),hold on,plot(x31,'k','LineWidth',2),hold on,plot(x32,'r','LineWidth',2)
figure(4),plot(tdm21x,'b','LineWidth',2),hold on,plot(tdm31x,'k','LineWidth',2),hold on,plot(tdm32x,'r','LineWidth',2)



pdm1y=load('fort.1041.outpdm1y.dat');
pdm10y=interp1(r,pdm1y,r11,'spline');
save('pdm1y.txt','pdm10y','-ascii')

pdm2y=load('fort.1041.outpdm2y.dat');
pdm20y=interp1(r,pdm2y,r11,'spline');
save('pdm2y.txt','pdm20y','-ascii')

pdm3y=load('fort.1041.outpdm3y.dat');
pdm30y=interp1(r,pdm3y,r11,'spline');
save('pdm3y.txt','pdm30y','-ascii')

y21=load('fort.1041.outtdm21y.dat');
for i=1:length(r)-0
    r01(i)=r(i+0);
    y211(i)=y21(i+0);
end
tdm21y=interp1(r01,y211,r11,'pchip');
save('tdm21y.txt','tdm21y','-ascii')

y31=load('fort.1041.outtdm31y.dat');
for i=1:length(r)-4
    r02(i)=r(i+4);
    y311(i)=y31(i+4);
end
tdm31y=interp1(r02,y311,r11,'pchip');
save('tdm31y.txt','tdm31y','-ascii')

y32=load('fort.1041.outtdm32y.dat');
for i=1:length(r)
    r03(i)=r(i);
    y321(i)=y32(i);
end
tdm32y=interp1(r03,y321,r11,'pchip');
save('tdm32y.txt','tdm32y','-ascii')

%figure(3),plot(y21,'b','LineWidth',2),hold on,plot(y31,'k','LineWidth',2),hold on,plot(y32,'r','LineWidth',2)
figure(5),plot(tdm21y,'b','LineWidth',2),hold on,plot(tdm31y,'k','LineWidth',2),hold on,plot(tdm32y,'r','LineWidth',2)


















