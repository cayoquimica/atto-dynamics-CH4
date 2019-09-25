clear all,close all
matrix=load('amp-matrix.txt');
q1=matrix(:,1);
for i=1:(length(matrix(1,:)))/2-1
    real(:,i)=matrix(:,2*i+1);ima(:,i)=matrix(:,2*i+2);
end

figure(1),plot((real)) 
figure(2),plot(ima)
figure(3),plot(real.^2 + ima.^2) 

hpsi=load('hpsi.txt');
psic=load('psic.txt');

soma1=0;soma2=0;
for i=1:61
    xxx(i,1)= (psic(i,1)*hpsi(i,1)) - (psic(i,2)*hpsi(i,2));
    xxx(i,2)= (psic(i,1)*hpsi(i,2)) + (psic(i,2)*hpsi(i,1))
    soma1=soma1+xxx(i,1)
    soma2=soma1+xxx(i,2)
end

%break

x=load('test.txt');
figure(4),plot(x)
break
for i=1:length(real)
plot(real(i),ima(i)),hold on%,pause
end