clear all; close all;format long
mredu=12766.336337294689656d0;
Nq1=61;
e0=0.005d0;
%k0= sqrt(2.d0*mredu*e0);
x0=2.09970623d0;
stepX=0.01d0;
im=sqrt(-1.0d0);
kf=4.d0*mredu*e0^2.d0;
w=sqrt(kf/mredu);
c0 = ((mredu*w)/pi)^(1.d0/4.d0);
c1 = (4.d0/pi*(mredu*w)^3.d0)^(1.d0/4.d0);

for i=1:Nq1
    ch=(i-(Nq1-1.0d0)/2.0d0-1.0d0)*stepX;
    x=x0+ch;
    expo = exp(-(x-x0)^2.d0*mredu*w/2.d0);
    vec0(i) = c0 * expo;
    vec1(i) = c1 * (x-x0) * expo;
end
%figure(1),plot((real(matrix(:,tt)).^2 +imag(matrix(:,tt)).^2)),hold on
%figure(1),plot(imag(matrix(:,tt))),hold on
%figure(2),plot(matrix2(:,tt)),hold on,pause
soma=0.0d0;
for i=1:Nq1
    soma=soma+conj(vec0(i))*vec0(i)*stepX;
end
sqrt(soma)
