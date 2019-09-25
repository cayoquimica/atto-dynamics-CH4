clear all; close all;format long
Mtotal=29244.18598791975d0;
a=4d0;
e0=03.674930360070d-1;
k0= sqrt(2.d0*Mtotal*e0);
const = ( (2.d0*a^2.d0)/pi )^(1.d0/4.d0);
im=sqrt(-1.d0);
time=100.d0;
step=4d-3;% time step
for j=0:time
    tt=j+1d0;
    j=j*step
    teta=( atan( (2.0d0*j)/(a.^2.0d0*Mtotal) ) ) /2.0d0;
    phi=-teta-(k0^2.0d0/(2.0d0*Mtotal))*j;
    term1 = exp(im*phi)/(a^4.0d0+(4.0d0*j/Mtotal^2.0d0))^(1.0d0/4.0d0);
    for i=1:61
        term2 = exp(im*k0*(i-16.0d0));
        term3= exp( - ((i-16.0d0) - (k0*j/Mtotal))^2.0d0 / (a^2.0d0 + (2.0d0*im*j/Mtotal)) );
        vec1(i) = const * term1 * term2 *term3;
        matrix(i,tt) = vec1(i);
        %matrix2(i,tt)= sqrt(2/(pi*a^2))  *  1 / (sqrt(1+(4*j^2/(Mtotal^2*a^4))))  *  exp(-2*a^2 * ((i-16)-k0*j/Mtotal)^2/(a^4+(4*j^2/Mtotal^2)));
        %matrix2(i,tt) = real(matrix2(i,tt))^2 + imag(matrix2(i,tt))^2;
        const2 = -k0^2.d0 - (2/a^2.d0) - (4*im*k0*(i-16.d0)/a^2.d0) + (4*(i-16)^2.d0/a^4.d0);
        deriv(i) = -1.d0/(2.d0*Mtotal) * const2 * exp(im*k0*(i-16) - ((i-16)^2.d0/a^4.d0));
        %kine(tt)=0.0d0;
        %kine(tt)=kine(tt) + conj(matrix(i,tt))'* deriv(i);
    end
    figure(1),plot((real(matrix(:,tt)).^2 +imag(matrix(:,tt)).^2)),hold on
    %figure(2),plot(matrix2(:,tt)),hold on,pause
    norm(matrix(:,tt))
    matrix
    %kine(tt)= conj(matrix(:,tt))'* deriv(:);
    %kine(tt)
    pause
end

