clear all; close all;format long
Mtotal=12766.3363372946896561d0;
a=0.2d0; % 10 times the step in x
Nq1=61;
e0=0.00534042d0;
k0= sqrt(2.d0*Mtotal*e0);
const = ( (2.d0*a^2.d0)/pi )^(1.d0/4.d0);
x0=2.09970623d0;
stepX=0.01d0;
im=sqrt(-1.0d0);
time=100.d0;
step=1.0d0;% time step

%fac=1.0d0/12.0d0/stepX^2.0d0
%y31=(-1.0d0/(2.0d0*Mtotal))*((-fac*complex(2.6402517221245114,-0.6280713790801997))+(16.0d0*fac*complex(2.7775342129801497,-0.3258183282476641))+(-30.0d0*fac*complex(2.8246850472175167,-0.0d0))+(16.0d0*fac*complex(2.7775342129801497,-0.3258183282476641))+(-fac*complex(2.6402517221245114,-0.6280713790801997)))
%break

for j=0:time
    tt=j+1d0;
    j=j*step;
    teta=( atan( (2.0d0*j)/(a.^2.0d0*Mtotal) ) ) /2.0d0;
    phi=-teta-(k0^2.0d0*j/(2.0d0*Mtotal));
    term1 = exp(im*phi)/(a^4.0d0+(4.0d0*j/Mtotal^2.0d0))^(1.0d0/4.0d0);
    for i=1:Nq1
        ch=(i-(Nq1-1)/2-1)*stepX;
        x=x0+ch;
        term2 = exp(im*k0*(x-x0));
        term3= exp( - ((x-x0) - (k0*j/Mtotal))^2.0d0 / (a^2.0d0 + (2.0d0*im*j/Mtotal)) );
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
    %figure(1),plot(imag(matrix(:,tt))),hold on
    %figure(2),plot(matrix2(:,tt)),hold on,pause
    soma=0.0d0;
    for i=1:Nq1
        soma=soma+conj(matrix(i,tt))*matrix(i,tt)*stepX;
    end
    sqrt(soma)
    
    %matrix
    %kine(tt)= conj(matrix(:,tt))'* deriv(:);
    %kine(tt)
    %pause
end

