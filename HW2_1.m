clear;

N = 1024;

%% signals
[sig1, noisy1] = wnoise('blocks',10,3);
[sig2, noisy2] = wnoise('bumps',10,3);
[sig3, noisy3] = wnoise('heavy sine',10,3);
[sig4, noisy4] = wnoise('doppler',10,3);

figure(1)
subplot(2,2,1)
plot(sig1);
title('Blocks')
subplot(2,2,2)
plot(sig2);
title('Bumps')
subplot(2,2,3)
plot(sig3);
title('Heavy Sine')
subplot(2,2,4)
plot(sig4);
title('Doppler')

figure(2)
subplot(2,2,1)
plot(noisy1);
title('Blocks contaminated')
subplot(2,2,2)
plot(noisy2);
title('Bumps contaminated')
subplot(2,2,3)
plot(noisy3);
title('Heavy Sine contaminated')
subplot(2,2,4)
plot(noisy4);
title('Doppler contaminated')


%% wavelets

[phi1,psi1,xval1] = wavefun('haar',5);
[phi2,psi2,xval2] = wavefun('db10',5);
[phi3,psi3,xval3] = wavefun('sym5',5);
[phi4,psi4,xval4] = wavefun('coif5',5);

figure(3)
subplot(2,2,1)
plot(xval1,psi1);
title('Haar')
subplot(2,2,2)
plot(xval2,psi2);
title('Daub5')
subplot(2,2,3)
plot(xval3,psi3);
title('Symlet5')
subplot(2,2,4)
plot(xval4,psi4);
title('Coiflet5')

%% A) orthogonality
ort1 = dot(phi1.', psi1);
ort2 = dot(phi2.', psi2);
ort3 = dot(phi3.', psi3);
ort4 = dot(phi4.', psi4);

%% B) 5-level MRA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HAAR WT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%a1 = zeros(5, 512);
%b1 = zeros(5, 512);
%d1 = zeros(5, 512);
%fr1 = zeros(1, 1024);
for i=1:1:N/2
    a1(1,i)=(sig4(2*i-1)+sig4(2*i))/sqrt(2);
    d1(1,i)=(sig4(2*i-1)-sig4(2*i))/sqrt(2);
end

for j=2:1:5
    for i=1:1:N/(2^j)
        a1(j,i)=(a1(j-1,2*i-1)+a1(j-1,2*i))/sqrt(2);
        d1(j,i)=(a1(j-1,2*i-1)-a1(j-1,2*i))/sqrt(2);
    end
end

b1(1:5,1:N/2)=0;
b1(5,1)=a1(5,1);

for j=5:-1:2
    for i=1:1:N/(2^j)
        b1(j-1,2*i-1)=(b1(j,i)+d1(j,i))/sqrt(2);
        b1(j-1,2*i)=(b1(j,i)-d1(j,i))/sqrt(2);
    end
end

for i=1:1:N/2
    fr4(2*i-1)=(b1(1,i)+d1(1,i))/sqrt(2);
    fr4(2*i)=(b1(1,i)-d1(1,i))/sqrt(2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Daub5 WT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% coeff of scaling function
al1=0.22641898/sqrt(2);
al2=0.85394354/sqrt(2);
al3=1.02432694/sqrt(2);
al4=0.19576696/sqrt(2);
al5=-0.34365671/sqrt(2);
al6=-0.04560113/sqrt(2);
al7=0.10970265/sqrt(2);
al8=-0.00882680/sqrt(2);
al9=-0.01779187/sqrt(2);
al10=0.0047174/sqrt(2);
%%coeff of wavelet function
be1=al10;
be2=-al9;
be3=al8;
be4=-al7;
be5=al6;
be6=-al5;
be7=al4;
be8=-al3;
be9=al2;
be10=-al1;
fv1=[al1 al2 al3 al4 al5 al6 al7 al8 al9 al10]; %scaling function
fw1=[be1 be2 be3 be4 be5 be6 be7 be8 be9 be10]; %wavelet function

Aa=[];
Da=[];
f = sig1;
n=length(f);

for i=1:n/2-4;
   a11 = fv1(1)*f(2*i-1) + fv1(2)*f(2*i) + fv1(3)*f(2*i+1) + fv1(4)*f(2*i+2) + fv1(5)*f(2*i+3) + fv1(6)*f(2*i+4) + fv1(7)*f(2*i+5) + fv1(8)*f(2*i+6) + fv1(9)*f(2*i+7) + fv1(10)*f(2*i+8);
   Aa=[Aa a11];
   d11 = fw1(1)*f(2*i-1) + fw1(2)*f(2*i) + fw1(3)*f(2*i+1) + fw1(4)*f(2*i+2) + fw1(5)*f(2*i+3) + fw1(6)*f(2*i+4) + fw1(7)*f(2*i+5) + fw1(8)*f(2*i+6) + fw1(9)*f(2*i+7) + fw1(10)*f(2*i+8);
   Da=[Da d11];
end 
   n4=n/2-3;
   a15=fv1(9)*f(1) + fv1(10)*f(2) + fv1(1)*f(2*n4-1) + fv1(2)*f(2*n4) + fv1(3)*f(2*n4+1) + fv1(4)*f(2*n4+2) + fv1(5)*f(2*n4+3) + fv1(6)*f(2*n4+4) + fv1(7)*f(2*n4+5) + fv1(8)*f(2*n4+6);
   d15=fw1(9)*f(1) + fw1(10)*f(2) + fw1(1)*f(2*n4-1) + fw1(2)*f(2*n4) + fw1(3)*f(2*n4+1) + fw1(4)*f(2*n4+2) + fw1(5)*f(2*n4+3) + fw1(6)*f(2*n4+4) + fw1(7)*f(2*n4+5) + fw1(8)*f(2*n4+6);
   n3=n/2-2;
   a14=fv1(7)*f(1) + fv1(8)*f(2) + fv1(9)*f(3) + fv1(10)*f(4) + fv1(1)*f(2*n3-1) + fv1(2)*f(2*n3) + fv1(3)*f(2*n3+1) + fv1(4)*f(2*n3+2) + fv1(5)*f(2*n3+3) + fv1(6)*f(2*n3+4);
   d14=fw1(7)*f(1) + fw1(8)*f(2) + fw1(9)*f(3) + fw1(10)*f(4) + fw1(1)*f(2*n3-1) + fw1(2)*f(2*n3) + fw1(3)*f(2*n3+1) + fw1(4)*f(2*n3+2) + fw1(5)*f(2*n3+3) + fw1(6)*f(2*n3+4);
   
   n1=n/2-1;
   a13=fv1(5)*f(1) + fv1(6)*f(2) + fv1(7)*f(3) + fv1(8)*f(4) + fv1(9)*f(5) + fv1(10)*f(6) + fv1(1)*f(2*n1-1) + fv1(2)*f(2*n1) + fv1(3)*f(2*n1+1) + fv1(4)*f(2*n1+2);
   d13=fw1(5)*f(1) + fw1(6)*f(2) + fw1(7)*f(3) + fw1(8)*f(4) + fw1(9)*f(5) + fw1(10)*f(6) + fw1(1)*f(2*n1-1) + fw1(2)*f(2*n1) + fw1(3)*f(2*n1+1) + fw1(4)*f(2*n1+2);
   
   n2=n/2;
   a12 = fv1(1)*f(2*n2-1) + fv1(2)*f(2*n2) + fv1(3)*f(1) + fv1(4)*f(2) + fv1(5)*f(3) + fv1(6)*f(4) + fv1(7)*f(5) + fv1(8)*f(6) + fv1(9)*f(7) + fv1(10)*f(8);
   d12 = fw1(1)*f(2*n2-1) + fw1(2)*f(2*n2) + fw1(3)*f(1) + fw1(4)*f(2) + fw1(5)*f(3) + fw1(6)*f(5) + fw1(7)*f(5) + fw1(8)*f(6) + fw1(9)*f(7) + fw1(10)*f(8);
   
   a1=[Aa a15 a14 a13 a12];
   d1=[Da d15 d14 d13 d12];
%%


figure(4)
subplot(4,2,1)
plot(noisy1,'r');
title('Blocks Original')
subplot(4,2,2)
plot(fr1,'b');
title('Blocks Reconstructed')
subplot(4,2,3)
plot(noisy2,'r');
title('Bumps Original')
subplot(4,2,4)
plot(fr2, 'b');
title('Bumps Recontructed')
subplot(4,2,5)
plot(noisy3,'r');
title('Heavy Sine Original')
subplot(4,2,6)
plot(fr3,'b');
title('Heavy Sine Reconstructed')
subplot(4,2,7)
plot(noisy4,'r');
title('Doppler Original')
subplot(4,2,8)
plot(fr4, 'b');
title('Doppler Recontructed')

%%

[c1, l1] = wavedec(noisy1, 5, 'Db5');
fr1 = waverec(c1, l1, 'Db5');
 
[c2, l2] = wavedec(noisy2, 5, 'Db5');
fr2 = waverec(c2, l2, 'Db5');
 
 
[c3, l3] = wavedec(noisy3, 5, 'Db5');
 fr3 = waverec(c3, l3, 'Db5');
 
[c4, l4] = wavedec(noisy4, 5, 'Db5');
fr4 = waverec(c4, l4, 'Db5');

figure(4)
subplot(4,2,1)
plot(noisy1,'r');
title('Blocks Original')
subplot(4,2,2)
plot(fr1,'b');
title('Blocks Reconstructed')
subplot(4,2,3)
plot(noisy2,'r');
title('Bumps Original')
subplot(4,2,4)
plot(fr2, 'b');
title('Bumps Recontructed')
subplot(4,2,5)
plot(noisy3,'r');
title('Heavy Sine Original')
subplot(4,2,6)
plot(fr3,'b');
title('Heavy Sine Reconstructed')
subplot(4,2,7)
plot(noisy4,'r');
title('Doppler Original')
subplot(4,2,8)
plot(fr4, 'b');
title('Doppler Recontructed')