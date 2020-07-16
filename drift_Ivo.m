%	---- High-pass filtering ----
% Yn=C1(Xn-X(n-1))+C2Y(n-1)
% C1 = 1/[1+tan(Fc*pi*T)]
% C2 = [1-tan(Fc*pi*T)]/[1+tan(Fc*pi*T)]
% Inputs:	X - signal;
%				Fc - cut-off frequency [Hz]
%				T -  sampling period [s]
% Output:	Y - filtered signal

function [Y]=drift_Ivo(X,Fc,T)

% make X as row
s=size(X);
if s(1)==1;
else X=X'; 
end

    
% Normalization of the input signal to the first point
X=X-X(1);

% Filter coefficients
c1 = 1/[1+tan(Fc*pi*T)];
c2 = [1-tan(Fc*pi*T)]/[1+tan(Fc*pi*T)];

b=[c1 -c1];	a=[1 -c2];
	Y1=filter(b,a,X);

Yf=fliplr(Y1);
Yf=Yf-Yf(1);
	Y=filter(b,a,Yf);
Y=fliplr(Y);

