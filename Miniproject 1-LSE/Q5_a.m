clc;
clear all;
close all;

%generating random bits
data = randi([0 1],1024,1,'int8');
%modulating as QPPSK symbols to generate X
sym = nrSymbolModulate(data,'QPSK','OutputDataType','single');
sym = sym * sqrt(2);
X = diag(sym);

%generating F using meshgrid
L = 32;
s = 1:512;
t = 1:L;
[I, J] = meshgrid(s,t);
const = 2j*pi/512;
F = exp(const*(I-1).*(J-1));
F = F.';

%noise generation
n=zeros(512,1);
sigma = 0.1;
for ii = (1:512)
    n(ii) = sigma*randn + 1i*sigma*randn;
end

%generating h (multipath Rayleigh fading channel vector)
k = (1:L).';
lambda = 0.2;
p = exp(-1*lambda*(k-1));
a = normrnd(0, 0.5, [L,1]);
b = normrnd(0, 0.5, [L,1]);
h = (1/norm(p)) * (a + 1i*b).*p;


% randomly generating 26 zero locations from the 32 available
z=randperm(L);      
tap=sort(z(27:end));  % 6 non-zero locations
z=sort(z(1:26));      % 26 zero locations

% Now new h vector contains only 6 non-zero taps
h(z)=0;

y = X*F*h + n;        %observations
H = X*F;


%algorithm for finding non zero tap locations
k0=6;   %number of non zero taps
S=[];   %for position of non zero taps
A=X*F;  
r=y;    %resdual
t=0;
for k=1:k0
    temp=-Inf;       %used to compute argmax
    for j=1:32
        val=norm((A(:,j))'*r); %argmax over j
        if val>temp
            t=j;
            temp=val; 
        end
    end
    S(k)=t;    %appending tk to S
    B=A(:,S);
    P=B*pinv(B);  
    r=(eye(512)-P)*y;
end

S=sort(S);
display(S)

total=1:32;
zfound=setdiff(total,S);    %subtracting non zero taps from all taps to get zero taps

% For constraint A*h_cest = b
A=zeros(26,L);  
b=zeros(26,1);

for i=1:26
    A(i,zfound(i))=1; % initialising  constraint matrix 
end

% Ordinary Least Squares estimation

h_est = inv(H'*H)*H'*y;


% Estimating using Lagrange Multiplier Technique
lambda2=inv(A*(inv(H'*H))*A')*((A*h_est) - b);

h_cest= h_est - (inv(H'*H)*A'*lambda2);

%finding error in estimation
error= (h - h_cest)'*(h - h_cest);

x= 1:32;

%plotting h vs hcest
figure
plot(x,abs(h.'),'.r',x,abs(h_cest.'),'.b','MarkerSize',10)
title('h vs h_{cest}')
xlabel('h_{index}')
ylabel('estimate')
xlim([0,L])
xticks(0:1:L)
legend('h_{actual}','h_{cest}')
