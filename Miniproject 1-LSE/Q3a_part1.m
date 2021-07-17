clc;
clear all;
close all;

%generating random bits
data = randi([0 1],1024,1,'int8');
%modulating as QPPSK symbols to generate X
sym = nrSymbolModulate(data,'QPSK','OutputDataType','single');
sym = sym * sqrt(2);
X = diag(sym);

% This is the only change from Q1 i.e introducing guard band 
%We set first and last 180 elements to 0
for i=1:180
    X(i,i)=0;
end

for i=(512-180 +1):512
    X(i,i)=0;
end

%generating F using meshgrid
L = 32;
s = 1:512;
t = 1:L;
[I, J] = meshgrid(s,t);
const = 2j*pi/512;
F = exp(const*(I-1).*(J-1));
F = F.';


%generating h (multipath Rayleigh fading channel vector)
k = (1:L).';
lambda = 0.2;
p = exp(-1*lambda*(k-1));
a = normrnd(0, 0.5, [L,1]);
b = normrnd(0, 0.5, [L,1]);
h = (1/norm(p)) * (a + 1i*b).*p;


sigmas=[0.1 sqrt(0.1)];

H_est=[];              %it stores the h_est for both sigmas for future reference

for sigma= sigmas
    
    %noise generation
    n=zeros(512,1);
    for ii = (1:512)
        n(ii) = sigma*randn + 1i*sigma*randn;
    end

    y = X*F*h + n;

    H = X*F;
    h_est = inv(H'*H)*H'*y;
    
    H_est=[H_est, h_est];


    x= 1:32;

    %plotting h vs hest
    figure
    set(gcf, 'Position',  [200,150, 800, 600])
    stem(x, abs(h.'),'MarkerFaceColor', 'red', 'MarkerSize',3,'Color','red',...
        'LineWidth',1)
    hold on
    plot(x,abs(h_est.'),'.b','MarkerSize',12)
    % plot(x,abs(h.'),'.r',x,abs(h_cest.'),'.b','MarkerSize',10)
    title(['h vs h_{est} for sigma =  ', num2str(sigma)] ,'FontSize',13)
    xlabel('h_{index}','FontSize',13)
    ylabel('estimate','FontSize',13)
    xlim([0,L])
    xticks(0:1:L)
    legend('h_{actual}','h_{cest}')
end

