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

% For constraint A*h_cest = b

A=zeros(26,L);  
b=zeros(26,1);

for i=1:26
    A(i,z(i))=1; % initialising  constraint matrix 
end

% Now new h vector contains only 6 non-zero taps
h(z)=0;
sigmas=[0.1 sqrt(0.1)];

H_cest=[];              %it stores the h_cest for both sigmas for future reference

for sigma= sigmas
    
    %noise generation
    n=zeros(512,1);
    for ii = (1:512)
        n(ii) = sigma*randn + 1i*sigma*randn;
    end

    y = X*F*h + n;

    H = X*F;
    h_est = inv(H'*H)*H'*y;


    % Estimating using Lagrange Multiplier Technique
    lambda2=inv(A*(inv(H'*H))*A')*((A*h_est) - b);

    h_cest= h_est - (inv(H'*H)*A'*lambda2);
    
    H_cest=[H_cest, h_cest];
%     %finding serror in estimation (h)
%     error= (h - h_cest)'*(h - h_cest);

    x= 1:32;

    %plotting h vs hcest
    figure
    set(gcf, 'Position',  [200,150, 800, 600])
    stem(x, abs(h.'),'MarkerFaceColor', 'red', 'MarkerSize',3,'Color','red',...
        'LineWidth',1)
    hold on
    plot(x,abs(h_cest.'),'.b','MarkerSize',12)
    % plot(x,abs(h.'),'.r',x,abs(h_cest.'),'.b','MarkerSize',10)
    title(['h vs h_{cest} for sigma =  ', num2str(sigma),' and non-zero taps  ', num2str(tap)] ,'FontSize',13)
    xlabel('h_{index}','FontSize',13)
    ylabel('estimate','FontSize',13)
    xlim([0,L])
    xticks(0:1:L)
    legend('h_{actual}','h_{cest}')
end