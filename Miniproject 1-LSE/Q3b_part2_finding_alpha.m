clc;
clear all;
close all;

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

% For constraint A*h_cest = b

A=zeros(26,L);  
b=zeros(26,1);

taps=5:5:30;
z=[];
for i=1:32
    if (~ismember(i,taps))
        z= [z,i];
    end
end


for i=1:26
    A(i,z(i))=1; % initialising  constraint matrix 
end

% Now new h vector contains only 6 non-zero taps
h(z)=0;

% % randomly generating 26 zero locations from the 32 available
% z=randperm(L);      
% tap=sort(z(27:end));  % 6 non-zero locations
% z=sort(z(1:26));      % 26 zero locations
% 
% 
% for i=1:26
%     A(i,z(i))=1; % initialising  constraint matrix 
% end
% 
% % Now new h vector contains only 6 non-zero taps
% h(z)=0;

% %Generating aplphas
alphas = 0.1:0.05:10;

MSE_list = [];

trials=1000;

    for i=1:trials
        
        se_a=[];


        %generating X
        data = randi([0 1],1024,1,'int8');

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
        
        %noise generation
        sigma=0.1;
        n=zeros(512,1);
        for ii = 1:512
            n(ii) = sigma*randn + 1i*sigma*randn;
        end

         for alpha = alphas
    
            % Ordinary Least Squares estimation
            y = X*F*h + n;
            Iden = eye(32);
            H = X*F;
            h_est = inv(H'*H + alpha*Iden)*H'*y;


            % Estimating using Lagrange Multiplier Technique
            lambda2=inv(A*(inv(H'*H + alpha*Iden))*A')*((A*h_est) - b);

            h_cest= h_est - (inv(H'*H + alpha*Iden)*A'*lambda2);


            %appending square of error to array
            
            se_a= [se_a, (h - h_cest)'*(h - h_cest)];
        


         end
         
    
    MSE_list= [MSE_list; se_a];
    end

MSE=mean(MSE_list);
plot(alphas, MSE);
title(['\alpha vs MSE (h_{est}))  for  sigma =  ', num2str(sigma)],'FontSize',13)
xlabel('\alpha','FontSize',13)
ylabel(' MSE (h_{est})','FontSize',13)
grid on;