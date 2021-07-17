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


% %Generating aplphas
alphas = 0.1:0.5:20;

MSE_list = [];

trials=1000;

    for i=1:trials
        
        se_a=[];
  
%         hf_est=zeros(32,trials);


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
    
            % Ordinary Least Squares Estimation

            y = X*F*h + n; %observations
            Iden = eye(32);

            H = X*F;
            h_est = inv(H'*H + alpha*Iden)*H'*y;

            %appending square of error to array
            se_a= [se_a, (h - h_est)'*(h - h_est)];
        

%         hf_est(:,i)=h_est;

         end
         
    %finding simulated error in estimation (h)
    %     mse_h=mean(se_h);
    %     fprintf('\nThe simulated mean square error in h (input vector) averaged over %d trials for aplha = %f is %f\n', trials,alpha,mse_h); 

    %     %finding theoretical error in estimation (h)
    %     mse_theoretical=mean(mse);
    %     fprintf('\nThe theoretical mean square error in h (input vector) for alpha = %f is %f\n',alpha,mse_theoretical); 
    
    MSE_list= [MSE_list; se_a];
    end

MSE=mean(MSE_list);
plot(alphas, MSE);
title(['\alpha vs MSE (h_{est}))  for  sigma =  ', num2str(sigma)],'FontSize',13)
xlabel('\alpha','FontSize',13)
ylabel(' MSE (h_{est})','FontSize',13)
grid on;