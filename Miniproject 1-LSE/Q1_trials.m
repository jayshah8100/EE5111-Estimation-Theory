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


sigmas=[0.1 sqrt(0.1)];


for sigma=sigmas
    
    se_y=[];
    se_h=[];
    mse=[];
    trials=1000;
    hf_est=zeros(32,trials);



    for i=1:trials

        %generating X
        data = randi([0 1],1024,1,'int8');

        sym = nrSymbolModulate(data,'QPSK','OutputDataType','single');
        sym = sym * sqrt(2);
        X = diag(sym);  


        %noise generation
        n=zeros(512,1);
        for ii = 1:512
            n(ii) = sigma*randn + 1i*sigma*randn;
        end


        % Ordinary Least Squares Estimation

        y = X*F*h + n; %observations

        H = X*F;
        h_est = inv(H'*H)*H'*y;


        %appending square of error to array
        se_y= [se_y, (y - H*h_est)'*(y - H*h_est)];
        se_h= [se_h, (h - h_est)'*(h - h_est)];
        
      
        mse_theory_h= (sigma.^2)*(inv(H'*H));
        mse=[mse, mse_theory_h(1,1)];
        hf_est(:,i)=h_est;

    end
    %finding simualted error in estimation (y)
    mse_y=mean(se_y);
    fprintf('\nThe simulated mean square error in y (observations) averaged over %d trials for sigma = %f is %f\n', trials,sigma,mse_y); 
 
    %finding simulated error in estimation (h)
    mse_h=mean(se_h);
    fprintf('\nThe simulated mean square error in h (input vector) averaged over %d trials for sigma = %f is %f\n', trials,sigma,mse_h); 
    
    %finding theoretical error in estimation (h)
    mse_theoretical=mean(mse);
    fprintf('\nThe theoretical mean square error in h (input vector) for sigma = %f is %f\n',sigma,mse_theoretical); 

    expected_h_cest=mean(hf_est,2);
    x= 1:32;

    %plotting h vs E(hcest)
    figure
    set(gcf, 'Position',  [200,150, 800, 600])
    plot(x,abs(h.'),'*r',x,abs(expected_h_cest.'),'.b','MarkerSize',10)
    title(['h vs E(h_{est})  for  sigma =  ', num2str(sigma)],'FontSize',13)
    xlabel('h_{index}','FontSize',13)
    ylabel('estimate','FontSize',13)
    xlim([0,L])
    xticks(0:1:L)
    legend('h_{actual}','E(h_{est})')


    fprintf("\nFrom the plot it can be seen that expected h_cest is close to h. Hence validating that h_cest is an unbiased estimator\n\n\n")    


end
