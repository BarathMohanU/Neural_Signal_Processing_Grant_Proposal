%%%%%%%%%%%%%%%%%%%%% Barath Mohan U & Sreekar Gunda %%%%%%%%%%%%%%%%%%%%%%

%Getting unipolar data
oNum=0; refType = 'unipolar';
[data,layout] = getData(oNum, refType); %function in EEGDataResources folder

%Frequency analysis
cfg                 = [];
cfg.output          = 'pow';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.tapsmofrq       = 4;  
cfg.foilim          = [0 50];
cfg.keeptrials      = 'yes';
cfg.pad             = 'nextpow2'; 
data_freq           = ft_freqanalysis(cfg, data);

x=data_freq.powspctrm;
x=reshape(x,[188,103*64]);
x=ica(x); %preprocessing - ICA
x=x-mean(x,1);
x=pca(x); %feature selection - PCA

%Creating a test set
test=x([28,29,53,54,77,78,99,100,125,126,146,147,168,169,187,188],:);
x=x([1:27,30:52,55:76,79:98,101:124,127:145,148:167,170:186],:);

%Creating y values
y=ones(172,1);
y(55:76)=2;
y(101:124)=3;
y(148:167)=4;

n=size(x,2);
m=size(x,1);
num=2;

Y=zeros(m,num);
for i=1:m
    if y(i)==1
        Y(i,1)=0;
        Y(i,2)=0;
    elseif y(i)==2
        Y(i,1)=1;
        Y(i,2)=0;
    elseif y(i)==3
        Y(i,1)=0;
        Y(i,2)=1;
    elseif y(i)==4
        Y(i,1)=1;
        Y(i,2)=1;
    end
end

x1=[ones(m,1) x];

%randomly initialising theta
epi=0.5;
Theta1=rand(n,n+1)*(2*epi)-epi;
Theta2=rand(n,n+1)*(2*epi)-epi;
Theta3=rand(num,n+1)*(2*epi)-epi;
Theta=[Theta1(:) ; Theta2(:) ; Theta3(:)];
lambda=0.05;%regularization parameter
alpha=0.9;%learning rate

%gradient descent
for i=1:30000
    [J,grad] = costFunction(m,n,Theta,x1,Y,lambda,num);
    Theta=Theta-alpha*grad;
end

%prediction on test set
p=predict(16,n,Theta,test,num);

function [z] = ica(x)
mu = mean(x,2);
T = sqrtm(inv(cov(x')));
Zcw = T * bsxfun(@minus,x,mu);
[W, ~, ~] = svd(bsxfun(@times,sum(Zcw.^2,1),Zcw) * Zcw');
z= W * Zcw;
end

function [z] = pca(x)
n=size(x,2);
m=size(x,1);
sigma=(1/m)*(x'*x);
[u,s,~]=svd(sigma);
for k=1:n
    if (sum(diag(s(1:k,1:k)))/sum(diag(s)))>=0.99
        break;
    end
end
u=u(:,1:k);
z=x*u;
end

function [p] = predict(m,n,Theta,test,num)

Theta1=reshape(Theta(1:n*(n+1)),n,n+1);
Theta2=reshape(Theta(1+n*(n+1):2*n*(n+1)),n,n+1);
Theta3=reshape(Theta(1+2*n*(n+1):end),num,n+1);
test=[ones(m,1) test];
p=sig([ones(m,1) sig([ones(m,1) sig(test*Theta1')]*Theta2')]*Theta3');

end

function [J,grad] = costFunction(m,n,Theta,x1,Y,lambda,num) 

Theta1=reshape(Theta(1:n*(n+1)),n,n+1);
Theta2=reshape(Theta(1+n*(n+1):2*n*(n+1)),n,n+1);
Theta3=reshape(Theta(1+2*n*(n+1):end),num,n+1);

h=sig([ones(m,1) sig([ones(m,1) sig(x1*Theta1')]*Theta2')]*Theta3');

J=((-1/m)*sum(sum((Y.*log(h))+((1-Y).*log(1-h)))))+(lambda/(2*m))*(sum(sum(Theta1(:,2:n+1).^2))+sum(sum(Theta2(:,2:n+1).^2))+sum(sum(Theta3(:,2:n+1).^2)));

delta4=h-Y;
delta3=(delta4*Theta3).*[ones(m,1) sig([ones(m,1) sig(x1*Theta1')]*Theta2')].*(1-[ones(m,1) sig([ones(m,1) sig(x1*Theta1')]*Theta2')]);
delta3=delta3(:,2:n+1);
delta2=(delta3*Theta2).*[ones(m,1) sig(x1*Theta1')].*(1-[ones(m,1) sig(x1*Theta1')]);
delta2=delta2(:,2:n+1);
theta1=Theta1;
theta2=Theta2;
theta3=Theta3;
theta1(:,1)=0;
theta2(:,1)=0;
theta3(:,1)=0;
Theta1_grad=(1/m)*(delta2'*x1+lambda*theta1);
Theta2_grad=(1/m)*(delta3'*[ones(m,1) sig(x1*Theta1')]+lambda*theta2);
Theta3_grad=(1/m)*(delta4'*[ones(m,1) sig([ones(m,1) sig(x1*Theta1')]*Theta2')]+lambda*theta3);
grad=[Theta1_grad(:) ; Theta2_grad(:) ; Theta3_grad(:)];
end

function g = sig(z)
g = 1.0 ./ (1.0 + exp(-z));
end

