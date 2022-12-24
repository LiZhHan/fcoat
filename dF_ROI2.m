%%% Jia et al_2011_In vivo two-photon imaging of 
% sensory-evoked dendritic calcium signals in cortical neurons.pdf

dbstop if error
clear 
close all

ca_data = load('Data path');
% Time window parameter setting
[N,M] = size(ca_data); % Capture saved data A through imageJ ROI
E = 90/(N-1);
e = E;
tao0 = 10*e;  
tao1 = 10*e; % neuron1
tao2 = 100*e;

t = 0:E:E*(N-1);
Ft = ca_data(:,4)-200;%
Ft = Ft/max(Ft);

F = 1/N*sum(Ft(:)); %ROI 
Y = zeros(N,1);
Fxmean = zeros(N,1);
for i = 1:N
    tao = (t(i)-tao1/2):E:(t(i)+tao1/2);
    tao(tao<0)=[];
    Tao = size(tao,2);
    for k = 1:Tao
        if isempty(tao)
           n = [];
        else
           [~,n] = min(abs(t-tao(k))); 
        end
        if isempty(n)
        else
           Y(k) = Ft(n);
        end
     end
        Y(Y == 0)=[];
        S = size(Y,1);
        Fxmean(i) = sum(Y(:))/S;
end
F0 = min(Fxmean);
Rt = (Ft-F0)./F0;
F1 = zeros(N,1);
F2 = ones(N,1);
f1 = zeros(N,1);
f2 = zeros(N,1);
for i = 1:N
    for j = 1:i
        if i-j==0
            f1(j) = 0;
        else
        f1(j) = Rt(i-j).*exp(-abs(t(j))/tao0);
        f2(j) = exp(-abs(t(j))/tao0);
        end
    end
    F1(i) = sum(f1(:));
    F2(i) = sum(f2(:));
end

dF_F = F1./F2;
figure(1)
plot(dF_F);
xlabel('time');
ylabel('deta_F / F')

%% Baseline leveling
x1 = 1:N;

%% Select the starting point of calcium signal
x = [4 54 104 154 204 254 304 354 404 454 504 554 590];
y = dF_F(x);

%%

y1=interp1(x,y,x1, 'PCHIP');
error_value = dF_F'-y1;
error_value = error_value/max(max(error_value));
[count, x]=imhist(error_value); 

x_num = find(error_value<0.05); 
figure(3)
plot(Rt);
hold on
plot(y1);
hold on
plot(x_num,dF_F(x_num),'*');

y_num = dF_F(x_num);

y2=interp1(x_num, y_num,x1, 'pchip');
figure(4)
plot(dF_F);
hold on
plot(y2);

Rt_co = Rt'-y2;
gs_75_122 = zeros(N,2);
gs_75_122(:,1) = t';
gs_75_122(:,2) = Rt_co';
figure(5),plot(t,Rt_co);
xlabel('time');
ylabel('deta_F / F')
csvwrite('Data path 1',Rt_co);






