clear all
clc

rng('default')
%% the truth line (represented as a cosin line)
omiga_1 = 0.002 * pi; phi = 0.2; T_1 = 2 * pi ./ omiga_1;
omiga_2 = 0.05 * pi; phi=0; T_2 = 2 * pi ./ omiga_2;
len = 350;

for i=1:1:len
    if i <= 50
        true_line(i) = cos(omiga_1 * i + phi * pi) + 1 ; % cos(wt+b) + 1
    elseif i > 50  && i <= 90
        true_line(i) = cos(omiga_1 * i + phi * pi) + 1   + 1* cos(omiga_2 * i + pi); %2
    elseif i > 90 && i <= 150
        true_line(i) = cos(omiga_1 * i + phi * pi) + 1;
    elseif i > 150 && i <= 170
        true_line(i) = cos(omiga_1 * i + phi * pi) + 1   + 1*cos(2*omiga_2 * i + 0.5*pi); %1.5
    elseif i> 170
        true_line(i) = cos(omiga_1 * i + phi * pi) + 1   + 0.75* (300/(4*i)) *cos(8*omiga_2 * i + 0.5*pi);
    end
end

mutiply = 5;
true_1 = mutiply.*true_line;

t_svar = 2; % if the noise grows bigger, the effect is more clear.
t_noise = normrnd(0, t_svar, [1, len]); 

true_2 = mutiply.*true_line + t_noise; % add the noise

%% the line of sub-models
n_std = 1;
% 1. lag (i.e. lay 10 points)
for i=1:1:len
    if i <= 50
        m_1(i) = cos(omiga_1 * i + phi * pi) + 1 ; % cos(wt+b) + 1
    elseif i > 50 && i <= 90
        m_1(i) = cos(omiga_1 * i + phi * pi) + 1 +   1* cos(omiga_2 * i + pi);
    elseif i > 90 && i <= 150
        m_1(i) = cos(omiga_1 * i + phi * pi) + 1;
    elseif i > 150 && i <= 170
        m_1(i) = cos(omiga_1 * i + phi * pi) + 1 +   1*cos(2*omiga_2 * i + 0.5*pi);
    elseif i> 170
        m_1(i) = cos(omiga_1 * i + phi * pi) + 1;
    end
end

n_1 = normrnd(0, n_std, [1, len]);
m_1 = mutiply .* [m_1(1:10) m_1];
m_1 = m_1(1:len) + n_1;

% 2. random error (twice the noise)
for i=1:1:len
    if i <= 50
        m_2(i) = cos(omiga_1 * i + phi * pi) + 1 ; % cos(wt+b) + 1
    elseif i > 50 && i <= 90
        m_2(i) = cos(omiga_1 * i + phi * pi) + 1 +   1* cos(omiga_2 * i + pi);
    elseif i > 90 && i <= 150
        m_2(i) = cos(omiga_1 * i + phi * pi) + 1;
    elseif i > 150 && i <= 170
        m_2(i) = cos(omiga_1 * i + phi * pi) + 1 +   1*cos(2*omiga_2 * i + 0.5*pi);
    elseif i> 170
        m_2(i) = cos(omiga_1 * i + phi * pi) + 1;
    end
end

n_2 = 2.*normrnd(0, n_std, [1, len]);
m_2 = mutiply .*m_2 + n_2;

% 3. smoothing sub-model
for i=1:1:len
    if i <= 50
        m_3(i) = cos(omiga_1 * i + phi * pi) + 1 ; % cos(wt+b) + 1
    elseif i > 50 && i <= 90
        m_3(i) = cos(omiga_1 * i + phi * pi) + 1 +   0.25* 1*cos(omiga_2 * i + pi);
    elseif i > 90 && i <= 150
        m_3(i) = cos(omiga_1 * i + phi * pi) + 1;
    elseif i > 150 && i <= 170
        m_3(i) = cos(omiga_1 * i + phi * pi) + 1 +   0.25 * 1*cos(2*omiga_2 * i + 0.5*pi);
    elseif i> 170
        m_3(i) = cos(omiga_1 * i + phi * pi) + 1 ;
    end
end
n_3 = 0.5 .* normrnd(0, n_std, [1, len]);
m_3 = mutiply .*m_3 + n_3;

% static error sub-model
for i=1:1:len
    if i <= 50
        m_4(i) = cos(omiga_1 * i + phi * pi) + 1 ; % cos(wt+b) + 1
    elseif i > 50 && i <= 90
        m_4(i) = cos(omiga_1 * i + phi * pi) + 1 +   1* cos(omiga_2 * i + pi);
    elseif i > 90 && i <= 150
        m_4(i) = cos(omiga_1 * i + phi * pi) + 1;
    elseif i > 150 && i <= 170
        m_4(i) = cos(omiga_1 * i + phi * pi) + 1 +   1*cos(2*omiga_2 * i + 0.5*pi);
    elseif i> 170
        m_4(i) = cos(omiga_1 * i + phi * pi) + 1 ;
    end
end
m_4 = mutiply .*m_4 + 1 + normrnd(0, n_std, [1, len]);


%% static estimation for the weighting vector (based on the measurement equation)
window_size = 10; member_model = [m_1; m_2; m_3; m_4]; 
w_s = static_est(member_model, true_2, window_size);

% test the weighting vector used for model combining
t_static(1) = true_1(1);
for i = 2:len
    t_static(i) = w_s(:, i-1)' * member_model(: ,i);
end

%% dynamic estimation for the weighting vector
r_svar = 0.1; s_svar= 0.1;
w_d = pf_dynamic_est(w_s, r_svar, s_svar);
t_dynamic(1) = true_1(1);
for i = 2:len
    t_dynamic(i) = w_d(:, i-1)' * member_model(: ,i);
end

%% result plot
figure()
plot(true_1, '-' ,'LineWidth', 2, 'color', [255, 0, 0]./255)
hold on 
plot(m_1, '-.',  'LineWidth', 1.5, 'color',  [0, 0, 255]./255)
set(gca,'FontSize',18);
xlabel('time k')
ylabel('value')
ylim([0, 20])
xlim([1, 350])
legend('true value', 'lag error')
set(gca,'FontSize',24);
set(gcf,'position',[0 0 1600 600])
print(gcf,'submodel1','-dpng','-r300')

figure()
plot(true_1, '-' ,'LineWidth', 2, 'color', [255, 0, 0]./255)
hold on 
plot(m_2, '-.',  'LineWidth', 1.5, 'color',  [0, 0, 255]./255)
set(gca,'FontSize',18);
xlabel('time k')
ylabel('value')
ylim([0, 20])
xlim([1, 350])
legend('true value', 'random error')
set(gca,'FontSize',24);
set(gcf,'position',[0 0 1600 600])
print(gcf,'submodel2','-dpng','-r300')

figure()
plot(true_1, '-' ,'LineWidth', 2, 'color', [255, 0, 0]./255)
hold on 
plot(m_3, '-.',  'LineWidth', 1.5, 'color',  [0, 0, 255]./255)
set(gca,'FontSize',18);
xlabel('time k')
ylabel('value')
ylim([0, 20])
xlim([1, 350])
legend('true value', 'smoothing error')
set(gca,'FontSize',24);
set(gcf,'position',[0 0 1600 600])
print(gcf,'submodel3','-dpng','-r300')

figure()
plot(true_1, '-' ,'LineWidth', 2, 'color', [255, 0, 0]./255)
hold on 
plot(m_4, '-.',  'LineWidth', 1.5, 'color',  [0, 0, 255]./255)
set(gca,'FontSize',18);
xlabel('time k')
ylabel('value')
ylim([0, 20])
xlim([1, 350])
legend('true value', 'static error')
set(gca,'FontSize',24);
set(gcf,'position',[0 0 1600 600])
print(gcf,'submodel4','-dpng','-r300')

figure()
plot(true_1, 'r', 'LineWidth', 2)
hold on 
plot(t_static(1,:), 'go-', 'LineWidth', 1.5)
plot(t_dynamic(1,:), 'bs-', 'LineWidth', 1.5)
ylim([0, 20])
xlim([1, 350])
xlabel('k')
ylabel('value')
legend('true value','least square method', 'particle filter')
set(gca,'FontSize',24);
set(gcf,'position',[0 0 1600 600])
print(gcf,'performance','-dpng','-r300')

%% error indices
[m1_rmse, m1_mae, m1_mape] = error_indices(m_1, true_1)
[m2_rmse, m2_mae, m2_mape] = error_indices(m_2, true_1)
[m3_rmse, m3_mae, m3_mape] = error_indices(m_3, true_1)
[m4_rmse, m4_mae, m4_mape] = error_indices(m_4, true_1)

[s_rmse, s_mae, s_mape] = error_indices(t_static, true_1)
[d_rmse, d_mae, d_mape] = error_indices(t_dynamic, true_1)

% csv
datacolumns = {'m1','m2','m3','m4', 'lsm', 'pf'};
data = table([m1_rmse, m1_mae, m1_mape]', [m2_rmse, m2_mae, m2_mape]', [m3_rmse, m3_mae, m3_mape]', [m4_rmse, m4_mae, m4_mape]', ...
    [s_rmse, s_mae, s_mape]', [d_rmse, d_mae, d_mape]', 'VariableNames', datacolumns);
writetable(data, 'forecasting performance.csv')

%% related functions
function w_v = static_est(member_model, true_line, window_size)
% least square method
% window_size is the number of points used to estimate
% the initial window_size-1 points of weighting vector is 1/numl(member models)
m1 = member_model(1,:);
m2 = member_model(2,:);
m3 = member_model(3,:);
m4 = member_model(4,:);
tr = true_line;

w_v = zeros(size(member_model));

w_v(:, 1:(window_size-1)) = 1./size(member_model, 1);

for i =  window_size: size(member_model, 2)
    y = tr(i-window_size+1:i)';
    x = member_model(: , i-window_size+1: i)';
    w_v (:, i) = inv(x' * x) * x' * y;
end
end

function w = pf_dynamic_est(z, r_svar, s_svar)
% particle fiter 
% w is the weighting vector for estimation
% r and s is the state and measurement noise covariance
particles = 200;
R = r_svar^2;

w(:,1) = z(:, 1); 

% generate particles 
pr = zeros(size(z,1), particles);

for i = 1 : particles
    pr(:,i) = z(:, 1) + normrnd(0, r_svar^2, [size(z, 1), 1]);
end

for i = 2 : size(z, 2)
    % sample (based on state equation)
    for j = 1 : particles
        xpr_(:, j) = pr(:, j) + normrnd(0, r_svar^3, [size(z, 1), 1]);
    end
    
    % compute weight (based on measurement equation)
    for j = 1 : particles
        zpr_(:, j) = xpr_(:, j) + normrnd(0, s_svar, [size(z, 1), 1]);
        weight(i, j) = sqrt((2*pi)^size(z,1))*sqrt(det(R))*...,
            exp(-0.5*(z(:,i)-zpr_(:,j))'*inv(R)*(z(:,i)-zpr_(:,j))) + 1e-99;
    end
    
    % normalize
    weight(i,:) = weight(i,:)./sum(weight(i,:));
    
    % resample
    outIndex = randomR(weight(i,:));
    pr = xpr_(:, outIndex);
    
    w(:,i) = mean(pr,2);
     
end    
end


function outIndex = randomR(weight)
% random resampling
L=length(weight);
outIndex=zeros(1,L);
u=unifrnd(0,1,1,L);
u=sort(u);
cdf=cumsum(weight);
i=1;
for j=1:L
    while (i<=L) & (u(i)<=cdf(j))
        outIndex(i)=j;
        i=i+1;
    end
end
end

function [rmse, mae, mape] = error_indices(pred, true)
% error indices
rmse = sqrt(sum((true - pred).^2)./size(true, 2));
mae = sum(abs(true-pred))./size(true, 2);
mape = sum(abs((true - pred)./true))/size(true, 2);
end