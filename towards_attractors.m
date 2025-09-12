% Корреляционный анализ параметров пристеночной плазмы из
% экспериментальной кампании весны 2014 года на Т-10. 
% 
% by Mark 06.09.2025
% 
%% Data import
clear
clc
file_num = 7;
    
addpath(genpath('C:\Users\MSI\Desktop\Курчатовский институт\Т-10 Анализ рядов')) % вспомогательные функции и все такое
addpath(genpath('C:\Users\MSI\Documents\GitHub\T-10-Nonlinear-analysis\Wolf Lyapunov exp')) % алгоритм Вольфа для показателей Ляпунова
addpath(genpath('C:\Users\MSI\Documents\GitHub\T-10-Nonlinear-analysis\RecurrencePlot_ToolBox')) % FNN,Reccurence plot, mutual information)
addpath(genpath('C:\Users\MSI\Documents\GitHub\T-10-Nonlinear-analysis\alpha_stable_distribution')) % for Levi process test
time_series = readmatrix(['T10_66131_Lpf' num2str(file_num) '.txt']);
cd(['C:\Users\MSI\Documents\GitHub\T-10-Nonlinear-analysis']);
addpath(genpath('C:\Users\MSI\Documents\GitHub\T-10-Nonlinear-analysis'));
time = time_series(:,1);
fp = time_series(:,2);



%% Visualization 1
figure;

plot(time,fp,'k')
title("$Floating\ potential\ as\ time\ series$",Interpreter="latex")
xlabel("$Time,\ ms$",Interpreter="latex")
ylabel("$Floating\ potential,\ V$",Interpreter="latex")
set(gca, 'FontSize', 16, 'LineWidth', 2)
set(gcf, 'Color', 'white')




%% Analysis interval 

% Выбираем для анализа интервал по времени 
time_min = 400;
time_max = 600;
% Создаем маску для индексов
mask = (time >= time_min) & (time <= time_max);
% Извлекаем данные
time = time(mask);
fp = fp(mask);

%% Visualization 2
figure;

plot(time,fp,'k')
title("$Floating\ potential\ as\ time\ series$",Interpreter="latex")
xlabel("$Time,\ ms$",Interpreter="latex")
ylabel("$Floating\ potential,\ V$",Interpreter="latex")
set(gca, 'FontSize', 16, 'LineWidth', 2)
set(gcf, 'Color', 'white')

%% PDF

% Данные (нормированные)
x = (fp - mean(fp)) / std(fp);  

% Оценка плотности с помощью ядра
[f_pdf, xi] = ksdensity(x,NumPoints=length(x)/1000);


figure; hold on;


plot(xi, f_pdf, '.', 'MarkerSize', 12, 'Color', 'k');

% Гауссиан
x_gauss = linspace(min(xi), max(xi), 500);
y_gauss = (1/sqrt(2*pi)) * exp(-0.5 * x_gauss.^2);
plot(x_gauss, y_gauss, 'b', 'LineWidth', 2);

xlabel('$Normalized\ floating\ potential\ (V)$', 'Interpreter','latex');
ylabel('$Probability\ density$', 'Interpreter','latex');
title('$PDF\ of\ normalized\ floating\ potential\ vs\ standard\ Gaussian$', 'Interpreter','latex');

set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');
legend({'$Experimental\ PDF$ ','$Standard\ Gaussian$'}, 'Interpreter','latex');
grid on


%% Fourier Spectrum

dt = mean(diff(time));
Fs = 1/dt;                % частота дискретизации
N = length(fp);           % число точек
FP_fft = fft(fp);            % прямое преобразование Фурье
FP_mag = abs(FP_fft)/N;      % амплитуда (нормируем)
f = (0:N-1)*(Fs/N);          % частоты от 0 до Fs
figure;
loglog(f(1:floor(N/2)), 2*FP_mag(1:floor(N/2)), 'k', 'LineWidth', 1.5); 
labels = [title('$Fourier\ spectrum\ of\ signal\ (log-log\ scale)$', 'FontSize', 14), ...
          xlabel('$f, Hz$', 'FontSize', 14), ...
          ylabel('$s(f), rel.units$', 'FontSize', 14)];
          
set(labels, 'Interpreter', 'latex');
set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');

%% ACF


fp0 = fp - mean(fp);

[acf, lags] = xcorr(fp0, 'coeff');  % нормировка на 0-й лаг
lags_time = lags * dt;


figure;
plot(lags_time, acf, 'k', 'LineWidth', 1.5);

labels = [title('$ACF\ of\ the\ signal$', 'FontSize', 14), ...
          xlabel('$\tau, \mu s$','FontSize',14), ...
          ylabel('$R(\tau)$','FontSize',14);];

set(labels, 'Interpreter', 'latex');
set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');
grid on;

xlim([-5 5]) % +- 5 мкс ось

%% Correlation integral 

L = 9;                            % embedding depth

p_values = 30:10:90;  % произвольные лаги ( в индексах)

M = 200;                           % число вложенных векторов

x = fp(:);                        % временной ряд

colors = lines(length(p_values));  % набор цветов

figure; hold on;

for k = 1:length(p_values)
    p = p_values(k);
    
   
    [r, C] = correlation_integral(x, L, p, M);  
    
   
    loglog(r, C,'o' , 'Color', colors(k,:), 'MarkerSize', 5, 'LineWidth', 1.5);
end

xlabel('$log\ r$', 'Interpreter','latex');
ylabel('$log\ C(r)$', 'Interpreter','latex');
title(['Correlation integral, L = ' num2str(L) ', M = ' num2str(M)], 'Interpreter','latex');
grid on;

legendStrings = arrayfun(@(p) ['p = ' num2str(p)], p_values, 'UniformOutput', false);
legend(legendStrings, 'Interpreter','latex', 'Location','best');

set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');

%% XY Attractor
p = 10;
M = 50;
Xplot = XY_Attractor(x,9,p,1000000); % XY-projection of trajectory in phase-space

figure
plot(Xplot(:,1), Xplot(:,2), 'k.', 'MarkerSize',10);
xlabel('$X$', 'Interpreter','latex');
ylabel('$Y$', 'Interpreter','latex');
title('XY projection of attractor', 'Interpreter','latex');


xlim_vals = xlim;
ylim_vals = ylim;
xpos = xlim_vals(2) - 0.2*diff(xlim_vals);
ypos = ylim_vals(2) - 0.08*diff(ylim_vals);

txt = sprintf('L(dimensions) = %d \n p(lag) = %d \n M(j-iterator) = %d', L, p, M);
text(xpos, ypos, txt, 'FontSize', 18, 'Color', 'k', 'Interpreter','latex');


set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');
grid on;

%% PCA 2D

% --- PCA для снижения размерности до 2 ---
[~, score, ~] = pca(Xplot);  % score - данные в главных компонентах
X2D = score(:,1:2);              % первые 2 главные компоненты

% --- 3D-проекция через PCA ---
figure;
plot(X2D(:,1), X2D(:,2), 'k.', 'MarkerSize', 10);
xlabel('PC1', 'Interpreter','latex');
ylabel('PC2', 'Interpreter','latex');
title(['2D PCA projection of attractor, L = ' num2str(L) ', p = ' num2str(p)], 'Interpreter','latex');
grid on;

set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');

%% XYZ attractor

% --- 3D-проекция первых трёх координат ---
figure;
plot3(Xplot(:,1), Xplot(:,2), Xplot(:,3), 'k.', 'MarkerSize', 10);
xlabel('X', 'Interpreter','latex');
ylabel('Y', 'Interpreter','latex');
zlabel('Z', 'Interpreter','latex');
title(['3D projection of attractor, L = ' num2str(L) ', p = ' num2str(p) ', M = ' num2str(M)], 'Interpreter','latex');
grid on;

view(45,30); % угол обзора
set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');

%% PCA 3D

% --- PCA для снижения размерности до 3 ---
[coeff, score, ~] = pca(Xplot);  % score - данные в главных компонентах
X2D = score(:,1:3);              % первые три главные компоненты

% --- 3D-проекция через PCA ---
figure;
plot3(X2D(:,1), X2D(:,2), X2D(:,3), 'k.', 'MarkerSize', 10);
xlabel('PC1', 'Interpreter','latex');
ylabel('PC2', 'Interpreter','latex');
zlabel('PC3', 'Interpreter','latex');
title(['3D PCA projection of attractor, L = ' num2str(L) ', p = ' num2str(p)], 'Interpreter','latex');
grid on;
view(45,30);
set(gca,'FontSize',16,'LineWidth',2);
set(gcf,'Color','white');

%%  MATLAB functions test 


eRange = 20;
lyapunovExponent(fp,Fs,10,L,'ExpansionRange',eRange)

[~,eLag,eDim] = phaseSpaceReconstruction(fp);
phaseSpaceReconstruction(fp,10,9);
% медленно работают, энтропию посчитать не тянут. Показатели ляпунова тоже
% не могут. Поэтому подтягиваем сторонние тулбоксы 

%% Data prep for Lyapunov exponents by Wolf et al.

fname = 'data.txt';
writematrix(fp,fname)

%% Wolf's Algo for Lyapunov Exponents

datcnt = 16340;
tau = 10;
ndim = 3;
ires = 10;
maxbox = 6000;

db = basgen(fname, tau, ndim, ires, datcnt, maxbox);

% dt = 1e-6;
evolve = 20;
dismin = 0.001;
dismax = 0.3;
thmax = 30;

[out, SUM] = fet(db, dt, evolve, dismin, dismax, thmax);

makeplot(db, out, evolve, 'NorthWest') % в figure будет значение старшей экспоненты Ляпунова
% либо последняя строка последний столбик в fetout.txt


%% Levi process test
% Проверка временного ряда на процесс Леви


% 1. Приращения
dX = diff(fp);

% 2. Проверка независимости (автокорреляция приращений)
figure;
subplot(2,2,1)
autocorr(dX)
title('Autocorrelation of increments')

subplot(2,2,2)
parcorr(dX)
title('Partial autocorrelation of increments')

% 3. Проверка стационарности приращений (визуально)
subplot(2,2,3)
histogram(dX(1:floor(end/2)), 'Normalization','pdf')
hold on
histogram(dX(floor(end/2)+1:end), 'Normalization','pdf')
legend('1-я половина','2-я половина')
title('Сравнение распределений приращений')

% 4. Подгонка устойчивого распределения
% Для этого нужен Stable Distribution Toolbox с FileExchange:
% https://www.mathworks.com/matlabcentral/fileexchange/37514-stable-distribution
% (функции stabfit, stblrnd и т.п.)

try
    params = stblfit(dX);  % [alpha, beta, gamma, delta]
    alpha = params(1); beta = params(2);
    gamma = params(3); delta = params(4);

    fprintf('Оцененные параметры устойчивого распределения:\n');
    fprintf('alpha = %.3f, beta = %.3f, gamma = %.3f, delta = %.3f\n',...
        alpha, beta, gamma, delta);

    % 5. Проверка согласия (К-С тест)
    pd = makedist('Stable','alpha',alpha,'beta',beta,...
                           'gam',gamma,'delta',delta);
    [h,p] = kstest(dX,'CDF',pd);

    if h == 0
        fprintf('KS-тест: распределение приращений совместимо с устойчивым (p=%.4f)\n', p);
    else
        fprintf('KS-тест: распределение приращений НЕ совместимо с устойчивым (p=%.4f)\n', p);
    end

catch
    warning(['Не найден Stable Distribution Toolbox. ', ...
             'Скачайте с File Exchange для оценки параметров.']);
end

%% additional tests for levi process
[h, pValue] = lbqtest(dX);

histogram(dX,'Normalization','pdf');
set(gca,'YScale','log');
set(gca,"XScale",'log')