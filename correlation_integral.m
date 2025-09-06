function [r, C] = correlation_integral(x, L, p, M)
    % correlationIntegral - считает корреляционный интеграл C(r)
    % по временному ряду x с вложением L, лагом p и числом векторов M
    %
    % Вход:
    %   x - временной ряд (вектор)
    %   L - размерность вложения
    %   p - лаг
    %   M - число векторов (обычно N - (L-1)*p) 
    %
    % Выход:
    %   r - радиусы
    %   C - значения корреляционного интеграла

    x = x(:);
    N = length(x);

    if nargin < 4
        M = N - (L-1)*p; % если не задано, берём максимум
    end

    % --- Формируем L-мерные векторы ---
    X = zeros(M, L);
    for m = 0:L-1
        X(:, m+1) = x((1:M) + m*p);
    end

    % --- Диапазон радиусов ---
    r = linspace(0, max(x)-min(x), 200);
    C = zeros(size(r));

    % --- Вычисляем C(r) ---
    for k = 1:length(r)
        cnt = 0;
        for i = 1:M
            for j = 1:M
                d = sum(abs(X(i,:) - X(j,:)));  % норма L1
                if d <= r(k)
                    cnt = cnt + 1;
                end
            end
        end
        C(k) = cnt / (N * M);
    end

    % --- Визуализация ---
    % figure;
    % loglog(r, C, 'k.-', 'MarkerSize', 12, 'LineWidth', 1.5);
    % xlabel('$r$', 'Interpreter','latex');
    % ylabel('$C(r)$', 'Interpreter','latex');
    % title('Correlation integral (log-log)', 'Interpreter','latex');
    % grid on;
    % 
    % set(gca,'FontSize',16,'LineWidth',2);
    % set(gcf,'Color','white');
end
