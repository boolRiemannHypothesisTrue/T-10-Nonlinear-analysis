function Xplot = XY_Attractor(x, L, p,Npoints,M)
    % plotXYAttractor - строит XY проекцию аттрактора из временного ряда
    %
    % Вход:
    %   x       - временной ряд (вектор)
    %   L       - размерность вложения
    %   p       - лаг (delay)
    %   Npoints - число точек для отображения
    
    x = x(:);               % вектор-столбец
    N = length(x);
            % число векторов
    if nargin < 5
        M = N - (L-1)*p; % если не задано, берём максимум
    end
    if Npoints > M
        Npoints = M;        % ограничиваем, если больше доступного
    end
    
    % --- Формируем L-мерные векторы X_L(j) ---
    X = zeros(M, L);
    for m = 0:L-1
        X(:, m+1) = x((1:M) + m*p);
    end
    
    % --- Выбираем нужное количество точек ---
    Xplot = X(1:Npoints, :);
    
    % % --- XY проекция ---
    % figure;
    % plot(Xplot(:,1), Xplot(:,2), 'k.', 'MarkerSize',10);
    % xlabel('$X$', 'Interpreter','latex');
    % ylabel('$Y$', 'Interpreter','latex');
    % title('XY projection of attractor', 'Interpreter','latex');
    % 
    % % --- Добавляем текст с параметрами на рисунке ---
    % xlim_vals = xlim;
    % ylim_vals = ylim;
    % xpos = xlim_vals(2) - 0.2*diff(xlim_vals);
    % ypos = ylim_vals(2) - 0.08*diff(ylim_vals);
    % 
    % txt = sprintf('L(dimensions) = %d \n p(lag) = %d \n M(j-iterator) = %d', L, p, M);
    % text(xpos, ypos, txt, 'FontSize', 18, 'Color', 'k', 'Interpreter','latex');
    % 
    % 
    % set(gca,'FontSize',16,'LineWidth',2);
    % set(gcf,'Color','white');
    % grid on;
end
