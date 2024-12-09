p=[-0.0032480562835404563, 0.13387641022560062, -1.3976143856453171, 5.300165914337288, -5.501800261100386, -4.892603997626288, 4.046960362958583];
x = chebfun('x',[0.5,6],'minSamples',500);

%% kJ/mol to eV
V= 0.0103643*(p(1)*x.^6+p(2)*x.^5+p(3)*x.^4+p(4)*x.^3+p(5)*x.^2+p(6)*x+p(7));

% Find Vmax and its location
[Vmax, x_at_Vmax] = max(V);

% Define energy level (En) and calculate turning points
En = Vmax - 0.01; % Example energy level
turning_points = roots(V - En);

% Identify specific turning points (assuming they are real and sorted)
zeta_nL = turning_points(2);
zeta_m = x_at_Vmax;
zeta_nR = turning_points(end);

% Plot the potential energy barrier
figure;
plot(x, V, 'LineWidth', 2); hold on;

% Mark Vmax
plot(x_at_Vmax, Vmax, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(x_at_Vmax, Vmax, 'V_{max}', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

% Mark energy level En
yline(En, 'k--', 'LineWidth', 1.5, 'Label', 'E_n', 'LabelHorizontalAlignment', 'left');

% Mark turning points
plot(zeta_nL, En, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
plot(zeta_m, En, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
plot(zeta_nR, En, 'ko', 'MarkerSize', 8, 'LineWidth', 2);

% Annotate turning points
text(zeta_nL, En, '\zeta_{nL}', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(zeta_m, En, '\zeta_{m}', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(zeta_nR, En, '\zeta_{nR}', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Style plot
xlabel('\zeta (coordinate)');
ylabel('V(\zeta) (eV)');
title('Potential Energy Barrier with Turning Points');
grid off;

% Adjust y-axis scale
ylim([min(V) - 0.01, max(V) + 0.02]); % Increase the top y-axis scale by 0.04
hold off;