p=[-0.0032480562835404563, 0.13387641022560062, -1.3976143856453171, 5.300165914337288, -5.501800261100386, -4.892603997626288, 4.046960362958583];
x = chebfun('x',[0.5,6],'minSamples',500);

%% kJ/mol to eV
V= 0.0103643*(p(1)*x.^6+p(2)*x.^5+p(3)*x.^4+p(4)*x.^3+p(5)*x.^2+p(6)*x+p(7));

%% Calculate eigenvectors and eigenvalues

[xmin, xmax] = domain(V);
L = chebop([xmin,xmax]);
L.lbc = 0; L.rbc = 0;
L.op = @(x,u) (-5.363975e-5)*diff(u,2) + V.*u;
L.bc = 0;

neigs = 120;
[EV, D] = eigs(L, neigs, 'sr');
eval = diag(D); % Extract the eigenvalues    


%%
% Select the eigenvector index
Evec_ind = 28;

% Plot the potential and the selected eigenvalue
tt = 0.5:6;
figure; grid on; plot(V, 'LineWidth', 2); hold on;
plot(tt, eval(Evec_ind)*ones(size(tt)), '--', 'LineWidth', 2);
title('Potential V and Selected Eigenvalue');
xlabel('x'); ylabel('V(x) and E');

%%
eval(Evec_ind)
ind=find((V==eval(Evec_ind)));
xl=ind(2);
xr=ind(3);

figure,
Evec= plot(EV(:,Evec_ind));
xl_ind=find(abs(Evec.XData-xl)<0.001);
xr_ind=find(abs(Evec.XData-xr)<0.001);

% Normalize the eigenvector
Evec.YData = normr(Evec.YData);

% Calculate the probability integral
PnInteg=trapz(Evec.XData(xl_ind:xr_ind),Evec.YData(xl_ind:xr_ind).^2);
DwellTime=(44.904e-15/sqrt(eval(Evec_ind))) * PnInteg;

% Calculate tunneling time (entropic time)
interval=0.0001;
tmp=1./sqrt(V(xl+interval:interval:xr)-eval(Evec_ind));
TauInteg=44.904e-15*trapz(xl+interval:interval:xr,tmp);

EntropicTime=-TauInteg/(PnInteg*log(PnInteg));

% Convert time results to picoseconds
DwellTime_ps = DwellTime * 1e12;
EntropicTime_ps = EntropicTime * 1e12;

% Display results in picoseconds
disp(['Dwell Time: ', num2str(DwellTime_ps), ' ps']);
disp(['Entropic Time: ', num2str(EntropicTime_ps), ' ps']);
%%
% Plot the potential
tt = 7:0.01:11.5;
figure; grid on; plot(V, 'LineWidth', 2, 'DisplayName', 'Potential'); hold on;
title('Potential V and Selected Eigenvectors');
xlabel('\zeta'); ylabel('Energy (eV)');

% List of desired eigenvector indices
desired_indices = [28, 31, 32, 34, 37];

% Normalization factor based on the maximum potential energy barrier
V_max = max(V);

% Plot each desired eigenvector separately
for i = 1:length(desired_indices)
    Evec_ind = desired_indices(i);
    
    % Retrieve the eigenvector and corresponding eigenvalue
    Evec = EV(:, Evec_ind);
    eigenvalue = eval(Evec_ind);
    
    % Normalize the eigenvector for better visualization
    Evec_ind = normr(Evec_ind);

    % Scale the eigenvector by the maximum potential energy barrier
    Evec = Evec * V_max / 8;  % Adjust the scale factor as needed
    
    % Offset the eigenvector by its eigenvalue for better visibility
    plot(Evec + eigenvalue, 'LineWidth', 1.5, 'DisplayName', ['E_' num2str(i) ' = ' num2str(eigenvalue, '%.4f')]);

    % Annotate the plot with the eigenvalue
    mid_index = round(length(tt) / 2);  % Find the midpoint index for annotation
    text(mid_index, eigenvalue, num2str(eigenvalue, '%.2f'), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
end

% Show the legend with a custom position
lgd = legend('show', 'Location', 'southeast');
%%
% Find the maximum value of the potential energy function
V_max = max(V);

% Display the maximum value
disp(['Maximum value of the potential energy function V(x) is: ', num2str(V_max, '%.4f')]);

% Calculate the barrier height for the current eigenenergy
barrier_height = V_max - eigenvalue;
disp(['Barrier height is: ', num2str(barrier_height, '%.4f')]);
