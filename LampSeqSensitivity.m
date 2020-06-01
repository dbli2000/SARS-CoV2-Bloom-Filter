close all, clear all

%Input x- and y- values
Input = [0, 20, 40, 80, 160];
Detected = [0, 10, 15, 22, 24];

%Calculate some representative curves
x2 = linspace(0, 160);
Rep1 = (24*(1-0.95.^x2)); %1/20
Rep2 = (24*(1-0.9875.^x2)); %1/80

% Prepare for curve fitting
[x, y] = prepareCurveData( Input, Detected );

% Set up fittype and options
ft = fittype( '24*(1-(1-a)^x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = 0.03;

% Fit model to data.
[fitresult, gof] = fit(x, y, ft, opts);
disp(fitresult)

% Plot figure
figure('Name', 'Sensitivity Fit' );
tiledlayout(1, 1);
h = plot(fitresult, x, y);
hold on
plot(x2, Rep1)
plot(x2, Rep2)

% Label figure
legend('Data', 'Best Fit (\approx1/38)',  '1/20', '1/80', 'Location', 'SouthEast');
xlabel('Input number of RNA molecules');
ylabel( 'Reactions positive');
grid off
