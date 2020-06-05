close all, clear all

%Input lambda and number of positive results
Lambda = [20, 40, 80, 160];
Positive = [10, 15, 22, 24]; %Out off 24

%Prepare curves
x = linspace(0, 0.06, 400);
probs = likelihood(x);

% Plot figure
figure(1);
tiledlayout(1, 1);
plot(x, probs)
hold on
xline(1/35.6)

% Label figure
legend("Likelihood Function", "1/35.6")
xlabel('Probability an individual RNA molecule is amplified');
ylabel( 'Relative Likelihood');
set(gca,'ytick',[])
xticks([0 0.01, 0.02, 0.03, 0.04, 0.05, 0.06])
xticklabels({'0', '1/100','1/50','1/33.3','1/25','1/20','1/16.7'})d
grid off

%Define likelihood function
function prob = likelihood(p)
    prob = nchoosek(24, 10).*exp(-20.*p).^14.*(1-exp(-20.*p)).^10.*...
        nchoosek(24, 15).*exp(-40.*p).^9.*(1-exp(-40.*p)).^15.*...
        nchoosek(24, 22).*exp(-80.*p).^2.*(1-exp(-80.*p)).^22.*...
        nchoosek(24, 24).*exp(-160.*p).^0.*(1-exp(-160.*p)).^24;
end


