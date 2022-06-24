

freq = fmin:df:fmax;
figure
load inertszilard
plot(freq,20*log10(abs(inert)),'LineWidth',2)

hold on
load inert_kirchoff_2020
plot(freq,20*log10(abs(inert)),'LineWidth',2)

hold on
load inert_kirchoff_2022
plot(freq,20*log10(abs(inert)),'LineWidth',2)

legend('Szilard','Kirccchoff 2020','Kirccchoff 2022')