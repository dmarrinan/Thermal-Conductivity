clear all
close all
clc

jan24 = xlsread('Consolidated_Jan242013_thermaldata.xlsx');
jan29 = xlsread('Consolidated_Jan292013_thermaldata.xlsx');
jan31 = xlsread('Consolidated_Jan312013_thermaldata.xlsx');
feb05 = xlsread('Consolidated_Feb052013_thermaldata.xlsx');
x = 0:3:30; %position of TCs along bar in inches
x = x * 2.54; %convert inch to cm (2.54 cm / inch);
L = x(length(x)); %length of bar in cm
x0 = [x 31.5*2.54]; %add extra point for longer bar approximation
L0 = x0(length(x0));
xplot = 0:2.54/3:L; %higher resolution for x to plot theoretical temp
x0plot = 0:2.54/3:L0;
alpha = sqrt(61.16/60); %square root of thermal diffusivity for copper in cm^2/sec

%% Room temp calculations for each thermocouple new
muroomtempnew = zeros(13,1);
stdevroomtemp = zeros(13,1);
for i=1:13
    muroomtempnew(i) = mean(jan24(i,2:12));
    stdevroomtemp(i) = sqrt(1/10*sum((jan24(i,2:12)-muroomtempnew(i,1)).^2));
end

figure
plot(jan24(1:13,1),stdevroomtemp,'x','LineWidth',1)
ylim([0 1.1*max(stdevroomtemp)]);
xlabel('Time (s)','Fontsize',24)
ylabel('Standard Deviation (^{\circ}C)','Fontsize',24)
title('Standard Deviation for Room Temperature Thermocouple Measurements','Fontsize',24)
set(gca,'Fontsize',24)
saveas(gcf,'variation_room_temp.jpg')

sigmaT = max(stdevroomtemp)*ones(1,11);
%sigmaT = 2*ones(1,11);

%% 45 degrees C 
%Data Taken 1/29/13
t0_29 = 309; %last point before boundary conditions imposed
[mintemp29, tf_29] = min(jan29(:,2));

%Temp vs. Time for all TCs
figure
hold all
for i=2:12
    if i<9
        plot(jan29(:,1),jan29(:,i),'o')
    else
        plot(jan29(:,1),jan29(:,i),'+')
    end
end
hold off
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')
xlabel('Time (s)')
ylabel('Temperature (^{\circ}C)')
title('Temperature vs. Time for all Thermocouples 45^{\circ}C')

%Temp vs. Position for t=0
figure
errorbar(x,jan29(t0_29,2:12),sigmaT,'x')
title('Temperature Distribution at t=0 for T=45^{\circ}C')
xlabel('Position along Bar (cm)')
ylabel('Temperature (^{\circ}C)')
linfit29 = polyfit(x,jan29(t0_29,2:12),1);
linval29 = polyval(linfit29,x);
K29 = linfit29(1);
hold all
plot(x,linval29,'r')
hold off
legendentry = sprintf('Fit: K = %f',K29);
legend('Data',legendentry,'Location','Northwest')

%Temp Difference vs. Time for all Positions
U29 = zeros(length(jan29),11);
for i = 2:12
    U29(:,i-1) = jan29(:,i) - jan29(:,2);
end
figure
plot(U29)

%Fourier Series to fit theoretical values
t29 = jan29(t0_29:length(jan29),1)-jan29(t0_29,1);
Utheo29 = zeros(length(t29),length(xplot));
Ttheo29 = zeros(size(Utheo29));
for n = 1:20
    for i = 1:length(t29)
        Utheo29(i,:) = Utheo29(i,:) - 2*K29*L/pi*(-1)^n/n*sin(n*pi*xplot/L)*exp(-n^2*pi^2*alpha^2/L^2*t29(i));
    end
end
for i = 1:length(xplot)
    Ttheo29(:,i) = Utheo29(:,i) + jan29(t0_29:length(jan29),2);
end

figure
hold all
for i = 1:9:length(xplot)
    plot(Ttheo29(:,i))
end
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')

figure
hold all
plot(x,linval29,'LineWidth',2)
for i = 4:9:46
    plot(xplot,Ttheo29(i,:),'LineWidth',2)
end
plot(xplot,Ttheo29(180,:),'LineWidth',2)
errorbar(x,jan29(t0_29,2:12),sigmaT,'x','LineWidth',2)
for i = 3:9:45
    errorbar(x,jan29(t0_29+i,2:12),sigmaT,'x','LineWidth',2)
end
errorbar(x,jan29(t0_29+179,2:12),sigmaT,'x','LineWidth',2)
hold off
title('a','Fontsize',24)
% xlabel('Position Along Bar (cm)','Fontsize',24)
% ylabel('Temperature (^{\circ}C)','Fontsize',24)
ylim([0 1.10*45])
xlim([0 85])
set(gca,'YTick',0:5:1.05*45,'Fontsize',24,'XTickLabel','')
hleg = legend('t = 0 min','t = 1 min','t = 4 min','t = 7 min','t = 10 min','t = 13 min','t = 60 min','Location','Northwest');
set(hleg,'Fontsize',20)
saveas(gcf,'uncorrected_45.jpg')
chi_exp_29 = [Ttheo29(1,1:9:91); Ttheo29(4:9:46,1:9:91); Ttheo29(180,1:9:91)];
chi_obs_29 = [jan29(t0_29,2:12); jan29((t0_29+3):9:(t0_29+45),2:12); jan29(t0_29+179,2:12)];
chi_sq_29 = zeros(1,7);
for i = 1:7
    chi_sq_29(i) = sum((chi_exp_29(i,:)-chi_obs_29(i,:)).^2./chi_exp_29(i,:));
end

%Fourier Series to fit theoretical values with extra bar approximation
linval29_0 = [linval29 jan29(tf_29,2)];

t29 = jan29(t0_29:length(jan29),1)-jan29(t0_29,1);
Utheo29 = zeros(length(t29),length(x0plot));
Ttheo29 = zeros(size(Utheo29));
for n = 1:20
    for i = 1:length(t29)
        Utheo29(i,:) = Utheo29(i,:) - 2*K29*L0/pi*(-1)^n/n*sin(n*pi*x0plot/L0)*exp(-n^2*pi^2*alpha^2/L0^2*t29(i));
    end
end
for i = 1:length(x0plot)
    Ttheo29(:,i) = Utheo29(:,i) + jan29(t0_29:length(jan29),2);
end

figure
hold all
for i = 1:9:length(x0plot)
    plot(Ttheo29(:,i))
end
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')

figure
hold all
plot(x0,linval29_0,'LineWidth',2)
for i = 4:9:46
    plot(x0plot,Ttheo29(i,:),'LineWidth',2)
end
plot(x0plot,Ttheo29(180,:),'LineWidth',2)
errorbar(x0,[jan29(t0_29,2:12) linval29_0(12)],[sigmaT 0],'x','LineWidth',2)
for i = 3:9:45
    errorbar(x0,[jan29(t0_29+i,2:12) linval29_0(12)],[sigmaT 0],'x','LineWidth',2)
end
errorbar(x0,[jan29(t0_29+179,2:12) linval29_0(12)],[sigmaT 0],'x','LineWidth',2)
hold off
title('b','Fontsize',24)
%xlabel('Position Along Bar (cm)','Fontsize',24)
%ylabel('Temperature (^{\circ}C)','Fontsize',24)
ylim([0 1.10*45])
set(gca,'YTick',0:5:1.05*45,'Fontsize',24,'XTickLabel','')
xlim([0 85])
hleg = legend('t = 0 min','t = 1 min','t = 4 min','t = 7 min','t = 10 min','t = 13 min','t = 60 min','Location','Northwest');
set(hleg,'Fontsize',20)
saveas(gcf,'corrected_45.jpg')

chi_exp_29 = [Ttheo29(1,1:9:91); Ttheo29(4:9:46,1:9:91); Ttheo29(180,1:9:91)];
chi_obs_29 = [jan29(t0_29,2:12); jan29((t0_29+3):9:(t0_29+45),2:12); jan29(t0_29+179,2:12)];
chi_sq_29_cor = zeros(1,7);
for i = 1:7
    chi_sq_29_cor(i) = sum((chi_exp_29(i,:)-chi_obs_29(i,:)).^2./chi_exp_29(i,:));
end
figure
hold all
plot(chi_sq_29,'x')
plot(chi_sq_29_cor,'ro')
legend('{\chi}^2','Corrected {\chi}^2');
saveas(gcf,'chi_sq_29.jpg')

%% 50 degrees C 
%Data from 1/24/13
t0_24 = 186;
[mintemp24, tf_24] = min(jan24(:,2));

%Temp vs. Time for all TCs
figure
hold all
for i=2:12
    if i<9
        plot(jan24(:,1),jan24(:,i),'o')
    else
        plot(jan24(:,1),jan24(:,i),'+')
    end
end
hold off
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')
xlabel('Time (s)')
ylabel('Temperature (^{\circ}C)')
title('Temperature vs. Time for all Thermocouples 50^{\circ}C')

%Temp vs. Position for t=0
figure
errorbar(x,jan24(t0_24,2:12),sigmaT,'x')
title('Temperature Distribution at t=0 for T=50^{\circ}C')
xlabel('Position Along Bar (cm)')
ylabel('Temperature (^{\circ}C)')
linfit24 = polyfit(x,jan24(t0_24,2:12),1);
linval24 = polyval(linfit24,x);
K24 = linfit24(1);
hold all
plot(x,linval24,'r')
hold off
legendentry = sprintf('Fit: K = %f',linfit24(1));
legend('Data',legendentry,'Location','Northwest')

%Temp Difference vs. Time for all TCs
U24 = zeros(length(jan24),11);
for i = 2:12
    U24(:,i-1) = jan24(:,i) - jan24(:,2);
end
figure
plot(U24)

%Fourier Series to fit theoretical values
t24 = jan24(t0_24:length(jan24),1)-jan24(t0_24,1);
Utheo24 = zeros(length(t24),length(xplot));
Ttheo24 = zeros(size(Utheo24));
for n = 1:20
    for i = 1:length(t24)
        Utheo24(i,:) = Utheo24(i,:) - 2*K24*L/pi*(-1)^n/n*sin(n*pi*xplot/L)*exp(-n^2*pi^2*alpha^2/L^2*t24(i));
    end
end
for i = 1:length(xplot)
    Ttheo24(:,i) = Utheo24(:,i) + jan24(t0_24:length(jan24),2);
end
figure
hold all
for i = 1:9:length(xplot)
    plot(Ttheo24(:,i))
end
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')
figure
hold all
plot(x,linval24,'LineWidth',2)
for i = 4:9:46
    plot(xplot,Ttheo24(i,:),'LineWidth',2)
end
plot(xplot,Ttheo24(180,:),'LineWidth',2)
errorbar(x,jan24(t0_24,2:12),sigmaT,'x','LineWidth',2)
for i = 3:9:45
    errorbar(x,jan24(t0_24+i,2:12),sigmaT,'x','LineWidth',2)
end
errorbar(x,jan24(t0_24+179,2:12),sigmaT,'x','LineWidth',2)
hold off
%title('Temperature vs. Position for Max T = 50^{\circ}C','Fontsize',24)
%xlabel('Position Along Bar (cm)','Fontsize',24)
%ylabel('Temperature (^{\circ}C)','Fontsize',24)
ylim([0 1.10*50])
set(gca,'YTick',0:5:1.05*50,'Fontsize',24,'XTickLabel','')
xlim([0 85])
hleg = legend('t = 0 min','t = 1 min','t = 4 min','t = 7 min','t = 10 min','t = 13 min','t = 60 min','Location','Northwest');
set(hleg,'Fontsize',20)
saveas(gcf,'uncorrected_50.jpg')
chi_exp_24 = [Ttheo24(1,1:9:91); Ttheo24(4:9:46,1:9:91); Ttheo24(180,1:9:91)];
chi_obs_24 = [jan24(t0_24,2:12); jan24((t0_24+3):9:(t0_24+45),2:12); jan24(t0_24+179,2:12)];
chi_sq_24 = zeros(1,7);
for i = 1:7
    chi_sq_24(i) = sum((chi_exp_24(i,:)-chi_obs_24(i,:)).^2./chi_exp_24(i,:));
end

%Fourier Series to fit theoretical values with extra bar approximation
linval24_0 = [linval24 jan24(tf_24,2)];

t24 = jan24(t0_24:length(jan24),1)-jan24(t0_24,1);
Utheo24 = zeros(length(t24),length(x0plot));
Ttheo24 = zeros(size(Utheo24));
for n = 1:20
    for i = 1:length(t24)
        Utheo24(i,:) = Utheo24(i,:) - 2*K24*L0/pi*(-1)^n/n*sin(n*pi*x0plot/L0)*exp(-n^2*pi^2*alpha^2/L0^2*t24(i));
    end
end
for i = 1:length(x0plot)
    Ttheo24(:,i) = Utheo24(:,i) + jan24(t0_24:length(jan24),2);
end

figure
hold all
for i = 1:9:length(x0plot)
    plot(Ttheo24(:,i))
end
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')

figure
hold all
plot(x0,linval24_0,'LineWidth',2)
for i = 4:9:46
    plot(x0plot,Ttheo24(i,:),'LineWidth',2)
end
plot(x0plot,Ttheo24(180,:),'LineWidth',2)
errorbar(x0,[jan24(t0_24,2:12) linval24_0(12)],[sigmaT 0],'x','LineWidth',2)
for i = 3:9:45
    errorbar(x0,[jan24(t0_24+i,2:12) linval24_0(12)],[sigmaT 0],'x','LineWidth',2)
end
errorbar(x0,[jan24(t0_24+179,2:12) linval24_0(12)],[sigmaT 0],'x','LineWidth',2)
hold off
%title('Corrected Temperature vs. Position for Max T = 50^{\circ}C','Fontsize',24)
%xlabel('Position Along Bar (cm)','Fontsize',24)
%ylabel('Temperature (^{\circ}C)','Fontsize',24)
ylim([0 1.10*50])
set(gca,'YTick',0:5:1.05*50,'Fontsize',24,'XTickLabel','')
xlim([0 85])
hleg = legend('t = 0 min','t = 1 min','t = 4 min','t = 7 min','t = 10 min','t = 13 min','t = 60 min','Location','Northwest');
set(hleg,'Fontsize',20)
saveas(gcf,'corrected_50.jpg')

chi_exp_24 = [Ttheo24(1,1:9:91); Ttheo24(4:9:46,1:9:91); Ttheo24(180,1:9:91)];
chi_obs_24 = [jan24(t0_24,2:12); jan24((t0_24+3):9:(t0_24+45),2:12); jan24(t0_24+179,2:12)];
chi_sq_24_cor = zeros(1,7);
for i = 1:7
    chi_sq_24_cor(i) = sum((chi_exp_24(i,:)-chi_obs_24(i,:)).^2./chi_exp_24(i,:));
end

figure
hold all
plot(chi_sq_24,'x')
plot(chi_sq_24_cor,'ro')
legend('{\chi}^2','Corrected {\chi}^2');
saveas(gcf,'chi_sq_24.jpg')
%% 55 degrees C 
%Data from 1/31/13
t0_31 = 247;
[mintemp31, tf_31] = min(jan31(:,2));

%Temp vs. Time for all TCs
figure
hold all
for i=2:12
    if i<9
        plot(jan31(:,1),jan31(:,i),'o')
    else
        plot(jan31(:,1),jan31(:,i),'+')
    end
end
hold off
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')
xlabel('Time (s)','Fontsize',24)
ylabel('Temperature (^{\circ}C','Fontsize',24)
title('Temperature vs. Time for all Thermocouples 55^{\circ}C','Fontsize',24)
set(gca,'Fontsize',24)
saveas(gcf,'tc_data_55.jpg');

%Temp vs. Position for t=0
figure
errorbar(x,jan31(t0_31,2:12),sigmaT,'x')
title('Temperature Distribution at t=0 for T=55^{\circ}C')
xlabel('Thermocouple')
ylabel('Temperature (^{\circ}C)')
linfit31 = polyfit(x,jan31(t0_31,2:12),1);
linval31 = polyval(linfit31,x);
K31 = linfit31(1);
hold all
plot(x,linval31,'r')
hold off
legendentry = sprintf('Fit: K = %f',linfit31(1));
legend('Data',legendentry,'Location','Northwest')

%Temp Difference vs. Time for all TCs
U31 = zeros(length(jan31),11);
for i = 2:12
    U31(:,i-1) = jan31(:,i) - jan31(:,2);
end
figure
plot(U31)

%Fourier Series to fit theoretical values
t31 = jan31(t0_31:length(jan31),1)-jan31(t0_31,1);
Utheo31 = zeros(length(t31),length(xplot));
Ttheo31 = zeros(size(Utheo31));
for n = 1:20
    for i = 1:length(t31)
        Utheo31(i,:) = Utheo31(i,:) - 2*K31*L/pi*(-1)^n/n*sin(n*pi*xplot/L)*exp(-n^2*pi^2*alpha^2/L^2*t31(i));
    end
end
for i = 1:length(xplot)
    Ttheo31(:,i) = Utheo31(:,i) + jan31(t0_31:length(jan31),2);
end
figure
hold all
for i = 1:9:length(xplot)
    plot(Ttheo31(:,i))
end
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')
figure
hold all
plot(x,linval31,'LineWidth',2)
for i = 4:9:46
    plot(xplot,Ttheo31(i,:),'LineWidth',2)
end
plot(xplot,Ttheo31(180,:),'LineWidth',2)
errorbar(x,jan31(t0_31,2:12),sigmaT,'x','LineWidth',2)
for i = 3:9:45
    errorbar(x,jan31(t0_31+i,2:12),sigmaT,'x','LineWidth',2)
end
errorbar(x,jan31(t0_31+179,2:12),sigmaT,'x','LineWidth',2)
hold off
%title('Temperature vs. Position for Max T = 55^{\circ}C','Fontsize',24)
%xlabel('Position Along Bar (cm)','Fontsize',24)
%ylabel('Temperature (^{\circ}C)','Fontsize',24)
ylim([0 1.10*55])
set(gca,'YTick',0:5:1.05*55,'Fontsize',24,'XTickLabel','')
xlim([0 85])
hleg = legend('t = 0 min','t = 1 min','t = 4 min','t = 7 min','t = 10 min','t = 13 min','t = 60 min','Location','Northwest');
set(hleg,'Fontsize',20)
saveas(gcf,'uncorrected_55.jpg')
chi_exp_31 = [Ttheo31(1,1:9:91); Ttheo31(4:9:46,1:9:91); Ttheo31(180,1:9:91)];
chi_obs_31 = [jan31(t0_31,2:12); jan31((t0_31+3):9:(t0_31+45),2:12); jan31(t0_31+179,2:12)];
chi_sq_31 = zeros(1,7);
for i = 1:7
    chi_sq_31(i) = sum((chi_exp_31(i,:)-chi_obs_31(i,:)).^2./chi_exp_31(i,:));
end

%Fourier Series to fit theoretical values with extra bar approximation
linval31_0 = [linval31 jan31(tf_31,2)];

t31 = jan31(t0_31:length(jan31),1)-jan31(t0_31,1);
Utheo31 = zeros(length(t31),length(x0plot));
Ttheo31 = zeros(size(Utheo31));
for n = 1:20
    for i = 1:length(t31)
        Utheo31(i,:) = Utheo31(i,:) - 2*K31*L0/pi*(-1)^n/n*sin(n*pi*x0plot/L0)*exp(-n^2*pi^2*alpha^2/L0^2*t31(i));
    end
end
for i = 1:length(x0plot)
    Ttheo31(:,i) = Utheo31(:,i) + jan31(t0_31:length(jan31),2);
end

figure
hold all
for i = 1:9:length(x0plot)
    plot(Ttheo31(:,i))
end
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')

figure
hold all
plot(x0,linval31_0,'LineWidth',2)
for i = 4:9:46
    plot(x0plot,Ttheo31(i,:),'LineWidth',2)
end
plot(x0plot,Ttheo31(180,:),'LineWidth',2)
errorbar(x0,[jan31(t0_31,2:12) linval31_0(12)],[sigmaT 0],'x','LineWidth',2)
for i = 3:9:45
    errorbar(x0,[jan31(t0_31+i,2:12) linval31_0(12)],[sigmaT 0],'x','LineWidth',2)
end
errorbar(x0,[jan31(t0_31+179,2:12) linval31_0(12)],[sigmaT 0],'x','LineWidth',2)
hold off
%title('Corrected Temperature vs. Position for Max T = 55^{\circ}C','Fontsize',24)
%xlabel('Position Along Bar (cm)','Fontsize',24)
%ylabel('Temperature (^{\circ}C)','Fontsize',24)
ylim([0 1.10*55])
set(gca,'YTick',0:5:1.05*55,'Fontsize',24,'XTickLabel','')
xlim([0 85])
hleg = legend('t = 0 min','t = 1 min','t = 4 min','t = 7 min','t = 10 min','t = 13 min','t = 60 min','Location','Northwest');
set(hleg,'Fontsize',20)
saveas(gcf,'corrected_55.jpg')
chi_exp_31 = [Ttheo31(1,1:9:91); Ttheo31(4:9:46,1:9:91); Ttheo31(180,1:9:91)];
chi_obs_31 = [jan31(t0_31,2:12); jan31((t0_31+3):9:(t0_31+45),2:12); jan31(t0_31+179,2:12)];
chi_sq_31_cor = zeros(1,7);
for i = 1:7
    chi_sq_31_cor(i) = sum((chi_exp_31(i,:)-chi_obs_31(i,:)).^2./chi_exp_31(i,:));
end

figure
hold all
plot(chi_sq_31,'x')
plot(chi_sq_31_cor,'ro')
legend('{\chi}^2','Corrected {\chi}^2');
saveas(gcf,'chi_sq_31.jpg')
%% 75 degrees C
%Data from 2/5/13
t0_05 = 225;
[mintemp05, tf_05] = min(feb05(:,2));

%Temp vs. Time for all TCs
figure
hold all
for i=2:12
    if i<9
        plot(feb05(:,1),feb05(:,i),'o')
    else
        plot(feb05(:,1),feb05(:,i),'+')
    end
end
hold off
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')
xlabel('Time (s)')
ylabel('Temperature (^{\circ}C)')
title('Temperature vs. Time for all Thermocouples 75^{\circ}C')

%Temp vs. Position for t=0
figure
errorbar(x,feb05(t0_05,2:12),sigmaT,'x')
title('Temperature Distribution at t=0 for T=75^{\circ}C')
xlabel('Thermocouple')
ylabel('Temperature (^{\circ}C)')
linfit05 = polyfit(x,feb05(t0_05,2:12),1);
linval05 = polyval(linfit05,x);
K05 = linfit05(1);
hold all
plot(x,linval05,'r')
hold off
legend('Data','Fit','Location','Northwest')
legendentry = sprintf('Fit: K = %f',linfit05(1));
legend('Data',legendentry,'Location','Northwest')

%Temp Difference vs. Time for all TCs
U05 = zeros(length(feb05),11);
for i = 2:12
    U05(:,i-1) = feb05(:,i) - feb05(:,2);
end
figure
plot(U05)

%Fourier Series to fit theoretical values
t05 = feb05(t0_05:length(feb05),1)-feb05(t0_05,1);
Utheo05 = zeros(length(t05),length(xplot));
Ttheo05 = zeros(size(Utheo05));
for n = 1:20
    for i = 1:length(t05)
        Utheo05(i,:) = Utheo05(i,:) - 2*K05*L/pi*(-1)^n/n*sin(n*pi*xplot/L)*exp(-n^2*pi^2*alpha^2/L^2*t05(i));
    end
end
for i = 1:length(xplot)
    Ttheo05(:,i) = Utheo05(:,i) + feb05(t0_05:length(feb05),2);
end
figure
hold all
for i = 1:9:length(xplot)
    plot(Ttheo05(:,i))
end
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')
figure
hold all
plot(x,linval05,'LineWidth',2)
for i = 4:9:46
    plot(xplot,Ttheo05(i,:),'LineWidth',2)
end
plot(xplot,Ttheo05(180,:),'LineWidth',2)
errorbar(x,feb05(t0_05,2:12),sigmaT,'x','LineWidth',2)
for i = 3:9:45
    errorbar(x,feb05(t0_05+i,2:12),sigmaT,'x','LineWidth',2)
end
errorbar(x,feb05(t0_05+179,2:12),sigmaT,'x','LineWidth',2)
hold off
%title('Temperature vs. Position for Max T = 75^{\circ}C','Fontsize',24)
xlabel('Position Along Bar (cm)','Fontsize',24)
%ylabel('Temperature (^{\circ}C)','Fontsize',24)
ylim([0 1.10*75])
set(gca,'YTick',0:5:1.05*75,'Fontsize',24)
xlim([0 85])
hleg = legend('t = 0 min','t = 1 min','t = 4 min','t = 7 min','t = 10 min','t = 13 min','t = 60 min','Location','Northwest');
set(hleg,'Fontsize',20)
saveas(gcf,'uncorrected_75.jpg')
chi_exp_05 = [Ttheo05(1,1:9:91); Ttheo05(4:9:46,1:9:91); Ttheo05(180,1:9:91)];
chi_obs_05 = [feb05(t0_05,2:12); feb05((t0_05+3):9:(t0_05+45),2:12); feb05(t0_05+179,2:12)];
chi_sq_05 = zeros(1,7);
for i = 1:7
    chi_sq_05(i) = sum((chi_exp_05(i,:)-chi_obs_05(i,:)).^2./chi_exp_05(i,:));
end

%Fourier Series to fit theoretical values with extra bar approximation
linval05_0 = [linval05 feb05(tf_05,2)];

t05 = feb05(t0_05:length(feb05),1)-feb05(t0_05,1);
Utheo05 = zeros(length(t05),length(x0plot));
Ttheo05 = zeros(size(Utheo05));
for n = 1:20
    for i = 1:length(t05)
        Utheo05(i,:) = Utheo05(i,:) - 2*K05*L0/pi*(-1)^n/n*sin(n*pi*x0plot/L0)*exp(-n^2*pi^2*alpha^2/L0^2*t05(i));
    end
end
for i = 1:length(x0plot)
    Ttheo05(:,i) = Utheo05(:,i) + feb05(t0_05:length(feb05),2);
end

figure
hold all
for i = 1:9:length(x0plot)
    plot(Ttheo05(:,i))
end
legend('TC 1','TC 2','TC 3','TC 4','TC 5','TC 6','TC 7','TC 8','TC 9','TC 10','TC 11')

figure
hold all
plot(x0,linval05_0,'LineWidth',2)
for i = 4:9:46
    plot(x0plot,Ttheo05(i,:),'LineWidth',2)
end
plot(x0plot,Ttheo05(180,:),'LineWidth',2)
errorbar(x0,[feb05(t0_05,2:12) linval05_0(12)],[sigmaT 0],'x','LineWidth',2)
for i = 3:9:45
    %plot(x,feb05(t0_05+i,2:12),'x')
    errorbar(x0,[feb05(t0_05+i,2:12) linval05_0(12)],[sigmaT 0],'x','LineWidth',2)
end
errorbar(x0,[feb05(t0_05+179,2:12) linval05_0(12)],[sigmaT 0],'x','LineWidth',2)
hold off
%title('Corrected Temperature vs. Position for Max T = 75^{\circ}C','Fontsize',24)
xlabel('Position Along Bar (cm)','Fontsize',24)
%ylabel('Temperature (^{\circ}C)','Fontsize',24)
ylim([0 1.10*75])
set(gca,'YTick',0:5:1.05*75,'Fontsize',24)
xlim([0 85])
hleg = legend('t = 0 min','t = 1 min','t = 4 min','t = 7 min','t = 10 min','t = 13 min','t = 60 min','Location','Northwest');
set(hleg,'Fontsize',20)
saveas(gcf,'corrected_75.jpg')
chi_exp_05 = [Ttheo05(1,1:9:91); Ttheo05(4:9:46,1:9:91); Ttheo05(180,1:9:91)];
chi_obs_05 = [feb05(t0_05,2:12); feb05((t0_05+3):9:(t0_05+45),2:12); feb05(t0_05+179,2:12)];
chi_sq_05_cor = zeros(1,7);
for i = 1:7
    chi_sq_05_cor(i) = sum((chi_exp_05(i,:)-chi_obs_05(i,:)).^2./chi_exp_05(i,:));
end

figure
hold all
plot(chi_sq_05,'x')
plot(chi_sq_05_cor,'ro')
legend('{\chi}^2','Corrected {\chi}^2');
saveas(gcf,'chi_sq_05.jpg')
%% Room temp calculations for each thermocouple old
muroomtemp=zeros(1,11);
munumerator = 0;
mudenominator = 0;
sigmathermocouple = 1; %error on thermocouple +- 1 degree celsius
sigmadaq = 1; %error on DAQ +- 1 degree celsius
sigma = sigmathermocouple+sigmadaq;

for i = 1:11
    for j = 1:13
        munumerator = munumerator + jan24(j,i+1)/sigma^2;
        mudenominator = mudenominator + 1/sigma^2;
    end
    muroomtemp(i) = munumerator/mudenominator;
    munumerator = 0;
    mudenominator = 0;
end

figure
plot(muroomtemp,'x')
xlabel('Thermocouple')
ylabel('Temperature (^{\circ}C)')
title('Average Room Temperature Measurement by Thermocouple')
xlim([0 12])


