% Odontocete In Vitro sound production

% This script plots data distributions of porpoise phonic lip opening and closing relative 
% to hydrophone signals in individual porpoises

% Author: Coen P.H. Elemans
% Date: 22/10/10


clc
close all
clear all

% Plot style
FontName = 'Arial';
FontSize = 18;
FontWeight = 'Normal';
FontColor = 'k';
FontName = 'Arial';
LabelFontSizeMultiplier=1.5;
colmap = lines(7);

%%
cd('/Users/coen/Dropbox/Scientific Publications/In Prep MS/Odontocetes sound prod/Data/In vitro/Combined final data')
IDs=[1:4];
for j =IDs
    switch j
        case 4
            load P24152_R_2_Photron_FinalData.mat
        case 3
            load P24152_R_MotionPro_FinalData.mat  
        case 2
            load P19928_R_MotionPro_FinalData.mat  
        case 1
            load P22291_R_MotionPro_FinalData.mat  
    end
    
IPI = Data.time_Click_HI1(2:end)-Data.time_Click_HI1(1:end-1);
ClickRate=1./IPI;

Data.IPI=IPI; Data.ClickRate=ClickRate;
Data.Open =  (Data.Bursa_Time_Cl-Data.Bursa_Time_Op); % [s]
Data.OQ= Data.Open(2:end)./Data.IPI;

Ind= ~isnan(Data.time_Click_HI1); % Find indices with HI1 clicks

figure (1)% histogram
subplot 211
hold on
pd = fitdist(1000*(Data.Bursa_Time_Cl-Data.time_Click_HI1)','Normal');
if j==4
    tim=(1000*(Data.Bursa_Time_Cl-Data.time_Click_HI1));
    yyaxis right
    h=histogram(tim,'BinWidth',.05, 'Normalization', 'count') ;
    h.FaceColor = colmap(j,:);
    pd = fitdist(tim','Normal');
end
yyaxis left
x_values = -5:.01:2;
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2, 'Color', colmap(j,:))
xlim ([-3 1])
title (['Phonic Lips closing' ])

mu_bursa_cl(j)=pd.mu

subplot 212
hold on
if j==4
    yyaxis right
    h=histogram(1000*(Data.Bursa_Time_Op-Data.time_Click_HI1),'BinWidth',.05, 'Normalization', 'count');
    h.FaceColor = colmap(j,:);
end
yyaxis left
pd = fitdist(1000*(Data.Bursa_Time_Op-Data.time_Click_HI1)','Normal');
x_values = -5:.01:2;
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2, 'Color', colmap(j,:))
xlabel('Time around click emission (ms)', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
title (['Phonic Lips opening' ])
xlim ([-3 1])
ax = gca;ax.FontName = FontName; ax.FontSize; ax.LabelFontSizeMultiplier=LabelFontSizeMultiplier;

mu_bursa_op(j)=pd.mu

end
%%
figure (2)
plot([1 1 1 1],mu_bursa_cl, 'o')
hold on
plot([2 2 2 2],mu_bursa_op, 'o')

['click radiated ' num2str(mean(mu_bursa_cl(IDs))) '+/- ' num2str(std(mu_bursa_cl(IDs)))  'ms after PL closing']
['click radiated ' num2str(-mean(mu_bursa_op(IDs))) '+/- ' num2str(std(mu_bursa_op(IDs)))  'ms after PL opening']
errorbar([1 2],[mean(mu_bursa_cl(IDs)) mean(mu_bursa_op)],[std(mu_bursa_cl(IDs)) std(mu_bursa_op(IDs))])
xlim([.5 2.5])
ylim([-3 .5])
xticks([1 2])
xticklabels({'Closing','Opening'})

ylabel('Time around click emission (ms)', 'FontName', FontName, 'FontSize', FontSize,  'FontWeight', FontWeight, 'Color', FontColor);
ax = gca;ax.FontName = FontName; ax.FontSize; ax.LabelFontSizeMultiplier=LabelFontSizeMultiplier;

['PL closed ' num2str(mean(1000*(Data.SI1_diff_min_time(Ind)-Data.Bursa_Time_Cl(Ind)))) '+/- ' num2str(std(1000*(Data.SI1_diff_min_time(Ind)-Data.Bursa_Time_Cl(Ind))))  'ms after closing signature in SI1']
['PL opened ' num2str(-mean(1000*(Data.time_Op_SI1(Ind)-Data.Bursa_Time_Op(Ind)))) '+/- ' num2str(std(1000*(Data.time_Op_SI1(Ind)-Data.Bursa_Time_Op(Ind))))  'ms after opening signature in SI1']
