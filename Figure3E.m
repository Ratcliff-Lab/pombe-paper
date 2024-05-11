%change directory to correct folder

%import data
data = importdata("SS_assay_data.txt");

%calculate cluster # for each strain and dilution
t334_4000 = data(:,1) - data(:,2);
t60_4000 = data(:,2);

t334_8000 = data(:,3) - data(:,4);
t60_8000 = data(:,4);

t334_2000_ss = data(:,5) - data(:,6);
t60_2000_ss = data(:,6);

t334_4000_ss = data(:,7) - data(:,8);
t60_4000_ss = data(:,8);

%calculate CFU/mL for each strain and dilution
cfu_t334_4000 = (t334_4000*4000)/.1;
cfu_t60_4000 = (t60_4000*4000)/.1;

cfu_t334_8000 = (t334_8000*8000)/.1;
cfu_t60_8000 = (t60_8000*8000)/.1;

cfu_t334_2000_ss = (t334_2000_ss*2000)/.1;
cfu_t60_2000_ss = (t60_2000_ss*2000)/.1;

cfu_t334_4000_ss = (t334_4000_ss*4000)/.1;
cfu_t60_4000_ss = (t60_4000_ss*4000)/.1;

%calculate frequency of t334 before and after ss for each dilution
t334_freq_4000 = cfu_t334_4000./(cfu_t334_4000+cfu_t60_4000);

t334_freq_8000 = cfu_t334_8000./(cfu_t334_8000+cfu_t60_8000);

t334_freq_2000_ss = cfu_t334_2000_ss./(cfu_t334_2000_ss+cfu_t60_2000_ss);

t334_freq_4000_ss = cfu_t334_4000_ss./(cfu_t334_4000_ss+cfu_t60_4000_ss);

%take the mean of three replicates (3 per column) with std 10%

mean_t334_freq_4000 = [] ;
std_t334_freq_4000 = [] ;

for i=[1 4 7]
    mean_t334_freq_4000 = [mean_t334_freq_4000 mean(t334_freq_4000(i:i+2), 'omitnan')] ;
    std_t334_freq_4000 = [std_t334_freq_4000 std(t334_freq_4000(i:i+2))] ;
end

mean_t334_freq_8000 = [] ;
std_t334_freq_8000 = [] ;

for i=[1 4 7]
    mean_t334_freq_8000 = [mean_t334_freq_8000 mean(t334_freq_8000(i:i+2))] ;
    std_t334_freq_8000 = [std_t334_freq_8000 std(t334_freq_8000(i:i+2))] ;
end

mean_t334_freq_2000_ss = [] ;
std_t334_freq_2000_ss = [] ;

for i=[1 4 7]
    mean_t334_freq_2000_ss = [mean_t334_freq_2000_ss mean(t334_freq_2000_ss(i:i+2))] ;
    std_t334_freq_2000_ss = [std_t334_freq_2000_ss std(t334_freq_2000_ss(i:i+2))] ;
end

mean_t334_freq_4000_ss = [] ;
std_t334_freq_4000_ss = [] ;

for i=[1 4 7]
    mean_t334_freq_4000_ss = [mean_t334_freq_4000_ss mean(t334_freq_4000_ss(i:i+2))] ;
    std_t334_freq_4000_ss = [std_t334_freq_4000_ss std(t334_freq_4000_ss(i:i+2))] ;
end

% mean of dilutions before selection and after selection
merged = [mean_t334_freq_4000;mean_t334_freq_8000];
mean_t334_freq = mean(merged, 1);
merged_ss = [mean_t334_freq_2000_ss;mean_t334_freq_4000_ss];
mean_t334_freq_ss = mean(merged_ss, 1);

smushed = [std_t334_freq_4000;std_t334_freq_8000];
std_t334_freq = mean(smushed, 1, "omitnan");
smushed_ss = [std_t334_freq_2000_ss;std_t334_freq_4000_ss];
std_t334_freq_ss = mean(smushed_ss, 1);

%% make plot
hold on
plot(0:1, 0:1, '--', 'color', '#464646');
errorbar(mean_t334_freq, mean_t334_freq_ss, std_t334_freq, "horizontal",'.', 'color', "#464646", 'MarkerSize', 30);
errorbar(mean_t334_freq, mean_t334_freq_ss, std_t334_freq_ss,'.', 'color', "#464646");
plot(mean_t334_freq, mean_t334_freq_ss, 'b.', 'MarkerSize', 35, 'Color', "#756bb1");

xlim([0.4 0.9]);
ylim([0.4 1]);

%title('2-D Line Plot')
xlabel('t334 freq. before selection','FontSize',20);
ylabel('t334 freq. after selection','FontSize',20);
shg
set(gca,'Fontsize',20,'FontName','Times New Roman');

%% Calculate selective advantage
mean(mean_t334_freq_ss./mean_t334_freq) % 1.29
std(mean_t334_freq_ss./mean_t334_freq) % 0.15

%% Calculate relative selection rate r
% I choose to calculate the relative selection rate r in this case, because
% we are looking for relative fitness during selection, in which both
% competitors decrease in frequency. 
% r = (ln(A1/A0)-ln(B1/B0))/day
% This is when looking at growth (log and time)
% In our case, the calculation of r can be: A1/A0 - B1/B0

% Lower dilution before selection
% cfu_t334_4000 ad cfu_t60_4000 
% Higher dilution before selection
% cfu_t334_8000 and cfu_t60_8000 

% Lower dilution after selection
% cfu_t334_2000_ss and cfu_t60_2000_ss 
% Higher dilution after selection
%cfu_t334_4000_ss and cfu_t60_4000_ss 


r_dilution1 = cfu_t334_2000_ss./cfu_t334_4000 - cfu_t60_2000_ss./cfu_t60_4000 ;
r_dilution2 = cfu_t334_4000_ss./cfu_t334_8000 - cfu_t60_4000_ss./cfu_t60_8000 ;


% Take means before calculating the selection rates
r_dilution2(6) = NaN ; % outlier value
r = mean([r_dilution1,r_dilution2]','omitnan') ;

% three replicates, start at line 1, 4, 7
mean_r = [] ;
std_r = [] ;
stde_r = [] ;

% Calculate standard error
% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1255808/)

for i=[1 4 7]
    mean_r = [mean_r mean(r(i:i+2), 'omitnan')] ;
    std_r = [std_r std(r(i:i+2), 'omitnan')] ;
    stde_r = [stde_r, std(r(i:i+2), 'omitnan') / sqrt(length(r(i:i+2)))] ;
end

%% make plot
figure
hold on
plot([0,4], [0,0], '--', 'Color', "#756bb1");
errorbar(1:3, mean_r, stde_r,'.', 'color', "#464646");
plot(1:3, mean_r, '.','MarkerSize', 35, 'Color', "#756bb1");
xticks(1:3)
xticklabels(1:3)

xlim([0 4]);
ylim([-0.2 0.2]);

xlabel('Replicates','FontSize',20);
ylabel('Selection rate','FontSize',20);
shg
set(gca,'Fontsize',20,'FontName','Times New Roman');




%% Save figure

% set(gcf, 'color', 'none');    
% orient(gcf,'landscape')
% set(gcf,'Position',[100 100 300 300])
% exportgraphics(gcf,'selection_rate.pdf',...   
%     'ContentType','vector',...
%     'BackgroundColor','none')

%% Statistical test on selection rates
r_dataset = [ r(1:3)' ,  r(4:6)' ,  r(7:9)'] ;


aov = anova(r_dataset) ;
m = multcompare(aov) % "tukey-kramer"  is the default test here








