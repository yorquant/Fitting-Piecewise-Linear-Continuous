% RNA_Polymerase_II_transcription_WANG_Yukun_LIU_Shixin20240513.m
% Copyright 2025 The Rockefeller University
% Author: Joel E. Cohen cohen@rockefeller.edu 20250925
% JEC  West Shokan, NY 20231123 24 1206 07 08 09 10 12 18 22 25
% 20240131 0204 05 0319 23 24 0422 0512 13 0612 13 14 15

% User must set:
% printdiary = 0 (no print diary, save figures) or = 1 (yes print diary, save figures)
% readmat = 0 means read from Excel spreadsheet; do NOT read from .mat file
%   readmat = 1 means read from .mat file; do NOT read from Excel spreadsheet
% Nsheets = number of conditions (sheets of Excel spreadsheet) 
%   Nsheets = 14 for all presently available conditions
% globalmaxchangepoints = 10 limits maximum number of changepoints (NCP)
% mintimegap = 3 sets minimal time (s) between changepoints;
%   if two changepoints are initially estimated within 'mintimegap',
%   they are replaced by one changepoint at the midpoint.

% TERMINOLOGY: duration measures time (s); length measures Bp (base pairs)

tic % Start timer to measure total duration

format compact
clc
clear all
close all
d = char([datetime('now','Format','yyyyMMddHHmmss')])
% addpath('C:\Users\Joel E. Cohen\Dropbox (Dropbox @RU)\JOEL\Mss\LiuShixinWangYukun2023\Rbeast')
addpath('.\Rbeast')
% https://cran.r-project.org/web/packages/Rbeast/index.html
% Rbeast: Bayesian Change-Point Detection and Time Series Decomposition
% Please use the canonical form https://CRAN.R-project.org/package=Rbeast to link to this page.

%% Set up printing options: 
% printdiary = 0 means do NOT save diary and figures
% printdiary = 1 means DO save diary and figures

printdiary = 1
if printdiary % Only if printdiary == 1
    diaryname = ['traces-14-conditions',d,'.dry']
    diary(diaryname)
end

%% Set up data reading options: 
% readmat = 0 means read from Excel spreadsheet; do NOT read from .mat file
% readmat = 1 means read from .mat file; do NOT read from Excel spreadsheet

readmat = 0


if readmat

    % Recover data saved in .mat file
    D = dir('traces*.mat');
    load([D.folder,'\',D.name])

    % By default, number of sheets or conditions Nsheets =  8. 
    % Nsheets can be reset below.
    
else

    % Describe spreadsheet

    datafilename = "traces-14-conditions.xlsx"
    
    sheetcoordinates{1} = 'A:HX';       % All-factors
    sheetcoordinates{2} = 'A:BP';       % ΔTFIIS
    sheetcoordinates{3} = 'A:CL';       % ΔIWS1
    sheetcoordinates{4} = 'A:BD';       % ΔDISF
    sheetcoordinates{5} = 'A:CJ';       % DSIFp-
    sheetcoordinates{6} = 'A:AR';       % ΔP-TEFb
    sheetcoordinates{7} = 'A:AH';       % Pol II p-
    sheetcoordinates{8} = 'A:AT';       % ΔPAF1C
    sheetcoordinates{9} = 'A:AX';       % ΔSPT6
    sheetcoordinates{10} = 'A:AP';      % ΔRTF1    
    sheetcoordinates{11} = 'A:AL';      % ΔRTF1&DSIF
    sheetcoordinates{12} = 'A:AH';      % ΔRTF1&PAF1C
    sheetcoordinates{13} = 'A:BF';      % Pol IIp+&DSIFp+
    sheetcoordinates{14} = 'A:CN';      % All-factors-16-kb

    sheetcoordinates

    Nsheets = length(sheetcoordinates)

    sheetcoordinatestxt{1} = 'A1:HX1';
    sheetcoordinatestxt{2} = 'A1:BP1';
    sheetcoordinatestxt{3} = 'A1:CL1';
    sheetcoordinatestxt{4} = 'A1:BD1';
    sheetcoordinatestxt{5} = 'A1:CJ1';
    sheetcoordinatestxt{6} = 'A1:AR1';
    sheetcoordinatestxt{7} = 'A1:AH1';
    sheetcoordinatestxt{8} = 'A1:AT1';
    sheetcoordinatestxt{9} = 'A1:AX1';
    sheetcoordinatestxt{10} = 'A1:AP1';   
    sheetcoordinatestxt{11} = 'A1:AL1';   
    sheetcoordinatestxt{12} = 'A1:AH1';   
    sheetcoordinatestxt{13} = 'A1:BF1';      
    sheetcoordinatestxt{14} = 'A1:CN1';    
    sheetcoordinatestxt

    Ntrajectories = zeros(1,Nsheets);
    sheetnames = {"All-factors", "ΔTFIIS", "ΔIWS1", "ΔDISF", ...
        "DSIFp-", "ΔP-TEFb", "Pol II p-", "ΔPAF1C", "ΔSPT6", "ΔRTF1","ΔRTF1&DSIF","ΔRTF1&PAF1C","Pol IIp+&DSIFp+","All-factors-16-kb"}  
    stringsheetnames=string(sheetnames)

    % Read spreadsheet data 
 
    for sheet = 1:Nsheets
        outnum{sheet} = readmatrix(datafilename,"FileType","spreadsheet", "Sheet" ,sheet, "Range",sheetcoordinates{sheet} );
        outtxt{sheet} = readcell(datafilename,"FileType","spreadsheet", "Sheet" ,sheet, "Range",sheetcoordinatestxt{sheet} );
        Ntrajectories(sheet) = length(outtxt{sheet})/2;
        trajectorynames{sheet} = outtxt{sheet}(1:2:end);
        for trace = 1:Ntrajectories(sheet)
            time{sheet} = outnum{sheet}(:,(1:2:(2*Ntrajectories(sheet))));
            bp{sheet} = outnum{sheet}(:,(2:2:(2*Ntrajectories(sheet))));
        end 
    end

    % Recycle OLD .mat file(s); save current values as NEW .mat file
    recycle("on")
    delete('*.mat')
    save(['traces-14-conditions_',d,'.mat'],'bp','Nsheets',...
        'Ntrajectories','outnum','outtxt','time','trajectorynames','sheetnames','stringsheetnames')
    recycle("off")
end
%% Initial values of arrays saved in .mat file and required for execution

bp
Nsheets
sheetnames
Ntrajectories
outnum
outtxt
sheetnames
stringsheetnames
time
trajectorynames

totalNtrajectories = sum(Ntrajectories)

%% Find piecewise linear changepoints using ischange
% Then find continuous piecewise linear approximation using
% fitPiecewiseLinearFunctionGolovchecko.m with these changepoints

globalmaxchangepoints = 15
mintimegap = 3  % Minimal time between change points
figno = 0

% Nsheets = 1   %Temporary value for debugging

for sheet =  1:Nsheets
        
    trajectorylengths{sheet} = zeros(1,length(trajectorynames{sheet}));
    % q Temporary container for Bayesian estimates of probabilities
    q = zeros(globalmaxchangepoints+1,Ntrajectories(sheet));
    % MaxCP Maximal number of changepoints required by data 
    MaxCP{sheet} = zeros(1,Ntrajectories(sheet)); 
    % SegmentStats has 1 row for each segment, 1 column for each attribute:
    % sheet,trace,SegmentNumber,StartTime,,LengthBp,slopeBpPerSecBetwCPs,slopeClass.
    % A fresh version of SegmentStats starts for each sheet and cumulates.
    SegmentStats = [];  

    % For each sheet, plot each trace, up to panels per figure
    for trace =  1:Ntrajectories(sheet)
        
        if 1 == mod(trace, 16)
            figno = figno + 1
            figure('WindowState','maximized')
            tiledlayout(4,4)
        end
    
        % Clean up data: suppress NaN, perturb slightly to eliminate duplicate times

        disp(['sheet ',num2str(sheet),', trajectory ',num2str(trace)])
        trajectorylengths{sheet}(1,trace) = sum(isfinite(bp{sheet}(:,trace)));
        bpreal{sheet,trace} = bp{sheet}(isfinite(bp{sheet}(:,trace)),trace);
        timereal{sheet,trace} = sort(time{sheet}(isfinite(bp{sheet}(:,trace)),trace)+0.000001*(rand(size(bpreal{sheet,trace}))-1/2));

       diary off  
       % Omit Bayesian estimation of changepoints from diary, if diary is ON.
        % https://github.com/zhaokg/Rbeast  Identify change points with cumulative probability >= 95%
        % Version edited by JEC beast_irregJEC20240614 currently has NO changes from original.
        o =  beast_irreg(bpreal{sheet,trace}, 'time',timereal{sheet,trace},...
            'tcp.minmax',[0,globalmaxchangepoints], ...
            'deltat',1,'season','none','print.options',0);
        q(:,trace) = o.trend.ncpPr;
        MaxCP{sheet}(1,trace) = sum(cumsum(q(:,trace))<0.99);
        % Equivalent to more concise version: MaxCP{sheet} = sum(cumsum(q)<0.95)
        if printdiary, diary on, end % Restart diary if printdiary == 1.

        % ischange  

        [tf,slopes,intercepts] = ischange(bpreal{sheet,trace},'linear',...
            'MaxNumChanges',MaxCP{sheet}(1,trace),'SamplePoints',timereal{sheet,trace});
        % tf, slopes, intercepts are all column vectors of the same size 
        %   as timereal{sheet,trace}.
        %   tf is vector of 0 with 1 at each position where slope changes.
        % slopes is slope of each linear segment; 
        %   slope changes to new value at position where tf == 1.
        % intercepts shows the intercept of each linear segment; 
        %   intercepts changes to new value at position where tf == 1.
        
        % Plot data and piecewise linear approximation
        
        nexttile
        plot(timereal{sheet,trace},bpreal{sheet,trace},'x')
        hold on                         
        % To plot piecewise linear but NOT continuous approximation from
        % "ischange", UNcomment the following 2 or 5 lines:
        % plot(timereal{sheet,trace},slopes.*timereal{sheet,trace} + intercepts,'.r') % linear regime
        % hold on
        % xlabel('time (s)')
        % ylabel('length (Bp)')
        % legend('data','linear','linear continuous','change','Location','southoutside')
        
    
        changepointtimes{sheet,trace}  =  timereal{sheet,trace}(tf);
        diffchangepointtimes = diff(changepointtimes{sheet,trace});
        NchangepointsBefore{sheet,trace} = length(changepointtimes{sheet,trace}(:));
                
        %Combine changepointtimes separated by less than mintimegap
        
        while min(diffchangepointtimes) < mintimegap
        
            for j = 1:length(diffchangepointtimes)
                if diffchangepointtimes(j) < mintimegap
                    changepointtimes{sheet,trace}(j) = ...
                        0.5*(changepointtimes{sheet,trace}(j) ...
                        + changepointtimes{sheet,trace}(j+1));
                    changepointtimes{sheet,trace}(j+1) = -1;
                end
            end
            changepointtimes{sheet,trace} = ...
                changepointtimes{sheet,trace}(changepointtimes{sheet,trace} >= 0 );
            diffchangepointtimes = diff(changepointtimes{sheet,trace});
        end

        NchangepointsAfter{sheet,trace}  = length(changepointtimes{sheet,trace}(:));
        
        
        % Find continuous piecewise linear approximation using fitPiecewiseLinearFunctionGolovchenko.m
    
        x0 = [min(timereal{sheet,trace});changepointtimes{sheet,trace};max(timereal{sheet,trace})];
        p = fitPiecewiseLinearFunctionGolovchenko(timereal{sheet,trace},bpreal{sheet,trace}, x0);
        %
        % function p = fitPiecewiseLinearFunction(x, y, x0)
        % Fit a piecewise continuous function f(x) to the pairs of data points (x,y)
        % such that the sum of squares of error is minimal.
        %% x0 - values of x that define ends of segments of function f(x)
        % p - end points of the segments p = f(x0)
        % See also: http://golovchenko.org/docs/ContinuousPiecewiseLinearFit.pdf
        % 4-May-2004 Nikolai Golovchenko.
        %
        durationSec{sheet,trace} = x0(end);
        lengthBp{sheet,trace} = p(end);
        line(x0',p,'LineWidth',2,'LineStyle','-','Color','k')
        if ~isempty(changepointtimes{sheet,trace})
            hold on
            plot([1;1]*changepointtimes{sheet,trace}',ylim,'g')            % changes in linear regime
        end
        hold off
        
        title(trajectorynames{sheet}(trace),'FontWeight','normal')
        axis tight
        sgtitle(cell2mat(['Pol2transcription ',d,' Sheet ',num2str(sheet),' ',sheetnames{sheet}]),'FontSize',10,'FontWeight','normal')
        set(gca,'FontSize',10)

    % Tabulate for EVERY SEGMENT: sheet, trace, segment number, start time, ...
    % segment duration seconds, segment length basepairs, slope, slope class 
    % slopeClass = 2 iff BpPerSecBetwCPs>10
    %   = 1 iff (slopeBpPerSecBetwCPs<=10) & (LengthBp>300)
    %   = 0 otherwise

    % SecBpAtCPs = [x0,p];
    SegmentNumber = (1:(1+NchangepointsAfter{sheet,trace}))';
    StartTime = x0(1:(end-1));
    DurationSec = diff(x0);
    LengthBp = diff(p);
    slopeBpPerSecBetwCPs = LengthBp./DurationSec;
    % DEFINE slopeClass:
    % slopeClass = 2 iff slopeBpPerSecBetwCPs>10
    %   = 1 iff (1<=slopeBpPerSecBetwCPs<=10) & (LengthBp>300)
    %   = 0 otherwise
    slopeClass = zeros(size(slopeBpPerSecBetwCPs)) ...
        + 2*(slopeBpPerSecBetwCPs>10) ...
        + ((slopeBpPerSecBetwCPs<=10) & (LengthBp>100) & (slopeBpPerSecBetwCPs>=1));
    SegmentStats = [SegmentStats;[sheet*ones(max(SegmentNumber),1),...
        trace*ones(max(SegmentNumber),1),SegmentNumber,StartTime,...
        DurationSec,LengthBp,slopeBpPerSecBetwCPs,slopeClass]];
        end    

    sgtitle(cell2mat(['Pol2transcription ',d,' Sheet ',num2str(sheet),' ',sheetnames{sheet}]),'FontSize',10,'FontWeight','normal')
    if printdiary
        print(['Pol2transcription',d,'_Sheet',num2str(sheet),' ',char(sheetnames{sheet})],'-dpng')
        savefig('compact')
    end
    % Output table SegmentStatistics separately for each sheet. Traj = trace number
    SegmentStatistics{sheet} = array2table(SegmentStats,...
    'VariableNames',{'Sheet','Traj','Segment','StartTime','DurationSec',...
        'LengthBp','BpPerSecBetwCPs','slopeClass'})
    SegmentStatistics{sheet}

    % Plot scattergram of segment statistics 'StartTime','DurationSec',
        % 'LengthBp','BpPerSecBetwCPs' versus each other 
        % plus histograms of each variable.
    figno = figno + 1
    figure('WindowState','maximized')
    [S,AX,BigAx,H,HAx] = plotmatrix(SegmentStats(:,4:7));
    title(BigAx,['Segments in ',char(sheetnames{sheet})])
% https://www.mathworks.com/matlabcentral/answers/183203-add-label-to-sub-axes-in-plotmatrix
    xlabel(AX(4,1),'Start Time (s)','FontSize',10)
    xlabel(AX(4,2),'Duration (s)','FontSize',10)
    xlabel(AX(4,3),'Length (Bp)','FontSize',10)
    xlabel(AX(4,4),'Slope Bp/s','FontSize',10)
    ylabel(AX(1,1),'Start Time (s)','FontSize',10)
    ylabel(AX(2,1),'Duration (s)','FontSize',10)
    ylabel(AX(3,1),'Length (Bp)','FontSize',10)
    ylabel(AX(4,1),'Slope Bp/s','FontSize',10)
    set(BigAx,'FontSize',10)
end

% End of level one loop 2==================================================
%% Save SegmentStatistics{sheet} to Excel spreadsheet ['SegmentStatistics',d,'.xlsx']
for sheet =  1:Nsheets
    writetable(SegmentStatistics{sheet},['SegmentStatistics',d,'.xlsx'],'Sheet',sheetnames{sheet}) 
end

%% Tabulate changepoints before and after combining changepoints 
% separated in time by less than mintimegap 

jj = 0;
for sheet = 1:Nsheets
    for trace = 1:Ntrajectories(sheet)
        jj = jj + 1;
        Sheet(jj,1) = sheet;
        Traj(jj,1) = trace;
        MaxChangepoints(jj,1) = MaxCP{sheet}(trace);
        % NCP = Number of Change Points
        NCPBefore(jj,1) = NchangepointsBefore{sheet,trace};
        NCPAfter(jj,1) = NchangepointsAfter{sheet,trace};
        deltaNCP(jj,1) = NCPAfter(jj,1) - NCPBefore(jj,1);
        DurationSec(jj,1) = durationSec{sheet,trace};
        CPperSec(jj,1) = NchangepointsAfter{sheet,trace}/x0(end);
        LengthBp(jj,1) = lengthBp{sheet,trace};
        CPperBp(jj,1) = NchangepointsAfter{sheet,trace}/p(end);
        BpPerSec(jj,1) = lengthBp{sheet,trace}/durationSec{sheet,trace};
    end
end

mintimegap

TableCP = table(Sheet,Traj,NCPBefore,NCPAfter,deltaNCP,...
    DurationSec,CPperSec,LengthBp,CPperBp,BpPerSec)
% sumTableCP = sum(TableCP) for 4 variables: NCPAfter,deltaNCP,DurationSec,LengthBp
sumTableCP = sum(table(NCPAfter,deltaNCP,DurationSec,LengthBp))

%% slopeClass frequency (absolute and relative) by condition (sheet)
slopeClassFreq = zeros(Nsheets+1,7);
for sheet = 1:Nsheets
    for col = 1:3
        slopeClassFreq(sheet, col) = ...
            sum(SegmentStatistics{sheet}.slopeClass == (col-1));
    end
    
    slopeClassFreq(sheet, 6) = slopeClassFreq(sheet, 2)/...
        (slopeClassFreq(sheet, 2)+slopeClassFreq(sheet, 3));
    slopeClassFreq(sheet, 7) = slopeClassFreq(sheet, 3)/...
        (slopeClassFreq(sheet, 2)+slopeClassFreq(sheet, 3));
end
slopeClassFreq(Nsheets+1,1:3) = sum(slopeClassFreq(:,1:3));
slopeClassFreq(:,4) = sum(slopeClassFreq(:,1:3)');
slopeClassFreq(:,5) = slopeClassFreq(:, 1)./slopeClassFreq(:, 4);
slopeClassFreq(Nsheets+1, 6) = slopeClassFreq(Nsheets+1, 2)/...
        (slopeClassFreq(Nsheets+1, 2)+slopeClassFreq(Nsheets+1, 3));
    slopeClassFreq(Nsheets+1, 7) = slopeClassFreq(Nsheets+1, 3)/...
        (slopeClassFreq(Nsheets+1, 2)+slopeClassFreq(Nsheets+1, 3));

TslopeClassFreq = array2table(slopeClassFreq,...
    'VariableNames',{'Pause','Slow','Fast','Sum','Pause/Sum',...
    'Slow/(S+F)','Fast/(S+F)'},...
    'RowNames',[stringsheetnames,"Total"])

%% slopeClass frequencies for all sheets (conditions)
slopeClassSummary = zeros(Nsheets,8);
IndependenceModel = zeros(Nsheets,3);
for sheet = 1:Nsheets
    sumSlopeClass = sum(histcounts(SegmentStatistics{sheet}.slopeClass));
    slopeClassSummary(sheet,:) = [sheet, ...
        histcounts(SegmentStatistics{sheet}.slopeClass),...
        sumSlopeClass,...
        histcounts(SegmentStatistics{sheet}.slopeClass)/sumSlopeClass];
end

TslopeClass = array2table(slopeClassSummary,...
    'VariableNames',{'Sheet','Freq.0','Freq.1','Freq.2','Sum','Prop.0','Prop.1','Prop.2'})

colsums = sum(slopeClassSummary(:,2:4));
rowsums = sum(slopeClassSummary(:,2:4)')';
IndependenceModel = rowsums*colsums/sum(colsums);
TableIndependenceModel = array2table(IndependenceModel,...
    'VariableNames',{'Freq.0','Freq.1','Freq.2'},...
    'RowNames',stringsheetnames)
chisq = slopeClassSummary(:,2:4)-IndependenceModel;
chisq = sum(chisq(:).^2./IndependenceModel(:))
df = (Nsheets-1)*2
pALLslopeClass = chi2cdf(chisq,df,'upper')
%% Pairwise comparisons of frequencies of slopeClass
pPairwiseslopeClass = zeros(Nsheets);
for sheet1 = 1:(Nsheets-1)
    for sheet2 = (sheet1+1):Nsheets
        sheet1, sheet2
        twosheetsslopeClass = [slopeClassSummary(sheet1,2:4);slopeClassSummary(sheet2,2:4)]
        colsums = sum(twosheetsslopeClass);
        rowsums = sum(twosheetsslopeClass')';
        IndependenceModel = rowsums*colsums/sum(colsums)
        chisq = twosheetsslopeClass-IndependenceModel;
        chisq = sum(chisq(:).^2./IndependenceModel(:));
        df = 2;
        pPairwiseslopeClass(sheet1,sheet2) = chi2cdf(chisq,df,'upper');
    end
end
TablepPairwiseslopeClass = array2table(pPairwiseslopeClass,'VariableNames',stringsheetnames,'RowNames',stringsheetnames)

%% Mean and standard deviation for each Sheet (condition)
% of trajectory NCP, DurationSec, LengthBp, BpPerSec

[GSheets,ID] = findgroups(Sheet);
Ntraj = Ntrajectories(ID)';
meanNCP = splitapply(@mean, NCPAfter,GSheets);
stdNCP = splitapply(@std, NCPAfter,GSheets);
mnDurationSec = splitapply(@mean, DurationSec,GSheets);
stdDurationSec = splitapply(@std, DurationSec,GSheets);
mnLengthBp = splitapply(@mean, LengthBp,GSheets);
stdLengthBp = splitapply(@std, LengthBp,GSheets);
mnBpPerSec = splitapply(@mean, BpPerSec,GSheets);
stdBpPerSec = splitapply(@std, BpPerSec,GSheets);

TMeanStdSheet = table(Ntraj,meanNCP,stdNCP,mnDurationSec,...
    stdDurationSec,mnLengthBp,stdLengthBp,mnBpPerSec,stdBpPerSec,...
    'RowNames',stringsheetnames)
if printdiary
    save(string(['Pol2transcription',d,'TableCP']), "TableCP", '-mat')
end

% End of level one loop 3==================================================

%% For each combination of sheet (condition) and slope class, find
% Mean and standard deviation of segments'
% DurationSec, LengthBp, BpPerSec
% For each variable, do Welch anova1 by sheet (condition), that is,
% for each slope class, do the means differ by sheet?

TSegmentStatistics = [];
for sheet = 1:Nsheets
    TSegmentStatistics = [TSegmentStatistics;SegmentStatistics{sheet}];
end

condition = reshape(repmat(stringsheetnames,3,1),(3*Nsheets),1);

[G,ID] = findgroups(TSegmentStatistics.Sheet,TSegmentStatistics.slopeClass);
SlopeClass = repmat([0;1;2],max(ID),1);
figno = figno + 1
figure('WindowState','maximized')
disp('WARNING: In log-log plots, means < 1 are replaced by 1, so log mean = 0.')

disp(' ')
disp('DurationSec All slope classes')
[n,m,v] = welchanovaJEC20240613([TSegmentStatistics.DurationSec, G], 0.01);
% function [n,m,v] = welchanovaJEC20240613(x,alpha)
%WELCHANOVA Welch ANOVA Test for Unequal Variances.
%The ANOVA F-test to compare the means of k normally distributed
%  populations is not applicable when the variances are unknown, and not
%  known to be equal. A spacial case, k=2, is the famous Behrens-Fisher 
%  problem (Behrens, 1929; Fisher, 1935). Welch (1951) test was proposed to
%  fill this void, a generalization to his 1947 previous paper (Welch, 
%  1947). 
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2012). welchanova: Welch 
%     ANOVA Test for Unequal Variances. [WWW document]. URL http://
%     www.mathworks.com/matlabcentral/fileexchange/37121-welchanova
% OUTPUT FORMAT MODIFIED 2024-06-13 BY JOEL E. COHEN; NO CHANGE TO
% CALCULATIONS.
TSegmentDurationSec = table(condition,SlopeClass,n,m,v,...
    'VariableNames',{'Condition','SlopeClass','NSegments','MnDurationSec','VrDurationSec'})


for sC= 0:1:2
    disp(' ')
    disp(['DurationSec slopeClass ',num2str(sC)])
    [GsC,IDsC] = findgroups(TSegmentStatistics(TSegmentStatistics.slopeClass==sC,1));
    [n,m,v] = welchanovaJEC20240613([TSegmentStatistics(TSegmentStatistics.slopeClass==sC,:).DurationSec, GsC], 0.01);
end

nexttile
loglog(max(1,m),v,'x','MarkerSize',15)
title('DurationSec')
xlabel('mean')
ylabel('variance')
set(gca,'FontSize',10)

disp(' ')
disp('LengthBp All slope classes')
[n,m,v] = welchanovaJEC20240613([TSegmentStatistics.LengthBp, G], 0.01);
TSegmentLengthBp = table(condition,SlopeClass,n,m,v,...
    VariableNames={'Condition','SlopeClass','NSegments','MnLengthBp','VrLengthBp'})

for sC= 0:1:2
    disp(['slopeClass ',num2str(sC)])
    [GsC,IDsC] = findgroups(TSegmentStatistics(TSegmentStatistics.slopeClass==sC,1));
    [n,m,v] = welchanovaJEC20240613([TSegmentStatistics(TSegmentStatistics.slopeClass==sC,:).LengthBp, GsC], 0.01);
end


nexttile
loglog(max(1,m),v,'x','MarkerSize',15)
title('LengthBp')
xlabel('mean')
ylabel('variance')
set(gca,'FontSize',10)

disp(' ')
disp('BpPerSecBetwCPs All slope classes')
[n,m,v] = welchanovaJEC20240613([TSegmentStatistics.BpPerSecBetwCPs, G], 0.01);
TBpPerSecBetwCPs = table(condition,SlopeClass,n,m,v,...
    VariableNames={'Condition','SlopeClass','NSegments','MnBpPerSec','VrBpPerSec'})

for sC= 0:1:2
    disp(['slopeClass ',num2str(sC)])
    [GsC,IDsC] = findgroups(TSegmentStatistics(TSegmentStatistics.slopeClass==sC,1));
    [n,m,v] = welchanovaJEC20240613([TSegmentStatistics(TSegmentStatistics.slopeClass==sC,:).BpPerSecBetwCPs, GsC], 0.01);
end

nexttile
loglog(max(1,m),v,'x','MarkerSize',15)
title('BpPerSec')
xlabel('mean')
ylabel('variance')
set(gca,'FontSize',10)


%% Histograms of trajectory features by sheet

VariableNames = {"NchangesAfter","DurationSec","CPperSec","LengthBp","CPperBp"}
disp('Histograms of trajectory features by sheet')
% More variables can be added here



% FOR Nsheets = 8, NsheetsRows = 2, NsheetsCols = 4; FOR Nsheets = 9, 3 3;
% FOR Nsheets = 14, 4 4
NsheetsRows = 4
NsheetsCols = 4

% NchangesAfter
figno = figno + 1
figure('WindowState','maximized')

for sheet=1:Nsheets
   
    subplot(NsheetsRows,NsheetsCols,sheet)
    histogram((NCPAfter(Sheet == sheet)), 'BinMethod','integers')
    title(sheetnames{sheet},'FontWeight','normal')
    if 1 == mod(sheet,3), ylabel('frequency'),end
    set(gca,'FontSize',10)
end
sgtitle(VariableNames{1},'FontSize',10,'FontWeight','bold')

% DurationSec
figno = figno + 1
figure('WindowState','maximized')

for sheet=1:Nsheets
   
    subplot(NsheetsRows,NsheetsCols,sheet)
    histogram((DurationSec(Sheet == sheet))) %, 'BinMethod','integers')
    title(sheetnames{sheet},'FontWeight','normal')
    if 1 == mod(sheet,3), ylabel('frequency'),end
    set(gca,'FontSize',10)
end
sgtitle(VariableNames{2},'FontSize',10,'FontWeight','bold')

% CPperSec
figno = figno + 1
figure('WindowState','maximized')

for sheet=1:Nsheets
   
    subplot(NsheetsRows,NsheetsCols,sheet)
    histogram((CPperSec(Sheet == sheet))) %, 'BinMethod','integers')
    title(sheetnames{sheet},'FontWeight','normal')
    if 1 == mod(sheet,3), ylabel('frequency'),end
    set(gca,'FontSize',10)
end
sgtitle(VariableNames{3},'FontSize',10,'FontWeight','bold')

% LengthBp
figno = figno + 1
figure('WindowState','maximized')

for sheet=1:Nsheets
   
    subplot(NsheetsRows,NsheetsCols,sheet)
    histogram((LengthBp(Sheet == sheet))) %, 'BinMethod','integers')
    title(sheetnames{sheet},'FontWeight','normal')
    if 1 == mod(sheet,3), ylabel('frequency'),end
    set(gca,'FontSize',10)
end
sgtitle(VariableNames{4},'FontSize',10,'FontWeight','bold')

% CPperBp
figno = figno + 1
figure('WindowState','maximized')

for sheet=1:Nsheets
   
    subplot(NsheetsRows,NsheetsCols,sheet)
    histogram((CPperBp(Sheet == sheet))) %, 'BinMethod','integers')
    title(sheetnames{sheet},'FontWeight','normal')
    if 1 == mod(sheet,3), ylabel('frequency'),end
    set(gca,'FontSize',10)
end
sgtitle(VariableNames{5},'FontSize',10,'FontWeight','bold')

% End of histograms of variables plotted by sheet (condition)
%% 
%% 
%% 
%% 

% %% anova1  one-way analysis of variance, comparing conditions for each variable 
% G = findgroups(Sheet);
% 
% % Compare conditions by total number of CP per trajectory
% 
% figno = figno + 1
% figure('WindowState','maximized')
% [~,anovaCP,statsCP]  = anova1(NCPAfter, G,'off');
% anovaCP
% [cCP, mCP] = multcompare(statsCP);   %,"Display","off"
% title(VariableNames{1},'FontSize',20,'FontWeight','bold')
% set(gca,'FontSize',18)
% tblNCPAfter = [table(sheetnames(1:Nsheets)','VariableNames',{'Condition'}), ...
%     array2table(mCP,'VariableNames',{'Mean','Standard Error'})]
% tCompareNCPAfter = array2table(cCP,"VariableNames", ...
%     ["Sheet","Comparison","Lower Limit","Difference","Upper Limit","P-value"])
% 
% % Compare conditions by DurationSec
% 
% figno = figno + 1
% figure('WindowState','maximized')
% 
% [~,anovaDurationSec,statsDurationSec]  = anova1(DurationSec, G,'off');
% anovaDurationSec
% [cDurationSec, mDurationSec] = multcompare(statsDurationSec);   %,"Display","off"
% title(VariableNames{2},'FontSize',20,'FontWeight','bold')
% set(gca,'FontSize',18)
% tblDurationSec = [table(sheetnames(1:Nsheets)','VariableNames',{'Condition'}), ...
%     array2table(mDurationSec,'VariableNames',{'Mean','Standard Error'})]
% tCompareDurationSec = array2table(cDurationSec,"VariableNames", ...
%     ["Sheet","Comparison","Lower Limit","Difference","Upper Limit","P-value"])
% 
% % Compare conditions by number of CP per second
% 
% figno = figno + 1
% figure('WindowState','maximized')
% 
% [~,anovaCPperSec,statsCPperSec]  = anova1(CPperSec, G,'off');
% anovaCPperSec
% [cCPperSec, mCPperSec] = multcompare(statsCPperSec);   %,"Display","off"
% title(VariableNames{3},'FontSize',20,'FontWeight','bold')
% set(gca,'FontSize',18)
% tblCPperSec = [table(sheetnames(1:Nsheets)','VariableNames',{'Condition'}), ...
%     array2table(mCPperSec,'VariableNames',{'Mean','Standard Error'})]
% tCompareCPperSec = array2table(cCPperSec,"VariableNames", ...
%     ["Sheet","Comparison","Lower Limit","Difference","Upper Limit","P-value"])
% 
% 
% 
% % Compare conditions by LengthBp
% 
% figno = figno + 1
% figure('WindowState','maximized')
% 
% [~,anovaLengthBp,statsLengthBp]  = anova1(LengthBp, G,'off');
% anovaLengthBp
% [cLengthBp, mLengthBp] = multcompare(statsLengthBp);   %,"Display","off"
% title(VariableNames{4},'FontSize',20,'FontWeight','bold')
% set(gca,'FontSize',18)
% tblLengthBp = [table(sheetnames(1:Nsheets)','VariableNames',{'Condition'}), ...
%     array2table(mLengthBp,'VariableNames',{'Mean','Standard Error'})]
% tLengthBp = array2table(cLengthBp,"VariableNames", ...
%     ["Sheet","Comparison","Lower Limit","Difference","Upper Limit","P-value"])
% 
% % Compare conditions by number of CP per base pair
% 
% figno = figno + 1
% figure('WindowState','maximized')
% 
% [~,anovaCPperBp,statsCPperBp]  = anova1(CPperBp, G,'off');
% anovaCPperBp
% [cCPperBp, mCPperBp] = multcompare(statsCPperBp);   %,"Display","off"
% title(VariableNames{5},'FontSize',20,'FontWeight','bold')
% set(gca,'FontSize',18)
% tblCPperBp = [table(sheetnames(1:Nsheets)','VariableNames',{'Condition'}), ...
%     array2table(mCPperBp,'VariableNames',{'Mean','Standard Error'})]
% tCompareCPperBp = array2table(cCPperBp,"VariableNames", ...
%     ["Sheet","Comparison","Lower Limit","Difference","Upper Limit","P-value"])


%% Slope as a function of duration (s) for slow and fast segments
% and average slope of pooled slow and fast segments weighted by
% transcript length (Bp)
figno = figno + 1
figure('WindowState','maximized')

DurationSec1 = cell(Nsheets,1);
DurationSec2 = cell(Nsheets,1);
slope1 = cell(Nsheets,1);
slope2 = cell(Nsheets,1);
Slope12WtedAve = zeros(Nsheets,1);

for sheet = 1:Nsheets
    DurationSec1{sheet} = SegmentStatistics{sheet}.DurationSec(SegmentStatistics{sheet}.slopeClass == 1);
    slope1{sheet} = SegmentStatistics{sheet}.BpPerSecBetwCPs(SegmentStatistics{sheet}.slopeClass == 1);
    DurationSec2{sheet} = SegmentStatistics{sheet}.DurationSec(SegmentStatistics{sheet}.slopeClass == 2);
    slope2{sheet} = SegmentStatistics{sheet}.BpPerSecBetwCPs(SegmentStatistics{sheet}.slopeClass == 2);
    
    LengthBpNotPaused = SegmentStatistics{sheet}.LengthBp...
        ((SegmentStatistics{sheet}.slopeClass == 1) | (SegmentStatistics{sheet}.slopeClass == 2));
    slope12 = SegmentStatistics{sheet}.BpPerSecBetwCPs...
        ((SegmentStatistics{sheet}.slopeClass == 1) | (SegmentStatistics{sheet}.slopeClass == 2));
    Slope12WtedAve(sheet,1) = sum(LengthBpNotPaused.*slope12)/sum(LengthBpNotPaused);



    subplot(NsheetsRows,NsheetsCols,sheet)
    if length(slope1) > 0
        plot(DurationSec1{sheet},slope1{sheet},'.k','MarkerSize',9)
        % loglog(DurationSec1{sheet},slope1{sheet},'.k','MarkerSize',15)
        hold on
    end
    if length(slope1) > 0
        plot(DurationSec2{sheet},slope2{sheet},'.k','MarkerSize',9)
        % loglog(DurationSec2{sheet},slope2{sheet},'xk','MarkerSize',15)
        hold off
    end
    xlabel('Duration (s)')
    ylabel('Slope (base pairs/sec)')
    axis([0 300 0 100])
    set(gca,'FontSize',10)
    title(cell2mat(['Sheet ',num2str(sheet),' ',sheetnames{sheet}]),'FontSize',11,'FontWeight','normal')
end
sgtitle(['Pol2transcription ',d, ' fast and slow segments, not paused'],'FontSize',14,'FontWeight','normal')   

disp(' ')
disp('Mean slope (Bp/s) of fast or slow segments (not paused)')
disp('weighted by length (Bp) of fast or slow segments')
Slope12WtedAve

% Slope as a function of length (bp) for slow and fast segments
figno = figno + 1
figure('WindowState','maximized')

LengthBp1 = cell(Nsheets,1);
LengthBp2 = cell(Nsheets,1);
% slope1 = cell(Nsheets,1);
% slope2 = cell(Nsheets,1);
Slope12WtedAve = zeros(Nsheets,1);

for sheet = 1:Nsheets
    LengthBp1{sheet} = SegmentStatistics{sheet}.LengthBp(SegmentStatistics{sheet}.slopeClass == 1);
    slope1{sheet} = SegmentStatistics{sheet}.BpPerSecBetwCPs(SegmentStatistics{sheet}.slopeClass == 1);
    LengthBp2{sheet} = SegmentStatistics{sheet}.LengthBp(SegmentStatistics{sheet}.slopeClass == 2);
    slope2{sheet} = SegmentStatistics{sheet}.BpPerSecBetwCPs(SegmentStatistics{sheet}.slopeClass == 2);
    
    % LengthBpNotPaused = SegmentStatistics{sheet}.LengthBp...
    %     ((SegmentStatistics{sheet}.slopeClass == 1) | (SegmentStatistics{sheet}.slopeClass == 2));
    % slope12 = SegmentStatistics{sheet}.BpPerSecBetwCPs...
    %     ((SegmentStatistics{sheet}.slopeClass == 1) | (SegmentStatistics{sheet}.slopeClass == 2));
    % Slope12WtedAve(sheet,1) = sum(LengthBpNotPaused.*slope12)/sum(LengthBpNotPaused);



    subplot(NsheetsRows,NsheetsCols,sheet)
    if length(slope1) > 0
        plot(LengthBp1{sheet},slope1{sheet},'.k','MarkerSize',9)
        hold on
    end
    if length(slope1) > 0
        plot(LengthBp2{sheet},slope2{sheet},'.k','MarkerSize',9)
        hold off
    end
    legend('slow','fast','Location','northeast')
    xlabel('Length (Bp)')
    ylabel('Slope (base pairs/sec)')
    axis([0 4000 0 100])
    set(gca,'FontSize',12)
    title(cell2mat(['Sheet ',num2str(sheet),' ',sheetnames{sheet}]),'FontSize',11,'FontWeight','normal')
end
sgtitle(['Pol2transcription ',d, ' fast and slow segments, not paused'],'FontSize',14,'FontWeight','normal')   

% disp(' ')
% disp('Mean slope (Bp/s) of fast or slow segments (not paused)')
% disp('weighted by length (Bp) of fast or slow segments')
% Slope12WtedAve

%% Slope as a function of duration (s) for paused segments 
% (excluding last paused segment of each trajectory including the last
% trajectory)

figno = figno + 1
figure('WindowState','maximized')

NlastSegmentsNotPaused = zeros(Nsheets,1);
Sum0EXcludingLastSegments = zeros(Nsheets,1);
FractionTimeSlopeClass0 = zeros(Nsheets,1);
DurationSec0 = cell(Nsheets,1);
DurationSec0EXcludingLastSegments = cell(Nsheets,1);

for sheet = 1:Nsheets
    % DurationSec0 = Duration (s) of each paused segment
    DurationSec0{sheet} = SegmentStatistics{sheet}.DurationSec(SegmentStatistics{sheet}.slopeClass == 0);
    % NotLastSegment = 1 for NOT last segment of each trajectory; 
    %                = 0 for last segment of each trajectory
    NotLastSegment = (diff([SegmentStatistics{sheet}.Segment;0])>0);
    NotfirstSegment =zeros(length(SegmentStatistics{sheet}.Segment),1);
    for iii=1:length(SegmentStatistics{sheet}.Segment)
        if SegmentStatistics{sheet}.Segment(iii)==1
            NotfirstSegment(iii)=0;
        else
            NotfirstSegment(iii)=1;
        end
    end

    % DurationSec0EXcludingLastSegments = Duration (s) of every NOT last paused segment 
    DurationSec0EXcludingLastSegments{sheet} = ...
        SegmentStatistics{sheet}.DurationSec((SegmentStatistics{sheet}.slopeClass == 0) & NotLastSegment & NotfirstSegment);
    % For this sheet or condition, number of not paused last segments 
    NlastSegmentsNotPaused(sheet,1) = sum(SegmentStatistics{sheet}.slopeClass(((~NotLastSegment)>0)|(~NotfirstSegment)>0));
    SumDurationSec0EXcludingLastSegments(sheet,1) = sum(DurationSec0EXcludingLastSegments{sheet});
    FractionTimeSlopeClass0(sheet,1) = SumDurationSec0EXcludingLastSegments(sheet,1)/...
        sum(SegmentStatistics{sheet}.DurationSec(NotLastSegment));
    
    
    binedgesT=[0:50:500];
    subplot(NsheetsRows,NsheetsCols,sheet)
    histogram(DurationSec0EXcludingLastSegments{sheet},'BinEdges',binedgesT)
    % legend('slow','fast','Location','northeast')
    xlabel('Duration (s)')
    ylabel('Number of segments')
    set(gca,'FontSize',18)
    title(cell2mat(['Sheet ',num2str(sheet),' ',sheetnames{sheet}]),'FontSize',14,'FontWeight','normal')
end
sgtitle(['Pol2transcription ',d, ' paused segments, excluding last segment'],...
    'FontSize',10,'FontWeight','normal') 
Tpaused = table(NlastSegmentsNotPaused,SumDurationSec0EXcludingLastSegments,...
    FractionTimeSlopeClass0,'RowNames',stringsheetnames)
%% Stop diary and quit

if printdiary
    % Save all figures to .png files
    for fign = 1:figno
        print(['-f',num2str(fign)],['./FiguresToShare/',d,'_figure',num2str(fign)],'-dpng')
    end
    diary off
end
%% Report total duration

toc  

%% Exploration of data prior to fitting
% for sheet = 1:Nsheets   
%     disp(['sheet = ',num2str(sheet)])
%     time{sheet}(2,:), range(time{sheet}(2,:))       % First times on sheet & their range
%     time{sheet}(3,:), range(time{sheet}(3,:))       % Second times on sheet & their range
%     deltat{sheet} = time{sheet}(3,:) - time{sheet}(2,:);    % Initial time steps
%     deltat{sheet}
% end

% return
