%% Vicon Nexus Treadmill Gait Data Analysis %%
% Originating Author: Andrew Hagen
% Last Revised: 3/20/2023

% This script uses data from Vicon Nexus during a treadmill walking trial to calculate a variety of gait parameters including Phase Coordination Index (PCI)...
% other temporal gait measures, step length asymmetry (SLA), ground reaction forces, joint angles, joint moments, and joint powers

% This script was used in a study in people with multiple sclerosis, and labeled each limb as affected or unaffected

% At the end of the script, the combined data can be exported as an xlsx file, and figures of PCI, SLA, and limb excursion can be created

%% Script initalization

% Must have a Nexus trail open on this machine
    % This trail must have gait events (heel strike and toe off) labeled and the dynamic gait model processed
    % Using the MATLAB integration in the Nexus pipeline makes this easy
   
% Table variable abreviations
    % UnA = Unaffected Limb
    % A = Affected limb
    % FH = First Half of Trial
    % SH = Second Half of Trial

clc 
clear 
close all

%% Set up connection with Vicon Nexus and set path
    %Make sure Nexus is open with the desired trial 
vicon = ViconNexus;

SubjectName = vicon.GetSubjectNames;
subnum = char(SubjectName);
Rate = vicon.GetFrameRate;
Time = 1/Rate;
Trial = inputdlg('Trial Name: ');
Alimb = inputdlg('Affected Limb? (L or R): '); 
TrialName = char(Trial);

 %PATH = pwd;
% RootPath = uigetdir(PATH); % Self select path
 RootPath = ['R:\SBNT Study\Pilot Testing\' subnum '\'];
cd(RootPath)

%% Calculate PCI Equation

    %     phi_CV is the coefficient of variation of the mean of phi values
    %         for each subject.   
    %         phi_CV = sigma_phi/mu_phi (STD of phi/ Mean of phi)
    %     Pphi_ABS is the conversion of phi_ABS to a percentage.
    %         Pphi_ABS = 100*phi_ABS/180
    %     Finally,  
    %         PCI = phi_CV + Pphi_ABS


L_HeelStrikes = vicon.GetEvents(subnum, 'Left', 'Foot Strike');
R_HeelStrikes = vicon.GetEvents(subnum, 'Right', 'Foot Strike');
L_ToeOffs = vicon.GetEvents(subnum, 'Left', 'Foot Off');
R_ToeOffs = vicon.GetEvents(subnum, 'Right', 'Foot Off');


%% Gait Cycle Duration [s] : - Heel strike to ipsilateral Heel strike 

L_StrideT = (diff([L_HeelStrikes(:,end) L_HeelStrikes],[],2));
R_StrideT = (diff([R_HeelStrikes(:,end) R_HeelStrikes],[],2));

L_StrideT(:,1) = []; % Crop first columns (always 0)
R_StrideT(:,1) = []; 

% Confirm stride count is the same so we can work with matrix
    if length(L_StrideT) ~= length(R_StrideT)
        if length(L_StrideT) > length(R_StrideT)
            L_StrideT(length(L_StrideT)) = [];
        else R_StrideT(length(R_StrideT)) = [];  
        end
    end

% Convert from frames to seconds 
L_StrideT = cast(L_StrideT, "double");
L_StrideT_sec = L_StrideT/Rate;
R_StrideT = cast(R_StrideT, "double");
R_StrideT_sec = R_StrideT/Rate;

% Crop 5 gait cycles
L_StrideT_Crop = L_StrideT_sec;
L_StrideT_Crop(:,(1:5)) = [];
R_StrideT_Crop = R_StrideT_sec;
R_StrideT_Crop(:,(1:5)) = [];

% First and second half of trial
StrideT_Half = length(L_StrideT_Crop)/2; %Halfway point of trial (will be same for L and R)
L_StrideT_FH = L_StrideT_Crop(:,1:fix(StrideT_Half)); % First half (FH) of trial
L_StrideT_SH = L_StrideT_Crop(:,(fix(StrideT_Half))+1:end); % Second half (SH) of trial
R_StrideT_FH = R_StrideT_Crop(:,1:fix(StrideT_Half));
R_StrideT_SH = R_StrideT_Crop(:,(fix(StrideT_Half))+1:end);

% Means
if strcmp(Alimb,'L')
    A_GCD_Avg = mean(L_StrideT_Crop);
    A_GCD_Avg_FH = mean(L_StrideT_FH);
    A_GCD_Avg_SH = mean(L_StrideT_SH);
    UnA_GCD_Avg = mean(R_StrideT_Crop);
    UnA_GCD_Avg_FH = mean(R_StrideT_FH);
    UnA_GCD_Avg_SH = mean(R_StrideT_SH);
elseif strcmp(Alimb,'R')
    UnA_GCD_Avg = mean(L_StrideT_Crop);
    UnA_GCD_Avg_FH = mean(L_StrideT_FH);
    UnA_GCD_Avg_SH = mean(L_StrideT_SH);
    A_GCD_Avg = mean(R_StrideT_Crop);
    A_GCD_Avg_FH = mean(R_StrideT_FH);
    A_GCD_Avg_SH = mean(R_StrideT_SH);
end
%% Step Duration [s] - Heel strike to contralateral heel strike

%START W/ LEFT STEP
if R_HeelStrikes(1) > L_HeelStrikes(1) % Starting w/ a left step (i.e. left heel strike)
    L_StepDur= zeros(1,length(L_HeelStrikes)-1);
    R_StepDur= zeros(1,length(L_HeelStrikes)-1);
            for i = 1:(length(L_HeelStrikes)-1) % The # of heelstrikes will always be equal or +1 for L_Heelstrikes depending on last step
                L_StepDur(i) = R_HeelStrikes(i) - L_HeelStrikes(i);
                R_StepDur(i) = L_HeelStrikes(i+1) - R_HeelStrikes(i);
            end

%START w/ RIGHT STEP
elseif L_HeelStrikes(1) > R_HeelStrikes(1) % Starting w/ a right step (i.e. right heel strike)
        L_StepDur= zeros(1,length(R_HeelStrikes)-1);
        R_StepDur= zeros(1,length(R_HeelStrikes)-1);
            for i = 1:(length(R_HeelStrikes)-1) % The # of heelstrikes will always be equal or +1 for R_Heelstrikes depending on last step
                R_StepDur(i) = L_HeelStrikes(i) - R_HeelStrikes(i);
                L_StepDur(i) = R_HeelStrikes(i+1) - L_HeelStrikes(i);
            end
end

% Confirm step count is the same so we can work with matrix
    if length(L_StepDur) ~= length(R_StepDur)
        if length(L_StepDur) > length(R_StepDur)
            L_StepDur(length(L_StepDur)) = [];
        else R_StepDur(length(R_StepDur)) = [];    
        end
    end

% Convert from frames to seconds 
L_StepDur = cast(L_StepDur, "double");
L_StepDur = L_StepDur/Rate;
R_StepDur = cast(R_StepDur, "double");
R_StepDur = R_StepDur/Rate;

% Crop 5 gait cycles
L_StepDur_Crop = L_StepDur;
L_StepDur_Crop(:,(1:5)) = [];
R_StepDur_Crop = R_StepDur;
R_StepDur_Crop(:,(1:5)) = [];

% First and second half of trial
StepDur_Half = length(L_StepDur_Crop)/2; % Halfway point of trial (will be same for L and R)
L_StepDur_FH = L_StepDur_Crop(:,1:fix(StepDur_Half)); % First half (FH) of trial
L_StepDur_SH = L_StepDur_Crop(:,(fix(StepDur_Half))+1:end); % Second half (SH) of trial
R_StepDur_FH = R_StepDur_Crop(:,1:fix(StepDur_Half));
R_StepDur_SH = R_StepDur_Crop(:,(fix(StepDur_Half))+1:end);

% Means
if strcmp(Alimb,'L')
    A_StepDur_Avg = mean(L_StepDur_Crop);
    A_StepDur_Avg_FH = mean(L_StepDur_FH);
    A_StepDur_Avg_SH = mean(L_StepDur_SH);
    UnA_StepDur_Avg = mean(R_StepDur_Crop);
    UnA_StepDur_Avg_FH = mean(R_StepDur_FH);
    UnA_StepDur_Avg_SH = mean(R_StepDur_SH);
elseif strcmp(Alimb,'R')
    UnA_StepDur_Avg = mean(L_StepDur_Crop);
    UnA_StepDur_Avg_FH = mean(L_StepDur_FH);
    UnA_StepDur_Avg_SH = mean(L_StepDur_SH);
    A_StepDur_Avg = mean(R_StepDur_Crop);
    A_StepDur_Avg_FH = mean(R_StepDur_FH);
    A_StepDur_Avg_SH = mean(R_StepDur_SH);
end

%% Swing Time [s] - Toe off to ipsilateral heel strike

if length(L_HeelStrikes) == length(L_ToeOffs) % Always will start with a toe off so heel strikes <= toe offs. Therefore do n or n-1 of L_HeelStrikes iterations.
    L_SwingT= zeros(1,length(L_HeelStrikes)-1);
    for i = 1:(length(L_HeelStrikes)-1) 
    L_SwingT(i) = L_HeelStrikes(i) - L_ToeOffs(i); % Since we always start with a toe off we can subtract from same index.
    end
else
    L_SwingT= zeros(1,length(L_HeelStrikes));
    for i = 1:(length(L_HeelStrikes))
    L_SwingT(i) = L_HeelStrikes(i) - L_ToeOffs(i); % Since we always start with a toe off we can subtract from same index.
    end
end

if length(R_HeelStrikes) == length(R_ToeOffs)
    R_SwingT= zeros(1,length(R_HeelStrikes)-1);% Always will start with a toe off so heel strikes > or = toe offs. Therefore do n or n-1 of R_HeelStrikes iterations.
    for i = 1:(length(R_HeelStrikes)-1) 
    R_SwingT(i) = R_HeelStrikes(i) - R_ToeOffs(i); % Since we always start with a toe off we can subtract from same index.
    end
else
    R_SwingT= zeros(1,length(R_HeelStrikes));
    for i = 1:length(R_HeelStrikes) 
    R_SwingT(i) = R_HeelStrikes(i) - R_ToeOffs(i); % Since we always start with a toe off we can subtract from same index.
    end
end

% Confirm swing count is the same so we can work with matrix
    if length(L_SwingT) ~= length(R_SwingT)
        if length(L_SwingT) > length(R_SwingT)
            L_SwingT(length(L_SwingT)) = [];

        else R_SwingT(length(R_SwingT)) = []; 
        end
    end

% Convert from frames to seconds
L_SwingT = cast(L_SwingT, "double");
L_SwingT_sec = L_SwingT/Rate;
R_SwingT = cast(R_SwingT, "double");
R_SwingT_sec = R_SwingT/Rate;

% Crop 5 gait cycles
L_SwingT_Crop = L_SwingT_sec;
L_SwingT_Crop(:,(1:5)) = [];
R_SwingT_Crop = R_SwingT_sec;
R_SwingT_Crop(:,(1:5)) = [];

% First and second half of trial
SwingT_Half = length(L_SwingT_Crop)/2; % Halfway point of trial (will be same for L and R)
L_SwingT_FH = L_SwingT_Crop(:,1:fix(SwingT_Half)); % First half (FH) of trial
L_SwingT_SH = L_SwingT_Crop(:,(fix(SwingT_Half))+1:end); % Second half (SH) of trial
R_SwingT_FH = R_SwingT_Crop(:,1:fix(SwingT_Half));
R_SwingT_SH = R_SwingT_Crop(:,(fix(SwingT_Half))+1:end);

% Means
if strcmp(Alimb,'L')
    A_SwingT_Avg = mean(L_SwingT_Crop);
    A_SwingT_Avg_FH = mean(L_SwingT_FH);
    A_SwingT_Avg_SH = mean(L_SwingT_SH);
    UnA_SwingT_Avg = mean(R_SwingT_Crop);
    UnA_SwingT_Avg_FH = mean(R_SwingT_FH);
    UnA_SwingT_Avg_SH = mean(R_SwingT_SH);
elseif strcmp(Alimb,'R')
    UnA_SwingT_Avg = mean(L_SwingT_Crop);
    UnA_SwingT_Avg_FH = mean(L_SwingT_FH);
    UnA_SwingT_Avg_SH = mean(L_SwingT_SH);
    A_SwingT_Avg = mean(R_SwingT_Crop);
    A_SwingT_Avg_FH = mean(R_SwingT_FH);
    A_SwingT_Avg_SH = mean(R_SwingT_SH);
end


% Determining Short and Long Swing Times 
if mean(L_SwingT) < mean(R_SwingT)
        ShortTime = L_StepDur_Crop; LongTime = R_StrideT_Crop;
    else
        ShortTime = R_StepDur_Crop; LongTime = L_StrideT_Crop;
end

%% Calculating phi (degrees) and PCI:
    % Total phi and PCI
     Phi= zeros(1,(length(L_SwingT_Crop)-1));
    for i = 1:(length(L_SwingT_Crop)-1)
        Phi(i) = 360*ShortTime(i)/LongTime(i);      
    end 

    Phi_diff = abs(Phi - 180); % Absolute difference of phi values
    Phi_abs = mean(Phi_diff); % Mean value of absolute differences (degrees)
    Pphi_abs = 100*Phi_abs/180; % Percentage-converted phi_ABS (%)
    Phi_cv = (std(Phi)/mean(Phi))*100; % Coefficient of variation of mean of phi (%)
    PCI = Phi_cv + Pphi_abs; % Calculate PCI


    % First half (FH) phi and PCI
    Phi_FH = Phi(1:(fix(length(Phi)/2)));
    Phi_diff_FH = abs(Phi_FH - 180); % Absolute difference of phi values
    Phi_abs_FH = mean(Phi_diff_FH); % Mean value of absolute differences (degrees)
    Pphi_abs_FH = 100*Phi_abs_FH/180; % Percentage-converted phi_ABS (%)
    Phi_cv_FH = (std(Phi_FH)/mean(Phi_FH))*100; % Coefficient of variation of mean of phi (%)
    PCI_FH = Phi_cv_FH + Pphi_abs_FH; % Calculate PCI

    % Second half (SH) phi and PCI
    Phi_SH = Phi((fix(length(Phi)/2)):end);
    Phi_diff_SH = abs(Phi_SH - 180); % Absolute difference of phi values
    Phi_abs_SH = mean(Phi_diff_SH); % Mean value of absolute differences (degrees)
    Pphi_abs_SH = 100*Phi_abs_SH/180; % Percentage-converted phi_ABS (%)
    Phi_cv_SH = (std(Phi_SH)/mean(Phi_SH))*100; % Coefficient of variation of mean of phi (%)
    PCI_SH = Phi_cv_SH + Pphi_abs_SH; % Calculate PCI

%% Step Length Asymetry (total, positional contribution, timing contribution, and velocity contribution) , Step Postion Asymetry, Step Time Asymetry
    % Contribution analyis based on temporal and spatial decomposition of SLA detailed in:
        % Finley JM, Long A, Bastian AJ, Torres-Oviedo G. Spatial and Temporal Control Contribute to Step Length Asymmetry During Split-Belt Adaptation and Hemiparetic Gait. Neurorehabilitation and Neural Repair. 2015;29(8):786-795. doi:10.1177/1545968314567149

    % If you only set one variable = vicon.GetTrajectory it will output only the x direction. If you want all variables use "[X, Y, Z, Exists]" as the variable
    % Because of our calibration, our x axis goes the width of the treadmill and y axis goes the length of the treadmill (opposite of what is typical)
    % so we grab the y values from our data and label them as x values for this code. 
    % Positive contributions correspond to larger values for the unaffected limb and negative contributions correspond to larger values for the affected limb
    % This analysis uses ankle position relative to ASIS postition to avoid the confound of translation on the treadmill

[L_AnkY,L_Ankx] = vicon.GetTrajectory(subnum, 'LANK'); % Switched the x and y axis because our calibration was opposite
[R_AnkY,R_Ankx] = vicon.GetTrajectory(subnum, 'RANK');
[L_AsisY,L_Asisx] = vicon.GetTrajectory(subnum, 'LASI');
[R_AsisY,R_Asisx] = vicon.GetTrajectory(subnum, 'RASI');

% Find trajectory locations at each heel strike
    % Ankle - ankle at each HS and ankle - hip at each HS

% Determine Start Leg
if R_HeelStrikes(1)<L_HeelStrikes(1)
    StartLeg = 1;
else
    StartLeg = 2;
end

% Created a copy of L & R heel strike matricies as lhs and rhs so it is not affected for other analysis
Rhs = R_HeelStrikes;
Lhs = L_HeelStrikes;
Rto = R_ToeOffs;
Lto = L_ToeOffs;

% Check if the start leg matches the slow leg; if not then skip first event. Then check if slow leg has one more event than the fast leg:if not
% then remove final fast leg event because it will not be used in calculations.

% There should always be one more slow leg heel strike than fast leg heel strike and this is made sure of by below if statement.

if StartLeg==1 && strcmpi(Alimb,'L')==1 % Right is first and slow
    if length(Rhs)<length(Lhs)
        Lhs = Lhs(1:length(Rhs));
        Lto = Lto(1:length(Rto));
    end
    if length(Rhs)==length(Lhs)
        Lhs = Lhs(1:end-1);
        Lto = Lto(1:end-1);
    end
    X_f = L_Ankx-L_Asisx; 
    X_s = R_Ankx-R_Asisx; 

elseif StartLeg==1 && strcmpi(Alimb,'R')==1 % Right is first but fast
    Rhs = Rhs(2:end);
    if length(Lhs)<length(Rhs)
        Rhs = Rhs(1:length(Lhs));
        Rto = Rto(1:length(Lto));
    end
    if length(Rhs)==length(Lhs)
        Rhs = Rhs(1:end-1);
        Rto = Rto(1:end-1);
    end
    X_s = L_Ankx-L_Asisx; 
    X_f = R_Ankx-R_Asisx; 

elseif StartLeg==2 && strcmpi(Alimb,'R')==1 %Left is first and slow
    if length(Lhs)<length(Rhs)
        Rhs = Rhs(1:length(Lhs));
        Rto = Rto(1:length(Lto));
    end
    if length(Rhs)==length(Lhs)
        Rhs = Rhs(1:end-1);
        Rto = Rto(1:end-1);
    end
    X_s = L_Ankx-L_Asisx; 
    X_f = R_Ankx-R_Asisx; 

elseif StartLeg==2 && strcmpi(Alimb,'L')==1 % Left is first but fast
    Lhs = Lhs(2:end);
    if length(Rhs)<length(Lhs)
        Lhs = Lhs(1:length(Rhs));
        Lto = Lto(1:length(Rto));
    end
    if length(Rhs)==length(Lhs)
        Lhs = Lhs(1:end-1);
        Lto = Lto(1:end-1);
    end
    X_f = L_Ankx-L_Asisx; 
    X_s = R_Ankx-R_Asisx; 
end

% Combine L and R vairables 
Num_Strikes = [length(Rhs), length(Lhs)];
Combined_HS = {Rhs,Lhs};
Combined_TO = {Rto,Lto};
Combined_Ank = {R_Ankx,L_Ankx};
Combined_Hip = {R_Asisx,L_Asisx};
Fast_Leg = strcmpi(Alimb,{'R','L'});

% Initialize variables to be used in loops
Xs_SHS1 = zeros(1,Num_Strikes(~Fast_Leg)-1); % Slow leg ankle-hip at slow leg heelstrke 1
Xf_SHS1 = Xs_SHS1; % Fast leg ankle-hip at slow leg heelstrke 1
Xs_SHS2 = Xs_SHS1; % Slow leg ankle-hip at slow leg heelstrke 2 (SHS1 shifted by one position)
Xf_SHS2 = Xs_SHS1; % Fast leg ankle-hip at slow leg heelstrke 2 (SHS1 shifted by one position)
Xf_FTO1 = Xs_SHS1; % Fast leg ankle-hip at fast toe off 1
Xs_STO1 = Xs_SHS1; % Slow leg ankle-hip at slow toe off 1
Xs_FHS = zeros(1,Num_Strikes(Fast_Leg));
Xf_FHS = Xs_FHS;
T_s = zeros(1,Num_Strikes(~Fast_Leg)-1);
T_f = T_s;

% Ankle location relative to hip location at slow leg heel strikes
for i = 1:Num_Strikes(~Fast_Leg)-1
    Xs_SHS1(i) = Combined_Ank{~Fast_Leg}(Combined_HS{~Fast_Leg}(i))-Combined_Hip{~Fast_Leg}(Combined_HS{~Fast_Leg}(i));
    Xf_SHS1(i) = Combined_Ank{Fast_Leg}(Combined_HS{~Fast_Leg}(i))-Combined_Hip{Fast_Leg}(Combined_HS{~Fast_Leg}(i));
    Xs_SHS2(i) = Combined_Ank{~Fast_Leg}(Combined_HS{~Fast_Leg}(i+1))-Combined_Hip{~Fast_Leg}(Combined_HS{~Fast_Leg}(i+1));
    Xf_SHS2(i) = Combined_Ank{Fast_Leg}(Combined_HS{~Fast_Leg}(i+1))-Combined_Hip{Fast_Leg}(Combined_HS{~Fast_Leg}(i+1));
end

% Ankle location relative to hip location at slow leg toe off
for i = 1:Num_Strikes(~Fast_Leg)-1
    Xs_STO1(i) = Combined_Ank{~Fast_Leg}(Combined_TO{~Fast_Leg}(i))-Combined_Hip{~Fast_Leg}(Combined_TO{~Fast_Leg}(i));
end

% Ankle location relative to hip location at fast leg heel strike
for i = 1:Num_Strikes(Fast_Leg)
    Xs_FHS(i) = Combined_Ank{~Fast_Leg}(Combined_HS{Fast_Leg}(i))-Combined_Hip{~Fast_Leg}(Combined_HS{Fast_Leg}(i));
    Xf_FHS(i) = Combined_Ank{Fast_Leg}(Combined_HS{Fast_Leg}(i))-Combined_Hip{Fast_Leg}(Combined_HS{Fast_Leg}(i));
end

% Ankle location relative to hip location at fast leg toe off
for i = 1:Num_Strikes(Fast_Leg)
    Xf_FTO1(i) = Combined_Ank{Fast_Leg}(Combined_TO{Fast_Leg}(i))-Combined_Hip{Fast_Leg}(Combined_TO{Fast_Leg}(i));
end

% Limb Excursion - Modified stride length - Distance covered from toe off to heel stike on same leg. On a treadmill true stride length will be 0
LE_s = Xs_SHS1-Xs_STO1; 
LE_f = Xf_FHS-Xf_FTO1;

% Limb Excursion Asymmetry
LEA = zeros(1,length(LE_s)-1);
for i = 1:length(LE_s)-1
    LEA(i) = LE_f(i)-LE_s(i);
 
end

% Step Length
SL_s = Xs_SHS2-Xf_SHS2;%Step length slow leg (ankle markers position at slow leg heel strike)
SL_f = Xf_FHS-Xs_FHS;%Step length fast leg (ankle markers position at fast leg heel strike)

% Step Length Asym (difference not a percentage): SL_f-SL_s
SLA = zeros(1,length(SL_s)-1);
for i = 1:length(SL_s)-1
     SLA(2*i-1) = SL_f(i)-SL_s(i);
     SLA(2*i) = SL_f(i)-SL_s(i+1);
end

% Step Times
for i = 1:Num_Strikes(~Fast_Leg)-1
    T_s(i) = (Combined_HS{~Fast_Leg}(i+1)-Combined_HS{Fast_Leg}(i));%time of slow HS2 - time of fast HS
    T_s(i) = cast(T_s(i),"double")/Rate; %change vairlable type for precision
    T_f(i) = (Combined_HS{Fast_Leg}(i)-Combined_HS{~Fast_Leg}(i)); %time of fast HS - time of slow HS1
    T_f(i) = cast(T_f(i),"double")/Rate; %change vairlable type for precision
end

% Step Timing Asymmetry
Step_Time = zeros(1,length(T_s)-1);
for i = 1:length(T_s)-1
    Step_Time(2*i-1) = (T_f(i)-T_s(i));
    Step_Time(2*i) = (T_f(i)-T_s(i+1));
end

% Step Position Asymmetry (also step position contribution to SlA)
% Positive step position contributions indicate that the unaffected limb progresses further in front of the trunk than the affected limb.
% Finley 2015 - "Individuals who have more paretic limb progression may prefer to walk with the paretic limb in increased flexion throughout the gait
% cycle while individuals with more non-paretic flexion may rely more heavily on non-paretic hip flexion for forward progress."

Alpha_s = Xs_SHS2-Xf_FHS;
Alpha_f = Xf_FHS-Xs_SHS1;
Alpha = zeros(1,length(Xf_FHS)-1);
for i = 1:length(Xf_FHS)-1
    Alpha(2*i-1) = Alpha_f(i)-Alpha_s(i);
    Alpha(2*i) = Alpha_f(i)-Alpha_s(i+1);
end

%Average foot velocity with respect to hips while during stance phase
V_s = (Xs_SHS1-Xs_FHS)./T_s;
V_f = (Xf_FHS-Xf_SHS2)./T_f;


%Step Velocity Asymmetry
Step_Vel = zeros(1,length(T_s)-1);
for i = 1:length(T_s)-1
    Step_Vel(2*i-1) = (V_f(i)-V_s(i));
    Step_Vel(2*i) = (V_f(i)-V_s(i+1));
end

% Step Timing Contribution to SLA 
% Negative step time contributions represent shorter stance-phase durations in the affected versus unaffected limb.
Step_Time_Contrib = zeros(1,length(T_s)-1);
for i = 1:length(T_s)-1
    Step_Time_Contrib(2*i-1) = ((V_s(i)+V_f(i))/2)*(T_s(i)-T_f(i));
    Step_Time_Contrib(2*i) = ((V_s(i+1)+V_f(i))/2)*(T_s(i+1)-T_f(i));
end

% Step Velocity Contribution to SLA 
% % Negative step velocity contributions are indicative of greater angular velocity of the unaffected limb, and this strategy would be expected from individuals with marked weakness of the paretic plantar flexor
Step_Vel_Contrib = zeros(1,length(T_s)-1);
for i = 1:length(T_s)-1
    Step_Vel_Contrib(2*i-1) = ((T_s(i)+T_f(i))/2)*(V_s(i)-V_f(i));
    Step_Vel_Contrib(2*i) = ((T_s(i+1)+T_f(i))/2)*(V_s(i+1)-V_f(i));
end
  
Step_Num = (1:length(SLA));

% SLA Table
SLA_Contrib_Results = [Step_Num', SLA', Alpha', Step_Time_Contrib', Step_Vel_Contrib'];
SLA_means = mean(SLA_Contrib_Results);


%% Find Peak and Mean Ground Reaction Force (z) During Stance (N/Kg)

L_GRFs = vicon.GetModelOutput(subnum, 'LGroundReactionForce');
R_GRFs= vicon.GetModelOutput(subnum, 'RGroundReactionForce');


% Create Matrix to store average and peak values of each stride for each leg
% Column 1 = L_Peak, Column 2 = L_Avg, Column 3 = R_Peak, Column 4 = R_Avg
GRF_StanceData = zeros(length(L_HeelStrikes),4);

   % Left Side
for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_GRFs(3,L_HeelStrikes(i):L_ToeOffs(i+1)); % Find window of GRF from each stance phase (heel strike to toe off), only want row 3 (z)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        GRF_StanceData(i,1) = 0;
        GRF_StanceData(i,2) = 0;
    else
    GRF_StanceData(i,1) = max(TempWindow);
    GRF_StanceData(i,2) = mean(TempWindow);
    end
end
  % Right Side
for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_GRFs(3,R_HeelStrikes(i):R_ToeOffs(i+1)); % Only row 3(z)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        GRF_StanceData(i,3) = 0;
        GRF_StanceData(i,4) = 0;
    else
    GRF_StanceData(i,3) = max(TempWindow);
    GRF_StanceData(i,4) = mean(TempWindow);
    end
end

if strcmp(Alimb,'L')

    % Means over all gait cycles in trial
    A_GRFPeakAvg = mean(nonzeros(GRF_StanceData(:,1)));
    A_GRFAvg = mean(nonzeros(GRF_StanceData(:,2)));
    UnA_GRFPeakAvg = mean(nonzeros(GRF_StanceData(:,3)));
    UnA_GRFAvg = mean(nonzeros(GRF_StanceData(:,4)));

    % First and second half means
    A_GRF_Half = length(GRF_StanceData(:,1))/2; % Halfway point of trial 
    A_GRFPeakAvg_FH = mean(nonzeros(GRF_StanceData(1:fix(A_GRF_Half),1))); % First half (FH) of trial
    A_GRFPeakAvg_SH = mean(nonzeros(GRF_StanceData((fix(A_GRF_Half)+1:end),1)));  % Second half (SH) of trial
    A_GRFAvg_FH = mean(nonzeros(GRF_StanceData(1:fix(A_GRF_Half),2))); 
    A_GRFAvg_SH = mean(nonzeros(GRF_StanceData((fix(A_GRF_Half)+1:end),2)));  

    UnA_GRF_Half = length(GRF_StanceData(:,3))/2; % Halfway point of trial 
    UnA_GRFPeakAvg_FH = mean(nonzeros(GRF_StanceData(1:fix(UnA_GRF_Half),3))); % First half (FH) of trial
    UnA_GRFPeakAvg_SH = mean(nonzeros(GRF_StanceData((fix(UnA_GRF_Half)+1:end),3)));  % Second half (SH) of trial
    UnA_GRFAvg_FH = mean(nonzeros(GRF_StanceData(1:fix(UnA_GRF_Half),4))); 
    UnA_GRFAvg_SH = mean(nonzeros(GRF_StanceData((fix(UnA_GRF_Half)+1:end),4))); 

elseif strcmp(Alimb,'R')

    % Means over all gait cycles in trial
    UnA_GRFPeakAvg = mean(nonzeros(GRF_StanceData(:,1)));
    UnA_GRFAvg = mean(nonzeros(GRF_StanceData(:,2)));
    A_GRFPeakAvg = mean(nonzeros(GRF_StanceData(:,3)));
    A_GRFAvg = mean(nonzeros(GRF_StanceData(:,4)));

    % First and second half means
    UnA_GRF_Half = length(GRF_StanceData(:,1))/2; % Halfway point of trial 
    UnA_GRFPeakAvg_FH = mean(nonzeros(GRF_StanceData(1:fix(UnA_GRF_Half),1))); % First half (FH) of trial
    UnA_GRFPeakAvg_SH = mean(nonzeros(GRF_StanceData((fix(UnA_GRF_Half)+1:end),1)));  % Second half (SH) of trial
    UnA_GRFAvg_FH = mean(nonzeros(GRF_StanceData(1:fix(UnA_GRF_Half),2))); 
    UnA_GRFAvg_SH = mean(nonzeros(GRF_StanceData((fix(UnA_GRF_Half)+1:end),2)));  

    A_GRF_Half = length(GRF_StanceData(:,3))/2; % Halfway point of trial 
    A_GRFPeakAvg_FH = mean(nonzeros(GRF_StanceData(1:fix(A_GRF_Half),3))); % First half (FH) of trial
    A_GRFPeakAvg_SH = mean(nonzeros(GRF_StanceData((fix(A_GRF_Half)+1:end),3)));  % Second half (SH) of trial
    A_GRFAvg_FH = mean(nonzeros(GRF_StanceData(1:fix(A_GRF_Half),4))); 
    A_GRFAvg_SH = mean(nonzeros(GRF_StanceData((fix(A_GRF_Half)+1:end),4))); 
 
end

%% Propulsion Forces - Same as above except using anterior GRFs and only Peak AVG
% Make sure we have the right direction (x, y or z) by grabbing the right row in L_GRFs and R_GRFs

% Create Matrix to store average and peak values of each stride for each leg
% Column 1 = L_Peak, Column 2 = R_Peak, 
AntGRF_StanceData = zeros(length(L_HeelStrikes),2);

   % Left Side
for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_GRFs(2,L_HeelStrikes(i):L_ToeOffs(i+1)); % Find window of GRF from each stance phase (heel strike to toe off), only want row 3 (z)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        AntGRF_StanceData(i,1) = 0;
    else
    AntGRF_StanceData(i,1) = max(TempWindow);
    end
end
  % Right Side
for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_GRFs(2,R_HeelStrikes(i):R_ToeOffs(i+1)); % Only row 3(z)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        AntGRF_StanceData(i,2) = 0;
    else
    AntGRF_StanceData(i,2) = max(TempWindow);
    end
end

if strcmp(Alimb,'L')

    % Means over all gait cycles in trial
    A_AntGRFPeakAvg = mean(nonzeros(AntGRF_StanceData(:,1)));
    UnA_AntGRFPeakAvg = mean(nonzeros(AntGRF_StanceData(:,2)));
   
    % First and second half means
    A_GRF_Half = length(AntGRF_StanceData(:,1))/2; % Halfway point of trial 
    A_AntGRFPeakAvg_FH = mean(nonzeros(AntGRF_StanceData(1:fix(A_GRF_Half),1))); % First half (FH) of trial
    A_AntGRFPeakAvg_SH = mean(nonzeros(AntGRF_StanceData((fix(A_GRF_Half)+1:end),1)));  % Second half (SH) of trial
   

    UnA_GRF_Half = length(AntGRF_StanceData(:,2))/2; % Halfway point of trial 
    UnA_AntGRFPeakAvg_FH = mean(nonzeros(AntGRF_StanceData(1:fix(UnA_GRF_Half),2))); % First half (FH) of trial
    UnA_AntGRFPeakAvg_SH = mean(nonzeros(AntGRF_StanceData((fix(UnA_GRF_Half)+1:end),2)));  % Second half (SH) of trial
  

elseif strcmp(Alimb,'R')

    % Means over all gait cycles in trial
    UnA_AntGRFPeakAvg = mean(nonzeros(AntGRF_StanceData(:,1)));
    A_AntGRFPeakAvg = mean(nonzeros(AntGRF_StanceData(:,2)));
    

    % First and second half means
    UnA_GRF_Half = length(AntGRF_StanceData(:,1))/2; % Halfway point of trial 
    UnA_AntGRFPeakAvg_FH = mean(nonzeros(AntGRF_StanceData(1:fix(UnA_GRF_Half),1))); % First half (FH) of trial
    UnA_AntGRFPeakAvg_SH = mean(nonzeros(AntGRF_StanceData((fix(UnA_GRF_Half)+1:end),1)));  % Second half (SH) of trial
   
    A_GRF_Half = length(AntGRF_StanceData(:,2))/2; % Halfway point of trial 
    A_AntGRFPeakAvg_FH = mean(nonzeros(AntGRF_StanceData(1:fix(A_GRF_Half),2))); % First half (FH) of trial
    A_AntGRFPeakAvg_SH = mean(nonzeros(AntGRF_StanceData((fix(A_GRF_Half)+1:end),2)));  % Second half (SH) of trial
 
end


%% Find Peak Ankle, Knee and Hip Flexion Angle During Swing

L_AnkleAng = vicon.GetModelOutput(subnum, 'LAbsAnkleAngle'); 
R_AnkleAng = vicon.GetModelOutput(subnum, 'RAbsAnkleAngle');
L_KneeAng = vicon.GetModelOutput(subnum, 'LKneeAngles');
R_KneeAng = vicon.GetModelOutput(subnum, 'RKneeAngles');
L_HipAng = vicon.GetModelOutput(subnum, 'LHipAngles');
R_HipAng = vicon.GetModelOutput(subnum, 'RHipAngles');

% Create Matrix to store peak values of each stride for each leg
% Column 1 = L_AnkPeak, Column 2 = RAnk_Peak, Column 3 = L_KneePeak, Column 4 = R_KneePeak, Column 5 = L_HipPeak, Column 6 = R_HipPeak
JointAng_Data = zeros(length(L_HeelStrikes),6);

   % Left Ankle Angles
for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_AnkleAng(1,L_ToeOffs(i):L_HeelStrikes(i)); % Find window of angles from each swing phase (toe off to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        JointAng_Data(i,1) = 0;
    else
        JointAng_Data(i,1) = max(TempWindow);
    end
end
    % Right Ankle Angles
for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_AnkleAng(1,R_ToeOffs(i):R_HeelStrikes(i)); % Find window of angles from each swing phase (toe off to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        JointAng_Data(i,2) = 0;
    else
        JointAng_Data(i,2) = max(TempWindow);
    end
end
  % Left Knee Angles
for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_KneeAng(1,L_ToeOffs(i):L_HeelStrikes(i)); % Only row 1 (x) for flexion/extension
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        JointAng_Data(i,3) = 0;
    else
        JointAng_Data(i,3) = max(TempWindow);
    end
end
    % Right Knee Angles
for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_KneeAng(1,R_ToeOffs(i):R_HeelStrikes(i)); % Only row 1 (x) for flexion/extension
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        JointAng_Data(i,4) = 0;
    else
        JointAng_Data(i,4) = max(TempWindow);
    end
end
    % Left Hip Angles
for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_HipAng(1,L_ToeOffs(i):L_HeelStrikes(i)); % Only row 1 (x) for flexion/extension
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        JointAng_Data(i,5) = 0;
    else
        JointAng_Data(i,5) = max(TempWindow);
    end
end    
    % Right Hip Angles
for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_HipAng(1,R_ToeOffs(i):R_HeelStrikes(i)); % Only row 1 (x) for flexion/extension
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        JointAng_Data(i,6) = 0;
    else
        JointAng_Data(i,6) = max(TempWindow);
    end
end

% Defining variables as affected (A) and unaffected (UnA)

if strcmp(Alimb,'L')

    % Means over all gait cycles in trial
    A_AnkleAng_PeakAvg = mean(nonzeros(JointAng_Data(:,1)));
    UnA_AnkleAng_PeakAvg = mean(nonzeros(JointAng_Data(:,2)));
    A_KneeAng_PeakAvg = mean(nonzeros(JointAng_Data(:,3)));
    UnA_KneeAng_PeakAvg = mean(nonzeros(JointAng_Data(:,4)));
    A_HipAng_PeakAvg = mean(nonzeros(JointAng_Data(:,5)));
    UnA_HipAng_PeakAvg = mean(nonzeros(JointAng_Data(:,6)));

    % First and second half means
    A_Angs_Half = length(JointAng_Data(:,1))/2; % Halfway point of trial 
    A_AnkleAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(A_Angs_Half),1))); % First half (FH) of trial
    A_AnkleAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(A_Angs_Half)+1:end),1)));  % Second half (SH) of trial
    A_KneeAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(A_Angs_Half),3))); 
    A_KneeAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(A_Angs_Half)+1:end),3)));  
    A_HipAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(A_Angs_Half),5))); 
    A_HipAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(A_Angs_Half)+1:end),5)));  

    UnA_Angs_Half = length(JointAng_Data(:,1))/2; % Halfway point of trial 
    UnA_AnkleAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(UnA_Angs_Half),2))); % First half (FH) of trial
    UnA_AnkleAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(UnA_Angs_Half)+1:end),2)));  % Second half (SH) of trial
    UnA_KneeAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(UnA_Angs_Half),4))); 
    UnA_KneeAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(UnA_Angs_Half)+1:end),4)));  
    UnA_HipAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(UnA_Angs_Half),6))); 
    UnA_HipAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(UnA_Angs_Half)+1:end),6)));  

elseif strcmp(Alimb,'R')

    % Means over all gait cycles in trial
    UnA_AnkleAng_PeakAvg = mean(nonzeros(JointAng_Data(:,1)));
    A_AnkleAng_PeakAvg = mean(nonzeros(JointAng_Data(:,2)));
    UnA_KneeAng_PeakAvg = mean(nonzeros(JointAng_Data(:,3)));
    A_KneeAng_PeakAvg = mean(nonzeros(JointAng_Data(:,4)));
    UnA_HipAng_PeakAvg = mean(nonzeros(JointAng_Data(:,5)));
    A_HipAng_PeakAvg = mean(nonzeros(JointAng_Data(:,6)));

    % First and second half means
    UnA_Angs_Half = length(JointAng_Data(:,1))/2; % Halfway point of trial 
    UnA_AnkleAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(UnA_Angs_Half),1))); % First half (FH) of trial
    UnA_AnkleAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(UnA_Angs_Half)+1:end),1)));  % Second half (SH) of trial
    UnA_KneeAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(UnA_Angs_Half),3))); 
    UnA_KneeAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(UnA_Angs_Half)+1:end),3)));  
    UnA_HipAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(UnA_Angs_Half),5))); 
    UnA_HipAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(UnA_Angs_Half)+1:end),5)));  

    A_Angs_Half = length(JointAng_Data(:,1))/2; % Halfway point of trial 
    A_AnkleAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(A_Angs_Half),2))); % First half (FH) of trial
    A_AnkleAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(A_Angs_Half)+1:end),2)));  % Second half (SH) of trial
    A_KneeAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(A_Angs_Half),4))); 
    A_KneeAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(A_Angs_Half)+1:end),4)));  
    A_HipAng_PeakAvg_FH = mean(nonzeros(JointAng_Data(1:fix(A_Angs_Half),6))); 
    A_HipAng_PeakAvg_SH = mean(nonzeros(JointAng_Data((fix(A_Angs_Half)+1:end),6)));  
end

%% Find Ankle Angle at Heel Strike
 % Vicon defines as standing (90 deg joint ang) as 0 deg
 % Does one frame give us an accurate assessment??
 % Not inluded in table for now - using peak ankle ang over gait cycle

  L_HeelStrikeAngs= zeros(1,length(L_HeelStrikes));
  R_HeelStrikeAngs= zeros(1,length(R_HeelStrikes));
 for i = 1:length(L_HeelStrikes) 
 L_HeelStrikeAngs(i) = (L_AnkleAng(1,L_HeelStrikes(i)));  % Finding ankle angle at each frame there is a heel strike
 end
 for i = 1:length(R_HeelStrikes) 
 R_HeelStrikeAngs(i) = (R_AnkleAng(1,R_HeelStrikes(i)));  % Finding ankle angle at each frame there is a heel strike
 end
 
%% Find Peak Ankle, Knee and Hip Moments and Powers During Each Gait Cycle

L_AnkleMom = vicon.GetModelOutput(subnum, 'LAnkleMoment'); 
R_AnkleMom = vicon.GetModelOutput(subnum, 'RAnkleMoment');
L_KneeMom = vicon.GetModelOutput(subnum, 'LKneeMoment');
R_KneeMom = vicon.GetModelOutput(subnum, 'RKneeMoment');
L_HipMom = vicon.GetModelOutput(subnum, 'LHipMoment');
R_HipMom = vicon.GetModelOutput(subnum, 'RHipMoment');
L_AnklePow = vicon.GetModelOutput(subnum, 'LAnklePower'); 
R_AnklePow = vicon.GetModelOutput(subnum, 'RAnklePower');
L_KneePow = vicon.GetModelOutput(subnum, 'LKneePower');
R_KneePow = vicon.GetModelOutput(subnum, 'RKneePower');
L_HipPow = vicon.GetModelOutput(subnum, 'LHipPower');
R_HipPow = vicon.GetModelOutput(subnum, 'RHipPower');

% Create Matrix to store peak values of each stride for each leg
% Column 1 = L_AnkMom, Column 2 = R_AnkMom, Column 3 = L_KneeMom, Column 4 = R_KneeMom, Column 5 = L_HipMom, Column 6 = R_HipMom, Column 7 = L_AnkPow, Column 8 = R_AnkPow, Column 9 = L_KneePow, Column 10 = R_KneePow, Column 11 = L_HipPow, Column 12 = R_HipPow
MomPow_Data = zeros(length(L_HeelStrikes),12);

    % Left Ankle Moments
    for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_AnkleMom(1,L_HeelStrikes(i):L_HeelStrikes(i+1)); % Find window of moments from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,1) = 0;
    else
        MomPow_Data(i,1) = max(TempWindow);
    end
    end

    % Right Ankle Moments
    for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_AnkleMom(1,R_HeelStrikes(i):R_HeelStrikes(i+1)); % Find window of moments from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,2) = 0;
    else
        MomPow_Data(i,2) = max(TempWindow);
    end
    end

    % Left Knee Moments
    for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_KneeMom(1,L_HeelStrikes(i):L_HeelStrikes(i+1)); % Find window of angles from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,3) = 0;
    else
        MomPow_Data(i,3) = max(TempWindow);
    end
    end

    % Right Knee Moments
    for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_KneeMom(1,R_HeelStrikes(i):R_HeelStrikes(i+1)); % Find window of moments from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,4) = 0;
    else
        MomPow_Data(i,4) = max(TempWindow);
    end
    end

    % Left Hip Moments
    for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_HipMom(1,L_HeelStrikes(i):L_HeelStrikes(i+1)); % Find window of angles from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,5) = 0;
    else
        MomPow_Data(i,5) = max(TempWindow);
    end
    end

    % Right Hip Moments
    for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_HipMom(1,R_HeelStrikes(i):R_HeelStrikes(i+1)); % Find window of moments from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,6) = 0;
    else
        MomPow_Data(i,6) = max(TempWindow);
    end
    end

    % Left Ankle Powers
    for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_AnklePow(3,L_HeelStrikes(i):L_HeelStrikes(i+1)); % Find window of angles from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,7) = 0;
    else
        MomPow_Data(i,7) = max(TempWindow);
    end
    end

    % Right Ankle Power
    for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_AnklePow(3,R_HeelStrikes(i):R_HeelStrikes(i+1)); % Find window of moments from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,8) = 0;
    else
        MomPow_Data(i,8) = max(TempWindow);
    end
    end

     % Left Knee Powers
    for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_KneePow(3,L_HeelStrikes(i):L_HeelStrikes(i+1)); % Find window of angles from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,9) = 0;
    else
        MomPow_Data(i,9) = max(TempWindow);
    end
    end

    % Right Knee Power
    for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_KneePow(3,R_HeelStrikes(i):R_HeelStrikes(i+1)); % Find window of moments from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,10) = 0;
    else
        MomPow_Data(i,10) = max(TempWindow);
    end
    end

     % Left Hip Powers
    for i = 1:(length(L_HeelStrikes)-1)
    TempWindow = L_HipPow(3,L_HeelStrikes(i):L_HeelStrikes(i+1)); % Find window of angles from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,11) = 0;
    else
        MomPow_Data(i,11) = max(TempWindow);
    end
    end

     % Right Hip Power
    for i = 1:(length(R_HeelStrikes)-1)
    TempWindow = R_HipPow(3,R_HeelStrikes(i):R_HeelStrikes(i+1)); % Find window of moments from each gait cycle (heel strike to heel strike)
    ZeroLocs = TempWindow(TempWindow == 0);
    HasZeros = isempty(ZeroLocs);
    if HasZeros == 0
        MomPow_Data(i,12) = 0;
    else
        MomPow_Data(i,12) = max(TempWindow);
    end
    end

if strcmp(Alimb,'L')

    % Means over all gait cycles in trial
    A_AnkleMom_PeakAvg = mean(nonzeros(MomPow_Data(:,1)));
    UnA_AnkleMom_PeakAvg = mean(nonzeros(MomPow_Data(:,2)));
    A_KneeMom_PeakAvg = mean(nonzeros(MomPow_Data(:,3)));
    UnA_KneeMom_PeakAvg = mean(nonzeros(MomPow_Data(:,4)));
    A_HipMom_PeakAvg = mean(nonzeros(MomPow_Data(:,5)));
    UnA_HipMom_PeakAvg = mean(nonzeros(MomPow_Data(:,6)));
    A_AnklePow_PeakAvg = mean(nonzeros(MomPow_Data(:,7)));
    UnA_AnklePow_PeakAvg = mean(nonzeros(MomPow_Data(:,8)));
    A_KneePow_PeakAvg = mean(nonzeros(MomPow_Data(:,9)));
    UnA_KneePow_PeakAvg = mean(nonzeros(MomPow_Data(:,10)));
    A_HipPow_PeakAvg = mean(nonzeros(MomPow_Data(:,11)));
    UnA_HipPow_PeakAvg = mean(nonzeros(MomPow_Data(:,12)));

    % First and second half means
    A_MomPow_Half = length(MomPow_Data(:,1))/2; % Halfway point of trial 
    A_AnkleMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),1))); % First half (FH) of trial
    A_AnkleMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),1)));  % Second half (SH) of trial
    A_KneeMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),3))); 
    A_KneeMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),3)));  
    A_HipMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),5))); 
    A_HipMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),5)));  
    A_AnklePow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),7))); % First half (FH) of trial
    A_AnklePow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),7)));  % Second half (SH) of trial
    A_KneePow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),9))); 
    A_KneePow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),9)));  
    A_HipPow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),11))); 
    A_HipPow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),11)));

    UnA_MomPow_Half = length(MomPow_Data(:,1))/2; % Halfway point of trial 
    UnA_AnkleMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),2))); % First half (FH) of trial
    UnA_AnkleMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),2)));  % Second half (SH) of trial
    UnA_KneeMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),4))); 
    UnA_KneeMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),4)));  
    UnA_HipMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),6))); 
    UnA_HipMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),6))); 
    UnA_AnklePow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),8))); % First half (FH) of trial
    UnA_AnklePow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),8)));  % Second half (SH) of trial
    UnA_KneePow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),10))); 
    UnA_KneePow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),10)));  
    UnA_HipPow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),12))); 
    UnA_HipPow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),12))); 

elseif strcmp(Alimb,'R')

    % Means over all gait cycles in trial
    UnA_AnkleMom_PeakAvg = mean(nonzeros(MomPow_Data(:,1)));
    A_AnkleMom_PeakAvg = mean(nonzeros(MomPow_Data(:,2)));
    UnA_KneeMom_PeakAvg = mean(nonzeros(MomPow_Data(:,3)));
    A_KneeMom_PeakAvg = mean(nonzeros(MomPow_Data(:,4)));
    UnA_HipMom_PeakAvg = mean(nonzeros(MomPow_Data(:,5)));
    A_HipMom_PeakAvg = mean(nonzeros(MomPow_Data(:,6)));
    UnA_AnklePow_PeakAvg = mean(nonzeros(MomPow_Data(:,7)));
    A_AnklePow_PeakAvg = mean(nonzeros(MomPow_Data(:,8)));
    UnA_KneePow_PeakAvg = mean(nonzeros(MomPow_Data(:,9)));
    A_KneePow_PeakAvg = mean(nonzeros(MomPow_Data(:,10)));
    UnA_HipPow_PeakAvg = mean(nonzeros(MomPow_Data(:,11)));
    A_HipPow_PeakAvg = mean(nonzeros(MomPow_Data(:,12)));

    % First and second half means
    UnA_MomPow_Half = length(MomPow_Data(:,1))/2; % Halfway point of trial 
    UnA_AnkleMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),1))); % First half (FH) of trial
    UnA_AnkleMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),1)));  % Second half (SH) of trial
    UnA_KneeMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),3))); 
    UnA_KneeMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),3)));  
    UnA_HipMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),5))); 
    UnA_HipMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),5)));  
    UnA_AnklePow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),7))); % First half (FH) of trial
    UnA_AnklePow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),7)));  % Second half (SH) of trial
    UnA_KneePow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),9))); 
    UnA_KneePow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),9)));  
    UnA_HipPow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(UnA_MomPow_Half),11))); 
    UnA_HipPow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(UnA_MomPow_Half)+1:end),11)));

    A_MomPow_Half = length(MomPow_Data(:,1))/2; % Halfway point of trial 
    A_AnkleMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),2))); % First half (FH) of trial
    A_AnkleMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),2)));  % Second half (SH) of trial
    A_KneeMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),4))); 
    A_KneeMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),4)));  
    A_HipMom_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),6))); 
    A_HipMom_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),6))); 
    A_AnklePow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),8))); % First half (FH) of trial
    A_AnklePow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),8)));  % Second half (SH) of trial
    A_KneePow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),10))); 
    A_KneePow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),10)));  
    A_HipPow_PeakAvg_FH = mean(nonzeros(MomPow_Data(1:fix(A_MomPow_Half),12))); 
    A_HipPow_PeakAvg_SH = mean(nonzeros(MomPow_Data((fix(A_MomPow_Half)+1:end),12))); 

end
%% Data Table

PCI_Cell = {subnum,Trial,Alimb,Pphi_abs,Phi_cv,PCI,Pphi_abs_FH, Phi_cv_FH, PCI_FH,Pphi_abs_SH,Phi_cv_SH, PCI_SH};
SLA_Cell = {mean(LE_f),mean(LE_f(1:(fix(length(LE_f)/2)))),mean(LE_f((fix(length(LE_f)/2)):end)),mean(LE_s),mean(LE_s(1:(fix(length(LE_s)/2)))),mean(LE_s((fix(length(LE_s)/2)):end)),mean(LEA),mean(LEA(1:(fix(length(LEA)/2)))),mean(LEA((fix(length(LEA)/2)):end)), mean(SL_f),mean(SL_f(1:(fix(length(SL_f)/2)))),mean(SL_f((length(SL_f/2)):end)),mean(SL_s),mean(SL_s(1:(fix(length(SL_s)/2)))),mean(SL_s((fix(length(SL_s)/2)):end)),mean(SLA),mean(SLA(1:(fix(length(SLA)/2)))),mean(SLA((fix(length(SLA)/2)):end)),mean(Alpha),mean(Alpha(1:(fix(length(Alpha)/2)))),mean(Alpha((fix(length(Alpha)/2)):end)), mean(Step_Time_Contrib),mean(Step_Time_Contrib(1:(fix(length(Step_Time_Contrib)/2)))),mean(Step_Time_Contrib((fix(length(Step_Time_Contrib)/2)):end)),mean(Step_Vel_Contrib),mean(Step_Vel_Contrib(1:(fix(length(Step_Vel_Contrib)/2)))),mean(Step_Vel_Contrib((fix(length(Step_Vel_Contrib)/2)):end)),mean(Step_Time),mean(Step_Time(1:(fix(length(Step_Time)/2)))),mean(Step_Time((fix(length(Step_Time)/2)):end)),mean(Step_Vel),mean(Step_Vel(1:(fix(length(Step_Vel)/2)))),mean(Step_Vel((fix(length(Step_Vel)/2)):end))};
Gait_Cell = {A_GCD_Avg,A_GCD_Avg_FH,A_GCD_Avg_SH,UnA_GCD_Avg,UnA_GCD_Avg_FH,UnA_GCD_Avg_SH,A_StepDur_Avg,A_StepDur_Avg_FH,A_StepDur_Avg_SH,UnA_StepDur_Avg,UnA_StepDur_Avg_FH,UnA_StepDur_Avg_SH,A_SwingT_Avg,A_SwingT_Avg_FH,A_SwingT_Avg_SH,UnA_SwingT_Avg,UnA_SwingT_Avg_FH,UnA_SwingT_Avg_SH};
GRF_Cell = {A_GRFPeakAvg,A_GRFPeakAvg_FH,A_GRFPeakAvg_SH,UnA_GRFPeakAvg,UnA_GRFPeakAvg_FH,UnA_GRFPeakAvg_SH,A_GRFAvg,A_GRFAvg_FH,A_GRFAvg_SH,UnA_GRFAvg,UnA_GRFAvg_FH,UnA_GRFAvg_SH,A_AntGRFPeakAvg,A_AntGRFPeakAvg_FH,A_AntGRFPeakAvg_SH,UnA_AntGRFPeakAvg,UnA_AntGRFPeakAvg_FH,UnA_AntGRFPeakAvg_SH};
JointAng_Cell = {A_AnkleAng_PeakAvg,A_AnkleAng_PeakAvg_FH,A_AnkleAng_PeakAvg_SH,UnA_AnkleAng_PeakAvg,UnA_AnkleAng_PeakAvg_FH,UnA_AnkleAng_PeakAvg_SH,A_KneeAng_PeakAvg,A_KneeAng_PeakAvg_FH,A_KneeAng_PeakAvg_SH,UnA_KneeAng_PeakAvg,UnA_KneeAng_PeakAvg_FH,UnA_KneeAng_PeakAvg_SH,A_HipAng_PeakAvg,A_HipAng_PeakAvg_FH,A_HipAng_PeakAvg_SH,UnA_HipAng_PeakAvg,UnA_HipAng_PeakAvg_FH,UnA_HipAng_PeakAvg_SH};
Moments_Cell = {A_AnkleMom_PeakAvg,A_AnkleMom_PeakAvg_FH,A_AnkleMom_PeakAvg_SH,UnA_AnkleMom_PeakAvg,UnA_AnkleMom_PeakAvg_FH,UnA_AnkleMom_PeakAvg_SH,A_KneeMom_PeakAvg,A_KneeMom_PeakAvg_FH,A_KneeMom_PeakAvg_SH,UnA_KneeMom_PeakAvg,UnA_KneeMom_PeakAvg_FH,UnA_KneeMom_PeakAvg_SH,A_HipMom_PeakAvg,A_HipMom_PeakAvg_FH,A_HipMom_PeakAvg_SH,UnA_HipMom_PeakAvg,UnA_HipMom_PeakAvg_FH,UnA_HipMom_PeakAvg_SH};
Powers_Cell = {A_AnklePow_PeakAvg,A_AnklePow_PeakAvg_FH,A_AnklePow_PeakAvg_SH,UnA_AnklePow_PeakAvg,UnA_AnklePow_PeakAvg_FH,UnA_AnklePow_PeakAvg_SH,A_KneePow_PeakAvg,A_KneePow_PeakAvg_FH,A_KneePow_PeakAvg_SH,UnA_KneePow_PeakAvg,UnA_KneePow_PeakAvg_FH,UnA_KneePow_PeakAvg_SH,A_HipPow_PeakAvg,A_HipPow_PeakAvg_FH,A_HipPow_PeakAvg_SH,UnA_HipPow_PeakAvg,UnA_HipPow_PeakAvg_FH,UnA_HipPow_PeakAvg_SH};
ViconBertec_DataTable = [PCI_Cell, SLA_Cell, Gait_Cell, GRF_Cell, JointAng_Cell, Moments_Cell, Powers_Cell];
ViconBertec_DataTable = cell2table(ViconBertec_DataTable,'VariableNames', {'Subject_Code', 'Gait_Trial', 'Affected_Limb', 'Phi_ABS', 'Phi_CV', 'PCI', 'Phi_ABS_FH', 'Phi_CV_FH', 'PCI_FH', 'Phi_ABS_SH', 'Phi_CV_SH', 'PCI_SH', 'LimbExcursion_A','LimbExcursion_A_FH','LimbExcursion_A_SH','LimbExcursion_UnA','LimbExcursion_UnA_FH','LimbExcursion_UnA_SH','LimbExcursion_Asymm','LimbExcursion_Asymm_FH','LimbExcursion_Asymm_SH','StepLength_A', 'StepLength_A_FH','StepLength_A_SH', 'StepLength_UnA', 'StepLength_UnA_FH','StepLength_UnA_SH', 'SLA', 'SLA_FH', 'SLA_SH', 'StepPositionAsym', 'StepPositionAsym_FH', 'StepPositionAsym_SH', 'StepTimeContrib', 'StepTimeContrib_FH', 'StepTimeContrib_SH', 'StepVeloctiyContrib', 'StepVelocityContrib_FH', 'StepVelocityContrib_SH', 'StepTimeAsym', 'StepTimeAsym_FH', 'StepTimeAsym_SH', 'StepVelocityAsym', 'StepVelocityAsym_FH', 'StepVelocityAsym_SH','GCD_A', 'GCD_A_FH', 'GCD_A_SH', 'GCD_UnA', 'GCD_UnA_FH', 'GCD_UnA_SH','StepDuration_A', 'StepDuration_A_FH', 'StepDuration_A_SH', 'StepDuration_UnA', 'StepDuration_UnA_FH', 'StepDuration_UnA_SH','SwingTime_A', 'SwingTime_A_FH', 'SwingTime_A_SH', 'SwingTime_UnA', 'SwingTime_UnA_FH', 'SwingTime_UnA_SH', 'GRFPeak_A', 'GRFPeak_A_FH', 'GRFPeak_A_SH', 'GRFPeak_UnA', 'GRFPeak_UnA_FH', 'GRFPeak_UnA_SH', 'GRFAvg_A', 'GRFAvg_A_FH', 'GRFAvg_A_SH', 'GRFAvg_UnA', 'GRFAvg_UnA_FH', 'GRFAvg_UnA_SH', 'AntGRFPeak_A','AntGRFPeak_A_FH','AntGRFPeak_A_SH','AntGRFPeak_UnA','AntGRFPeak_UnA_FH','AntGRFPeak_UnA_SH','AnkleAngPeak_A', 'AnkleAngPeak_A_FH', 'AnkleAngPeak_A_SH','AnkleAngPeak_UnA','AnkleAngPeak_UnA_FH',' AnkleAngPeak_UnA_SH','KneeAngPeak_A', 'KneeAngPeak_A_FH', 'KneeAngPeak_A_SH','KneeAngPeak_UnA', 'KneeAngPeak_UnA_FH','KneeAngPeak_UnA_SH','HipAngPeak_A', 'HipAngPeak_A_FH', 'HipAngPeak_A_SH','HipAngPeak_UnA','HipAngPeak_UnA_FH',' HipAngPeak_UnA_SH','AnkleMomentPeak_A', 'AnkleMomentPeak_A_FH', 'AnkleMomentPeak_A_SH','AnkleMomentPeak_UnA','AnkleMomentPeak_UnA_FH','AnkleMomentPeak_UnA_SH','KneeMomentPeak_A', 'KneeMomentPeak_A_FH', 'KneeMomentPeak_A_SH','KneeMomentPeak_UnA','KneeMomentPeak_UnA_FH','KneeMomentPeak_UnA_SH','HipMomentPeak_A', 'HipMomentPeak_A_FH', 'HipMomentPeak_A_SH','HipMomentPeak_UnA','HipMomentPeak_UnA_FH','HipMomentPeak_UnA_SH','AnklePowerPeak_A', 'AnklePowerPeak_A_FH', 'AnklePowerPeak_A_SH','AnklePowerPeak_UnA','AnklePowerPeak_UnA_FH','AnklePowerPeak_UnA_SH','KneePowerPeak_A', 'KneePowerPeak_A_FH', 'KneePowerPeak_A_SH','KneePowerPeak_UnA','KneePowerPeak_UnA_FH','KneePowerPeak_UnA_SH','HipPowerPeak_A', 'HipPowerPeak_A_FH', 'HipPowerPeakAvg_A_SH','HipPowerPeak_UnA','HipPowerPeak_UnA_FH','HipPowerPeak_UnA_SH'});

% Export Trial to Individual Data
ExportIndivData = inputdlg('Save and Export Data to Individual xlsx? Y/N: ');
 if strcmp(ExportIndivData,'Y')
    if exist([RootPath '\ViconBertec_Data.mat'], 'file') == 2
        load ([RootPath '\ViconBertec_Data.mat'])
        Indiv_Data = [Indiv_Data;ViconBertec_DataTable];
        save([RootPath '\ViconBertec_Data.mat'], 'Indiv_Data')
    else
         Indiv_Data = ViconBertec_DataTable;
        save([RootPath '\ViconBertec_Data.mat'], 'Indiv_Data')
    end
   writetable(Indiv_Data,[RootPath '\ViconBertec_Data.csv']);
 end

% Export Trial to Group Data
ExportGroupData = inputdlg('Save and Export to Group Data xlsx? Y/N: ');
 if strcmp(ExportGroupData,'Y')
      Outpath = 'R:\SBMS Study\Group_Data';
    if exist([Outpath '\ViconBertec_Group_Data.mat'], 'file') == 2
        load ([Outpath '\ViconBertec_Group_Data.mat'])
        Group_Data = [Group_Data;ViconBertec_DataTable];
        save("R:\SBMS Study\Group_Data\ViconBertec_Group_Data", 'Group_Data')
    else
         Group_Data = ViconBertec_DataTable;
        save("R:\SBMS Study\Group_Data\ViconBertec_Group_Data", 'Group_Data')
    end
    writetable(Group_Data,[Outpath '\ViconBertec_Group_Data.csv']);
 end

%% Create Figure of 3 PCI plots
PCIplot1 = figure (1);
        width=1500;
        height=1000;
        set(gcf,'position',[0,0,width,height])

    subplot(3,1,1);
    plot(Phi,'MarkerSize',8,'Marker','o','LineStyle','none');%Create Plot
    ylim([90 270]); %Set Y-axis
    title('PCI for Whole Trial','FontSize',12,'FontName','Arial');% Create title
    xlabel('Stride Number','FontSize',12,'FontName','Arial');% Create xlabel
    ylabel('Stepping Phase (°)','FontSize',12,'FontName','Arial');% Create ylabel
    hold on
    Xline = refline([0 180]);
    Xline.Color = 'r';
    Xline.LineWidth = 2;
    WPCI = num2str(PCI);
    text(2,100,['PCI = ' WPCI]);
    
    subplot(3,1,2);
    plot(Phi_FH,'MarkerSize',8,'Marker','o','LineStyle','none');%Create Plot
    ylim([90 270]); %Set Y-axis
    title('PCI for First Half of Trial','FontSize',12,'FontName','Arial');% Create title
    xlabel('Stride Number','FontSize',12,'FontName','Arial');% Create xlabel
    ylabel('Stepping Phase (°)','FontSize',12,'FontName','Arial');% Create ylabel
    hold on
    Xline = refline([0 180]);
    Xline.Color = 'r';
    Xline.LineWidth = 2;
    WPCI = num2str(PCI_FH);
    text(2,100,['PCI = ' WPCI]);

     subplot(3,1,3);
    plot(Phi_SH,'MarkerSize',8,'Marker','o','LineStyle','none');%Create Plot
    ylim([90 270]); %Set Y-axis
    title('PCI for Second Half of Trial','FontSize',12,'FontName','Arial');% Create title
    xlabel('Stride Number','FontSize',12,'FontName','Arial');% Create xlabel
    ylabel('Stepping Phase (°)','FontSize',12,'FontName','Arial');% Create ylabel
    hold on
    Xline = refline([0 180]);
    Xline.Color = 'r';
    Xline.LineWidth = 2;
    WPCI = num2str(PCI_SH);
    text(2,100,['PCI = ' WPCI]);

   ExportFig = inputdlg('Export Figures? Y/N: ');
    if strcmp(ExportFig,'Y')
        PlotPCI = 'Bertec PCI Figure' ;
        savefig(PCIplot1,PlotPCI);
        print([RootPath, '\', TrialName, ' Bertec PCI Figure'],'-dpng','-r0')
    else
    end

%% Plot SLA, aplha, step_time and step_velocity
SLA_Plot = figure(2);
width=1500;
height=1000;
set(gcf,'position',[0,0,width,height])
hold on
set(gca,'FontSize',14);
ylabel('Asymmetry (Absolute Difference in mm)','FontName','Arial')
xlabel('Step Number','FontName','Arial')
title('Step Length Asymmetry Contributions','FontName','Arial')

plot(SLA_Contrib_Results(:,1),SLA_Contrib_Results(:,2),'b','LineWidth',2)
plot(SLA_Contrib_Results(:,1),SLA_Contrib_Results(:,3),'r','LineWidth',1)
plot(SLA_Contrib_Results(:,1),SLA_Contrib_Results(:,4),'m','LineWidth',1)
plot(SLA_Contrib_Results(:,1),SLA_Contrib_Results(:,5),'g','LineWidth',1)
legend('Step Length Asymmetry','Spatial Component','Step Time Component','Step Velocity Component')


    if strcmp(ExportFig,'Y')
      PlotSLA = 'SLA Figure' ;
        savefig(SLA_Plot,PlotSLA);
        print([RootPath, '\', TrialName, ' SLA Figure'],'-dpng','-r0')
    else
    end

%% Plot LE and LEA
LEA_Plot = figure(3);
width=1500;
height=1000;
set(gcf,'position',[0,0,width,height])
hold on
set(gca,'FontSize',14);
ylabel('Excursion Distance (mm)','FontName','Arial')
xlabel('Gait Cycle Number','FontName','Arial')
title('Limb Excursion Asymmetry','FontName','Arial')

plot(LEA,'b','LineWidth',2)
plot(LE_s,'r','LineWidth',1)
plot(LE_f, 'g','LineWidth',1)
legend('Limb Excursion Asymmetry','Limb Excursion Unaffected','Limb Excursion Affected')


    if strcmp(ExportFig,'Y')
      PlotLEA = 'LEA Figure' ;
        savefig(LEA_Plot,PlotLEA);
        print([RootPath, '\', TrialName, ' LEA Figure'],'-dpng','-r0')
    else
    end
