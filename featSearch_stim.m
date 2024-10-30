% featSearch.m  -- 474 TRs, 355.5 sec, 5 min 55.5 secs
% 
% pilot script for motion/color feature guided search paradigm (DDT & TCS) 
% 
% Adapted from inCapture_2afc_bb
% 
% WHAT'S NEW?: -- added 'test' trials where there is a single item to see
%                 if there is anything wrong with reconstructions 
% 
% TODO 11/4/2021: -- 
% 
% 
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% subj info, filenames, etc  
clear;clc
% Screen('Preferences','SkipSyncTests',1); 
Screen('Preference','SyncTestSettings',.01,50,.1,5); % use this to avoid sync test annoyances
    % .01 = tolerance of deviation between actual and expected timing (no
        % larger than 10ms deviation
    % 50 = # of times it can fail before giving an error
    % .1 = the interval of time for each check
    % 5 = the # of times it checks an interval before counting it as a failure

p.root = pwd;

if ~isfolder([p.root '/data/'])
    mkdir([p.root '/data/']);
end

%get subject info
prompt =      {'SubjID','Session number', 'Run Number','Random seed', 'display (1 = fMRI, 2 = behavioral room)','Eye tracking','training?'};

%grab a random number and seed the generator (include this as a gui field
%in case we want to repeat an exact sequence)
s = round(sum(120*clock));
%put in some default answers
default_ans = {'sub999',    '1',              'X',                               num2str(s),     '1',                  '0',           '1'};
box = inputdlg(prompt,'Enter Subject Information...', 1, default_ans);


if length(box)==length(default_ans)
    p.subjID=char(box{1});
    p.session_num=str2double(box{2});
    p.run_num = str2double(box{3});
    p.rnd_seed=str2double(box{4});
    p.display=str2double(box{5});
    p.eye_tracking = str2double(box{6});
%     motion_help = str2double(box{7});
%     color_help = str2double(box{8}); 
    training = str2double(box{7}); 

    rng('default');
    rng(p.rnd_seed);  %actually seed the random number generator
else
    return
end

% FILENAMES
p.date = datestr(now,30);
p.fn = sprintf('%s/data/%s_sess%02.f_r%02.f_featSearch_%s.mat',p.root,p.subjID,p.session_num,p.run_num,p.date);

if p.eye_tracking == 1
    % for eyelink (8 char)
    p.et_fn = sprintf('s%s_T%02.f',p.subjID((end-2):end),p.run_num);
end

p.show_eye_pos = 0;

% font stuff
p.fontSize = 24;
p.fontName = 'ARIAL';
p.textColor = [100, 100, 100];

% difficulty adjustment stuff 
motion_help = 0; 
color_help = 1; 

%% general experiment structure
% load('color360.mat'); 
fullcolormatrix = hsv(360)*255; 


% trial structure setup 
p.numFeats = 8; % not sure if I will keep, number of features/locations

p.regTrials = 24; 
p.singleTarg = 6; 
p.targNdist = 6; 

p.ntrials = p.regTrials + p.singleTarg + p.targNdist; 


% conditions info by column: 1 target loc; 2 distractor loc; 3 remembered
% feature; 4 trial color; 5 trial motion; 6 distractor type; 7 target
% feature; 8 search target orientation; 9 single stim search
p.conditions = nan(p.ntrials,9); 

% making sure equal number of all trial types 
trial_type(:,1) = repmat([1;2],p.ntrials/2,1); 
trial_type(:,2) = repmat([1;2;3],p.ntrials/3,1); 

% should update column 3 and 6 so that they are flexible with num of trials
p.conditions(:,3) = trial_type(:,1);  % 1 = remember color; 2 = remember motion
% p.conditions(:,3) = repmat(1,1,p.ntrials)';
p.conditions(:,4) = randsample(linspace(0,315,p.numFeats),p.ntrials,'true');  % 8 possible bins for template feature
p.conditions(:,5) = randsample(linspace(0,315,p.numFeats),p.ntrials,'true');  % 8 possible color feature bins for search display 'other' items
p.conditions(:,6) = randsample(linspace(0,315,p.numFeats),p.ntrials,'true');  % 8 possible motion feature bins for search display 'other' items
p.conditions(:,7) = trial_type(:,2);  % 1 = no crit dist; 2 = non-feature dist; 3 = same-feature dist
% p.conditions(:,7) = repmat(3,1,p.ntrials)'; 
p.conditions(:,8) = randsample([1:2],p.ntrials,'true');  % 1 = horz target; 2 = vert target
p.conditions(p.ntrials-(p.singleTarg+p.targNdist)+1:end,9) = repmat([1 2],1,(p.singleTarg + p.targNdist)/2); % making only six trials single stim search and six two stim











% Daniel did something here to make all items appear on an array
p.conditions(:,9) = nan(p.ntrials,1);

% search target stimulus locations
p.conditions(:,1) = randsample([0:45:315]',p.ntrials,'true');

% search distractor stimulus locations
p.conditions(:,2) = randsample([0:45:315]',p.ntrials,'true');

p.conditions = Shuffle(p.conditions,2); 

% make sure distractor is not in the same position as target
for tt = 1:p.ntrials
    while p.conditions(tt,2) == p.conditions(tt,1)
        p.conditions(tt,2) = randsample([0:45:315]',1,'true');
    end 
end

% with target and distractor locations known, generate non-sing distractor locations
dist_locs = repmat([0:45:315],p.ntrials,1)';
for tt = 1:p.ntrials
    for dd = 1:p.numFeats
        if dist_locs(dd,tt) == p.conditions(tt,1) 
            dist_locs(dd,tt) = nan; 
        end 
        
        if dist_locs(dd,tt) == p.conditions(tt,2) 
            dist_locs(dd,tt) = nan; 
        end 
    end
end

dist_locs(isnan(dist_locs)) = [];
dist_locs = reshape(dist_locs,p.numFeats-2,p.ntrials);

% 8 feature bins, random value as close as one point away from the neighboring feature
% column 1 = memory sample jitter; column 2 = distractor jitter; column
% 3 = search others/target jitter. Also, spatial jitter.
for tt = 1:p.ntrials
    p.feature_jitter(tt,:) = repmat(randperm(360/p.numFeats,1),1,3); % currently three rows in case I want independent jitter for each stimulus 
end

p.spatial_jitter   = randperm(360/p.numFeats,p.ntrials);

% add spatial jitter to all search stimuli locations (target, distractor,
% and other) 
for tt = 1:p.ntrials 
    p.conditions(tt,1) = p.conditions(tt,1) + p.spatial_jitter(tt); 
    p.conditions(tt,2) = p.conditions(tt,2) + p.spatial_jitter(tt); 
    dist_locs(:,tt)    = dist_locs(:,tt)    + p.spatial_jitter(tt); 
end 

% add feature jitter to the color and motion directions of stimuli
for tt = 1:p.ntrials 
   p.conditions(tt,4) = p.conditions(tt,4) + p.feature_jitter(tt,1); 
   p.conditions(tt,5) = p.conditions(tt,5) + p.feature_jitter(tt,1); 
   p.conditions(tt,6) = p.conditions(tt,6) + p.feature_jitter(tt,1); 
end 

% jitter_range = 20; % features cannot be this similar

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% define general stimulus parameters: TODO: add values for the response circle
p.aperture_size  = 9.15; % radius from fixation (for full gray aperture)
p.fix_size       = 0.25;
p.fix_size_inner = 0.04;
p.lineApt_size   = 0.5; % size of aperture within moving stimulus to help perceive lines
p.stim_size      = 1.5;  % SHOULD MATCH SIZE OF MAPPING STIMULUS
p.stim_ecc       = 5;
p.mem_resp_size  = 8;

% stim loc in deg
p.distLocsDeg = [cosd(p.conditions(:,2)) sind(p.conditions(:,2))]*p.stim_ecc.*[1 -1];

p.dot_size    = 0.05 ; % was .025
p.dot_speed   = 60*0.075 + motion_help; % deg/sec (approx...rotational speed harder to define this way...) -- was 60*.075
p.dot_life    = 0.2; % seconds (remember to turn to frames!) -- was 0.05
p.dot_density = 15;  % dots/deg^2 (TODO: convert to # of dots!) -- was 40

% compute p.ndots using density, area
p.ndots = round( p.stim_size.^2 * pi * p.dot_density );

% COLORS
p.fix_color  = [180 180 180];
p.bg_color   = [0 0 0];
p.text_color = [180 180 180];

% 1 = motion cond colors; 2 = color cond colors
p.color_cond = round(linspace(1,360,p.numFeats)); 
p.dot_color = cell(length(unique(p.conditions(:,3))),1);

% testing some code to adjust saturation of colors 
fullcolormatrix = rgb2hsv(fullcolormatrix/255); 
fullcolormatrix = [fullcolormatrix(:,1) fullcolormatrix(:,2)/color_help fullcolormatrix(:,3)]; 
fullcolormatrix = hsv2rgb(fullcolormatrix)*255; 


p.dot_color{1} = fullcolormatrix; 
% p.dot_color{1} = fullcolormatrix(p.color_cond,:); 
p.dot_color{2} = [180]; % assumes uniform color inputs 

p.feature_cond = [0:45:315]; 
% p.feature_cond = repmat(90,8,1)'; 


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% define experiment timing parameters -- seconds

% beginning of experiment: wait after first trigger (after dummy scans)
p.start_wait = 3.0; % after first trigger, wait a bit

% end of experiment, after last trial, how long to keep scanning
p.end_wait   = 10.5; % was 10.5


% each trial we:
% - show attn cue at fixation
% - blank (?)
% - stimulus onset (TR-locked)
% - ITI

% NOTE: cue_dur, post_cue_dur, and ITI should sum to integer # of TRs
p.mem_dur        = 1.0;
p.postmem_dur    = 0.5;
p.search_dur     = 2.0; 

ITI_vals = [6 8 10] * 0.75; % was 8 10 12

if training == 1
    ITI_vals = [4 4 4] * 0.75;
end
% OR: generate a list
%ITI_vals = [8 9 9 9 10 10 10 11 11 13 13 14]
% ITI_vals = repmat(2,p.ntrials,1); % TODO: will have to update this before pilot
% p.ITIs = ITI_vals; 
p.ITIs = Shuffle(repmat(ITI_vals.',p.ntrials/length(ITI_vals),1));
% p.ITIs = p.ITIs(randperm(p.ntrials));


% estimated experiment duration:
p.expt_dur = p.start_wait + p.ntrials * (p.mem_dur + p.postmem_dur + p.search_dur) + sum(p.ITIs) + p.end_wait;
% p.expt_end - p.expt_start should be very close to p.expt_dur


% set up the variables we save timing info into
t.expt_start  = nan;
t.trial_start = nan(p.ntrials,1);
t.search_start  = nan(p.ntrials,1);
t.search_end    = nan(p.ntrials,1);
t.trial_end = nan(p.ntrials,1); 
t.expt_end    = nan;

% (other things...)


% font stuff
p.font_size = 13; % TODO: scale this by an expected ppd?
p.font_name = 'ARIAL';

%Response keys
if ismac
    p.escape = 41;
else
    p.escape = 27;
end


% (keyboard)

% TODO: also allow left/right arrow...
p.keys = [KbName('1!'),KbName('2@'),KbName('3#'),KbName('4$')]; % keys 1 through 4
p.space = KbName('space');
p.start = KbName('5%');


% print out total experiment duration
fprintf('TOTAL EXPERIMENT DURATION: %.02f s\n',p.expt_dur);


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% define condition-specific stimulus parameters
% including monitor info
if p.display == 1     % fMRI
    p.scr_height = 35.5; % cm
    p.viewing_distance = 110; % cm
elseif p.display == 2 % behavior
    p.scr_height = 30;
    p.viewing_distance = 52; % cm
end

% SCREEN COMES ALIVE
AssertOpenGL; % bail if current version of PTB does not use
s=max(Screen('Screens'));

HideCursor;
[w, p.s_rect] = Screen('OpenWindow', s, [0 0 0]);


% turn on eyetracker (hello eyetracker!)
if p.eye_tracking == 1
    
    
    el=EyelinkInitDefaults(w);
    
    el.backgroundcolour=p.bg_color(1);  % TODO: fix this?
    el.calibrationtargetcolour=p.fix_color(1);
    %el.calibrationtargetwidth=1;
    el.msgfontcolour=p.fix_color(1);
    p.foregroundcolour=p.fix_color(1);
    
    EyelinkUpdateDefaults(el);
    
    
    Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
    
    % SCANNER: right eye!!!!!! TODO: script this, record output...
    Eyelink('command','calibration_type=HV5'); % updating number of callibration dots
    Eyelink('command','link_sample_data=LEFT,RIGHT,GAZE,AREA');% (,GAZERES,HREF,PUPIL,STATUS,INPUT');
    Eyelink('command', 'sample_rate=500');
    Eyelink('command','screen_pixel_coords=%ld %ld %ld %ld', 0, 0, p.s_rect(3)-1,p.s_rect(4)-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, p.s_rect(3)-1,p.s_rect(4)-1);
    Eyelink('command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    Eyelink('command','file_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS');
    
    
    % make sure that we get gaze data from the Eyelink
    
    %------ calibrate the eye tracker --------
    EyelinkDoTrackerSetup(el);
    Eyelink('openfile',p.et_fn);
    
end

Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

p.ifi=Screen('GetFlipInterval', w);   % inter-frame-time
p.fps = round(1/p.ifi);
if p.fps==0                           % if fps does not register, then set the fps based on ifi
    p.fps=1/p.ifi;
end


% TODO: error handling here....

% compute ppd from viewing distance, height
p.ppd = 0.5*(p.s_rect(4)-p.s_rect(2)) / atan2d(p.scr_height/2, p.viewing_distance);

p.dot_life_frames = ceil(p.dot_life * p.fps);
p.dot_speed_frames = p.dot_speed/p.fps; % distance traveled per frame
p.motion_coh =1; 

p.memstim_frames = p.fps*p.mem_dur;
p.search_frames  = p.fps*p.search_dur;

% for both p.ds and p.dc, two features are initially generated for each
% trial (columns 1 and 2). Memory test values are taken from these. For the
% value that isn't used (e.g., motion on a remember color trial), this
% feature becomes the singleton feature on non-feature dimension distractor
% trials. within-feature dimension values generated using while loop.
% Distractor values are placed in column three. Finally, the other
% distractors in the search array are placed in column 4. 

% compute all dot sequences -- this is where separate dot sequences should
% be computed for each part of the trial (mem and search)
p.ds = cell(p.ntrials,3);
for tt = 1:p.ntrials
    % generating dots for the search template -- motion
    try
        p.ds{tt,1} = make_dot_seq_planar(p.ndots,p.conditions(tt,4),p.fps*p.search_dur,p.dot_speed_frames,p.motion_coh,p.dot_life_frames);
    catch
        p.ds{tt,1} = make_dot_seq_planar(p.ndots,p.conditions(tt,4),p.fps*p.search_dur,p.dot_speed_frames,p.motion_coh,p.dot_life_frames);
    end
    
    % generating dots for the search template -- color
    try
        p.dsc{tt} = make_dot_seq_planar(p.ndots,p.conditions(tt,4),p.fps*p.search_dur,0,p.motion_coh,p.dot_life_frames);
    catch
        p.dsc{tt} = make_dot_seq_planar(p.ndots,p.conditions(tt,4),p.fps*p.search_dur,0,p.motion_coh,p.dot_life_frames);
    end
    
    % shuffle the dots!
    p.ds{tt,1} = p.ds{tt,1}(:,randperm(size(p.ds{tt,1},2)),:);
end


% generating dot colors for each condition -- 1st column color; 2nd
% column gray
p.dc = cell(p.ntrials,3);
p.dcg{1} = repmat(repelem(p.dot_color{2},3,1)',p.ndots,1)';
for tt = 1:p.ntrials
    p.dc{tt,1} = repmat(p.dot_color{1}(p.conditions(tt,4),:),p.ndots,1)';
end

%%  computing dot motion and color for the target/other dot arrays in search
% and changing distractor stimulus if non-feature dimension singleton trial
for tt = 1:p.ntrials
    % search array 'other' items dimensions
    
    % color
    while p.conditions(tt,4) == p.conditions(tt,5)
        p.conditions(tt,5) = p.feature_cond(randperm(p.numFeats,1)) + p.feature_jitter(tt,1);
    end
    
    % index colors
    p.dc{tt,3} = repelem(p.dot_color{1}(p.conditions(tt,5),:)',1,p.ndots);
    
 
    % motion
    while p.conditions(tt,4) == p.conditions(tt,6)
        p.conditions(tt,6) = p.feature_cond(randperm(p.numFeats,1)) + p.feature_jitter(tt,1);
    end
    
    % generating dots for the search array
    try
        p.ds{tt,3} = make_dot_seq_planar(p.ndots,p.conditions(tt,6),p.fps*p.search_dur,p.dot_speed_frames,p.motion_coh,p.dot_life_frames);
    catch
        p.ds{tt,3} = make_dot_seq_planar(p.ndots,p.conditions(tt,6),p.fps*p.search_dur,p.dot_speed_frames,p.motion_coh,p.dot_life_frames);
    end
    
    p.ds{tt,3} = p.ds{tt,3}(:,randperm(size(p.ds{tt,3},2)),:);

    % adjusting feature of target item if depending of template feature
    if p.conditions(tt,3) == 1
        p.ds{tt,1} = p.ds{tt,3}; 
    elseif p.conditions(tt,3) == 2
        p.dc{tt,1} = p.dc{tt,3};
    end
    
    % singleton feature values for all singleton-present conditions
    if p.conditions(tt,4) > 180
        p.singleton_feat(tt) = p.conditions(tt,4) - 180;
    elseif p.conditions(tt,4) <= 180
        p.singleton_feat(tt) = p.conditions(tt,4) + 180;
    end
    
    % singleton dimensions when it is absent
    if p.conditions(tt,7) == 1
        p.dc{tt,2} = p.dc{tt,3};
        p.ds{tt,2} = p.ds{tt,3};
        
    elseif p.conditions(tt,3) == 1 && p.conditions(tt,7) == 2 % color template & motion singleton
        % generating dots for the search array singleton
        try
            p.ds{tt,2} = make_dot_seq_planar(p.ndots,p.singleton_feat(tt),p.fps*p.search_dur,p.dot_speed_frames,p.motion_coh,p.dot_life_frames);
        catch
            p.ds{tt,2} = make_dot_seq_planar(p.ndots,p.singleton_feat(tt),p.fps*p.search_dur,p.dot_speed_frames,p.motion_coh,p.dot_life_frames);
        end
        
        % shuffle the dots!
        p.ds{tt,2} = p.ds{tt,2}(:,randperm(size(p.ds{tt,2},2)),:);
        p.dc{tt,2} = p.dc{tt,3}; 
        
    elseif p.conditions(tt,3) == 1 && p.conditions(tt,7) == 3 % color template & color singleton
        p.dc{tt,2} = repelem(p.dot_color{1}(p.singleton_feat(tt),:)',1,p.ndots);
        p.ds{tt,2} = p.ds{tt,3}; 
        
    elseif p.conditions(tt,3) == 2 && p.conditions(tt,7) == 2 % motion template & color singleton
        p.dc{tt,2} = repelem(p.dot_color{1}(p.singleton_feat(tt),:)',1,p.ndots);
        p.ds{tt,2} = p.ds{tt,3}; 
        
    elseif p.conditions(tt,3) == 2 && p.conditions(tt,7) == 3 % motion template & motion singleton
         % generating dots for the search array singleton
        try
            p.ds{tt,2} = make_dot_seq_planar(p.ndots,p.singleton_feat(tt),p.fps*p.search_dur,p.dot_speed_frames,p.motion_coh,p.dot_life_frames);
        catch
            p.ds{tt,2} = make_dot_seq_planar(p.ndots,p.singleton_feat(tt),p.fps*p.search_dur,p.dot_speed_frames,p.motion_coh,p.dot_life_frames);
        end
        
        % shuffle the dots!
        p.ds{tt,2} = p.ds{tt,2}(:,randperm(size(p.ds{tt,2},2)),:);
        p.dc{tt,2} = p.dc{tt,3}; 
        
    end
end

% parameters for lines in the middle of search array stimuli -- in pixels
p.xLineTarg = 12; 
p.yLineTarg = 12; 
line_width  = 4; 


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pre-experiment setup 

% allocate response variables
p.resp     = nan(p.ntrials,1); % for search
p.RT       = nan(p.ntrials,1); % for search 
p.correct  = nan(p.ntrials,1); % for search 
p.mem_resp = nan(p.ntrials,1); % for mem test 
p.mem_RT   = nan(p.ntrials,1); 

% initialize variable to keep track of flips
t.frame_start = nan(p.ntrials,p.search_frames);

% beginning of experiment - welcome screen

% store center point
%center_point = [p.s_rect(3)/2 p.s_rect(4)/2]; % OR
[center_point(1),center_point(2)] = RectCenter(p.s_rect);

% draw aperture
aperture_rect = CenterRect([0 0 p.aperture_size*2 p.aperture_size*2] * p.ppd,p.s_rect);
outline_rect = CenterRect([0 0 p.aperture_size*1.99 p.aperture_size*1.99] * p.ppd,p.s_rect);
rim_rect = CenterRect([0 0 p.aperture_size*1.98 p.aperture_size*1.98] * p.ppd,p.s_rect);

Screen('FillRect',w,0);
Screen('FillOval',w,p.bg_color,aperture_rect);
Screen('FillOval',w,p.fix_color,outline_rect);
Screen('FillOval',w,p.bg_color,rim_rect);

% draw fixation (dimmed/big) (TODO: correctly compute line width...)
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color*0.8,center_point,3);
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color,     center_point,3);
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color*0.8,center_point,3);


Screen('Flip',w);


% WAIT TO START SCANNER
resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.start, p.escape);
    if resp == -1
        Screen('CloseAll'); ShowCursor;
        if p.eye_tracking==1
            Eyelink('ShutDown'); 
        end
        return;
    end
end
clear resp;
t.expt_start = GetSecs;

if p.eye_tracking == 1
    Eyelink('Message','xDAT %i', 101);
    Eyelink('StartRecording'); % make 1 big edf file (save time)
end


Screen('FillRect',w,0);
Screen('FillOval',w,p.bg_color,aperture_rect);
Screen('FillOval',w,p.fix_color,outline_rect);
Screen('FillOval',w,p.bg_color,rim_rect);

% draw fixation (dimmed/big) (TODO: correctly compute line width...)
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color,center_point,3);
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color, center_point,3);
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color,center_point,3);

Screen('Flip',w);

% start wait
resp = 0;
while (GetSecs-t.expt_start) < p.start_wait
    [resp, ~] = checkForResp([], p.escape);
    if resp == -1
        sca;
        fprintf('ESC pressesd during post-trigger wait time\n');
        ShowCursor;
        if p.eye_tracking==1
            Eyelink('ShutDown');
        end
        return;
    end
end


%% trial loop
start_check = GetSecs; 

% for each trial
for tt = 1:p.ntrials
    
    t.trial_start(tt) = GetSecs;
    
    % present memory stimulus
    for ff = 1:p.memstim_frames
        
        if p.eye_tracking == 1
            
            Eyelink('Message','xDAT %i',1);
            
            Eyelink('command', 'record_status_message "Trial %d of %d"', tt, p.ntrials);
            
        end
        
        % background
        Screen('FillRect',w,0);
        Screen('FillOval',w,p.bg_color,aperture_rect);
        Screen('FillOval',w,p.fix_color,outline_rect);
        Screen('FillOval',w,p.bg_color,rim_rect);
        
        % stim
        if p.conditions(tt,3) == 1
            Screen('DrawDots',w,p.ppd*p.stim_size*p.dsc{1}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,1},center_point,3);
        elseif p.conditions(tt,3) == 2
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,1}(:,:,ff),p.dot_size*p.ppd*2,p.dcg{1},center_point,3);
        end
        
        % draw fixation
        Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color,center_point,3);
        Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color, center_point,3);
        Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color,center_point,3);
        
        Screen('Flip',w);

        % quick screenshot test
%         if p.conditions(tt,3) == 1
%             if ff == 1 && p.conditions(tt,7) == 1
%                 myImageMemNone = Screen('GetImage',w);
% 
%             elseif ff == 1 && p.conditions(tt,7) == 2
%                 myImageMemMismatch = Screen('GetImage',w);
% 
%             elseif ff == 1 && p.conditions(tt,7) == 3
%                 myImageMemMatch = Screen('GetImage',w);
% 
%             end
%         end
%         
        % check for ACTUAL responses
        [resp, timeStamp] = checkForResp(p.keys, p.escape);
        if resp==-1
            sca; ShowCursor;
            if p.eye_tracking==1
                Eyelink('ShutDown');
            end
            return;
        end
    end
    
    if p.eye_tracking == 1
        
        Eyelink('Message','xDAT %i',1);
        
        Eyelink('command', 'record_status_message "Trial %d of %d"', tt, p.ntrials);
        
    end
    
    % background
    Screen('FillRect',w,0);
    Screen('FillOval',w,p.bg_color,aperture_rect);
    Screen('FillOval',w,p.fix_color,outline_rect);
    Screen('FillOval',w,p.bg_color,rim_rect);
    
    % draw fixation
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color,center_point,3);
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color, center_point,3);
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color,center_point,3);
    
    Screen('Flip',w);
    
    
    while GetSecs < (t.trial_start(tt) + p.mem_dur)
        [resp, ~] = checkForResp([], p.escape);
        if resp == -1
            sca;
            ShowCursor;
            if p.eye_tracking==1
                Eyelink('ShutDown');
            end
            return;
        end
    end
    
    % pre-search wait
    Screen('FillRect',w,0);
    Screen('FillOval',w,p.bg_color,aperture_rect);
    Screen('FillOval',w,p.fix_color,outline_rect);
    Screen('FillOval',w,p.bg_color,rim_rect);
    
    % draw fixation
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color,center_point,3);
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color, center_point,3);
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color,center_point,3);
    
    Screen('Flip',w);
    
    if p.eye_tracking == 1
        Eyelink('Message','xDAT %i',2);
    end
    
    while GetSecs < (t.trial_start(tt) + p.mem_dur + p.postmem_dur)
        [resp, ~] = checkForResp([], p.escape);
        if resp == -1
            sca;
            ShowCursor;
            if p.eye_tracking==1
                Eyelink('ShutDown');
            end
            
            return;
        end
    end

    
    % SEARCH STIMULUS
    
    % stimuli positions: (Cartesian coordinates; dva)
    
    % target
    p.targ_pos(tt,:) = [cosd(p.conditions(tt,1)) sind(p.conditions(tt,1))]*p.stim_ecc;
    stim_pos_pix(1,:) = center_point + p.ppd*p.targ_pos(tt,:).*[1 -1]; % center in pixels
    
    % distractor
    p.dist_pos(tt,:) = [cosd(p.conditions(tt,2)) sind(p.conditions(tt,2))]*p.stim_ecc; % repeated above
    stim_pos_pix(2,:) = center_point + p.ppd*p.dist_pos(tt,:).*[1 -1]; % center in pixels
    
    % others
    for dd = 1:p.numFeats-2
        p.other_pos(dd,:) = [cosd(dist_locs(dd,tt)) sind(dist_locs(dd,tt))]*p.stim_ecc;
        stim_pos_pix(dd+2,:) = center_point + p.ppd*p.other_pos(dd,:).*[1 -1]; % center in pixels
    end
    
    line_values = [repelem(stim_pos_pix(:,1),2)';repelem(stim_pos_pix(:,2),2)']; 
    line_ori    = randsample([1,2],8,'true'); 
    line_ori(1) = p.conditions(tt,8);
    
    for dd = 1:p.numFeats
        if line_ori(dd) == 1
            line_values(1,dd*2-1) = line_values(1,dd*2-1) - p.xLineTarg;
            line_values(1,dd*2) = line_values(1,dd*2) + p.xLineTarg;
        elseif line_ori(dd) == 2
            line_values(2,dd*2-1) = line_values(2,dd*2-1) + p.yLineTarg;
            line_values(2,dd*2) = line_values(2,dd*2) - p.yLineTarg;
        end
    end
    
    % might need to change this for more effecient timing later
    for ff = 1:p.search_frames
        
        % fixation point & aperture
        Screen('FillRect',w,0);
        Screen('FillOval',w,p.bg_color,aperture_rect);
        Screen('FillOval',w,p.fix_color,outline_rect);
        Screen('FillOval',w,p.bg_color,rim_rect);
        
        % draw fixation (dimmed/big) (TODO: correctly compute line width...)
        Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color,center_point,3);
        Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color, center_point,3);
        Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color,center_point,3);

        if p.conditions(tt,9) == 1 % single target cond

            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,1}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,1},stim_pos_pix(1,:),3); % target
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(1,:),3); % drawing circles around lines
            Screen('DrawLines',w,line_values(:,1:2),line_width,p.fix_color-40); % should make this into a unique color variable...

        elseif p.conditions(tt,9) == 2 % two stimulus (targ and dist) cond

            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,1}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,1},stim_pos_pix(1,:),3); % target
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,2}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,2},stim_pos_pix(2,:),3); % distractor

            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(1,:),3); % drawing circles around lines
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(2,:),3);

            Screen('DrawLines',w,line_values(:,1:4),line_width,p.fix_color-40); % should make this into a unique color variable...

        else % normal stimuli 

            % draw search stimuli!
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,1}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,1},stim_pos_pix(1,:),3); % target
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,2}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,2},stim_pos_pix(2,:),3); % distractor
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,3}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,3},stim_pos_pix(3,:),3); % others
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,3}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,3},stim_pos_pix(4,:),3);
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,3}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,3},stim_pos_pix(5,:),3);
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,3}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,3},stim_pos_pix(6,:),3);
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,3}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,3},stim_pos_pix(7,:),3);
            Screen('DrawDots',w,p.ppd*p.stim_size*p.ds{tt,3}(:,:,ff),p.dot_size*p.ppd*2,p.dc{tt,3},stim_pos_pix(8,:),3);

            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(1,:),3); % drawing circles around lines
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(2,:),3);
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(3,:),3);
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(4,:),3);
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(5,:),3);
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(6,:),3);
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(7,:),3);
            Screen('DrawDots',w,[0;0],p.ppd*p.lineApt_size,p.bg_color,stim_pos_pix(8,:),3);

            Screen('DrawLines',w,line_values,line_width,p.fix_color-40); % should make this into a unique color variable...
        end
        
        Screen('Flip',w);
        t.frame_start(tt,ff) = GetSecs;
        if ff == 1
            t.search_start(tt) = t.frame_start(tt,ff);
            
            if p.eye_tracking == 1
                Eyelink('Message','xDAT %i',3);
            end
        end


%         if p.conditions(tt,3) == 1
%             if ff == 1 && p.conditions(tt,7) == 1
%                 myImageSearchNone = Screen('GetImage',w);
% 
%             elseif ff == 1 && p.conditions(tt,7) == 2
%                 myImageSearchMismatch = Screen('GetImage',w);
% 
%             elseif ff == 1 && p.conditions(tt,7) == 3
%                 myImageSearchMatch = Screen('GetImage',w);
% 
%             end
%         end
        
        % check for ACTUAL responses
        [resp, timeStamp] = checkForResp(p.keys, p.escape);
        if resp==-1
            sca; ShowCursor;
            if p.eye_tracking==1
                Eyelink('ShutDown');
            end
            return;
        end
        
        if resp && find(p.keys==resp) && isnan(p.resp(tt))
            p.resp(tt) = find(p.keys==resp); % (FIRST response is that recorded)
            p.RT(tt) = timeStamp-t.search_start(tt);
        end
        
    end
    
    % aperture
    Screen('FillRect',w,0);
    Screen('FillOval',w,p.bg_color,aperture_rect);
    Screen('FillOval',w,p.fix_color,outline_rect);
    Screen('FillOval',w,p.bg_color,rim_rect);
    
    % draw fixation
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color,center_point,3);
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color, center_point,3);
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color,center_point,3);
    
    Screen('Flip',w);
    
    % ITI
    % fixation point & aperture
    Screen('FillRect',w,0);
    Screen('FillOval',w,p.bg_color,aperture_rect);
    Screen('FillOval',w,p.fix_color,outline_rect);
    Screen('FillOval',w,p.bg_color,rim_rect);
    
    % draw fixation (dimmed/big) (TODO: correctly compute line width...)
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color,center_point,3);
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color, center_point,3);
    Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color,center_point,3);
    
    Screen('Flip',w);

    
    if p.eye_tracking == 1
        Eyelink('Message','xDAT %i',4);
    end
    
    while GetSecs < (t.trial_start(tt) + p.mem_dur + p.postmem_dur + p.search_dur + p.ITIs(tt))
        [resp, timeStamp] = checkForResp(p.keys, p.escape);
        if resp == -1
            sca;
            ShowCursor;
            return;
        end
    end
    
    t.trial_end(tt) = GetSecs;
    
end
 
end_check  = GetSecs;  

Screen('FillRect',w,0);
Screen('FillOval',w,p.bg_color,aperture_rect);
Screen('FillOval',w,p.fix_color,outline_rect);
Screen('FillOval',w,p.bg_color,rim_rect);

% draw fixation (dimmed/big) (TODO: correctly compute line width...)
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color,center_point,3);
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color, center_point,3);
Screen('DrawDots',w,[0;0],p .ppd*p.fix_size_inner*2, p.fix_color,center_point,3);

Screen('Flip',w);

if p.eye_tracking == 1
    
    Eyelink('Message','xDAT %i',0);
    
    Eyelink('command', 'record_status_message "Waiting for scanner to stop"');
    
end

% end-experiment wait
resp = 0;
waitstart = GetSecs;

while GetSecs < (t.trial_start(tt) + p.mem_dur + p.postmem_dur + p.search_dur + p.ITIs(tt)) + p.end_wait
    [resp, ~] = checkForResp([], p.escape);
    if resp == -1
        sca;
        ShowCursor;
        if p.eye_tracking==1
            Eyelink('ShutDown');
        end
        return;
    end
end

t.expt_end = GetSecs; 

for tt = 1:p.ntrials
    p.correct(tt) = p.resp(tt) == p.conditions(tt,8); 
end 

p.feedback = mean(p.correct); 
% quick save
save(p.fn,'p','t');

% provide behavior feedback
txt = sprintf('Search Accuracy: %.02f',p.feedback); % should include missed responses eventually

% draw aperture
Screen('FillRect',w,0);
Screen('FillOval',w,p.bg_color,aperture_rect);
Screen('FillOval',w,p.fix_color,outline_rect);
Screen('FillOval',w,p.bg_color,rim_rect);

% draw fixation (dimmed/big) (TODO: correctly compute line width...)
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2,       p.fix_color*0.8,center_point,3);
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size*2*0.9,   p.bg_color,    center_point,3);
Screen('DrawDots',w,[0;0],p.ppd*p.fix_size_inner*2, p.fix_color*0.8,center_point,3);
 
DrawFormattedText(w,txt,'center',center_point(2)-3*p.ppd,p.fix_color);
 
Screen('Flip',w);

if p.eye_tracking == 1
    Eyelink('StopRecording');
    Eyelink('ReceiveFile',[p.et_fn '.edf'],[p.et_fn '.edf']);

      p.eyedatafile_renamed = [p.fn(1:(end-3)) 'edf'];
    movefile([p.et_fn '.edf'],p.eyedatafile_renamed);

    Eyelink('ShutDown');
end

resp = 0;
fprintf('waiting for final space\n');
while resp == 0
    [resp, ~] = checkForResp(KbName('space'),  p.escape);
    if resp == -1
        sca;
        fprintf('ESC pressesd during feedback time\n');
        ShowCursor;
        return;
    end
end
clear resp ;

% quick RT and ACC analyses 
% nanmean(p.RT(p.conditions(:,7) == 1));
% nanmean(p.RT(p.conditions(:,7) == 2));
% nanmean(p.RT(p.conditions(:,7) == 3));
% nanmean(p.RT(p.conditions(:,7) == 4));

% save one final time
save(p.fn,'p','t');

Screen('CloseAll');
ShowCursor;
