function wmDistort_pilot1(subj,run)
Screen('Preference', 'SkipSyncTests', 1);

try
p.expt_name = 'wmDistort_pilot1';

p.do_et = 0;
p.subj = subj;
p.run = run;

% data directory 2 above the current directory
p.filename = sprintf('../../data/wmDistort_data/%s_r%02.f_%s_%s.mat',p.subj,p.run,p.expt_name,datestr(now,30));
if p.do_et == 1
    p.eyedatafile = sprintf('%s_r%02.f',p.subj(1:min(length(p.subj),3)),p.run);
end

p.rng_seed = cputime*1000;
rng(p.rng_seed);0

% ------ size of relevant stim features, etc ------ %
p.wm_ecc = 8;     % deg [in behavioral room, max ecc of circle 15 deg]
p.cue_size = 0.55; % deg
p.wm_size = 0.65;  % deg, size of WM dots
p.sep_ang = 30;    % deg polar angle, maximum stim separation distance
p.aperture_size = 15; % [or max ecc?]

p.fix_size_in  = 0.075; % radius, deg
p.fix_size_out = 0.30; % radius, deg
p.fix_pen = 1.5;
% now we draw out (frameoval), cue (if end of trial), in (filloval) for fix

% ------- keyboard stuff --------------------------- %
KbName('UnifyKeyNames');
if ismac == 1
    p.esc_key = KbName('escape'); % press this key to abort
else
    p.esc_key = KbName('escape');
end
%p.start_key = [KbName('5%') KbName('5')];  % should be lower-case t at prisma? (or %5, which is top-row 5, or 5, which is numpad 5)
p.space = KbName('space');


% ------ color of relevant stim features ----------- %
p.bg_color  = 20*[1 1 1];
p.fix_color = 75*[1 1 1];%[150 150 150];        % during trial/delay/etc
% p.wm_colors = [200 0   0;         % red
%                  0 0 255          % blue
%                  180 0 180;       % purple
%                  130 130 0];      % yellow
p.wm_colors = 75*[1 1 1];
                 

p.choose_color = 130*[1 1 1];%[255 255 255]; % when subj should choose, color of fix



% ------ conditions ------ %
p.r_cond = 1; % 1: R1
p.repetitions = 64; 
p.ntrials = length(p.r_cond)*p.repetitions;

p.conditions = nan(p.ntrials,1);
cnt = 1;
for cc = 1:length(p.r_cond)
    for rr = 1:p.repetitions
        p.conditions(cnt) = p.r_cond(cc);
        cnt = cnt+1;
    end
end

p.rnd_idx = randperm(p.ntrials);
p.conditions = p.conditions(p.rnd_idx,:);


% ------ timing of trial events --------- %
p.black_dur = 0;
p.blank_dur = 2;
p.targ_dur = 0.5;
p.delay_dur = 2.5;
p.cue_dur = 0.8; % a bit longer than usual
p.feedback_dur = 0.8; 
p.iti_range = [2 4]; % randomly choose between those

p.itis = linspace(p.iti_range(1),p.iti_range(2),p.ntrials);
p.itis = p.itis(randperm(p.ntrials));
% p.itis = p.itis(randperm(32));



% ------- things to save -------- %
p.targ_coords = cell(2,1);
p.targ_coords{1} = nan(p.ntrials,2);
p.targ_coords{2} = nan(p.ntrials,2);

p.targ_colors = cell(2,1);
p.targ_colors{1} = nan(p.ntrials,3);
p.targ_colors{2} = nan(p.ntrials,3);

p.targ_angs = nan(p.ntrials,2);


% ------- Screen setup, optics --------- %

p.screen_height = 30; % cm, in the experiment room
p.viewing_distance = 56; % cm, in the experiment room (inside lab)

% open a screen, to get the resolution
s = max(Screen('Screens'));
HideCursor;
[w, p.scr_rect] = Screen('OpenWindow',s,[0 0 0]);
% [w, p.scr_rect] = Screen('OpenWindow',max(Screen('Screens')),[0 0 0]); HideCursor;
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);

p.ifi = Screen('GetFlipInterval',w);

p.center = p.scr_rect([3 4])/2;

p.ppd = p.scr_rect(4)/(2*atan2d(p.screen_height/2,p.viewing_distance));


p.aperture_rect = CenterRectOnPoint([0 0 2 2]*p.ppd*p.aperture_size,p.center(1),p.center(2));
p.fix_rect_out  = CenterRectOnPoint([0 0 2 2] * p.ppd  * p.fix_size_out,p.center(1),p.center(2));
p.fix_rect_in   = CenterRectOnPoint([0 0 2 2] * p.ppd  * p.fix_size_in, p.center(1),p.center(2));

p.screen_width = p.scr_rect(3); % Screen width in pixels
p.screen_height = p.scr_rect(4); % Screen height in pixels
p.longest_side = p.screen_height; % The height of the monitor as the longest side
p.diameter_length = (p.screen_height + p.screen_height*0.618)/2

% Set aperture size based on the screen height
aperture_diameter = p.diameter_length; % For the circle
aperture_wide_rect = [0 0 p.screen_height * 0.618 p.screen_height];
aperture_tall_rect = [0 0 p.screen_height p.screen_height * 0.618];
aperture_square = [0 0 p.diameter_length p.diameter_length];


% Center the apertures
p.aperture_circle = CenterRectOnPoint([0 0 aperture_diameter aperture_diameter], p.center(1), p.center(2));
p.aperture_wide_rect = CenterRectOnPoint(aperture_wide_rect, p.center(1), p.center(2));
p.aperture_tall_rect = CenterRectOnPoint(aperture_tall_rect, p.center(1), p.center(2));
p.aperture_square = CenterRectOnPoint(aperture_square, p.center(1), p.center(2));



% --------- eyetracking ----------- %
if p.do_et == 1
    
    el=EyelinkInitDefaults(w);
    
    el.backgroundcolour=p.bg_color(1);  % TODO: fix this?
    el.calibrationtargetcolour=p.fix_color(1);
    
    el.msgfontcolour=p.fix_color(1);
    p.foregroundcolour=p.fix_color(1);
    
    EyelinkUpdateDefaults(el);

    
    Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
   
    Eyelink('command','calibration_type=HV13'); % updating number of callibration dots
    s=Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA');% (,GAZERES,HREF,PUPIL,STATUS,INPUT');
    s=Eyelink('command', 'sample_rate = 1000');
    s=Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    s=Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, p.scr_rect(3)-1,p.scr_rect(4)-1);
    
    
    
    % make sure that we get gaze data from the Eyelink
    % if a widescreen display:
    if p.scr_rect(3)==2560 && p.scr_rect(4)==1440
        %         Eyelink('command', 'generate_default_targets = NO');
        %
        %         % 13 samples
        %         Eyelink('command','calibration_samples = 13');
        %         Eyelink('command','validation_samples = 13');
        %
        %         % set up two random orders:
        %         calib_str = sprintf('%d,', randperm(13)-1);
        %         valid_str = sprintf('%d,', randperm(13)-1);
        %
        %         Eyelink('command',sprintf('calibration_sequence = %s',calib_str(1:end-1)));
        %
        %
        %         Eyelink('command',sprintf('validation_sequence = %s',valid_str(1:end-1)));
        
        Eyelink('command', 'calibration_area_proportion 0.59 0.83');
        %Eyelink('command', 'calibration_area_proportion 0.59 0.83');
        
    end
    
    
    %------ calibrate the eye tracker --------
    EyelinkDoTrackerSetup(el);
    if s~=0
        error('link_sample_data error, status: ',s)
    end
    Eyelink('openfile',p.eyedatafile);

end

shapes = {'circle', 'wide_rect', 'tall_rect', 'square'};
% Repeat the shape list enough times to cover all trials
shape_list = repmat(shapes, 1, p.ntrials / length(shapes));
% Shuffle the shape order for randomness
randomized_shapes = shape_list(randperm(p.ntrials));

% Setting the number of conditions and degrees of increment and the
% starting wm
p.num_shapes = length(shapes);
p.num_conditions = p.ntrials/4;
p.incre = 22.5;
p.start_deg = 11.25;

% Generate degrees
degrees = p.start_deg + (0:(p.num_conditions-1)) * p.incre; % 11.25, 33.75, ..., 348.75

% Create a list of unique conditions (shape, degree)
conditions = cell(p.num_shapes * p.num_conditions, 2); % Preallocate cell array for conditions
index = 1;
for shape_index = 1:p.num_shapes
    for degree_index = 1:p.num_conditions
        conditions{index, 1} = shapes{shape_index};  % Shape
        conditions{index, 2} = degrees(degree_index); % Degree
        index = index + 1;
    end
end

% Shuffle the conditions for randomness
% randomized_conditions = conditions(randperm(size(conditions, 1)), :);

% START OF EXPERIMENT

draw_aperture = @() Screen('FillOval',w,p.bg_color,p.aperture_circle);

Screen('FillRect',w,[0 0 0]);
% draw_aperture();

% TODO: instructions/greeting
% instructions
txt = 'Remember dot position precisely';
Screen('TextSize', w, 30);
DrawFormattedText(w,txt,'center',p.center(2)-4*p.ppd,p.fix_color);
%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
%Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
% Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2);
% Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);
% Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2);
%Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
Screen('Flip',w);


% check for esc, space.... 

resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.space, p.esc_key);
    if resp == -1
        Screen('CloseAll'); ShowCursor;
        Eyelink('ShutDown'); 
        return;
    end
end
clear resp;
p.start_expt = GetSecs;


% blank screen
Screen('FillRect',w,[0 0 0]);
% draw_aperture();
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);

Screen('Flip',w);

if p.do_et == 1
    Eyelink('Message','xDAT %i', 10);
    Eyelink('StartRecording'); % make 1 big edf file (save time)
end


WaitSecs(1.5); % wait a bit before first trial

% ------ generate file based on run number ------- %
filename = '../../data/wmDistort_data/temp_randomized_conditions.mat';  % Define filename

if mod(p.run, 2) == 1
    % If run number is odd, generate a new .mat file with 1 to 32 trials
    p.ntrials = 32;
    temp_randomized_conditions = conditions(randperm(size(conditions, 1)), :);
    save(filename, 'temp_randomized_conditions');
    randomized_conditions = temp_randomized_conditions(1:32, :);  % Use trials 33-64
else
    % If run number is even, load the file and use trials 33 to 64
    if isfile(filename)
        load(filename, 'temp_randomized_conditions');
        p.ntrials = 32;
        randomized_conditions = temp_randomized_conditions(33:64, :);  % Use trials 33-64
    else
        error('No existing .mat file found for this even run.');
    end
end

for tt = 1:p.ntrials
    disp("-----------------------------------------------------------")
    disp(tt)
    % disp(randomized_shapes(tt))
    % current_shape = [randomized_shapes{tt}];
    current_shape = randomized_conditions{tt, 1}; % Get the current shape
    current_degree = randomized_conditions{tt, 2}; % Get the current degree
    disp(['Shape: ', current_shape, ', Degree: ', num2str(current_degree)]);
    
    % this trial's position(s)
    this_ang = nan(1,2);
    this_ang(1) = current_degree;
    % this_ang(1) = 360*rand(1);
    disp(this_ang)

    
    p.targ_coords{1}(tt,:) = p.wm_ecc * [cosd(this_ang(1)) sind(this_ang(1))];
    
    if p.conditions(tt,1) ~= 1
        this_ang(2) = this_ang(1)+p.sep_ang + (360-2*p.sep_ang)*rand(1); % compute second position if necessary
        p.targ_coords{2}(tt,:) = p.wm_ecc * [cosd(this_ang(2)) sind(this_ang(2))];
    end
    
    p.targ_angs(tt,:) = this_ang; % save these for convenience
    
    
    % this trial's colors
    tmp_color_idx = randperm(size(p.wm_colors,1));
    
    p.targ_colors{1}(tt,:) = p.wm_colors(tmp_color_idx(1),:);
    
    if p.conditions(tt,1) ~= 1
        p.targ_colors{2}(tt,:) = p.wm_colors(tmp_color_idx(2),:);
    end
    
    
    
    % targets (XDAT 1) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    trial_start = GetSecs;
    
    
    if p.do_et == 1
        %Eyelink('Message','TarX %s', num2str(p.targ_coords{1}(tt,1)));
        %Eyelink('Message','TarY %s', num2str(p.targ_coords{1}(tt,2)));
        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
        
        Eyelink('Message','xDAT %i',1);
        
        Eyelink('command', 'record_status_message "TRIAL %d of %d"', tt, p.ntrials);
        
    end

    while GetSecs < trial_start + p.black_dur
        Screen('FillRect',w,[0 0 0]);
        Screen('DrawDots',w,[0;0], p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0], p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0], p.fix_size_in*p.ppd*2, p.fix_color,p.center,2); 
        Screen('Flip',w);
    end
    

    while GetSecs < trial_start + p.black_dur + p.blank_dur

        % aperture
        Screen('FillRect',w,[0 0 0]);
        % draw_aperture();
        switch current_shape
            case {'circle'}
                Screen('FillOval', w, p.bg_color, p.aperture_circle);
            case 'wide_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
            case 'tall_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
            case 'square'
                Screen('FillRect', w, p.bg_color, p.aperture_square);
            otherwise
                error('Unknown shape: %s', shape);
        end
        Screen('DrawDots',w,[0;0], p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0], p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0], p.fix_size_in*p.ppd*2, p.fix_color,p.center,2); 
        Screen('Flip',w);
    end

    while GetSecs < trial_start + p.black_dur + p.blank_dur + p.targ_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        %draw_aperture();
        switch current_shape
            case {'circle'}
                Screen('FillOval', w, p.bg_color, p.aperture_circle);
            case 'wide_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
            case 'tall_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
            case 'square'
                Screen('FillRect', w, p.bg_color, p.aperture_square);
            otherwise
                error('Unknown shape: %s', shape);
        end
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.targ_coords{1}(tt,:).', p.wm_size*p.ppd, p.targ_colors{1}(tt,:), p.center, 2);
        
        if p.conditions(tt,1)~=1
            
            % if necessary, target 2
            Screen('DrawDots',w,p.ppd*[1;-1].*p.targ_coords{2}(tt,:).', p.wm_size*p.ppd, p.targ_colors{2}(tt,:), p.center, 2);
        end
        
        % fixation
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
        
    end
    
    % delay (XDAT 2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    if p.do_et == 1
        Eyelink('Message','xDAT %i',2);    
    end
    
    
    while GetSecs < trial_start + p.black_dur + p.blank_dur + p.targ_dur + p.delay_dur
        % aperture
        Screen('FillRect',w,[0 0 0]);
        %draw_aperture();
        switch current_shape
            case {'circle'}
                Screen('FillOval', w, p.bg_color, p.aperture_circle);
            case 'wide_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
            case 'tall_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
            case 'square'
                Screen('FillRect', w, p.bg_color, p.aperture_square);
            otherwise
                error('Unknown shape: %s', shape);
        end
        
%         if p.conditions(tt,1) ~= 3
%             Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.choose_color,p.center,2);
%         else
%             Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.targ_colors{1}(tt,:),p.center,2);
%         end
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
        
    end
    
    % go cue (XDAT 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','xDAT %i',3);
    end
    
    while GetSecs < trial_start + p.black_dur + p.blank_dur + p.targ_dur + p.delay_dur + p.cue_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        %draw_aperture();
        switch current_shape
            case {'circle'}
                Screen('FillOval', w, p.bg_color, p.aperture_circle);
            case 'wide_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
            case 'tall_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
            case 'square'
                Screen('FillRect', w, p.bg_color, p.aperture_square);
            otherwise
                error('Unknown shape: %s', shape);
        end
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        if p.conditions(tt,1) == 3
            Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.choose_color,p.center,2);
        else
            Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.targ_colors{1}(tt,:),p.center,2);
        end
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
    end
    
    % feedback (XDAT 4, tarx, tary) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(p.targ_coords{1}(tt,1)));
        Eyelink('Message','TarY %s', num2str(p.targ_coords{1}(tt,2)));
        % NOTE: incorrect on 1/3 trials
        Eyelink('Message','xDAT %i',4);
        
    end
    
    while GetSecs < trial_start + p.black_dur + p.blank_dur + p.targ_dur + p.delay_dur + p.cue_dur + p.feedback_dur
        
        % aperture
        Screen('FillRect',w,[0 0 0]);
        %draw_aperture();
        switch current_shape
            case {'circle'}
                Screen('FillOval', w, p.bg_color, p.aperture_circle);
            case 'wide_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
            case 'tall_rect'
                Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
            case 'square'
                Screen('FillRect', w, p.bg_color, p.aperture_square);
            otherwise
                error('Unknown shape: %s', shape);
        end
        
        % target 1
        Screen('DrawDots',w,p.ppd*[1;-1].*p.targ_coords{1}(tt,:).', p.wm_size*p.ppd, p.targ_colors{1}(tt,:), p.center, 2);
        
        if p.conditions(tt,1)~=1
            
            % if necessary, target 2
            Screen('DrawDots',w,p.ppd*[1;-1].*p.targ_coords{2}(tt,:).', p.wm_size*p.ppd, p.targ_colors{2}(tt,:), p.center, 2);
        end
        
        % fixation
%        Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        
         %       Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        if p.conditions(tt,1) == 3
            Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.choose_color,p.center,2);
        else
            Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.targ_colors{1}(tt,:),p.center,2);
        end
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');
            save(p.filename,'p');
            return;
        end
 
    end
    
    % ITI (XDAT 5) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
        Eyelink('Message','xDAT %i',5);
    end

    while GetSecs < trial_start + p.black_dur + p.blank_dur + p.targ_dur + p.delay_dur + p.cue_dur + p.feedback_dur + p.itis(tt)
        
        % aperture
        % Screen('FillRect',w,[0 0 0]);
        %draw_aperture();
        % switch current_shape
        %     case {'circle'}
        %         Screen('FillOval', w, p.bg_color, p.aperture_circle);
        %     case 'wide_rect'
        %         Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
        %     case 'tall_rect'
        %         Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
        %     case 'square'
        %         Screen('FillRect', w, p.bg_color, p.aperture_square);
        %     otherwise
        %         error('Unknown shape: %s', shape);
        % end
        
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            Eyelink('StopRecording');
            Eyelink('ShutDown');

            save(p.filename,'p');
            return;
        end
        
    end
    

    disp("end of one tt")
    % save [note: in scanner, do this at beginning of ITI after first flip]
    save(p.filename,'p');
    
end

% END OF EXPERIMENT - TEXT
if p.do_et == 1
    Eyelink('Message','xDAT %i',11);
end

Screen('FillRect',w,[0 0 0]);
draw_aperture();
txt = sprintf('End of run %i',p.run);
DrawFormattedText(w,txt,'center',p.center(2)-4*p.ppd,p.fix_color);
%Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);

Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

Screen('Flip',w);

resp = 0;
while resp == 0
    [resp, timeStamp] = checkForResp(p.space, p.esc_key);
end
clear resp;

if p.do_et == 1
    Eyelink('StopRecording');
    Eyelink('ReceiveFile',[p.eyedatafile '.edf'],[p.eyedatafile '.edf']);
    
    p.eyedatafile_renamed = [p.filename(1:(end-3)) 'edf'];
    movefile([p.eyedatafile '.edf'],p.eyedatafile_renamed);
    
    Eyelink('ShutDown');
end

Screen('CloseAll');
ShowCursor;
catch
    Screen('CloseAll');
    ShowCursor;
    le = lasterror;
    rethrow(le);
end


return