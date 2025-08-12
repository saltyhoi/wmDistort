function wmDistort_pilot3()
Screen('Preference', 'SkipSyncTests', 1);

% Create input dialog for subject number, run number, and eye tracker condition
prompt = {'Enter Subject Number (e.g., sub001):                                                            .', 
    'Enter Run Number (e.g., 1):', 
    'Eye Tracker Condition (0 = Off, 1 = On):',
    'Response Mode (0 = Eye, 1 = Mouse):'};
dlg_title = 'Subject & Run Info Input';
num_lines = 1; 
default_input = {'sub001', '1', '0', '1'};  % Default: no eye tracker, mouse response


% Get user input via dialog box
user_input = inputdlg(prompt, dlg_title, num_lines, default_input);

% Extract the inputs from the user input dialog
subj = strtrim(user_input{1});       % Subject number
run  = str2double(strtrim(user_input{2}));  % Run number  , convert to double
et   = str2double(strtrim(user_input{3}));   % Eye tracker condition (0 or 1)

% Check if eye tracker condition is valid
if et ~= 0 && et ~= 1
    error('Invalid input for Eye Tracker condition. Must be 0 or 1.');
end

try
p.expt_name = 'wmDistort_pilot3';
p.subj = subj;
p.run = run; 
p.do_et = et;
p.response_mode = str2double(user_input{4});
 
% Display the experiment settings for confirmation
disp(['Running experiment for  Subject: ', p.subj, ', Run: ', num2str(p.run)]);
disp(['Eye Tracker Condition: ', num2str(p.do_et)]);  % Shows 0 or 1

% data directory 2 above the current directory
p.filename = sprintf('../../data/wmDistort_data/%s_r%02.f_%s_%s.mat',p.subj,p.run,p.expt_name,datestr(now,30));
if p.do_et == 1
    p.eyedatafile = sprintf('%s_r%02.f',p.subj(1:min(length(p.subj),3)),p.run);
end

p.rng_seed = cputime*1000;
rng(p.rng_seed);

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
p.fix_color = 75*[1 1 1];
p.wm_colors = 75*[1 1 1];                
p.choose_color = 130*[1 1 1];


% ------ conditions: aperture setting & angles setting ------ %
p.shapes = {'circle', 'square', 'wide_rect', 'tall_rect'};
p.repetitions = 16; % number of target angles per each aperture condition
p.r_cond = length(p.shapes); % numer of aperture condition
p.overall_ntrials = p.r_cond * p.repetitions;
p.ntrials = p.overall_ntrials/2; % run 1 full cycle into 2 runs

% ------- mouse stuff --------------------------- %

p.mouse_response_dur = 3; % seconds, or however long you want
p.mouse_rt = nan(p.ntrials, 1);
p.this_mouse_traj = cell(p.ntrials, 1);

% ------- mouse stuff --------------------------- %

% Generate degrees
p.incre = 360/p.repetitions; % the angle should be 22.5 degrees here
p.start_deg = 11.25;
p.degrees = p.start_deg + (0:(p.repetitions-1)) * p.incre; % 11.25, 33.75, ..., 348.75

% ------ generate condition file ------- %
filename = sprintf('../../data/wmDistort_data/fixed_conditions.mat');  % Define filename for the fixed file

% Check if the file already exists
if ~exist(filename, 'file')
    conditions = cell(p.overall_ntrials, 3);
    index = 1;
    for shape_index = 1:p.r_cond
        for degree_index = 1:p.repetitions
            conditions{index, 1} = p.shapes{shape_index};  % Shape
            conditions{index, 2} = p.degrees(degree_index); % Degree
            conditions{index, 3} = fix((index-1)/p.repetitions) + 1; % Numeric values used to distinguish aperture conditions 1:circle, 2:square, 3:wide_rect, 4:tall_rect
            index = index + 1;
        end
    end

    % Save the file
    save(filename, 'conditions');
    fprintf('Generated and saved fixed_conditions.mat\n');
else
    fprintf('File fixed_conditions.mat already exists. No action taken.\n');
end


% ------ generate random permutation or load from previous .mat file ------- %
if mod(p.run, 2) == 1   % Odd-numbered run, use trials 1 to 32
    load(filename, 'conditions');
    random_permutation = randperm(size(conditions, 1));  % Generate random permutation
    randomized_conditions = conditions(random_permutation(1:p.ntrials), :);  % Use first half (1:32)
    % Save the random permutation into the run-specific file
    p.random_permutation = random_permutation;
    save(p.filename, 'p');  % Save random_permutation along with p
else                    % Even-numbered run, use trials 33 to 64
    % Generate the base part of the filename without the timestamp
    prev_filename_pattern = sprintf('../../data/wmDistort_data/%s_r%02.f_%s_*.mat', p.subj, p.run-1, p.expt_name);
    
    % Find all matching files with the same base name (ignore the timestamp)
    files = dir(prev_filename_pattern);
    if ~isempty(files)
        % Assume the most recent file is the correct one (based on timestamp)
        [~, idx] = max([files.datenum]);
        prev_filename = fullfile(files(idx).folder, files(idx).name);
        loaded_data = load(prev_filename, 'p');
        random_permutation = loaded_data.p.random_permutation;  % Retrieve the stored permutation
        load(filename, 'conditions');  % Ensure conditions are loaded
        randomized_conditions = conditions(random_permutation(p.ntrials+1:end), :);  % Use second half (33:64)
        p.random_permutation = random_permutation;  % Save this permutation for the current run
        save(p.filename, 'p');  % Save the `p` struct and the permutation for this run
        fprintf('Second run complete, using second half of the permutation.\n');
    else
        error('Previous run file not found. Make sure you run an odd-numbered run first.');
    end
end

% finally define conditions
p.conditions = cell2mat(randomized_conditions(:, 3));



% ------ timing of trial events --------- %
p.aperture_dur = 2;
p.targ_dur = 0.5;
p.delay_dur = 2.5;
p.cue_dur = 3; % a bit longer than usual 0.8
p.feedback_dur = 0.8; 
p.iti_range = [2 4]; % randomly choose between those

p.itis = linspace(p.iti_range(1),p.iti_range(2),p.ntrials);
p.itis = p.itis(randperm(p.ntrials));

% ------ timing of trial events --------- %
% Simulate by setting durations to a small number, e.g., 0.01 or just 0
% p.aperture_dur = 0.01;  % small duration
% p.targ_dur = 0.01;      % small duration
% p.delay_dur = 0.01;     % small duration
% p.cue_dur = 0.01;       % small duration
% p.feedback_dur = 0.01;  % small duration
% p.iti_range = [0.01 0.01]; % small ITI range for simulation
% 
% % Randomized ITI for simulation
% p.itis = linspace(p.iti_range(1), p.iti_range(2), p.ntrials);
% p.itis = p.itis(randperm(p.ntrials));


% ------- things to save -------- %
p.targ_coords = cell(1,1);
p.targ_coords{1} = nan(p.ntrials,2);
p.targ_colors = cell(1,1);
p.targ_colors{1} = nan(p.ntrials,3);
p.targ_angs = nan(p.ntrials,1);
p.targ_apertures = cell(p.ntrials,1);


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
p.diameter_length = (p.screen_height + p.screen_height*0.618)/2;

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


% START OF EXPERIMENT

draw_aperture = @() Screen('FillOval',w,p.bg_color,p.aperture_circle);

Screen('FillRect',w,[0 0 0]);
% draw_aperture();

% TODO: instructions/greeting
% instructions
txt = 'Remember dot position precisely';
Screen('TextSize', w, 30);
DrawFormattedText(w,txt,'center',p.center(2)-4*p.ppd,p.fix_color);
% Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2);
Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
Screen('Flip',w);


% check for esc, space.... 

resp = 0;
while resp == 0
    [resp, ~] = checkForResp(p.space, p.esc_key);
    if resp == -1
        Screen('CloseAll'); ShowCursor;
        if p.do_et == 1
            Eyelink('ShutDown'); 
        end
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

for tt = 1:p.ntrials
    %disp("-----------------------------------------------------------")
    %disp(tt)
    current_shape = randomized_conditions{tt, 1}; % Get the current shape
    current_degree = randomized_conditions{tt, 2}; % Get the current degree
    %disp(['Shape: ', current_shape, ', Degree: ', num2str(current_degree)]);
    

    % this trial's position(s)
    this_ang = nan(1);
    this_ang(1) = current_degree;
    this_shape = strings(1); % Create a 1x1 string array
    this_shape(1) = string(current_shape); % Convert character array to string
    %disp(this_ang)

    
    p.targ_coords{1}(tt,:) = p.wm_ecc * [cosd(this_ang(1)) sind(this_ang(1))];
    p.targ_angs(tt) = this_ang; % save these for convenience
    p.targ_apertures(tt) = {this_shape}; % save this aperture shape
    
    % this trial's colors
    tmp_color_idx = randperm(size(p.wm_colors,1));
    
    p.targ_colors{1}(tt,:) = p.wm_colors(tmp_color_idx(1),:);
    
    
    % trial starts/aperture display (XDAT 1) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    trial_start = GetSecs;
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        Eyelink('Message','xDAT %i',1);
        Eyelink('command', 'record_status_message "TRIAL %d of %d"', tt, p.ntrials);
    end
    
    Screen('DrawDots',w,[0;0], p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0], p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
    Screen('DrawDots',w,[0;0], p.fix_size_in*p.ppd*2, p.fix_color,p.center,2); 
    Screen('Flip',w);


    while GetSecs < trial_start + p.aperture_dur

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

    % target  (XDAT 2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if p.do_et == 1
        Eyelink('Message','xDAT %i',2);
    end


    while GetSecs < trial_start + p.aperture_dur + p.targ_dur
        
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

        
        % fixation
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    % delay (XDAT 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    if p.do_et == 1
        Eyelink('Message','xDAT %i',3);    
    end
    
    
    while GetSecs < trial_start + p.aperture_dur + p.targ_dur + p.delay_dur
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
        
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;  
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end

    % Mouse response period (XDAT 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.response_mode == 1
    
        if p.do_et == 1
            Eyelink('Message','xDAT %i',4);
        end

        SetMouse(p.center(1), p.center(2), w);  % move mouse to center
        ShowCursor(1);                          % show mouse cursor

        this_frame = 1;  % initialize frame counter
        max_frames = round(p.mouse_response_dur * 120); % assuming 120 Hz
        p.this_mouse_traj{tt} = nan(max_frames, 2);     % preallocate trajectory

        no_resp = 1;     % response flag
        tStart = GetSecs;
    
        while GetSecs < trial_start + p.aperture_dur + p.targ_dur + p.delay_dur + p.cue_dur && no_resp
    
            while (GetSecs < tStart + p.mouse_response_dur) && no_resp
    
                % fill aperture
                Screen('FillRect', w, [0 0 0]);
                switch current_shape
                    case 'circle'
                        Screen('FillOval', w, p.bg_color, p.aperture_circle);
                    case 'wide_rect'
                        Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
                    case 'tall_rect'
                        Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
                    case 'square'
                        Screen('FillRect', w, p.bg_color, p.aperture_square);
                    otherwise
                        error('Unknown shape: %s', current_shape);
                end
    
                % draw fixation
    
                Screen('DrawDots', w, [0;0], p.fix_size_out*p.ppd*2 + p.fix_pen, p.fix_color, p.center, 2);
                Screen('DrawDots', w, [0;0], p.fix_size_out*p.ppd*2 - p.fix_pen, p.bg_color, p.center, 2);
                % Screen('DrawDots', w, [0;0], p.cue_size*p.ppd, p.targ_colors{1}(tt,:), p.center, 2); % center grey dot
                Screen('DrawDots', w, [0;0], p.fix_size_in*p.ppd*2, p.fix_color, p.center, 2);
    
                Screen('Flip', w);

                % % Screenshotting the cursor -----------------%
                % [x, y] = GetMouse(w);
                % % Draw your own cursor (e.g., a small circle)
                % Screen('DrawDots', w, [x; y], 10, [255 0 0], [], 2);
                % 
                % imageArray = Screen('GetImage', w);
                % imwrite(imageArray, 'screenshot_with_cursor.png'); 
                % % Screenshotting the cursor -----------------%


                % get mouse position
                [x, y, buttons] = GetMouse(w);
                if this_frame <= max_frames
                    % p.this_mouse_traj{tt}(this_frame, :) = [x, y] / p.ppd;        % Alison's code save the format this way
                    p.this_mouse_traj{tt}(this_frame, :) = ([1, -1] .* ([x, y] - p.center)) / p.ppd;
                    this_frame = this_frame + 1;
                end
    
                % check for click
                if any(buttons)
                    no_resp = 0;
                    tEnd = GetSecs;
                    p.mouse_rt(tt) = tEnd - tStart;
                    mouse_pos_deg = ([1, -1] .* ([x, y] - p.center)) / p.ppd;
                    p.mouse_click_pos(tt, :) = mouse_pos_deg;
                   % p.mouse_click_pos(tt, :) = [x, y] / p.ppd;                     % Alison's code save the format this way
                
                    % === Draw feedback frame ===
                    Screen('FillRect', w, [0 0 0]);
                    switch current_shape
                        case 'circle'
                            Screen('FillOval', w, p.bg_color, p.aperture_circle);
                        case 'wide_rect'
                            Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
                        case 'tall_rect'
                            Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
                        case 'square'
                            Screen('FillRect', w, p.bg_color, p.aperture_square);
                    end
                
                    % Draw central fixation and cue
                    Screen('DrawDots', w, [0;0], p.fix_size_out*p.ppd*2 + p.fix_pen, p.fix_color, p.center, 2);
                    Screen('DrawDots', w, [0;0], p.fix_size_out*p.ppd*2 - p.fix_pen, p.bg_color, p.center, 2);
                    Screen('DrawDots', w, [0;0], p.cue_size*p.ppd, p.targ_colors{1}(tt,:), p.center, 2);
                    Screen('DrawDots', w, [0;0], p.fix_size_in*p.ppd*2, p.fix_color, p.center, 2);
                
                    % === Feedback marker: red dot where user clicked ===
                    Screen('DrawDots', w, [x; y], 10, [255 0 0], [], 2);  % red dot
                
                    Screen('Flip', w);  % Show feedback frame
                    WaitSecs(0.5);      % Brief pause to let participant see feedback
                end
    
                % optional escape check
                [resp] = checkForResp([], p.esc_key);
                if resp == -1
                    Screen('CloseAll'); ShowCursor;
                    if p.do_et == 1
                        Eyelink('StopRecording');
                        Eyelink('ShutDown');
                    end
                    save(p.filename,'p');
                    return;
                end
    
            end % end of inner while mouse response
    
            HideCursor;
    
        end % end of outer while timing loop
    
    % go cue (XDAT 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
    
        if p.do_et == 1
            Eyelink('Message','xDAT %i',4);
        end
    
        while GetSecs < trial_start + p.aperture_dur + p.targ_dur + p.delay_dur + p.cue_dur
    
            % aperture
            Screen('FillRect', w, [0 0 0]);
            switch current_shape
                case 'circle'
                    Screen('FillOval', w, p.bg_color, p.aperture_circle);
                case 'wide_rect'
                    Screen('FillRect', w, p.bg_color, p.aperture_wide_rect);
                case 'tall_rect'
                    Screen('FillRect', w, p.bg_color, p.aperture_tall_rect);
                case 'square'
                    Screen('FillRect', w, p.bg_color, p.aperture_square);
                otherwise
                    error('Unknown shape: %s', current_shape);
            end
    
            Screen('DrawDots', w, [0;0], p.fix_size_out*p.ppd*2 + p.fix_pen, p.fix_color, p.center, 2);
            Screen('DrawDots', w, [0;0], p.fix_size_out*p.ppd*2 - p.fix_pen, p.bg_color, p.center, 2);
            Screen('DrawDots', w, [0;0], p.cue_size*p.ppd, p.targ_colors{1}(tt,:), p.center, 2);
            Screen('DrawDots', w, [0;0], p.fix_size_in*p.ppd*2, p.fix_color, p.center, 2);
    
            Screen('Flip', w);
    
            % check for escape key
            [resp] = checkForResp([], p.esc_key);
            if resp == -1
                Screen('CloseAll'); ShowCursor;
                if p.do_et == 1
                    Eyelink('StopRecording');
                    Eyelink('ShutDown');
                end
                save(p.filename, 'p');
                return;
            end
    
        end % end while go cue timing loop
    
    end % end if p.response_mode == 1

    % feedback (XDAT 5, tarx, tary) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(p.targ_coords{1}(tt,1)));
        Eyelink('Message','TarY %s', num2str(p.targ_coords{1}(tt,2)));
        % NOTE: incorrect on 1/3 trials
        Eyelink('Message','xDAT %i',5);
        
    end
    % Determine feedback start time based on actual RT (adjust cue duration)
    if ~isempty(p.mouse_rt(tt)) && p.mouse_rt(tt) < p.cue_dur
        % Participant responded early, so shift feedback earlier
        feedback_start = trial_start + p.aperture_dur + p.targ_dur + p.delay_dur + p.mouse_rt(tt);
    else
        % No early response, use full cue duration
        feedback_start = trial_start + p.aperture_dur + p.targ_dur + p.delay_dur + p.cue_dur;
    end


    while GetSecs < feedback_start + p.feedback_dur 
        
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
        
        
        % fixation
%        Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 

        
        %Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        
         %       Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 


        Screen('DrawDots',w,[0;0],p.cue_size*p.ppd,p.targ_colors{1}(tt,:),p.center,2);
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 

        
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
 
    end
    
    % ITI (XDAT 6) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if p.do_et == 1
        Eyelink('Message','TarX %s', num2str(0));
        Eyelink('Message','TarY %s', num2str(0));
        
        Eyelink('Message','xDAT %i',6);
    end

    while GetSecs < feedback_start + p.feedback_dur + p.itis(tt)
        Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2+p.fix_pen,p.fix_color,p.center,2); Screen('DrawDots',w,[0;0],p.fix_size_out*p.ppd*2-p.fix_pen,p.bg_color,p.center,2); 
        Screen('DrawDots',w,[0;0],p.fix_size_in*p.ppd*2,p.fix_color,p.center,2); 
        
        %Screen('DrawDots',w,[0;0],p.fix_size*p.ppd,p.fix_color,p.center,2);
        
        Screen('Flip',w);
        
        % check for esc.... 
        [resp] = checkForResp([], p.esc_key); % TODO: maybe turn on/off gaze indicator?
        if resp == -1
            Screen('CloseAll'); ShowCursor;
            if p.do_et == 1
                Eyelink('StopRecording');
                Eyelink('ShutDown');
            end
            save(p.filename,'p');
            return;
        end
        
    end
    
    % disp("end of one tt")
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
% catch
%     Screen('CloseAll');
%     ShowCursor;
%     le = lasterror; 
%     rethrow(le); 
    
catch ME
    sca;  
    ShowCursor;  
    rethrow(ME); 
end



return