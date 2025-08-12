function wmDistort_pilot4_search()
Screen('Preference', 'SkipSyncTests', 1);

p.expt_name = 'wmDistort_pilot4';
p.subj = 123;
p.run = 1; 
p.do_et = 0;

% data directory 2 above the current directory
p.filename = sprintf('../../data/wmDistort_data/%s_r%02.f_%s_%s.mat',p.subj,p.run,p.expt_name,datestr(now,30));
if p.do_et == 1
    p.eyedatafile = sprintf('%s_r%02.f',p.subj(1:min(length(p.subj),3)),p.run);
end

p.rng_seed = cputime*1000;
rng(p.rng_seed);

KbName('UnifyKeyNames');

% Basic experiment params
p.ntrials = 12;
searchSetSize = 90; % We want the RT to be around 3-10 seconds
p.shapes = {'circle', 'square', 'wide_rect', 'tall_rect'};
shapeMap = containers.Map({'circle', 'square', 'tall_rect', 'wide_rect'}, [1, 2, 3, 4]);
s = max(Screen('Screens'));
[w, p.scr_rect] = Screen('OpenWindow',s,[0 0 0]); 

p.fix_size_in  = 0.075; % radius, deg
p.fix_size_out = 0.30; % radius, deg
p.fix_pen = 1.5;
p.fix_color = [130 130 130]; % red fixation color, change if you want
p.bg_color = [20 20 20]; % background color matching your screen background



% ------- Screen setup, optics --------- %
p.screen_height = 30; % cm, in the experiment room
p.viewing_distance = 56; % cm, in the experiment room (inside lab)
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);
p.ifi = Screen('GetFlipInterval',w);
p.center = p.scr_rect([3 4])/2;
p.ppd = p.scr_rect(4)/(2*atan2d(p.screen_height/2,p.viewing_distance));
p.letter_size_px = round(0.74 * p.ppd);
p.halfLetter = round(p.letter_size_px / 2); %variable required to make buffer distance from aperture bounary




% % ------- things to save -------- %
p.targ_coords = cell(1,1);
p.targ_coords{1} = nan(p.ntrials,2);
p.targ_angs = nan(p.ntrials,1);
p.targ_rot = nan(p.ntrials,1);
p.targ_rot_dir = cell(p.ntrials, 1);
p.targ_apertures = cell(p.ntrials,1);
p.conditions = cell(p.ntrials, 1);
p.shapes_used = cell(p.ntrials, 1);
p.clickPos = nan(p.ntrials, 2);     % mouse click X, Y
p.RT = nan(p.ntrials, 1);           % response time
p.correct = nan(p.ntrials, 2);  % 2 columns: [correctLoc, correctRot] % correct (1) or incorrect (0)


% ------ timing of trial events --------- %
p.fixation_dur = 2;
p.aperture_dur = 2;
p.task_dur = 0.5;
p.response_dur = 1.5;
p.feedback_dur = 0.8; 
p.iti_range = [2 4]; % randomly choose between those
p.itis = linspace(p.iti_range(1),p.iti_range(2),p.ntrials);
p.itis = p.itis(randperm(p.ntrials));

% Setup screen & parameters
Screen('Preference', 'SkipSyncTests', 1);
[win, p.scr_rect] = Screen('OpenWindow', max(Screen('Screens')), [20 20 20]);
HideCursor();
p.center = p.scr_rect([3 4]) / 2;
p.screen_height = p.scr_rect(4);
p.screen_width = p.scr_rect(3);

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
        Eyelink('command', 'calibration_area_proportion 0.59 0.83');        
    end
    
    
    %------ calibrate the eye tracker --------
    EyelinkDoTrackerSetup(el);
    if s~=0
        error('link_sample_data error, status: ',s)
    end
    Eyelink('openfile',p.eyedatafile);

end

% Show instruction screen
showInstructions(win, p);
% Aperture sizes from your info
p.diameter_length = (p.screen_height + p.screen_height*0.618)/2;
aperture_diameter = p.diameter_length;
aperture_tall_rect = [0 0 p.screen_height * 0.618 p.screen_height]; 
aperture_wide_rect = [0 0 p.screen_height p.screen_height * 0.618]; 
aperture_square = [0 0 p.diameter_length p.diameter_length];

% Centered apertures
p.aperture_circle = CenterRectOnPoint([0 0 aperture_diameter aperture_diameter], p.center(1), p.center(2));
p.aperture_wide_rect = CenterRectOnPoint(aperture_wide_rect, p.center(1), p.center(2));
p.aperture_tall_rect = CenterRectOnPoint(aperture_tall_rect, p.center(1), p.center(2));
p.aperture_square = CenterRectOnPoint(aperture_square, p.center(1), p.center(2));


% % Data to store
% results = struct('trial', {}, 'shape', {}, 'clickPos', {}, 'RT', {}, 'correct', {});

Screen('Flip',win);

if p.do_et == 1
    Eyelink('Message','xDAT %i', 10);
    Eyelink('StartRecording'); % make 1 big edf file (save time)
end

WaitSecs(1.5); % wait a bit before first trial

% % --------- fMRI code implementing (timing of trial events) --------- % %

% % WAIT TO START SCANNER
% resp = 0;
% while resp == 0
%     [resp, ~] = checkForResp(p.start, p.escape);
%     if resp == -1
%         Screen('CloseAll'); ShowCursor;
%         if p.eye_tracking==1
%             Eyelink('ShutDown'); 
%         end
%         return;
%     end
% end
% clear resp;
% t.expt_start = GetSecs;

try
    for t = 1:p.ntrials
        % Randomize shape for this trial
        shapeName = p.shapes{randi(numel(p.shapes))};

        % Get aperture rect for shape (returns [left top right bottom])
        apertureRect = getApertureRect(shapeName, p);

        % Create aperture mask (logical matrix, size of screen)
        apertureMask = createApertureMask(shapeName, p);

        % Store aperture rectangle in cell (as vector)
        p.targ_apertures{t} = apertureRect;
        % Store condition code as integer from shapeMap
        p.conditions{t} = shapeMap(shapeName);

        % Fixation-only stage (XDAT 1) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if p.do_et == 1
            Eyelink('Message','TarX %s', num2str(0));
            Eyelink('Message','TarY %s', num2str(0));
            Eyelink('Message','xDAT %i',1);
            Eyelink('command', 'record_status_message "TRIAL %d of %d"', tt, p.ntrials);
        end

        % Show fixation cross
        drawFixation(win, p);
        Screen('Flip', win);
        timedWaitWithEscape(p.fixation_dur, p.do_et);

        % Aperture-Display stage (XDAT 2) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if p.do_et == 1
            Eyelink('Message','xDAT %i',2);
        end

        % Show aperture + fixation stage (pause for ~1 second)
        showApertureWithFixation(win, shapeName, apertureRect, p);
        timedWaitWithEscape(p.aperture_dur, p.do_et);

        % Clear screen with background color
        Screen('FillRect', win, p.bg_color);
        % Draw the aperture filled shape
        drawApertureFilled(win, shapeName, apertureRect);
        % Now draw search items inside aperture
        % [itemList, itemPositions, itemAngles] = generateSearchItems(searchSetSize, apertureMask, p.center);
        gridSize = 50;
        candidatePositions = generateCandidatePositions(apertureMask, p, gridSize);
        [itemList, itemPositions, itemAngles] = generateSearchItemsGrid(searchSetSize, candidatePositions, p.center);


        targetIdx = find(itemList == 'L');  % find index of target 'L'
        targetPos = itemPositions(:, targetIdx);  % << ADD THIS LINE

        if t > 1
            if prevRT < 2  % RT from previous trial
                Freeviewing = GetSecs;
                flickerSearchItems(win, itemList, itemPositions, itemAngles, p.do_et, p.bg_color, shapeName, apertureRect, p.letter_size_px);
            else
                Freeviewing = GetSecs;
                drawSearchItems(win, itemList, itemPositions, itemAngles, p.letter_size_px);
            end
        else
            % First trial — default to static display
            Freeviewing = GetSecs;
            drawSearchItems(win, itemList, itemPositions, itemAngles, p.letter_size_px);
        end

        % ----------------------------- IMPORTANT ------------------------------------------------
        % drawSearchItems(win, itemList, itemPositions, itemAngles); 
        % ----------------------------- IMPORTANT ------------------------------------------------


        % Save target coordinates (x, y) as row vector
        p.targ_coords{1}(t, :) = itemPositions(:, targetIdx)';
        % Save target angle with sign (negative = left rotation, positive = right)
        p.targ_rot(t) = itemAngles(targetIdx);

        if p.targ_rot(t) < 0
            p.targ_rot_dir{t} = 'l';  % left rotation
        else
            p.targ_rot_dir{t} = 'r';  % right rotation (including zero or positive)
        end

        % Visual_Search-Task stage (XDAT 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if p.do_et == 1
            Eyelink('Message','xDAT %i',3);
        end
        % Flip to show stimuli and record flip time
        vbl = Screen('Flip', win);


        % Response stage (XDAT 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ XDAT inside the function 
        [clickPos, RT, correctLoc, correctRot] = collectClickResponse(apertureMask, targetPos, p.center, win, p.response_dur, p.do_et, p.targ_rot_dir{t}, Freeviewing);
        fprintf('Trial %d — Loc: %d | Rot: %d | Overall: %d\n', t, correctLoc, correctRot, correctLoc && correctRot);

        % Feedback stage (XDAT 5) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % Provide feedback for location
        if correctLoc 
            feedbackLoc = 'You correctly clicked the target!';
        else
            feedbackLoc = 'You missed the target.';
        end
        
        % Provide feedback for rotation
        if correctRot
            feedbackRot = 'You correctly identified the rotation!';
        else
            feedbackRot = 'You identified the rotation incorrectly.';
        end
        
        % Combine feedback messages (can be displayed together or separately)
        feedbackText = sprintf('%s\n\n\n%s', feedbackLoc, feedbackRot);

        DrawFormattedText(win, feedbackText, 'center', 'center', [255 255 255]);

        if p.do_et == 1
            Eyelink('Message','xDAT %i',5);
            Eyelink('Message','TarX %s', num2str(p.targ_coords{1}(t,1)));
            Eyelink('Message','TarY %s', num2str(p.targ_coords{1}(t,2)));
        end

        Screen('Flip', win);
        timedWaitWithEscape(p.feedback_dur, p.do_et);

        % Store trial results
        p.clickPos(t, :) = clickPos;   % e.g., [x y]
        p.RT(t) = RT;
        p.correct(t, 1) = correctLoc;
        p.correct(t, 2) = correctRot;
        p.shapes_used{t} = shapeName;

        % Inter-trial interval stage (XDAT 6) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if p.do_et == 1
            Eyelink('Message','xDAT %i',6);
        end
        timedWaitWithEscape(p.itis(t), p.do_et);
        disp(RT)
        prevRT = RT;
    end

    % End of experiment message
    DrawFormattedText(win, 'All done!\n\nThank you!', 'center', 'center', [255 255 255]);
    Screen('Flip', win);
    % KbWait;
    while true
        checkForEscape(p.do_et);  % allows termination during end screen
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            break;
        end
        WaitSecs(0.01);
    end

    save(p.filename, 'p');

catch ME
    % Save data on error or escape
    save(p.filename, 'p');
    Screen('CloseAll');
    ShowCursor;
    if ~strcmp(ME.message, 'Experiment terminated by user (ESC key).')
        rethrow(ME);
    end
    % rethrow(ME);
end

% END OF EXPERIMENT - TEXT
if p.do_et == 1
    Eyelink('Message','xDAT %i',11);
end

Screen('CloseAll');
ShowCursor;
end

function rect = getApertureRect(shapeName, p)
    switch shapeName
        case 'circle'
            rect = p.aperture_circle;
        case 'square'
            rect = p.aperture_square;
        case 'wide_rect'
            rect = p.aperture_wide_rect;
        case 'tall_rect'
            rect = p.aperture_tall_rect;
    end
end

function apertureMask = createApertureMask(shapeName, p)
    apertureMask = false(p.scr_rect(4), p.scr_rect(3));
    [xGrid, yGrid] = meshgrid(1:p.scr_rect(3), 1:p.scr_rect(4));
    cx = p.center(1);
    cy = p.center(2);
    
    switch shapeName
        case 'circle'
            radius = p.diameter_length / 2;
            apertureMask = ((xGrid - cx).^2 + (yGrid - cy).^2) <= radius^2;
        case 'square'
            half_side = p.diameter_length / 2;
            apertureMask = abs(xGrid - cx) <= half_side & abs(yGrid - cy) <= half_side;
        case 'wide_rect'
            half_width = p.screen_height / 2;
            half_height = (p.screen_height * 0.618) / 2;
            apertureMask = abs(xGrid - cx) <= half_width & abs(yGrid - cy) <= half_height;
        case 'tall_rect'
            half_width = (p.screen_height * 0.618) / 2;
            half_height = p.screen_height / 2;
            apertureMask = abs(xGrid - cx) <= half_width & abs(yGrid - cy) <= half_height;
        otherwise
            error('Unknown shape');
    end
end

function drawFixation(w, p)
    % Draw a multi-layer fixation point as described
    
    % Outer dot (thicker, fix_color)
    Screen('DrawDots', w, [0;0], p.fix_size_out * p.ppd * 2 + p.fix_pen, p.fix_color, p.center, 2);
    % Inner dot (background color, slightly smaller)
    Screen('DrawDots', w, [0;0], p.fix_size_out * p.ppd * 2 - p.fix_pen, p.bg_color, p.center, 2);
    % Innermost dot (fix_color, smaller)
    Screen('DrawDots', w, [0;0], p.fix_size_in * p.ppd * 2, p.fix_color, p.center, 2);
end


function drawApertureFilled(win, shapeName, apertureRect)
    grey = [130 130 130];
    switch shapeName
        case 'circle'
            Screen('FillOval', win, grey, apertureRect);  % Proper circle
        case {'square', 'wide_rect', 'tall_rect'}
            Screen('FillRect', win, grey, apertureRect);  % Rectangular shapes
    end
end


function [items, positions, angles] = generateSearchItemsGrid(nItems, candidatePositions, center)
    % candidatePositions: Nx2 array of [x,y] candidate points
    
    % Select target position randomly (could add near/far logic here)
    targetIdx = randi(size(candidatePositions,1));
    targetPos = candidatePositions(targetIdx, :)';
    
    % Remove target from candidates to avoid overlap
    remainingCandidates = candidatePositions;
    remainingCandidates(targetIdx, :) = [];
    
    % Randomly select distractors from remaining candidates
    distractorIdx = randsample(size(remainingCandidates,1), nItems - 1);
    distractorPos = remainingCandidates(distractorIdx, :)';
    
    % Combine positions
    positions = [targetPos, distractorPos];
    
    % Create items vector: target 'L' + distractors 'T'
    items = repmat('T', 1, nItems);
    items(1) = 'L';
    
    % Generate angles as before (excluding -5 to +5)
    angles1 = randi([-90, -6], 1, ceil(nItems/2));
    angles2 = randi([6, 90], 1, floor(nItems/2));
    angles = [angles1, angles2];
    angles = angles(randperm(length(angles)));
    
    % Shuffle the items and positions together
    shuffleOrder = randperm(nItems);
    positions = positions(:, shuffleOrder);
    items = items(shuffleOrder);
    angles = angles(shuffleOrder);
end


function drawSearchItems(win, items, positions, angles, letter_size_px)
    Screen('TextFont', win, 'Arial');
    Screen('TextSize', win, letter_size_px);
    col = [0 0 0];
    
    for i = 1:numel(items)
        % Set color based on item type
        if items(i) == 'L'
            col = [255 0 0]; % Red color for T (target)
        else
            col = [0 0 0];   % Black color for L (distractors)
        end
    
        Screen('glPushMatrix', win);
        Screen('glTranslate', win, positions(1,i), positions(2,i), 0);
        Screen('glRotate', win, angles(i), 0, 0, 1);
        Screen('DrawText', win, items(i), -8, -12, col); % adjust offset as needed
        Screen('glPopMatrix', win);
    end
end

function [clickPos, RT, correctLoc, correctRot] = collectClickResponse(apertureMask, targetPos, center, win, responseDur, do_et, targetRot, Freeviewing)
    % Initialize outputs with default values
    clickPos = [NaN, NaN];
    RT = NaN;
    correctLoc = false;
    correctRot = false;

    % ----------- PHASE 1: ROTATION RESPONSE (Free Viewing) -----------
    disp('Waiting for participant to indicate they found the target and its rotation...');
    
    responseRot = 'x';  % Default invalid
    while true
        [~, ~, buttons] = GetMouse(win);
        checkForEscape(do_et);

        if any(buttons)
            if buttons(1)
                responseRot = 'l';  % Left click
            elseif buttons(3)
                responseRot = 'r';  % Right click
            else
                responseRot = 'x';
            end
            break;  % First response registered
        end

        WaitSecs(0.01);
    end

    % ✅ Record rotation accuracy
    correctRot = strcmp(responseRot, targetRot);
    fprintf('Rotation Response: %s | Correct: %d\n', responseRot, correctRot);

    % Wait until button is released
    while true
        [~, ~, buttons] = GetMouse(win);
        if ~any(buttons)
            break;
        end
        WaitSecs(0.01);
    end

    % ----------- PHASE 2: LOCATION RESPONSE (Timed) -----------
    if do_et == 1
        Eyelink('Message','xDAT %i',4);
    end

    ShowCursor(win);
    SetMouse(center(1), center(2), win);

    disp('Cursor shown. Please click on the target location...');
    startTime = GetSecs;

    while (GetSecs - startTime) < responseDur
        [x, y, buttons] = GetMouse(win);
        checkForEscape(do_et);

        if any(buttons)
            clickPos = [x, y];
            % RT = GetSecs - startTime; % measuring from click to click
            RT = GetSecs - Freeviewing;

            % Location accuracy
            dist = sqrt((x - targetPos(1))^2 + (y - targetPos(2))^2);
            correctLoc = dist <= 30;

            fprintf('Location dist: %.2f | Click: (%d,%d)\n', dist, x, y);
            break;
        end

        WaitSecs(0.01);
    end

    HideCursor(win);
end

function inside = isInsideAperture(x, y, apertureMask)
    x = round(x);
    y = round(y);
    x = max(1, min(size(apertureMask, 2), x));
    y = max(1, min(size(apertureMask, 1), y));
    inside = apertureMask(y, x);
end

function fits = isLetterFullyInside(x, y, apertureMask, halfLetter)
    % Check all 4 corners of the letter bounding box
    fits = true;
    
    % Define corners relative to center (x,y)
    corners = [
        x - halfLetter, y - halfLetter;  % top-left
        x + halfLetter, y - halfLetter;  % top-right
        x - halfLetter, y + halfLetter;  % bottom-left
        x + halfLetter, y + halfLetter;  % bottom-right
    ];
    
    % Check each corner inside aperture
    for k = 1:size(corners,1)
        if ~isInsideAperture(corners(k,1), corners(k,2), apertureMask)
            fits = false;
            return
        end
    end
end


function showInstructions(win, p)
    % Draw fixation with your custom style
    Screen('DrawDots', win, [0;0], p.fix_size_out*p.ppd*2 + p.fix_pen, p.fix_color, p.center, 2);
    Screen('DrawDots', win, [0;0], p.fix_size_out*p.ppd*2 - p.fix_pen, p.bg_color, p.center, 2);
    Screen('DrawDots', win, [0;0], p.fix_size_in*p.ppd*2, p.fix_color, p.center, 2);
    
    % Draw instructions text
    instructions = ['Try to find the letter "L" among "T"s.\n\n' ...
                    'Click on the "L" as quickly and accurately as possible.\n\n' ...
                    'Press SPACE to start.'];
    DrawFormattedText(win, instructions, 'center', p.center(2) + 100, [255 255 255]);
    
    % Show it on screen
    Screen('Flip', win);
    
    % Wait for space key press
    while true
        checkForEscape(p.do_et);
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown && keyCode(KbName('space'))
            break; % Exit loop to start trials
        end
        WaitSecs(0.01);
    end
end

function showApertureWithFixation(win, shapeName, apertureRect, p)
    % Draw gray background and aperture shape
    Screen('FillRect', win, p.bg_color);
    drawApertureFilled(win, shapeName, apertureRect);
    % Draw fixation point on top
    % Outer dot (thicker, fix_color)
    % Screen('DrawDots', win, [0;0], p.fix_size_out * p.ppd * 2 + p.fix_pen, p.bg_color, p.center, 2);
    Screen('DrawDots', win, [0;0], p.fix_size_out * p.ppd * 2 + p.fix_pen, [50 50 50], p.center, 2);
    % Inner dot (background color, slightly smaller)
    Screen('DrawDots', win, [0;0], p.fix_size_out * p.ppd * 2 - p.fix_pen, p.fix_color, p.center, 2);
    % Innermost dot (fix_color, smaller)
    % Screen('DrawDots', win, [0;0], p.fix_size_in * p.ppd * 2, p.bg_color, p.center, 2);    % Flip to screen
    Screen('DrawDots', win, [0;0], p.fix_size_in * p.ppd * 2, [50 50 50], p.center, 2);    % Flip to screen
    Screen('Flip', win);
end

function checkForEscape(do_et)
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown && keyCode(KbName('ESCAPE'))
        Screen('CloseAll');
        ShowCursor;
        if do_et == 1
            Eyelink('ShutDown'); 
        end
        error('Experiment terminated by user (ESC key).');
    end
end

function timedWaitWithEscape(waitTime, do_et)
    startTime = GetSecs;
    while GetSecs - startTime < waitTime
        checkForEscape(do_et)
        WaitSecs(0.01);
    end
end


function candidatePositions = generateCandidatePositions(apertureMask, p, gridSize)
    % apertureMask: logical mask of the aperture (height x width)
    % p: struct with screen parameters including center
    % gridSize: number of grid cells along each dimension (e.g., 50)
    
    [height, width] = size(apertureMask);
    
    % Calculate grid cell size
    cellWidth = width / gridSize;
    cellHeight = height / gridSize;
    
    candidatePositions = [];
    
    for i = 1:gridSize
        for j = 1:gridSize
            % Calculate center pixel coordinate of each grid cell
            xCenter = round((j - 0.5) * cellWidth);
            yCenter = round((i - 0.5) * cellHeight);
            
            % Check boundary conditions (in case rounding goes outside)
            xCenter = min(max(xCenter, 1), width);
            yCenter = min(max(yCenter, 1), height);

            if isLetterFullyInside(xCenter, yCenter, apertureMask, p.halfLetter)
                candidatePositions = [candidatePositions; xCenter, yCenter];
            end

            % Check if this center point is inside the aperture mask
            % if apertureMask(yCenter, xCenter)
            %     candidatePositions = [candidatePositions; xCenter, yCenter];
            % end
        end
    end
end

function responseRot = flickerSearchItems(win, items, positions, angles, do_et, bg_color, shapeName, apertureRect, letter_size_px)
    ifi = Screen('GetFlipInterval', win);
    nItems = numel(items);

    flickerFreqs = 0.5 + 1.5 * rand(1, nItems);
    toggleIntervals = round(1 ./ (2 .* flickerFreqs) ./ ifi);

    visible = true(1, nItems);
    framesSinceToggle = zeros(1, nItems);

    responseRot = '';  % Initialize output

    vbl = Screen('Flip', win);

    while true
        for i = 1:nItems
            framesSinceToggle(i) = framesSinceToggle(i) + 1;
            if framesSinceToggle(i) >= toggleIntervals(i)
                visible(i) = ~visible(i);
                framesSinceToggle(i) = 0;
            end
        end
        Screen('FillRect', win, bg_color);
        drawApertureFilled(win, shapeName, apertureRect);

        % Draw all visible items at once
        drawSearchItems(win, items(visible), positions(:, visible), angles(visible), letter_size_px);

        vbl = Screen('Flip', win, vbl + 0.5 * ifi);

        checkForEscape(do_et);

        [~, ~, buttons] = GetMouse(win);
        if buttons(1)
            responseRot = 'l';
            break;
        elseif buttons(3)
            responseRot = 'r';
            break;
        end

        WaitSecs(0.001);
    end

    Screen('FillRect', win, bg_color);
    drawApertureFilled(win, shapeName, apertureRect);
    drawSearchItems(win, items, positions, angles, letter_size_px);

end