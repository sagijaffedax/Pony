% Random Dot Kinectmatograms, RDK
% Naiqi G. Xiao
% Princeton U
% March 2017
% nx@princeton.edu

% In this program gaze will start directional motion
% looking at the box will trigger directional motion
% once the motion starts, it won't change until .5 seconds after audio
% offset

%% Update history

% April 3rd,
% 1. Added a .25 s pause between the soundplay and motion change.
% This will remove the changes in direction of dots random movements, which could be a cue of
% directional motion.
% 2. trial order is not randomized based on the direction of motion. Check
% ParametersRDK.m
% 3. the vertical location of AttentionalGetter varies in each trial.

% April 15th,
% 1. Added a learning test trials to examine whether particiants formed an association
% 2. Saved the main RDK codes into separate m files, making the script
% cleaner and simpler

% April 17th,
% 1. Integrated eye tracking codes.

% April 18th,
% 1. Shorter session. Keep the duration of the first induced motion and
% reduce the length of the following motion within each familiarization
% block.
% 2. Enlarge the target objct by half.

% May 4th,
% 1. Added support for Tobii eye tracker
% 2. Changed Calibration to 3 points
% 3. Reduced the cumulative looking time criterion to trigger the motion in the Eyelink setting

% May 26th,
% 1. Add a movie as an attention getter to attract infants attention before
% 2. starting calibration.
% 3. Press "SPACE" key to end the movie. For Tobii, any key pressing will end movie. 

% June 8th,
% 1. Decrease the propability of ommision trials to 25%;
% 2. Both learning and test blocks contain 4 trials;
% 3. Only 60% Coherencr is used;
% 4. Add 4 blocks of expectation violated blocks. Each contains 1 out of 4
% motion direction reversed trials;
% 5. 10 Test blocks for each direction.

% June 19th,
% dots will stop moving when the object appears
% when the dots stops, the brightness of the dots decreased to highlight
% the object.

% June 20th,
% Change the attentional getter movement to vertical moving at the center
% of the screen.
%%
sca;
close all;
clearvars;

Screen('Preference', 'SkipSyncTests', 1)

[T, G] = ParametersRDK();

% [T, G] = ParametersRDKEXPT3();

KbName('UnifyKeyNames');
SpaceKey = KbName('SPACE');
DCKey = KbName('RightShift'); %use left shift to force recalibration

if ~exist('Results', 'dir')
    
    mkdir('Results')
    
end

prompt = 'Participant Numbr:';
Participant = input(prompt, 's');

symbols = ['a':'z' 'A':'Z' '0':'9'];
MAX_ST_LENGTH = 50;
stLength = randi(MAX_ST_LENGTH);
nums = randi(numel(symbols),[1 3]);
st = symbols(nums);
    
Participant = [Participant, '_', st];

try
    
    Eyetracking = G.Eyetracking; % whether use eye tracker
    
    nframes = G.nframes; % number of animation frames in loop
    mon_width = G.mon_width;  % horizontal dimension of viewable screen (cm)
    v_dist = G.v_dist;   % viewing distance (cm)
    
    dot_speed = G.dot_speed;    % dot speed (deg/sec)
    f_kill = G.f_kill; % fraction of dots to kill each frame (limited lifetime)
    
    ndots = G.ndots;  % number of dots
    dot_w = G.dot_w;  % width of dot (deg)
    
    ShowGaze = G.ShowGaze;
    
    differentcolors = G.differentcolors; % Use a different color for each point if == 1. Use common color white if == 0.
    differentsizes = G.differentsizes; % Use different sizes for each point if >= 1. Use one common size if == 0.
    waitframes = G.waitframes;     % Show new dot-images at each waitframes'th monitor refresh.
    
    screens=Screen('Screens');
    screenNumber=max(screens);
    [w, rect] = Screen('OpenWindow', screenNumber, [0 0 0]);
    
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [center(1), center(2)] = RectCenter(rect);
    
    fps=Screen('FrameRate',w);      % frames per second
    ifi=Screen('GetFlipInterval', w);
    
    if fps == 0
        
        fps = 1/ifi;
        
    end
    
    white = WhiteIndex(w);
    %     HideCursor; % Hide the mouse cursor
    Priority(MaxPriority(w));
    
    % Do initial flip...
    vbl = Screen('Flip', w);
    
    % ---------------------------------------
    % initialize dot positions and velocities
    % ---------------------------------------
    
    ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
    pfs = dot_speed * ppd / fps;                            % dot speed (pixels/frame)
    %     pfs = pfs * (rand(ndots,1) / 2 + .5);
    pfs = repmat(pfs, ndots, 1);
    s = dot_w * ppd;                                        % dot size (pixels)
    
    X = rect(3) * rand(ndots,1);
    Y = rect(4) * rand(ndots,1);
    
    xy = [X Y];   % dot positions in Cartesian coordinates (pixels from center)
    
    t = 2 * pi * rand(ndots,1);
    %     index = randperm(size(t, 1));
    %     index = index(:,1: round(size(t, 1) * Coherence));
    %     if length(index) > 0
    %         t(index) = angle;
    %     end
    cs = [cos(t), sin(t)];
    
    dxdy = [pfs pfs] .* cs;                    % change in x and y per frame (pixels)
    
    % Create a vector with different colors for each single dot, if
    % requested:
    if (differentcolors==1)
        colvect = uint8(round(rand(3,ndots)*255));
    else
        colvect=white;
    end
    
    % Create a vector with different point sizes for each single dot, if
    % requested:
    if (differentsizes > 0)
        s = (0.4 + rand(1, ndots) * 0.6) .* s;
    end
    
    % Clamp point sizes to range supported by graphics hardware:
    [minsmooth, maxsmooth] = Screen('DrawDots', w);
    s = min(max(s, minsmooth), maxsmooth);
    
    xymatrix = transpose(xy);
    
    % --------------
    % define the box area on the screen
    % --------------
    baseRect = [0 0 200 200];
    centeredRect = CenterRectOnPointd(baseRect, center(1), center(2));
    
    %% Audio preparation
    %---------------
    % Sound Setup
    %---------------
    
    % Initialize Sounddriver
    InitializePsychSound(0);
    
    devices = PsychPortAudio('GetDevices');
    
    psychlasterror('reset');
    
    % Start immediately (0 = immediately)
    startCue = 0;
    
    % Should we wait for the device to really start (1 = yes)
    % INFO: See help PsychPortAudio
    waitForDeviceStart = 1;
    
    % --------------
    % audio files
    % --------------
    cd('~/Documents/MATLAB/Topdown_Preception/AudioStimuli/Rising and Falling Pitch')
    
    SoundList = dir('*.wav');
    
    AudioStimuli = cell(size(SoundList, 1), 3);
    
    for i = 1:size(SoundList, 1)
        
        filename = SoundList(i).name;
        
        [y, freq] = psychwavread(filename);
        wavedata = y';
        
        nrchannels = size(wavedata, 1);
        
        AudioStimuli{i, 1} = wavedata;
        AudioStimuli{i, 2} = freq;
        AudioStimuli{i, 3} = size(wavedata, 1);
        AudioStimuli{i, 4} = filename;
        
    end
    
    cd('~/Documents/MATLAB/Topdown_Preception')
    
    Screen('TextSize', w, 96);
    
    %% Eye Tracker
    
    switch Eyetracking
        
        case 1 % Eye Link
            
            EyeLinkSetup;
            
            eye_used = -1;
            
        case 2 % Tobii
            
            TobiiSetup;
            
    end
    
    %%
    
    % Whether the motion is triggered by gaze or experimenter pressing space bar.
    GazeTrigger = zeros(sum(T.NofDisplay), 5);
    GazeTrigger(:, 4) = 1:sum(T.NofDisplay);
    % GazeTrigger(:, 5) = [reshape(T.Direction(1:8, :), [], 1); reshape(T.Direction(9:30, 1), [], 1)];
    
    %% Camera Recording
    
    if G.Recording
        
        vidobj = videoinput('macvideo', 1);
        
        set(vidobj, 'FramesPerTrigger', 4000)
        
        set(vidobj, 'FrameGrabInterval', 3)
        
        set(vidobj, 'ReturnedColorSpace', 'grayscale')
        
        set(vidobj, 'TriggerRepeat', 1)
        
        triggerconfig(vidobj, 'manual')
        
        vidobj.TriggerFrameDelay = 0;
        
        FramesFetch = 400;
        
        ExportStartsAt = 1;
        
        if ~exist('Recordings', 'dir')
            
            mkdir('Recordings')
            
        end
        
        start(vidobj)
        
        trigger(vidobj)
        tic
        start = clock;
        
        WaitSecs(.5);
        [a, b, c] = getdata(vidobj, 1);
        
        RecordTimingCorrection = c.AbsTime(6) - start(6);
       
    else
        
        tic;
        
    end
    
    if ~exist('Participant', 'var')
        
        Participant = '9999'; % test purpose
        
    end
    
    %% Experiment Starts
    for Trial = 1:size(T, 1) %looping the trials
        
        disp('***************************');
        
        disp(sprintf('Trial # %d', Trial));
        
        ii = sum(T.NofDisplay(1:(Trial - 1)));
        ExportStartsAt = ii + 1;
        
        if mod(Trial, 4) == 1
            
            DrawFormattedText(w, '+', 'center', 'center', white);
            Screen('Flip', w);
            
            StartTime = GetSecs;
            [keyIsDown, keyCode, secs] = KbPressWait([], StartTime + .5);
            
        end
        %% Drive correction if necessary
        % Eyelink only
        switch Eyetracking
            
            case 1
                % STEP 7.2
                % Do a drift correction at the beginning of each trial
                % Performing drift correction (checking) is optional for
                % EyeLink 1000 eye trackers.
                
                if keyCode(DCKey)
                    
                    EyelinkDoDriftCorrection(el);
                    
                end
                
        end
        
        %%
        Audio = T.SoundIndex(Trial);
        Shape = T.VisualShape(Trial, :);
        ShapeColor = T.ColorIndex(Trial, :);
        angle = T.Direction(Trial, :);
        angle = angle * pi;
        Coherence = T.Coherence(Trial, :);
        NofDisplay = T.NofDisplay(Trial);
        LearningTest = T.LearningTest(Trial);
        TrialDuration = T.TrialDuration(Trial, :);
        
        ISI = 0;
        BrightnessIndex = 1;
        
        pahandle = PsychPortAudio('Open', 1, [], 0, AudioStimuli{Audio, 2}, AudioStimuli{Audio, 3});
        PsychPortAudio('FillBuffer', pahandle, AudioStimuli{Audio, 1});
        
        AudioDuration = size(AudioStimuli{Audio, 1}, 2) ./ AudioStimuli{Audio, 2};
        
        GazeCount = 0;
        PreviousGaze = 0;
        Switch = 0;
        DirectionInit = 1;
        AudioIniti = 1;
        RandomIniti = 0;
        DirectionalCount = 0;
        
        AGcenter(1) = center(1);
        AGcenter(2) = center(2) + (rand() - .5) * center(2) / 2;
        
        baseRect(4) = rect(4)/4;
        centeredRect = CenterRectOnPointd(baseRect, AGcenter(1), AGcenter(2));
        
        [mx, my, buttons] = GetMouse(screenNumber);
        
        if any(buttons)
            
            break;
            
        end
        
        %% Eye Tracker
        
        switch Eyetracking
            
            case 1 % Eye Link
                
                % Sending a 'TRIALID' message to mark the start of a trial in Data
                % Viewer.  This is different than the start of recording message
                % START that is logged when the trial recording begins. The viewer
                % will not parse any messages, events, or samples, that exist in
                % the data file prior to this message.
                Eyelink('Message', 'TRIALID %d', Trial);
                % This supplies the title at the bottom of the eyetracker display
                % Eyelink('command', 'record_status_message "TRIAL %d/%d  %s"', i,3);
                % Before recording, we place reference graphics on the host display
                % Must be offline to draw to EyeLink screen
                Eyelink('Command', 'set_idle_mode');
                % clear tracker display and draw box at center
                Eyelink('Command', 'clear_screen 0');
                
                % transfer the image to show it on the eye tracker computers so you
                % can watch the eye gaze relative to what they are actually looking
                % at.  but this is not necessary to present the stimuli on the
                % display PC
                
                % Eyelink('command', 'draw_box %d %d %d %d 15', width/2-50, height/2-50, width/2+50, height/2+50);
                
                %% STEP 7.2 Do a drift correction at the beginning of each trial
                % Performing drift correction (checking) is optional for
                % EyeLink 1000 eye trackers.
                % EyelinkDoDriftCorrection(el);  not mandatory
                
                %% STEP 7.3
                % start recording eye position (preceded by a short pause so that
                % the tracker can finish the mode transition)
                % The paramerters for the 'StartRecording' call controls the
                % file_samples, file_events, link_samples, link_events availability
                Eyelink('Command', 'set_idle_mode');  % definitely want to do before recording**
                WaitSecs(0.05); % more than enough time
                % Eyelink('StartRecording', 1, 1, 1, 1);
                Eyelink('StartRecording');
                % record a few samples before we actually start displaying
                % otherwise you may lose a few msec of data
                WaitSecs(0.1);  % build up a buffer of eye information over the link but would not be necessary if you are not doing gaze contingent stuff.
                
                Eyelink('Message', 'TrialStart');
                
            case 2 % Tobii
                
                talk2tobii('START_TRACKING')
                
                %check status of TETAPI
                cond_res = check_status(7, max_wait, tim_interv,1);
                
                tmp = find(cond_res==0);
                
                if( ~isempty(tmp) )
                    error('check_status has failed');
                end
                
                talk2tobii('EVENT', 'TrialStart', 0, 'Trial#', Trial);
                
                talk2tobii('RECORD');
                
        end
        %%
        
        % --------------
        % animation loop
        % --------------        
        imageF = cell(1, 1);
        
        for i = 1:nframes
            
            ISI = ISI + 1;
            
            if DirectionalCount >= NofDisplay % the number of directional motion displayed
                
                break;
                
            end
            
            if ISI >= 50 && ISI < 60
                
                BrightnessIndex = (60 - ISI) / 20 + .45;
                
            end
            
            % Draw nice dots:
            Screen('DrawDots', w, xymatrix, s, colvect .* BrightnessIndex, [0 0], 1);  % change 1 to 0 or 4 to draw square dots
            
            % Break out of animation loop if any key on keyboard or any button
            % on mouse is pressed:
            [mx, my, buttons] = GetMouse(screenNumber);
            
            if any(buttons)
                
                ii = ii + 1;
                break;
                
            end
            
            %% Eye Tracker
            
            switch Eyetracking
                
                case 1 % Eye Link
                    
                    if Eyelink('newfloatsampleavailable')
                        
                        % get the sample in the form of an event structure
                        evt = Eyelink('newestfloatsample');
                        
                        if eye_used ~= -1 % do we know which eye to use yet?
                            % if we do, get current gaze position from sample
                            gx = evt.gx(eye_used + 1); % +1 as we're accessing MATLAB array
                            gy = evt.gy(eye_used + 1);
                            
                            % if valid eye data
                            
                        else % if we don't, first find eye that's being tracked
                            
                            eye_used = Eyelink('eyeavailable'); % get eye that's tracked
                            
                            if eye_used == el.BINOCULAR % if both eyes are tracked
                                
                                eye_used = el.RIGHT_EYE; % use right eye
                                
                            end
                            
                        end
                        
                    end
                    
                case 2 % Tobii
                    
                    gaze = talk2tobii('GET_SAMPLE');
                    
                    gazeL_x = gaze(1) .* rect(3);
                    gazeL_y = gaze(2) .* rect(4);
                    gazeR_x = gaze(3) .* rect(3);
                    gazeR_y = gaze(4) .* rect(4);
                    
                    if (gazeL_x > 0) && (gazeR_x > 0)
                        
                        gx = (gazeL_x + gazeR_x) / 2;
                        
                    elseif (gazeL_x > 0)
                        
                        gx = gazeL_x;
                        
                    elseif (gazeR_x > 0)
                        
                        gx = gazeR_x;
                        
                    else
                        
                        gx = -1;
                        
                    end
                    
                    if (gazeL_y > 0) && (gazeR_y > 0)
                        
                        gy = (gazeL_y + gazeR_y) / 2;
                        
                    elseif (gazeL_y > 0)
                        
                        gy = gazeL_y;
                        
                    elseif (gazeR_y > 0)
                        
                        gy = gazeR_y;
                        
                    else
                        
                        gy = -1;
                        
                    end
                    
            end
            
            %%
            [KeyIsDown, KeyTime, KeyCode] = KbCheck;
            
            if i > 60 % initial period to display the RDK.
                
                if ~Switch
                    
                    switch Eyetracking
                        
                        case 1 % Eye Link
                            % if eye tracker is used
                            
                            if GazeCount <= 5
                                
                                if gx ~= el.MISSING_DATA && gy ~= el.MISSING_DATA && evt.pa(eye_used + 1) || KeyCode(SpaceKey)
                                    
                                    if IsInRect(gx, gy, centeredRect) || KeyCode(SpaceKey)
                                        
                                        if PreviousGaze
                                            
                                            GazeCount = GazeCount + 1;
                                            
                                        else
                                            
                                            GazeCount = 1;
                                            
                                        end
                                        
                                        PreviousGaze = 1;
                                        
                                    else
                                        
                                        PreviousGaze = 0;
                                        
                                    end
                                    
                                end
                                
                            else
                                
                                Switch = 1;
                                RandomIniti = 1;
                                
                                ii = ii + 1;
                                GazeTrigger(ii, 1) = ~KeyCode(SpaceKey);
                                GazeTrigger(ii, 2) = toc;
                                
                            end
                            
                        case 2 % Tobii
                            
                            if GazeCount <= 5
                                
                                if IsInRect(gx, gy, centeredRect) || KeyCode(SpaceKey)
                                    
                                    if PreviousGaze
                                        
                                        GazeCount = GazeCount + 1;
                                        
                                    else
                                        
                                        GazeCount = 1;
                                        
                                    end
                                    
                                    PreviousGaze = 1;
                                    
                                else
                                    
                                    PreviousGaze = 0;
                                    
                                end
                                
                            else
                                
                                Switch = 1;
                                RandomIniti = 1;
                                
                                ii = ii + 1;
                                GazeTrigger(ii, 1) = ~KeyCode(SpaceKey);
                                GazeTrigger(ii, 2) = toc;
                                
                            end
                            
                        case 0 % no eye tracker
                            
                            if GazeCount <= 5
                                
                                if IsInRect(mx, my, centeredRect) || KeyCode(SpaceKey)
                                    
                                    if PreviousGaze
                                        
                                        GazeCount = GazeCount + 1;
                                        
                                    else
                                        
                                        GazeCount = 1;
                                        
                                    end
                                    
                                    PreviousGaze = 1;
                                    
                                else
                                    
                                    PreviousGaze = 0;
                                    
                                end
                                
                            else
                                
                                Switch = 1;
                                RandomIniti = 1;
                                
                                ii = ii + 1;
                                GazeTrigger(ii, 1) = ~KeyCode(SpaceKey);
                                GazeTrigger(ii, 2) = toc;
                                
                            end
                            
                    end
                    
                end
                
            end
            
            %%
            if Switch
                
                AudioStatus = PsychPortAudio('GetStatus', pahandle);

                if (~AudioStatus.Active) && AudioIniti
                    
                    PsychPortAudio('Start', pahandle, 1, 0, 1);
                    
                    AudioStartTime = GetSecs();
                    
                    AudioIniti = 0;
                    
                    switch Eyetracking
                        
                        case 1 % Eye Link
                            
                            Eyelink('Message', 'SoundPlay');
                            
                        case 2 % Tobii
                            
                            talk2tobii('EVENT', 'SoundPlay', 0, 'Trial#', Trial);
                            
                    end
                    
                else
                    
                    if GetSecs() - AudioStartTime >= .5 && GetSecs() - AudioStartTime <= .6
                        
                        %                         Screen('FillRect', w, [0 0 0], rect);
                        WaitSecs(.1);
                        
                    end
                    
                    if GetSecs() - AudioStartTime >= .6 && GetSecs() - AudioStartTime <= (TrialDuration + AudioDuration) %#ok<*BDSCI>
                        
                        BrightnessIndex = 1;
                        
                        if LearningTest % whether this is a learning test trial
                            
                            RDK_LearningTest;
                            
                        else
                            
                            RDK_Familiarization;
                            
                        end
                        
                    else
                        
                        %                         RDK_Random;
                        RDK_Stop;
                        
                        % Terminate each trial within a block
                        % Reset parameters for the following trial
                        if GetSecs() - AudioStartTime > (.6 + AudioDuration)
                            
                            Switch = 0;
                            PreviousGaze = 0;
                            GazeCount = 0;
                            DirectionInit = 1;
                            AudioIniti = 1;
                            
                            DirectionalCount = DirectionalCount + 1;
                            
                            ISI = floor(rand() * 30); % jitter
                            
                        end
                        
                    end
                    
                end
                
            else
                
                if RandomIniti
                    
                    WaitSecs(.1);
                    
                    t = 2 * pi * rand(ndots,1);
                    cs = [cos(t), sin(t)];
                                        
                    dxdy = [pfs pfs] .* cs;
                    
                    RandomIniti = 0;
                    
                end
                %    
                % ==============
                % Draw a box to get attention
                % ==============
                
                if ISI > 60
                    
                    %                     Size = sin(i / 10) * 30 + 150;
                    %                     AG = CenterRectOnPointd([0 0 Size, Size], AGcenter(1), AGcenter(2));
                    
                    if strcmp(Shape, 'FillPoly')
                        
                        L = linspace(0,2.* pi, 6 + 1) + pi/6;
                        
                        Vloc = sin(i / 2 / 10) * rect(4)/10 + AGcenter(2);
                        
                        xv = cos(L)' .* 90 + AGcenter(1);
                        yv = sin(L)' .* 90 + Vloc;
                        
                        Screen(Shape, w, ShapeColor, [xv yv], 1);
                        
                    else
                        
                        Vloc = sin(i / 2 / 10) * rect(4)/10 + AGcenter(2);
                        AG = CenterRectOnPointd([0 0 180, 180], AGcenter(1), Vloc);
                        
                        Screen(Shape, w, ShapeColor, AG);
                        
                    end
                    
                    RDK_Stop;
                    
                else
                    
                    RDK_Random;
                    
                end
                
            end
            
            %% show the location of eye, if eye tracking is on; or the location of mouse, if eye tracking is off.
            
            if ShowGaze
                
                switch Eyetracking
                    
                    case 1 % Eye Link
                        
                        if exist('gx', 'var') == 1 && exist('gy', 'var') == 1
                            
                            cursorPos = CenterRectOnPointd([0 0 20 20], gx, gy);
                            Screen('FillOval', w, [255 255 0], cursorPos);
                            
                        end
                        
                    case 2 % Tobii
                        
                        cursorPos = CenterRectOnPointd([0 0 20 20], gx, gy);
                        Screen('FillOval', w, [255 255 0], cursorPos);
                        
                    case 0 % no eye tracker
                        
                        cursorPos = CenterRectOnPointd([0 0 20 20], mx, my);
                        Screen('FillOval', w, [255 255 0], cursorPos);
                        
                end
                
            end
            
            Screen('DrawingFinished', w);
            
            xy = xy + dxdy; % move dots
            size(dxdy(dxdy(:, 1) == 0, :), 1)
            xymatrix = transpose(xy);
            
            vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
            
%             imageF{i, 1} = Screen('GetImage', w);
%             
%             WaitSecs(20/60);
            
        end
        
        %% Eye Tracker
        switch Eyetracking
            
            case 1 % Eye Link
                
                %% STEP 7.6
                % Clear the display
                Eyelink('Message', 'BLANK_SCREEN');
                % adds 100 msec of data to catch final events
                WaitSecs(0.1);
                % stop the recording of eye-movements for the current trial
                Eyelink('StopRecording');
                
                %% STEP 7.7
                % Send out necessary integration messages for data analysis
                % Send out interest area information for the trial
                % See "Protocol for EyeLink Data to Viewer Integration-> Interest
                % Area Commands" section of the EyeLink Data Viewer User Manual
                % IMPORTANT! Don't send too many messages in a very short period of
                % time or the EyeLink tracker may not be able to write them all
                % to the EDF file.
                % Consider adding a short delay every few messages.
                
                % Send messages to report trial condition information
                % Each message may be a pair of trial condition variable and its
                % corresponding value follwing the '!V TRIAL_VAR' token message
                % See "Protocol for EyeLink Data to Viewer Integration-> Trial
                % Message Commands" section of the EyeLink Data Viewer User Manual
                WaitSecs(0.001); %%% in between every 5-10 messages wait at least 1ms because we could flood the channel and start to miss some messages
                Eyelink('Message', '!V TRIAL_VAR MotionDirection %d', rad2deg(angle(1)));% rad2deg(angle))
                Eyelink('Message', '!V TRIAL_VAR Coherence %d', Coherence(1) * 100);
                Eyelink('Message', '!V TRIAL_VAR LearningTest %d', LearningTest);
                % Eyelink('Message', '!V TRIAL_VAR imgfile %s', imgfile) NEED TO WORK ON THIS.
                % can send any number of variables here for groups for data analysis later (e.g, conditions etc)
                
                %% STEP 7.8
                % Sending a 'TRIAL_RESULT' message to mark the end of a trial in
                % Data Viewer. This is different than the end of recording message
                % END that is logged when the trial recording ends. The viewer will
                % not parse any messages, events, or samples that exist in the data
                % file after this message.
                Eyelink('Message', 'TRIAL_RESULT 0');
                
            case 2 % Tobii
                
                talk2tobii('EVENT', 'TrialEnd', 0, 'Trial#', Trial);
                
                talk2tobii('STOP_RECORD');
                
                talk2tobii('Stop_Tracking');
                
                talk2tobii('SAVE_DATA', ['Results/', Participant, '_eye_trackin_data.txt'],...
                    ['Results/', Participant, '_events.txt'], 'APPEND');
                
        end
        
        %% save recordings
        
        if G.Recording
            
            DrawFormattedText(w, '+', 'center', 'center', white);
            Screen('Flip', w);
            
            RecordingExport;
            
        end
        
        
    end
    
    if exist('vidobj', 'var')
        
        stop(vidobj)
        
    end
    
    %% Eye Tracker
    switch Eyetracking
        
        case 1 % Eye Link
            
            % Close EyeLink
            EyeLinkFinalize;
            
        case 2 % Tobii
            
            % Close Tobii
            TobiiFinalize;
            
    end
    
    Priority(0);
    ShowCursor;
    sca;
    
    save(['Results/', 'RDK_', Participant, '.mat'])
    
catch ME
    
    cd('~/Documents/MATLAB/Topdown_Preception')
    
    if exist('vidobj', 'var')
        
        stop(vidobj)
        
    end
    
    switch Eyetracking
        
        case 1 % Eye Link
            
            % Close EyeLink
            EyeLinkFinalize;
            
        case 2 % Tobii
            
            % Close Tobii
            TobiiFinalize;
            
    end
    
    Priority(0);
    ShowCursor;
    sca;
    
%     cd('~/Documents/MATLAB/Topdown_Preception/')
    save(['Results/', 'RDK_', Participant, '.mat'])
    
    rethrow(ME)
    
end

ShowCursor(0);
sca;