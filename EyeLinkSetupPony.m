%% STEP 1 Create EDF File

% Added a dialog box to set your own EDF file name before opening
% experiment graphics. Make sure the entered EDF file name is 1 to 8
% characters in length and only numbers or letters are allowed.
% prompt = {'Participant Numbr:'};
% dlg_title = 'Startup';
% num_lines = 1;
% defaultans = {''};
% Participant = inputdlg(prompt,dlg_title,num_lines,defaultans);
% %edfFile= 'DEMO.EDF'
% Participant = Participant{1, 1};


edfFile = [subName,'_',num2str(session),'.edf'];
fprintf('EDFFile: %s\n', edfFile);

%% STEP 2 Psychtoolbox Setup
% Screen Setup
% PsychDefaultSetup(2);

% Get the screen numbers
% screens = Screen('Screens');
% 
% % Select the external screen if it is present, else revert to the native
% % screen
% screenNumber = max(screens);
% 
% % Define black, white and grey
% black = BlackIndex(screenNumber);
% white = WhiteIndex(screenNumber);
% grey = 149;
% 
% % Open an on screen window and color it grey
% %windowRect = [0,0, 900,600];
% %window = PsychImaging('OpenWindow', screenNumber, grey, windowRect);
% 
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);
% %use above if i have external monitor
% 
% % Set the blend funciton for the screen
% Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% 
% % Get the size of the on screen window in pixels
% [screenXpixels, screenYpixels] = Screen('WindowSize', window);
% 
% % Query the frame duration
% ifi = Screen('GetFlipInterval', window);
% 
% % Get the centre coordinate of the window in pixels
% [xCenter, yCenter] = RectCenter(windowRect);
% 
% % Set the text size
% Screen('TextFont', window, 'Ariel');
% Screen('TextSize', window, 20);

%% Play a movie to attract infants' attention.
moviename = ['/Users/eyelink/Documents/MATLAB/Pony/Stimuli/Pimpa.mp4']; %access video (muppets, to play during set-up)

% Open movie file:
AttentionGetterMovie = Screen('OpenMovie', w, moviename);   
% Start playback engine:
Screen('PlayMovie', AttentionGetterMovie, 1);    
% Playback loop:play movie til end of movie or keypress:

[KeyIsDown, KeyTime, KeyCode] = KbCheck;

while KeyCode(nextKey) %play movie until button press
    % Wait for next movie frame, retrieve texture handle to it
    tex = Screen('GetMovieImage', w, AttentionGetterMovie);
    % Valid texture returned? A negative value means end of movie reached:
    if tex <= 0
        % We're done, break out of loop:
        break;
    end
    % Draw the new texture immediately to screen:
    Screen('DrawTexture', w, tex);   
    % Update display:
    Screen('Flip', w);
    % Release texture:
    Screen('Close', tex);
end  
% Stop playback:
Screen('PlayMovie', AttentionGetterMovie, 0);    
% Close movie:
Screen('CloseMovie', AttentionGetterMovie);    

%% STEP 3 Eyelink Setup
% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).
el = EyelinkInitDefaults(w);
% We are changing calibration to a black background with white targets,
% no sound and smaller targets
el.backgroundcolour = BlackIndex(el.window);
% el.backgroundcolour = GrayIndex(el.window);
el.foregroundcolour = BlackIndex(el.window); %this is deleted in AntiSaccade.m
% el.foregroundcolour = GrayIndex(el.window); %this is deleted in AntiSaccade.m
el.msgfontcolour  = WhiteIndex(el.window);
el.imgtitlecolour = WhiteIndex(el.window);
el.targetbeep = 0;
el.calibrationtargetcolour = WhiteIndex(el.window);
% for lower resolutions you might have to play around with these values
% a little. If you would like to draw larger targets on lower res
% settings please edit PsychEyelinkDispatchCallback.m and see comments
% in the EyelinkDrawCalibrationTarget function
el.calibrationtargetsize = 1;
el.calibrationtargetwidth = 0.5;
el.calanimationtargetfilename = '/Users/eyelink/Documents/MATLAB/Pony/Stimuli/AG5GS.mov';
% [PsychtoolboxRoot 'PsychDemos/MovieDemos/cymbals.mov' ]; %DualDiscs.mov

EyelinkUpdateDefaults(el);

%% STEP 4 Initialization of the connection with the Eyelink Gazetracker.
% exit program if this fails.
if ~EyelinkInit(0) % dummymode
    fprintf('Eyelink Init aborted.\n');
    cleanup;  % cleanup function
    return;
end

[v vs]=Eyelink('GetTrackerVersion');  % not necessary, we only use EL1000
fprintf('Running experiment on a ''%s'' tracker.\n', vs );

edfFile
% open file to record data to
i = Eyelink('Openfile', edfFile);
if i~=0
    fprintf('Cannot create EDF file ''%s'' ', edfFile);
%     cleanup;
    % Eyelink( 'Shutdown');
    return;
end

% creates a custom message at the top of the EDF file, that last
% section of text, 245 character limit
Eyelink('command', 'Pony');

[width, height]=Screen('WindowSize', screenNumber);

%% STEP 5 SET UP TRACKER CONFIGURATION
% Setting the proper recording resolution, proper calibration type,
% as well as the data file content;
% Help for Data Viewer: Protocol
Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1); % tells to do something
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1); %logs it
% set calibration type.
Eyelink('command', 'calibration_type = HV5');  % set to HV5 or 3 for the smaller calibration for babies
% set sampling rate
Eyelink('command', 'sample_rate = 500');

% you must send this command with value NO for custom calibration
% you must also reset it to YES for subsequent experiments
Eyelink('command', 'generate_default_targets = NO');

% STEP 5.1 modify calibration and validation target locations
Eyelink('command','calibration_samples = 5');
Eyelink('command','calibration_sequence = 0,1,2,3,4');
Eyelink('command','calibration_targets = %d,%d %d,%d %d,%d %d,%d',...
    round(width*.5),round(height*.25),round(width*.75),round(height*.5),...
    round(width*.5),round(height*.75),round(width*.25),round(height*.5));

Eyelink('command','validation_samples = 5');
Eyelink('command','validation_sequence = 0,1,2,3,4');
Eyelink('command','validation_targets = %d,%d %d,%d %d,%d %d,%d',...
    round(width*.3),round(height*.3), round(width*.7),round(height*.3),...
    round(width*.7),round(height*.7), round(width*.3),round(height*.7));

% set parser (conservative saccade thresholds)
Eyelink('command', 'saccade_velocity_threshold = 35');
Eyelink('command', 'saccade_acceleration_threshold = 9500');
% set EDF file contents
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,HTARGET'); % added HTARGET for eyelink 1000
% set link data (used for gaze cursor)
Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON'); % link things available over the link for gaze contingent stuff
Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,HTARGET');% added HTARGET for eyelink 1000
% allow to use the big button on the eyelink gamepad to accept the
% calibration/drift correction target
Eyelink('command', 'button_function 5 "accept_target_fixation"');

% throw in when you jump between m files to make sure that we're still connected
%if Eyelink('IsConnected')~=1
%    cleanup;
%    return;
%end;

%% Play Setup Movie
% DrawFormattedText(window, 'Press any button to Start!','center','center', white);
% Screen('Flip',window)
% % Wait for experimenter to press button
% while ~KbCheck
% end
% WaitSecs(1);

moviename = ['/Users/eyelink/Documents/MATLAB/Pony/Stimuli/Pimpa.mp4']; %access video (muppets, to play during set-up)

% Open movie file:
AttentionalGetter = Screen('OpenMovie', w, moviename);   
% Start playback engine:
Screen('PlayMovie', AttentionalGetter, 1);    

[keyIsDown,secs,keyCode] = KbCheck;

while ~keyCode(nextKey) %play movie until button press
    
    % Playback loop:play movie til end of movie or keypress:
    [keyIsDown,secs,keyCode] = KbCheck;
    
    % Wait for next movie frame, retrieve texture handle to it
    tex = Screen('GetMovieImage', w, AttentionalGetter);
    % Valid texture returned? A negative value means end of movie reached:
    if tex <= 0
        % We're done, break out of loop:
        break;
    end
    % Draw the new texture immediately to screen:
    Screen('DrawTexture', w, tex);   
    % Update display:
    Screen('Flip', w);
    % Release texture:
    Screen('Close', tex);
end  
% Stop playback:
Screen('PlayMovie', AttentionalGetter, 0);    
% Close movie:
Screen('CloseMovie', AttentionalGetter);    


%% STEP 6 Calibrate the eye tracker
% setup the proper calibration foreground and background colors
%el.backgroundcolour = 128;
%el.foregroundcolour = 0;
%EyelinkUpdateDefaults(el); % need to do this anytime you update the defaults

% Hide the mouse cursor;
Screen('HideCursorHelper', w);

EyelinkDoTrackerSetup(el); % the whole schzam (camera setup, calibration, validation) and no code will pass until you leave this section (press 'o')
