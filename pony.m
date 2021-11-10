clear
clc
rng('shuffle')

%% setup parameters
pinWheelSize = 4;% degrees of visual angle
screenSize = 32;% in degrees, assuming baby is 60cm from the 1280X1024 screen
screenWidth = 1280;% in pixels
screenHeight = 1024;% in pixels
locRadius = 8;% degrees
margin = 2;% degrees
nLocs = 8;
pinSpeed = 360;% deg/sec
fps = 30;
pre = .2;
isi = .5;% in sec
dur = .3;% in sec once baby is on it
post = .2;
lookOnTh = .3;% in sec for attention grabber
feedback = 1;% in sec
timeOut = 6;% in sec
minAlpha = .1;
maxAlpha = 1;
boxGray = .8;
nTrials = 105;

%% setup stimuli order
deltaAlpha = zeros(nTrials,1);
for ii = 1:nTrials
    if ii <= 10
        deltaRange = [2,4];
    elseif ii <= 20
        deltaRange = [1.75,3];
    elseif ii <= 30
        deltaRange = [1.5,2.5];
    elseif ii <= 40
        deltaRange = [1.25,2];
    elseif ii <= 50
        deltaRange = [1.2,1.75];
    else
        deltaRange = [1.1,1.5];
    end
    deltaAlpha(ii) = diff(deltaRange)*rand + deltaRange(1);
end
% deltaAlpha = [4*ones(floor(nTrials/4),1);2*ones(floor(nTrials/4),1);1.5*ones(ceil(nTrials/4),1);1.2*ones(ceil(nTrials/4),1)];
zeroInd = (1:20)';
while any(zeroInd(1:end-1) + 1 == zeroInd(2:end))
    zeroInd = (5:5:100)' + randi(5,20,1);
end
deltaAlpha(zeroInd) = 1;
alphas = zeros(nTrials,2);
for ii = 1:nTrials
    temp1 = 0;
    temp2 = 0;
    while temp1 < minAlpha || temp1 > maxAlpha || temp2 < minAlpha || temp2 > maxAlpha
        temp1 = minAlpha + (maxAlpha - minAlpha)*rand;
        if rand > .5
            temp2 = temp1*deltaAlpha(ii);
        else
            temp2 = temp1/deltaAlpha(ii);
        end
    end
    if rand > .5
        alphas(ii,1) = temp1;
        alphas(ii,2) = temp2;
    else
        alphas(ii,1) = temp2;
        alphas(ii,2) = temp1;
    end
end
[val,correctInd] = max(alphas,[],2);
correctInd(zeroInd) = 0;
chosenRect = nan(nTrials,1);
rt = nan(nTrials,1);

%% setup locations array
screenCenter = round([screenWidth,screenHeight]/2);
rectSize = round(pinWheelSize*screenWidth/screenSize);
locRadius = round(locRadius*screenWidth/screenSize);
locRects = nan(nLocs,4);
marginSize = round(margin*screenWidth/screenSize);
for ii = 1:nLocs
    theta = (ii - 1)*(2*pi/nLocs);
    rectCenter = round(screenCenter + locRadius*[cos(theta),sin(theta)]);
    locRects(ii,:) = round([rectCenter - rectSize/2,rectCenter + rectSize/2]);
end
centerRect = [screenCenter - rectSize/2,screenCenter + rectSize/2];

%% setup locations order
locInds = nan(nTrials,2);
locInds(1,1) = randi(nLocs);
temp = find(~ismember(1:nLocs,locInds(1,1) + [-1,0,1,nLocs - 1,1 - nLocs]));
temp = Shuffle(temp);
locInds(1,2) = temp(1);
for ii = 2:nTrials
    temp = find(~ismember(1:nLocs,locInds(ii - 1,:)));
    temp = Shuffle(temp);
    locInds(ii,1) = temp(1);
    temp = find(~ismember(1:nLocs,[locInds(ii - 1,:),locInds(ii,1) + [-1,0,1,nLocs - 1,1 - nLocs]]));
    temp = Shuffle(temp);
    locInds(ii,2) = temp(1);
end

%% get subject name and deal with files
subName = input('Subject name? ','s');
if ~exist(['../Data/RawData/',subName],'dir')
    mkdir(['../Data/RawData/',subName])
end
prevFileList = dir(['../Data/RawData/',subName,'/',subName,'_pony_*.mat']);
if strcmpi(subName,'XXXX9999')
    session = 0;
else
    session = length(prevFileList) + 1;
end
fileName = [subName,'_pony_',num2str(session),'.mat'];
addpath('EyeLinkFunctions')

%% initialize Psychtoolbox stuff
FlushEvents;
AssertOpenGL;
InitializePsychSound(1);
commandwindow;
KbName('UnifyKeyNames');
quitKey = KbName('ESCAPE');
nextKey = KbName('SPACE');
driftKey = KbName('d');
recalKey = KbName('c');

%% setup screen
screens=Screen('Screens');
screenNumber=max(screens);
Screen('Preference', 'SkipSyncTests', 1);

% Hide the mouse cursor and supress output in command window:
%     HideCursor;
ListenChar(2);

% Returns as default the mean gray value of screen:
% gray=GrayIndex(screenNumber);
gray=BlackIndex(screenNumber);
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
boxGray = uint8(boxGray*white);

%% Audio settings
fs = 44100;
nChannels = 2;
devices = PsychPortAudio('GetDevices');
deviceIndex = length(devices)-1;
pa = PsychPortAudio('Open',deviceIndex);
[origTestSound,origFs] = audioread('../Stimuli/sounds/Ao-Camptown Races-less24dB-plus10dB-mono.wav');
testSound = resample(origTestSound,fs,origFs);
PsychPortAudio('FillBuffer',pa,repmat(testSound',nChannels,1));
PsychPortAudio('Start',pa,0);% inf repetitions
str=sprintf('\nAdjust volume level.\nPress any key to continue...\n');
disp(str);
[KeyIsDown, KeyTime, KeyCode]=KbCheck;
while KeyCode(nextKey) == 0 && KeyCode(quitKey) == 0
    [KeyIsDown, KeyTime, KeyCode]=KbCheck;
end
PsychPortAudio('Stop',pa);
clear testSound origTestSound

%% fix the attention grabber background
fileList = dir('../Stimuli/attnPNGs/img*.png');
fileList = [fileList;flipud(fileList)];
attnIms = struct('im',cell(length(fileList),1));
for ii = 1:length(fileList)
    im = imread(['../Stimuli/attnPNGs/',fileList(ii).name]);
    imHSV = rgb2hsv(im);
    blackInd = imHSV(:,:,2) < .5;
    imR = im(:,:,1);
    imG = im(:,:,2);
    imB = im(:,:,3);
    imR(blackInd) = gray;
    imG(blackInd) = gray;
    imB(blackInd) = gray;
    sizeIm = min(size(imR));
    imR = imR(size(imR,1)/2 - sizeIm/2 + 1:size(imR,1)/2 + sizeIm/2,...
        size(imR,2)/2 - sizeIm/2 + 1:size(imR,2)/2 + sizeIm/2);
    imG = imG(size(imG,1)/2 - sizeIm/2 + 1:size(imG,1)/2 + sizeIm/2,...
        size(imG,2)/2 - sizeIm/2 + 1:size(imG,2)/2 + sizeIm/2);
    imB = imB(size(imB,1)/2 - sizeIm/2 + 1:size(imB,1)/2 + sizeIm/2,...
        size(imB,2)/2 - sizeIm/2 + 1:size(imB,2)/2 + sizeIm/2);
    fixIm = uint8(zeros(sizeIm,sizeIm,3));
    fixIm(:,:,1) = imR;
    fixIm(:,:,2) = imG;
    fixIm(:,:,3) = imB;
    resizeIm = imresize(fixIm,[rectSize,rectSize]);
    grayIm = rgb2gray(resizeIm);
    attnIms(ii).im = grayIm;
end

%% fix the pinwheel background and rotation speed
degPerFrame = pinSpeed/fps;
nFrames = ceil(360/degPerFrame);
pinwheelIms = struct('im',cell(nFrames,1));
im = imread('../Stimuli/pinwheelClipped.png');
for ii = 1:nFrames
    rotateIm = imrotate(im,ii*degPerFrame,'nearest','crop');
    imHSV = rgb2hsv(rotateIm);
    blackInd = imHSV(:,:,2) < .5;
    imR = rotateIm(:,:,1);
    imG = rotateIm(:,:,2);
    imB = rotateIm(:,:,3);
    imR(blackInd) = gray;
    imG(blackInd) = gray;
    imB(blackInd) = gray;
    sizeIm = min(size(imR));
    imR = imR(size(imR,1)/2 - sizeIm/2 + 1:size(imR,1)/2 + sizeIm/2,...
        size(imR,2)/2 - sizeIm/2 + 1:size(imR,2)/2 + sizeIm/2);
    imG = imG(size(imG,1)/2 - sizeIm/2 + 1:size(imG,1)/2 + sizeIm/2,...
        size(imG,2)/2 - sizeIm/2 + 1:size(imG,2)/2 + sizeIm/2);
    imB = imB(size(imB,1)/2 - sizeIm/2 + 1:size(imB,1)/2 + sizeIm/2,...
        size(imB,2)/2 - sizeIm/2 + 1:size(imB,2)/2 + sizeIm/2);
    fixIm = uint8(zeros(sizeIm,sizeIm,3));
    fixIm(:,:,1) = imR;
    fixIm(:,:,2) = imG;
    fixIm(:,:,3) = imB;
    resizeIm = imresize(fixIm,[rectSize,rectSize]);
    pinwheelIms(ii).im = resizeIm;
end

%% load sounds
soundFileNames = {'blueballzoom.wav';'boing1b.wav';'chimes.wav';...
    'chimesb.wav';'slidewhistle.wav';'windows nt 5 shu_1s.wav'};
soundMat = nan(length(soundFileNames),feedback*fs);
for ii = 1:length(soundFileNames)
    [origSound,origFs] = audioread(['../Stimuli/sounds/',soundFileNames{ii}]);
    while size(origSound,1) < feedback*origFs
        origSound = repmat(origSound,2,1);
    end
    croppedSound = origSound(1:feedback*origFs,1)';
    normSound = croppedSound./max([max(croppedSound),-min(croppedSound)]);
    soundMat(ii,:) = resample(normSound,fs,origFs);

end
clear origSound normSound croppedSound
%% enter experiment mode
try
    %% open screen
    [w,rect] = Screen('OpenWindow',screenNumber, gray);
    
    %% setup eyelink
    EyeLinkSetupPony
    eye_used = Eyelink('eyeavailable'); % get eye that's tracked
    if eye_used == el.BINOCULAR % if both eyes are tracked
        eye_used = el.RIGHT_EYE; % use right eye
    end
    %% trials loop
    for tt = 1:nTrials
        Eyelink('Message', 'TRIALID %d', tt);
        Eyelink('Command', 'set_idle_mode');
        Eyelink('Command', 'clear_screen 0');
        Eyelink('Command', 'set_idle_mode');  % definitely want to do before recording**
        WaitSecs(0.05); % more than enough time
        % Eyelink('StartRecording', 1, 1, 1, 1);
        Eyelink('StartRecording');
        WaitSecs(0.1);  % build up a buffer of eye information over the link but would not be necessary if you are not doing gaze contingent stuff.
        Eyelink('Message', 'TrialStart');
        
        %% load the feedback sound into the sound card buffer
        PsychPortAudio('FillBuffer',pa,repmat(soundMat(randi(size(soundMat,1)),:),nChannels,1));
        
        %% present attention grabber in center until baby is looking for onTime
        lookTime = 0;
        nextFrame = 1;
        lookOn = 0;
        lookWasOn = 0;
        frame = 1;
        % send trial number message to EyeLink
        Screen('FillRect',w,gray*ones(1,3));
        Screen('Flip',w);
        [KeyIsDown, KeyTime, KeyCode]=KbCheck;
        while KeyCode(quitKey)==0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                lookTime < lookOnTh
            if nextFrame
                tex = Screen('MakeTexture', w, attnIms(frame).im);
                Screen('DrawTexture', w, tex,[],centerRect);
                Screen('Flip', w);
                [KeyIsDown, frameTime, KeyCode]=KbCheck;
                Screen('Close', tex);
                frame = mod(frame,size(attnIms,1)) + 1;
                nextFrame = 0;
            end
            
            [KeyIsDown, KeyTime, KeyCode]=KbCheck;
%             [x,y] = GetMouse;% change to get the x,y from EyeLink
            if Eyelink('newfloatsampleavailable') 
                evt = Eyelink('newestfloatsample');% get the sample in the form of an event structure
                x = evt.gx(eye_used + 1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used + 1);
            end
            
            if x >= centerRect(1) - marginSize && x <= centerRect(1) + rectSize + marginSize &&...
                    y >= centerRect(2) - marginSize && y <= centerRect(2) + rectSize + marginSize
                lookOn = 1;
            else
                lookOn = 0;
            end
            
            if lookOn && ~lookWasOn
                lookOnset = KeyTime;
                lookWasOn = true;
            elseif ~lookOn && lookWasOn
                lookOffset = KeyTime;
                lookWasOn = false;
            end
            
            if lookOn
                lookTime = KeyTime - lookOnset;
            else
                lookTime = 0;
            end
            
            if KeyTime - frameTime >= 1/fps
                nextFrame = 1;
            end
        end
        
        %% pre blank
        Screen('FillRect',w,gray*ones(1,3));
        Screen('Flip',w);
        [KeyIsDown, KeyTime, KeyCode]=KbCheck;
        offTime = KeyTime;
        while KeyCode(quitKey) == 0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                KeyTime - offTime <= post
            [KeyIsDown, KeyTime, KeyCode]=KbCheck;
        end
        
        %% present first pinwheel until baby looks at it for dur
        lookTime = 0;
        nextFrame = 1;
        lookOn = 0;
        lookWasOn = 0;
        frame = 1;
        % send stim1 message to EyeLink
        currRect = locRects(locInds(tt,1),:);
        [KeyIsDown, KeyTime, KeyCode]=KbCheck;
        while KeyCode(quitKey)==0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                lookTime < dur
            if nextFrame
                im = pinwheelIms(frame).im;
                imHSV = rgb2hsv(im);
                imS = imHSV(:,:,2);
                imV = imHSV(:,:,3);
                nonGrayInd = imS > .1 & imV > .1;
                imS(nonGrayInd) = alphas(tt,1)*imS(nonGrayInd);
%                 imV(nonGrayInd) = alphas(tt,1)*imV(nonGrayInd);
                imHSV(:,:,2) = imS;
                imHSV(:,:,3) = imV;
                imRGB = uint8(255*hsv2rgb(imHSV));
                tex = Screen('MakeTexture', w, imRGB);
                Screen('DrawTexture', w, tex,[],currRect);
                Screen('Flip', w);
                [KeyIsDown, frameTime, KeyCode]=KbCheck;
                Screen('Close', tex);
                frame = mod(frame,size(pinwheelIms,1)) + 1;
                nextFrame = 0;
            end
            
            [KeyIsDown, KeyTime, KeyCode]=KbCheck;
%             [x,y] = GetMouse;% change to get the x,y from EyeLink
            if Eyelink('newfloatsampleavailable')
                evt = Eyelink('newestfloatsample');% get the sample in the form of an event structure
                x = evt.gx(eye_used + 1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used + 1);
            end

            if x >= currRect(1) - marginSize && x <= currRect(1) + rectSize + marginSize &&...
                    y >= currRect(2) - marginSize && y <= currRect(2) + rectSize + marginSize
                lookOn = 1;
            else
                lookOn = 0;
            end
            
            if lookOn && ~lookWasOn
                lookOnset = KeyTime;
                lookWasOn = true;
            elseif ~lookOn && lookWasOn
                lookOffset = KeyTime;
                lookWasOn = false;
            end
            
            if lookOn
                lookTime = KeyTime - lookOnset;
            else
                lookTime = 0;
            end
            
            if KeyTime - frameTime >= 1/fps
                nextFrame = 1;
            end
        end
        
        %% ISI
        Screen('FillRect',w,gray*ones(1,3));
        Screen('Flip',w);
        [KeyIsDown, KeyTime, KeyCode] = KbCheck;
        offTime = KeyTime;
        while KeyCode(quitKey) == 0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                KeyTime - offTime <= isi
            [KeyIsDown, KeyTime, KeyCode] = KbCheck;
        end
        
        %% present second pinwheel until baby looks at it for dur
        lookTime = 0;
        nextFrame = 1;
        lookOn = 0;
        lookWasOn = 0;
        frame = 1;
        % send stim2 message to EyeLink
        currRect = locRects(locInds(tt,2),:);
        [KeyIsDown, KeyTime, KeyCode] = KbCheck;
        while KeyCode(quitKey) == 0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                lookTime < dur
            if nextFrame
                im = pinwheelIms(frame).im;
                imHSV = rgb2hsv(im);
                imS = imHSV(:,:,2);
                imV = imHSV(:,:,3);
                nonGrayInd = imS > .1 & imV > .1;
                imS(nonGrayInd) = alphas(tt,2)*imS(nonGrayInd);
%                 imV(nonGrayInd) = alphas(tt,2)*imV(nonGrayInd);
                imHSV(:,:,2) = imS;
                imHSV(:,:,3) = imV;
                imRGB = uint8(255*hsv2rgb(imHSV));
                tex = Screen('MakeTexture', w, imRGB);
                Screen('DrawTexture', w, tex,[],currRect);
                Screen('Flip', w);
                [KeyIsDown, frameTime, KeyCode] = KbCheck;
                Screen('Close', tex);
                frame = mod(frame,size(pinwheelIms,1)) + 1;
                nextFrame = 0;
            end
            
            [KeyIsDown, KeyTime, KeyCode]=KbCheck;
%             [x,y] = GetMouse;% change to get the x,y from EyeLink
            if Eyelink('newfloatsampleavailable') 
                evt = Eyelink('newestfloatsample');% get the sample in the form of an event structure
                x = evt.gx(eye_used + 1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used + 1);
            end

            if x >= currRect(1) - marginSize && x <= currRect(1) + rectSize + marginSize &&...
                    y >= currRect(2) - marginSize && y <= currRect(2) + rectSize + marginSize
                lookOn = 1;
            else
                lookOn = 0;
            end
            
            if lookOn && ~lookWasOn
                lookOnset = KeyTime;
                lookWasOn = true;
            elseif ~lookOn && lookWasOn
                lookOffset = KeyTime;
                lookWasOn = false;
            end
            
            if lookOn
                lookTime = KeyTime - lookOnset;
            else
                lookTime = 0;
            end
            
            if KeyTime - frameTime >= 1/fps
                nextFrame = 1;
            end
        end
         
        %% post blank and wait for baby to disengage
        
        Screen('FillRect',w,gray*ones(1,3));
        Screen('Flip',w);
        [KeyIsDown, KeyTime, KeyCode] = KbCheck;
        offTime = KeyTime;
        while KeyCode(quitKey) == 0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                KeyTime - offTime <= post
            [KeyIsDown, KeyTime, KeyCode] = KbCheck;
        end
        offTime = KeyTime;
        
        lookTime = 0;
        nextFrame = 1;
        lookOn = 0;
        lookWasOn = 0;
        frame = 1;
        while KeyCode(quitKey)==0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                lookTime < post
            if nextFrame
                tex = Screen('MakeTexture', w, attnIms(frame).im);
                Screen('DrawTexture', w, tex,[],centerRect);
                Screen('Flip', w);
                [KeyIsDown, frameTime, KeyCode]=KbCheck;
                Screen('Close', tex);
                frame = mod(frame,size(attnIms,1)) + 1;
                nextFrame = 0;
            end
            
            [KeyIsDown, KeyTime, KeyCode]=KbCheck;
%             [x,y] = GetMouse;% change to get the x,y from EyeLink
            if Eyelink('newfloatsampleavailable') 
                evt = Eyelink('newestfloatsample');% get the sample in the form of an event structure
                x = evt.gx(eye_used + 1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used + 1);
            end
            
            if x >= centerRect(1) - marginSize && x <= centerRect(1) + rectSize + marginSize &&...
                    y >= centerRect(2) - marginSize && y <= centerRect(2) + rectSize + marginSize
                lookOn = 1;
            else
                lookOn = 0;
            end
            
            if lookOn && ~lookWasOn
                lookOnset = KeyTime;
                lookWasOn = true;
            elseif ~lookOn && lookWasOn
                lookOffset = KeyTime;
                lookWasOn = false;
            end
            
            if lookOn
                lookTime = KeyTime - lookOnset;
            else
                lookTime = 0;
            end
            
            if KeyTime - frameTime >= 1/fps
                nextFrame = 1;
            end
        end
        
        %% wait for baby to choose the correct box, or any box for zero trials
        Screen('FillRect',w,boxGray*uint8(ones(1,3)),[locRects(locInds(tt,1),:)',locRects(locInds(tt,2),:)']);
        lookTime = 0;
        lookOn = 0;
        lookWasOn = 0;
        if correctInd(tt) > 0
            corrRect = locRects(locInds(tt,correctInd(tt)),:);
            incoRect = locRects(locInds(tt,~ismember(1:2,correctInd(tt))),:);
        else
            corrRect = locRects(locInds(tt,1),:);
            incoRect = locRects(locInds(tt,2),:);
        end
        Screen('Flip',w);
        [KeyIsDown, KeyTime, KeyCode] = KbCheck;
        onset = KeyTime;
        while KeyCode(quitKey) == 0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                lookTime < lookOnTh && KeyTime - onset < timeOut
            
            [KeyIsDown, KeyTime, KeyCode] = KbCheck;
%             [x,y] = GetMouse;% change to get the x,y from EyeLink
            if Eyelink('newfloatsampleavailable')
                evt = Eyelink('newestfloatsample');% get the sample in the form of an event structure
                x = evt.gx(eye_used + 1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used + 1);
            end

            if correctInd(tt) > 0
                if x >= corrRect(1) - marginSize && x <= corrRect(1) + rectSize + marginSize &&...
                        y >= corrRect(2) - marginSize && y <= corrRect(2) + rectSize + marginSize
                    lookOn = 1;
                    if isnan(chosenRect(tt))
                        chosenRect(tt) = correctInd(tt);
                    end
                elseif x >= incoRect(1) - marginSize && x <= incoRect(1) + rectSize + marginSize &&...
                        y >= incoRect(2) - marginSize && y <= incoRect(2) + rectSize + marginSize
                    lookOn = 0;
                    if isnan(chosenRect(tt))
                        chosenRect(tt) = find(~ismember(1:2,correctInd(tt)));
                    end
                else
                    lookOn = 0;
                end
            else
                if (x >= corrRect(1) - marginSize && x <= corrRect(1) + rectSize + marginSize &&...
                        y >= corrRect(2) - marginSize && y <= corrRect(2) + rectSize + marginSize)
                    lookOn = 1;
                    chosenRect(tt) = 1;
                elseif (x >= incoRect(1) - marginSize && x <= incoRect(1) + rectSize + marginSize &&...
                        y >= incoRect(2) - marginSize && y <= incoRect(2) + rectSize + marginSize)
                    lookOn = 1;
                    chosenRect(tt) = 2;
                else
                    lookOn = 0;
                end
            end
            
            if lookOn && ~lookWasOn
                lookOnset = KeyTime;
                lookWasOn = true;
            elseif ~lookOn && lookWasOn
                lookOffset = KeyTime;
                lookWasOn = false;
            end
            
            if lookOn
                lookTime = KeyTime - lookOnset;
            else
                lookTime = 0;
            end
        end
        
        if KeyTime - onset < timeOut
            rt(tt) = KeyTime - onset - dur;
        else
            rt(tt) = nan;
        end
        % send response message to EyeLink

        %% present feedback
        if correctInd(tt) > 0
            currRect = locRects(locInds(tt,correctInd(tt)),:);
        elseif isnan(chosenRect(tt))
            currRect = locRects(locInds(tt,randi(2)),:);
        else
            currRect = locRects(locInds(tt,chosenRect(tt)),:);
        end
        
        lookTime = 0;
        nextFrame = 1;
        lookOn = 0;
        lookWasOn = 0;
        frame = 1;
        PsychPortAudio('Start',pa,0);% inf repetitions
        [KeyIsDown, KeyTime, KeyCode] = KbCheck;
        while KeyCode(quitKey) == 0 && KeyCode(nextKey) == 0 &&...
                KeyCode(driftKey) == 0 && KeyCode(recalKey) == 0 &&...
                lookTime < feedback
            if nextFrame
                im = pinwheelIms(frame).im;
                imHSV = rgb2hsv(im);
                imS = imHSV(:,:,2);
                imV = imHSV(:,:,3);
                nonGrayInd = imS > .1 & imV > .1;
                if correctInd(tt) > 0
                    imS(nonGrayInd) = alphas(tt,correctInd(tt))*imS(nonGrayInd);
                else
                    imS(nonGrayInd) = alphas(tt,1)*imS(nonGrayInd);
                end    
%                 imV(nonGrayInd) = alphas(tt,2)*imV(nonGrayInd);
                imHSV(:,:,2) = imS;
                imHSV(:,:,3) = imV;
                imRGB = uint8(255*hsv2rgb(imHSV));
                tex = Screen('MakeTexture', w, imRGB);
                Screen('DrawTexture', w, tex,[],currRect);
                Screen('Flip', w);
                [KeyIsDown, frameTime, KeyCode] = KbCheck;
                Screen('Close', tex);
                frame = mod(frame,size(pinwheelIms,1)) + 1;
                nextFrame = 0;
            end
            
            [KeyIsDown, KeyTime, KeyCode] = KbCheck;
%             [x,y] = GetMouse;% change to get the x,y from EyeLink
            if Eyelink('newfloatsampleavailable')
                evt = Eyelink('newestfloatsample');% get the sample in the form of an event structure
                x = evt.gx(eye_used + 1); % +1 as we're accessing MATLAB array
                y = evt.gy(eye_used + 1);
            end

%             if x >= currRect(1) - marginSize && x <= currRect(1) + rectSize + marginSize &&...
%                     y >= currRect(2) - marginSize && y <= currRect(2) + rectSize + marginSize
%                 lookOn = 1;
%             else
%                 lookOn = 0;
%             end

            lookOn = 1;

            if lookOn && ~lookWasOn
                lookOnset = KeyTime;
                lookWasOn = true;
            elseif ~lookOn && lookWasOn
                lookOffset = KeyTime;
                lookWasOn = false;
            end
            
            if lookOn
                lookTime = KeyTime - lookOnset;
            else
                lookTime = 0;
            end
            
            if KeyTime - frameTime >= 1/fps
                nextFrame = 1;
            end
        end
        PsychPortAudio('Stop',pa);

        if isnan(rt(tt))
            rt(tt) = KeyTime - onset - feedback;
        end
        
        save(['../Data/RawData/',subName,'/',fileName],'-regexp','^(?!(attnIms|pinwheelIms|*Im|im*|pinwheelRGB)$).')
        
        if KeyCode(nextKey)
            continue
        end
        
        if KeyCode(driftKey)            
            EyelinkDoDriftCorrection(el);
            eye_used = Eyelink('eyeavailable'); % get eye that's tracked
            if eye_used == el.BINOCULAR % if both eyes are tracked
                eye_used = el.RIGHT_EYE; % use right eye
            end
        end
        
        if KeyCode(recalKey)            
            EyelinkDoTrackerSetup(el,'c');
            eye_used = Eyelink('eyeavailable'); % get eye that's tracked
            if eye_used == el.BINOCULAR % if both eyes are tracked
                eye_used = el.RIGHT_EYE; % use right eye
            end
        end
        
        if KeyCode(quitKey)
            error('Aborted by user.');
        end
    end
    %% close EyeLink
    EyeLinkFinalize
    
    %% cleanup
    PsychPortAudio('Close');
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    ListenChar(0);
    FlushEvents;
    
catch ME
    sca;
    ShowCursor;
    ListenChar(0);
    save;
    EyeLinkFinalize
    rethrow(ME)
end
