% parameters for RDK experiment
function [T, G] = ParametersRDK()
rng('default')
%% General parameters

G.Eyetracking = 0;  %0: no eye tracker, 1: EyeLink, 2: Tobii
G.ShowGaze = 1;
G.Recording = 0;

G.nframes     = 10000; % number of animation frames in loop
G.mon_width   = 35;   % horizontal dimension of viewable screen (cm)
G.v_dist      = 60;   % viewing distance (cm)

G.dot_speed   = 8;    % dot speed (deg/sec)
G.f_kill      = 0.05; % fraction of dots to kill each frame (limited lifetime)

G.ndots       = 1200;  % number of dots
G.dot_w       = 0.2;  % width of dot (deg)

G.differentcolors = 1; % Use a different color for each point if == 1. Use common color white if == 0.
G.differentsizes = 0; % Use different sizes for each point if >= 1. Use one common size if == 0.
G.waitframes = 1;     % Show new dot-images at each waitframes'th monitor refresh.

LearningTest = 0; % whether add learning test phase

%% Stimuli parameters
rng('shuffle')

SoundIndex = randperm(4);
SoundIndex = SoundIndex(1:2);
SoundIndex = SoundIndex(randperm(2));
SoundIndex = SoundIndex';
SoundIndex(2) = SoundIndex(2) + 4;

ColorIndex = [[180 32 98]; [32 180 98]];
ColorIndex = ColorIndex(randperm(2), :);

Direction = [1; 0];
Direction = Direction(randperm(2));

VisualShape = ['FillRect'; 'FillOval'];
VisualShape = VisualShape(randperm(2), :);

Fami_TrialDuration = 2;

%% Familiarization phase
NofFamiBlockperDirection = 3;
Fami_NofDisplayPerBlock = 4;
Fami_Coherence = repmat(.60, 1, Fami_NofDisplayPerBlock);

Fami_Coherence = repmat(Fami_Coherence, 2 * NofFamiBlockperDirection, 1);
% Fami_Coherence = reshape(Fami_Coherence, [], 1);

Fami_Stimuli_Order = repmat([1 2; 2 1], NofFamiBlockperDirection, 1);
Fami_Stimuli_Order = Fami_Stimuli_Order(randperm(NofFamiBlockperDirection), :);
Fami_Stimuli_Order = reshape(Fami_Stimuli_Order', 1, []);

Fami_VisualShape = VisualShape(Fami_Stimuli_Order', :);

Fami_SoundIndex = SoundIndex(Fami_Stimuli_Order', :);

Fami_ColorIndex = ColorIndex(Fami_Stimuli_Order', :);

Fami_Direction = Direction(Fami_Stimuli_Order', :);
Fami_Direction = repmat(Fami_Direction, 1, 4);

Fami_NofDisplay = repmat(Fami_NofDisplayPerBlock, size(Fami_Coherence, 1), 1);

Fami_LearningTest = zeros(size(Fami_Coherence, 1), 1);

Fami_TrialDuration = repmat(Fami_TrialDuration, size(Fami_Coherence, 1), 1);

FamiliarizationPhase = table(Fami_Coherence, Fami_VisualShape, Fami_SoundIndex,...
    Fami_ColorIndex, Fami_Direction, Fami_NofDisplay, Fami_LearningTest, Fami_TrialDuration);

%% Learning test phase

if LearningTest
    
    LearningTest_Coherence = [.50; .50; .50; .50];
    
    LearningTest_Order = [1, 2, 2, 1];
    
    LearningTest_VisualShape = VisualShape(LearningTest_Order, :);
    
    LearningTest_SoundIndex = SoundIndex(LearningTest_Order, :);
    
    LearningTest_ColorIndex = ColorIndex(LearningTest_Order, :);
    
    LearningTest_Direction = Direction(LearningTest_Order, :);
    LearningTest_Direction = repmat(LearningTest_Direction, 1, 4);
    
    LearningTest_NofDisplay = ones(size(LearningTest_Order', 1), 1);
    
    LearningTest_LearningTest = ones(size(LearningTest_Order', 1), 1);
    
    LearningTest_TrialDuration = 5;
    
    LearningTest_TrialDuration = repmat(LearningTest_TrialDuration, size(LearningTest_Order', 1), 1);
    
    LearningTest_Phase = table(LearningTest_Coherence, LearningTest_VisualShape, LearningTest_SoundIndex,...
        LearningTest_ColorIndex, LearningTest_Direction, LearningTest_NofDisplay, LearningTest_LearningTest, LearningTest_TrialDuration);
    
end
%% Test phase

NofTestBlockperCondition = 10;

Test_Coherence = [.6, .6, .6, 0];
NofItemperBlock = length(Test_Coherence);

m = perms(1:NofItemperBlock);
m = m(randperm(size(m, 1)), :);

m = m(1:NofTestBlockperCondition * 2, :);

m = reshape(m, 1, 4 * 2 * NofTestBlockperCondition);

Test_Coherence = Test_Coherence(m);
Test_Coherence = reshape(Test_Coherence, 2 * NofTestBlockperCondition, 4);

Test_Order = repmat([1 2; 2 1], NofTestBlockperCondition, 1);
Test_Order = Test_Order(randperm(NofTestBlockperCondition), :);
Test_Order = reshape(Test_Order', 1, []);

Test_VisualShape = VisualShape(Test_Order, :);

Test_SoundIndex = SoundIndex(Test_Order, :);

Test_ColorIndex = ColorIndex(Test_Order, :);

Test_Direction = Direction(Test_Order, :);
Test_Direction = repmat(Test_Direction, 1, 4);

Test_NofDisplay = ones(size(Test_Order, 2), 1) * 4;

Test_LearningTest = zeros(size(Test_Order, 2), 1);

Test_TrialDuration = 2;

Test_TrialDuration = repmat(Test_TrialDuration, size(Test_Order, 2), 1);

TestPhase = table(Test_Coherence, Test_VisualShape, Test_SoundIndex,...
    Test_ColorIndex, Test_Direction, Test_NofDisplay, Test_LearningTest, Test_TrialDuration);
%% Reversed direction test phase

ReversedPhase = FamiliarizationPhase(1:4, :);

for i = 1:size(ReversedPhase, 1)
    
    j = randi(4);
    
    ReversedPhase.Fami_Direction(i, j) = 1 - ReversedPhase.Fami_Direction(i, j);
    
end
%%
FamiliarizationPhase.Properties.VariableNames = {'Coherence' 'VisualShape' 'SoundIndex' 'ColorIndex' 'Direction' 'NofDisplay' 'LearningTest', 'TrialDuration'};
TestPhase.Properties.VariableNames = {'Coherence' 'VisualShape' 'SoundIndex' 'ColorIndex' 'Direction' 'NofDisplay' 'LearningTest', 'TrialDuration'};
ReversedPhase.Properties.VariableNames = {'Coherence' 'VisualShape' 'SoundIndex' 'ColorIndex' 'Direction' 'NofDisplay' 'LearningTest', 'TrialDuration'};

if LearningTest
    
    LearningTest_Phase.Properties.VariableNames = {'Coherence' 'VisualShape' 'SoundIndex' 'ColorIndex' 'Direction' 'NofDisplay' 'LearningTest', 'TrialDuration'};
    
    T = [FamiliarizationPhase; LearningTest_Phase; TestPhase];
    
else
    
    T = [FamiliarizationPhase; TestPhase; ReversedPhase];
    
end

