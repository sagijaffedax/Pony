%% STEP 8 End of Experiment; close the file first
% close graphics window, close data file and shut down tracker
% reset so tracker uses defaults calibration for other experiemnts
Eyelink('command', 'generate_default_targets = YES')
Eyelink('Command', 'set_idle_mode');
WaitSecs(0.5);
Eyelink('CloseFile');

% download data file  % try to make sure that we get the entire file
try
    fprintf('Receiving data file ''%s''\n', edfFile );
    status=Eyelink('ReceiveFile');
    if status > 0
        fprintf('ReceiveFile status %d\n', status);
    end
    if 2==exist(edfFile, 'file')
        fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
    end
catch
    fprintf('Problem receiving data file ''%s''\n', edfFile );
end

%% STEP 9
% run cleanup function (close the eye tracker and window).
% cleanup;