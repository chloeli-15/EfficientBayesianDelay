function bayesDriftExp(expMat,block,subject,priorType)

%% General set-up
skipsync = 1; % set to zero before testing
use_mouse = 1;
debug_mode = 1;
use_knob = 0;
stimulus = 2; % 1 - gabor; 2 - line
blockLength = 50;

if debug_mode
    dbstop if error
    PsychDebugWindowConfiguration
end
commandwindow

if use_knob, hndl = PsychPowerMate('Open'); end % open response knob
if ~use_mouse, HideCursor; else ShowCursor('Arrow'); end

%% Clean up previous initializations
try
    PsyFile('Close')
    Eyelink('CloseLiveFileStream');
    PsychPortAudio('Close');
    Eyelink('StopLiveRecording');
catch
end

%% Save data
wd = cd;
[FILENAME,PATHNAME] = uiputfile([wd '\Data\sub' num2str(subject) '_block_' num2str(block) '.mat']); % will need \Data\ folder in cd
if FILENAME==0, return; end


%% Screen measurements
screen.dist = 60;
screen.cm = [60; 33.75]; % dark room
screen.size = 2 * atan(screen.cm*0.5/screen.dist)/pi*180;

%% Initialize Psychtoolbox
PsychDefaultSetup(2);
Screen('Preference', 'VisualDebugLevel', 3);
Screen('Preference', 'TextRenderer', 0);
Screen('Preference', 'SkipSyncTests', skipsync);

[win, windowRect] = PsychImaging('OpenWindow', 0, [0 0 0],[],[],[],[],4);
screen.res = windowRect(3:4)';
Screen('BlendFunction', win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
screen.hz = Screen('FrameRate', win);

[xCentre, yCentre] = RectCenterCustom(windowRect);
screen.centre = [xCentre; yCentre];

gammaTable_orig = Screen('ReadNormalizedGammaTable', win);

%% Experiment parameters
% trial timing (in ms)
fixPrepare = 500;  % fixation cross presented before fixation is regitstered
fixWait = 500; % fixation cross presented before stimuli are presented
stimExposure = 500; % exposure time

% stimulus settings (sizes in dva and px)
exp.ecc_deg = 5;                                                    % stimuli eccentricity
exp.stimEcc_pix  = PosToPixels(exp.ecc_deg, screen, 1);
exp.stimSize_deg = 2.8;                                             % stimuli size
exp.stimSize_pix = PosToPixels(exp.stimSize_deg, screen, 1);
exp.startPos_deg = [0; 0];                                          % screen center
exp.startPos_px = PosToPixels(exp.startPos_deg, screen);
exp.kLoc = 6;                                                       % number of possible locations
exp.annRad = (exp.stimSize_deg+.3)/2;                               % placeholder radius
exp.annRadCue = (exp.stimSize_deg+.7)/2;
exp.annWidth = .1;                                                  % placeholder width
exp.annCol = [.4 .4 .4];                                            % placeholder colour
exp.barWidth = .3;                                                  % bar width
exp.barHeight = exp.stimSize_deg;                                   % bar height
exp.barColour = [.2 .2 .2];                                         % bar colour

% preferences
pref.outerAnnRad = 0.25;
pref.outerAnnCol = [0.2 0.2 0.2];
pref.innerAnnRadPre = 0.15;
pref.innerAnnRadPost = 0.2;
pref.innerAnnCol = [0.5 0.5 0.5];

pref.fixDeviation = 2.0;
pref.fixDeviationTime = .3;

gabor.amp = .3; % 0.5;
gabor.wavel = .65;
gabor.sd = 0.5;
gabor.phase = 0;

%% Load display calibration
screen.gammaTableFile = 'gammaTableColour_asusSwift_darkroom_sRGB_2017_03_17.mat';
load(fullfile('helperF', screen.gammaTableFile), 'gammaTable');
screen.gammaTable = gammaTable;
Screen('LoadNormalizedGammaTable', win, screen.gammaTable);

% clear to grey
Screen('FillRect', win, [0.5 0.5 0.5]);

%% Initialize Eyelink
if ~use_mouse
    EyelinkInit;
    el = EyelinkInitDefaults(win);
    Eyelink('Command','link_sample_data = LEFT,RIGHT,GAZE,PUPIL,AREA');
    Eyelink('Command','screen_pixel_coords = 0 0 %ld %ld', screen.res(1), screen.res(2));
    Eyelink('Command','online_dcorr_refposn %ld %ld', screen.centre(1), screen.centre(2));
    Eyelink('Command','calibration_type = HV9'); % HV13
    Eyelink('Command','generate_default_targets = YES');
    EyelinkDoTrackerSetup(el);
    Eyelink('StartRecording', 0, 0, 1, 1);
end

% Eyelink extensions && New Live recording stream
[~, filename_noext, ~] = fileparts(FILENAME);
PsyFile('Create',[PATHNAME filename_noext '.eye']);
Eyelink('SetMotionAnalysis',  'Basic')
Eyelink('SetSaccadeAnalysis', 'BuiltIn')
Eyelink('StartLiveRecording');
Eyelink('WaitForSamples');
Eyelink('SetRunningFileStream',1);

if use_mouse
    Eyelink('SetDataSource', 'Mouse');  
    PsyFile('AddMessage', 'NB: Using Mouse');
else            
    Eyelink('SetDataSource', 'Eye');    
    PsyFile('AddMessage', 'NB: Using Eye');
end

%% trial structure
exp.conditions = expMat(1+(blockLength*(block-1)):blockLength*block,:);

clear trial
for i = 1:length(exp.conditions)
    trialS(i).trialInd         = i;
    trialS(i).nStim            = exp.conditions(i,1);
    trialS(i).delay            = exp.conditions(i,2);
    trialS(i).orientation      = exp.conditions(i,4:6);             % all stimuli
    trialS(i).probe            = randi(trialS(i).nStim);            % probe
    trialS(i).orientation(trialS(i).probe) = exp.conditions(i,3);   % replace one stimulus with the target orientation
end


%% Initialize temp variables
rest_every = 600;
trial = 1; stage = 0; intro = 1; initTrial = 1; total_score = 0;
dialPos = 0; fixating = 1; checkBlinkDur = 1; respOn = 0; fixTime = nan;

% fieldSize settings - for Gabors
fieldSizePix = PosToPixels(exp.stimSize_deg, screen, 1)';
if mod(fieldSizePix(1),2) == 0; fieldSize = [fieldSizePix(1), fieldSizePix(1)]; exp.stimSize_pix = [fieldSizePix(1), fieldSizePix(1)];
elseif mod(fieldSizePix(1),2) == 1; fieldSize = [fieldSizePix(1)+1, fieldSizePix(1)+1]; exp.stimSize_pix = [fieldSizePix(1)+1, fieldSizePix(1)+1]; end

%% main trial loop
if use_mouse, SetMouse(xCentre, yCentre, win); end
out.screen = screen;
while trial  <= length(trialS)
    if respOn == 0 && use_knob == 1
    PsychPowerMate('Close', hndl); end

    [~, ~, keyCode] = KbCheck;
    if keyCode(KbName('LeftShift'))
        if exist('gammaTable_orig'); Screen('LoadNormalizedGammaTable', win, gammaTable_orig); end
        Eyelink('StopLiveRecording');
        Eyelink('SetRunningFileStream',0);  % close filestream
        PsyFile('Close');
        sca
    end

    if keyCode(KbName('d'))
        Eyelink('command','online_dcorr_trigger');
        while keyCode(KbName('d')), [~, ~, keyCode] = KbCheck; end
    end

    if keyCode(KbName('c'))
        EyelinkDoTrackerSetup(el);
        Eyelink('StartRecording', 0, 0, 1, 1);
    end
    
    if initTrial == 1
        %if use_mouse, HideCursor; end 
        if debug_mode, ShowCursor('Arrow'); end
        if ~use_mouse
            HideCursor;
            eyeLinkClearScreen(0);
            eyeLinkDrawText(900,900,4,['Trial ' num2str(trial)]); end
        if use_mouse, SetMouse(xCentre, yCentre, win); end
        respOn = 0;
        eventTime = nan(10,1);
        
        % get locations
        locAngles = linspace(0,2*pi,exp.kLoc+1); locAngles = locAngles(1:end-1);
        [Xs1, Ys1] = pol2cart(locAngles, repelem(exp.stimEcc_pix(1),exp.kLoc));
        trialS(trial).initialRects = CenterRectOnPoint([0 0 exp.stimSize_pix(1), exp.stimSize_pix(2)], (Xs1(1,:)+xCentre)', (Ys1(1,:)+yCentre)')'; % offsets the rect to center it around an x and y position
        [trialS(trial).xPix, trialS(trial).yPix] = RectCenter(trialS(trial).initialRects); % returns the integer x,y point closest to the center of a rect
        locRandom = randperm(exp.kLoc); % randomize locations
        locPicked = locRandom(1:trialS(trial).nStim); % pick nStim locations
        trialS(trial).xyPositionAll = PixelsToPos([trialS(trial).xPix; trialS(trial).yPix], screen); % returns the position of the center of a rect in dva
        trialS(trial).xyPosition = trialS(trial).xyPositionAll(:,locPicked); % used locations
        trialS(trial).targetPix = [trialS(trial).xPix(locPicked(trialS(trial).probe)),trialS(trial).yPix(locPicked(trialS(trial).probe))];
       
        if dialPos ~= 0 && use_knob == 1, PsychPowerMate('Close', hndl); end
        initTrial = 0;
    end

    if trial == 1 && intro == 1
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPre, pref.innerAnnCol)
        Screen('Flip', win);
        pause(1); 
        intro = 0;
    end

    % get current eye position/velocity
    sample = Eyelink('GetLiveSample');
    eye_px = [sample.GazePosX; sample.GazePosY]; eye_pos = PixelsToPos(eye_px, screen);
    eye_vx = [sample.GazeVelX; sample.GazeVelY]; eye_vel = PixelsToPos(eye_vx, screen, 1);
    eye_blink = sample.IsBlinking; if use_mouse; eye_blink = 0; end
    
    % check if observers are blinking during delay
    if ~use_mouse
        fixating = norm(eye_pos-exp.startPos_deg)<pref.fixDeviation;
        if fixating
            checkBlinkDur = 1;
            fixTime = nan;
        elseif ~fixating && stage == 4 && checkBlinkDur
            checkBlinkDur = 0;
            fixTime = GetSecs;
        end
    end
    
    % draw fixation circle and annulus
    if stage == 0
        PsyFile('AddMessage', ['Trial ' num2str(trial)])        
%         if (trial > 1 && mod(trial,rest_every) == 1)
%             DrawCircle(win, screen, [0; 0], 2, [1 1 1]);
%             DrawText(win, screen, [-1.35; -0.3], 'CONTINUE', [0 0 0], 10);
%             ShowCursor('Arrow');
%             Screen('Flip', win);
%             buttons = 0;
%             while (~buttons | sum(([mx; my]-PosToPixels([0; 0], screen)).^2) > mean(PosToPixels([2; 2], screen, 1)).^2 )
%                 [mx my buttons] = GetMouse(win);
%             end
%             if ~use_mouse, HideCursor(); end
%         end

        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPre, pref.innerAnnCol) 
        Screen('Flip', win);

        if isnan(eventTime(1)), eventTime(1) = sample.ELTimestamp; end
        PsyFile('AddMessage','Initial fixation');
        stage = stage + 1;
    end
    
    % change the color of fixation circle to grey  for 500ms and prepare the stimuli
    if stage == 1 && (~eye_blink) && (norm(eye_pos-exp.startPos_deg) < pref.fixDeviation) && (sample.CurrentFixationStartAgeMs > fixPrepare) && (sample.ELTimestamp-eventTime(1) > fixPrepare)
        if isnan(eventTime(2)), eventTime(2) = sample.ELTimestamp; PsyFile('AddMessage','Fixation satisfied'); end
        
        % change fixation
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPost, pref.innerAnnCol)
        Screen('Flip', win);

        % prepate presentation phase
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPost, pref.innerAnnCol)
        if stimulus == 1
            DrawGabor(win, screen, trialS(trial).initialRects(:,locPicked), gabor.amp, gabor.wavel, trialS(trial).orientation, gabor.phase, gabor.sd, fieldSize(1))
        elseif stimulus == 2
            for lk= 1:trialS(trial).nStim
                DrawTextBar(win, screen, trialS(trial).xyPosition(:,lk), trialS(trial).orientation(lk), exp.barWidth, exp.barHeight, exp.barColour); end
        end
        for jj= 1:trialS(trial).nStim
            DrawAnnulus(win, screen, trialS(trial).xyPosition(:,jj), exp.annRad, exp.annWidth, exp.annCol);
        end
        stage = stage + 1;
    end
    
    % presentation phase
    if stage == 2 && (~eye_blink) && (norm(eye_pos-exp.startPos_deg) < pref.fixDeviation) && (sample.ELTimestamp-eventTime(2) > fixWait)
        if isnan(eventTime(3)), eventTime(3) = sample.ELTimestamp; PsyFile('AddMessage','Stimuli presentation');end
        
        Screen('Flip', win);
        
        % prepare delay phase
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPost, pref.innerAnnCol)
        for jj= 1:trialS(trial).nStim
            DrawAnnulus(win, screen, trialS(trial).xyPosition(:,jj), exp.annRad, exp.annWidth, exp.annCol);
        end

        probe = trialS(trial).probe;
        
        % initialize variables
        firstTouch = nan;
        firstTouchConf = nan;
        randOffset = rand*2*pi;
        dialPos = 0;

        stage = stage + 1;
    end
    
    % delay phase
    if stage == 3 && (~eye_blink) && (norm(eye_pos-exp.startPos_deg) < pref.fixDeviation) && sample.ELTimestamp-eventTime(3) > stimExposure
        if isnan(eventTime(4)), eventTime(4) = sample.ELTimestamp; PsyFile('AddMessage','Delay phase'); end
        if use_mouse, SetMouse(xCentre, yCentre, win); end
        
        Screen('Flip', win);

        % prepare response phase
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPost, pref.innerAnnCol)
        for jj= 1:trialS(trial).nStim
            DrawAnnulus(win, screen, trialS(trial).xyPosition(:,jj), exp.annRad, exp.annWidth, exp.annCol);
        end
        DrawAnnulus(win, screen, trialS(trial).xyPosition(:,probe), exp.annRadCue, exp.annWidth, exp.annCol);
        stage = stage + 1;
    end
    
    % response phase orientation
    % TO DO:
    % fix mouse drift
    if stage == 4 && sample.ELTimestamp-eventTime(4) > trialS(trial).delay *1000
        if isnan(eventTime(5)) && use_knob == 1, hndl = PsychPowerMate('Open'); elseif isnan(eventTime(5)) && use_knob == 0, ShowCursor('Arrow'); end
        if isnan(eventTime(5)), eventTime(5) = sample.ELTimestamp; Screen('Flip', win); SetMouse(xCentre, yCentre, win); respOn = 1;  PsyFile('AddMessage','Response phase'); end
        
        [mx, my, buttons] = GetMouse(win);
        mousePos = PixelsToPos([mx-trialS(trial).targetPix(1); my-trialS(trial).targetPix(2)], screen,1);
        [mouseTheta, ~] = cart2pol(mousePos(1), mousePos(2));
        if use_knob, [buttonKnob, dialPos] = PsychPowerMate('Get', hndl); end % 3.8298
        if isnan(firstTouch) && (dialPos ~= 0 || mx ~= xCentre)
            firstTouch = GetSecs - eventTime(5); end
        
        % draw response bar
        if ~isnan(firstTouch) && use_knob
            tempResp = randOffset + (mod(dialPos,360))*pi/180;
            if stimulus == 1
                DrawGabor(win, screen, trialS(trial).initialRects(:,locPicked(probe)), gabor.amp, gabor.wavel, tempResp+pi/2, gabor.phase, gabor.sd, fieldSize(1))
            elseif stimulus == 2
                DrawTextBar(win, screen, trialS(trial).xyPosition(:,probe), tempResp+pi/2, exp.barWidth, exp.barHeight, exp.barColour);
            end
        elseif  ~isnan(firstTouch) && ~use_knob
            tempResp = mouseTheta + pi/2;
            if stimulus == 1
                DrawGabor(win, screen, trialS(trial).initialRects(:,locPicked(probe)), gabor.amp, gabor.wavel, tempResp, gabor.phase, gabor.sd, fieldSize(1))
            elseif stimulus == 2
                DrawTextBar(win, screen, trialS(trial).xyPosition(:,probe), tempResp, exp.barWidth, exp.barHeight, exp.barColour);
            end
        end
        % draw fixation & placeholders
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPost, pref.innerAnnCol)
        for jj= 1:trialS(trial).nStim
            DrawAnnulus(win, screen, trialS(trial).xyPosition(:,jj), exp.annRad, exp.annWidth, exp.annCol);
        end
        DrawAnnulus(win, screen, trialS(trial).xyPosition(:,probe), exp.annRadCue, exp.annWidth, exp.annCol);

        if use_knob == 1
            if buttonKnob(1)
                latency_finish = sample.ELTimestamp - eventTime(5);
                trialS(trial).response =  tempResp+pi/2; %radians
                arc_rect = [trialS(trial).initialRects(1:2,locPicked(probe));trialS(trial).initialRects(3:4,locPicked(probe))];
                stage = stage + 1;
                %DrawTextBar(win, screen, trialS(trial).xyPosition(:,probe), trialS(trial).response, exp.barWidth, exp.barHeight, exp.barColour.*2);
                pause(1);
            end
        elseif use_knob == 0
            if buttons(1)
                latency_finish = sample.ELTimestamp - eventTime(5);
                trialS(trial).response =  tempResp;
                arc_rect = [trialS(trial).initialRects(1:2,locPicked(probe));trialS(trial).initialRects(3:4,locPicked(probe))];
                stage = stage + 1;
                %DrawTextBar(win, screen, trialS(trial).xyPosition(:,probe), trialS(trial).response, exp.barWidth, exp.barHeight, exp.barColour.*2);
                pause(1);
            end
        end
        Screen('Flip', win);
    end
    
    %confidence level
    if stage == 5
        if isnan(eventTime(6)), eventTime(6) = sample.ELTimestamp; ShowCursor('Arrow'); PsyFile('AddMessage','Confidence estimate'); end
        % firstTouchConf
        %probe mouse position
        [mx, my, buttons] = GetMouse(win);
        mousePos = PixelsToPos([mx-trialS(trial).targetPix(1); my-trialS(trial).targetPix(2)], screen,1);
        [mouseTheta, ~] = cart2pol(mousePos(1), mousePos(2));

        % draw fixation & placeholders
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPost, pref.innerAnnCol)
        for jj= 1:trialS(trial).nStim
            DrawAnnulus(win, screen, trialS(trial).xyPosition(:,jj), exp.annRad, exp.annWidth, exp.annCol);
        end
        DrawAnnulus(win, screen, trialS(trial).xyPosition(:,probe), exp.annRadCue, exp.annWidth, exp.annCol);

        if use_knob, [buttonKnob, dialPos] = PsychPowerMate('Get', hndl); end
        if isnan(firstTouchConf) && (dialPos ~= 0 || mx ~= xCentre)
            firstTouchConf = GetSecs - eventTime(6); end

        if isnan(firstTouchConf)
            DrawTextBar(win, screen, trialS(trial).xyPosition(:,probe), trialS(trial).response, exp.barWidth, exp.barHeight, exp.barColour.*2);
        else
    		%draw arc
    		distance = abs(wrap((mouseTheta + pi/2 - trialS(trial).response),pi/2)); %value coreesponds to half of the arc degree
            startArcAngle = rad2deg(wrap((trialS(trial).response - distance)));
            sizeArcAngle = rad2deg(distance*2);

    		Screen('FillArc', win, exp.barColour, arc_rect, startArcAngle, sizeArcAngle);
            Screen('FillArc', win, exp.barColour, arc_rect, 180+startArcAngle, sizeArcAngle);
        end
		
		if use_knob == 1
            if buttonKnob(1)
				stage = stage +1;
			end
		elseif use_knob == 0
			if buttons(1)
				stage = stage + 1;
			end
		end
        Screen('Flip', win);
	end

	%display the progress bar - true/ false, at which stages below
    displaystages = [6 6];
    display_progress_bar = true;
    progress_bar_col = [0 0 0];
    progress_bar_width = 5;
    if stage >= displaystages(1) && stage <=displaystages(2) && display_progress_bar
        percentage_done = trial / length(trialS);
        right = round(screen.res(1) * percentage_done);
        Screen('FillRect',win,progress_bar_col,[0,0,right,progress_bar_width]);
    end

	%score
	if stage == 6
		eulers = 2.7182818284590452353602874713527;
		if abs(wrap(trialS(trial).response*2 - trialS(trial).orientation(probe)*2) ) <= distance %if true orientation is within arc
			%score = round((180-sizeArcAngle)/180*100);
            score = round(100*eulers.^(-0.0295*sizeArcAngle)); %score range of 0-100 points inversely correlated to arc size
		else
			score = 0;
		end
		total_score = total_score + score;
		total_score_display = sprintf('Total Score: %d', total_score);
		score_display = sprintf('Score: %d', score);
		DrawText(win, screen, [-1.0; -2.0], score_display, [0 0 0], 15);
		DrawText(win, screen, [-1.0; -1.0], total_score_display, [0 0 0], 15);
		trialS(trial).score =  score;
        Screen('Flip', win);
		pause(1);

		stage = stage +1;
		
	end

    %saving the data
    if stage == 7
        if ~use_mouse, HideCursor; end
        if use_mouse, SetMouse(xCentre, yCentre, win); end
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPost, pref.innerAnnCol)
        Screen('Flip', win);
        eventTime(7) = sample.ELTimestamp;
        
        % save trial data
        out.trial(trial,:)              = trialS(trial).trialInd;           % trial num
        out.nStim(trial)                = trialS(trial).nStim;              % number of stim
        out.delay(trial)                = trialS(trial).delay;              % memory delay time
        out.orientation{trial,:}        = trialS(trial).orientation;        % presented directions
        out.probe(trial)                = probe;                            % probed direction
        out.response(trial)             = trialS(trial).response;           % response, last recorded movement
        out.target(trial)               = trialS(trial).orientation(probe);
        out.recallError(trial)          = wrap(trialS(trial).response*2 - trialS(trial).orientation(probe)*2);    % recall error when response is finished
        out.timeStart(trial)            = firstTouch;                       % response started, first touch detected
        out.timeFinish(trial)           = latency_finish;                   % time when response finished (no movement)
        out.confidence(trial)           = distance;
        out.confidenceArcLength(trial)  = sizeArcAngle;
        out.locAngles{trial}            = trialS(trial).xyPosition;         % locations of stimuli
        out.locPicked{trial}            = locPicked;
        out.eventTimes(trial,:)         = eventTime';                      % eye tracker times
        out.ExpMat{1}                   = expMat;
        out.block(1)                    = block;
        out.date(1)                     = datetime;
        out.priorType                   = priorType;
        out.score(trial)				= trialS(trial).score;

        [~, filename_noext, ~] = fileparts(FILENAME);
        save([PATHNAME filename_noext '_data.mat'],'trialS','out','screen','pref','exp');
        
        % inter-trial delay
        pause(0.5);
        stage = 0; 
        trial = trial+1;
        initTrial = 1;
	end
    
    if stage >= 1 && stage <=3 && ~respOn && eye_blink && norm(eye_pos-exp.startPos_deg) > pref.fixDeviation
        Screen('FillRect', win, [0.5 0.5 0.5]);
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPre, pref.innerAnnCol)
        DrawText(win, screen, [-5.55; -1.3], 'Don''t look away from the fixation dot', [0 0 0], 20);
        Screen('Flip', win);
        trialS(end+1) = trialS(trial);
        trialS(end).trialInd = length(trialS);
        pause(1);
        stage = 0;
        initTrial = 1;
    end 
    if stage == 4 && ~respOn && norm(eye_pos-exp.startPos_deg) > pref.fixDeviation && (GetSecs - fixTime) >= pref.fixDeviationTime
        Screen('FillRect', win, [0.5 0.5 0.5]);
        DrawCircle(win, screen,  exp.startPos_deg, pref.outerAnnRad, pref.outerAnnCol)
        DrawCircle(win, screen,  exp.startPos_deg, pref.innerAnnRadPre, pref.innerAnnCol)
        DrawText(win, screen, [-5.55; -1.3], 'Don''t look away from the fixation dot', [0 0 0], 20);
        Screen('Flip', win);
        trialS(end+1) = trialS(trial);
        trialS(end).trialInd = length(trialS);
        pause(1);
        stage = 0;
        initTrial = 1;
        trial = trial + 1;
    end
end

%% clean up on exit
Eyelink('StopLiveRecording');
PsyFile('AddMessage',['BLOCK END' ' [' num2str(GetSecs*1000) ']']);
Eyelink('SetRunningFileStream',0);  % close filestream
PsyFile('Close')
Screen('LoadNormalizedGammaTable', win, gammaTable_orig);

end

%%
function px = PosToPixels (pos, screen, uncentred)
if nargin<3, uncentred = 0; end
px = round((pos./repmat(screen.size,1,size(pos,2)) + 0.5*(uncentred==0)).*repmat(screen.res,1,size(pos,2)));
end

%%

function pos = PixelsToPos (px, screen, uncentred)
if nargin<3, uncentred = 0; end
pos = (px./repmat(screen.res,1,size(px,2)) - 0.5*(uncentred==0)).*repmat(screen.size,1,size(px,2));
end

%%

function DrawBar (window, screen, centre, w, h, angle, colour)
posVector = [-w/2 -h/2; w/2 -h/2; w/2 h/2; -w/2 h/2; -w/2 -h/2]';
posVector = rotmtx(-angle)*posVector;
posVector = posVector + repmat(centre,1,5);
Screen('FillPoly', window, colour, PosToPixels(posVector, screen)', 1);
end

%%
function DrawArrow(window, screen, centre, w, h, angle, colour)
posVector = [0 -h*.9; -w*1.5 -h/2; -w/2 -h/2; -w/2 h/2; w/2 h/2; w/2 -h/2; w*1.5 -h/2; 0 -h*.9]';
posVector = rotmtx(-angle)*posVector;
posVector = posVector + repmat(centre,1,8);
Screen('FillPoly', window, colour, PosToPixels(posVector, screen)', 1);
end

%%
function DrawCrosshair (window, screen, centre, width, size, colour)
DrawBar(window, screen, centre, width, size, 0, colour);
DrawBar(window, screen, centre, width, size, pi/2, colour);
end

%%

function DrawCircle (window, screen, centre, r, colour)
posVector = [centre-repmat(r,2,1) centre+repmat(r,2,1)];
Screen('FillOval', window, colour, PosToPixels(posVector, screen));
end

%%
function DrawAnnulus (window, screen, centre, r, width, colour)
posVector = [centre-repmat(r,2,1) centre+repmat(r,2,1)];
w = PosToPixels ([width; width], screen, 1);
Screen('FrameOval', window, colour, PosToPixels(posVector, screen), mean(w));
end

%%

function DrawText (window, screen, topleft, text, colour, textsize)
Screen('TextSize', window, textsize);
px = PosToPixels(topleft, screen);
Screen('DrawText', window, text, px(1), px(2), colour);
end
%%
function DrawTextBar(window, screen, position, angle, widthDeg, heightDeg, color)
px = PosToPixels(widthDeg, screen,1);
py = PosToPixels(heightDeg, screen,1);
b = [];
stimTexture = Screen('MakeTexture', window, cat(3, ones(px(1),py(1)).*color(1), ones(px(1),py(1)).*color(2), ones(px(1),py(1)).*color(3)));

posPx = PosToPixels(position, screen)';
destRect = [-ceil(px(1)/2), -ceil(py(2)/2), floor(px(1)/2), floor(py(1)/2)] + [posPx, posPx];

Screen('DrawTexture', window, stimTexture, [], destRect, rad2deg(angle));
end

%%
function DrawGabor (window, screen, centre, amp, wavel, angle, phase, sd, imSize)

for ii=1:length(angle)
    wavel_px = PosToPixels(wavel,screen,1); % sf = 1/wavel
    sd_px = PosToPixels(sd,screen,1);
    %imSize = sd_px*10;
    imSize = [imSize, imSize];
    [px py] = meshgrid(1:imSize(1),1:imSize(2));
    x = px - imSize(1)/2;
    y = py - imSize(2)/2;
    z = sin(2*pi*(x*cos(angle(ii))/wavel_px(1) + y*sin(angle(ii))/wavel_px(2) ) + phase);
    z = z*amp.*exp(-( x.^2 / (2*sd_px(1)^2) + y.^2 / (2*sd_px(2)^2)));
    gabortex(ii) = Screen('MakeTexture',window,(z+1)/2);
end
Screen('DrawTextures', window, gabortex,[], centre)
end

%%
function k = repelems(values,reps)

k = [];
for i = 1:reps
    k = [k values];
end

end
%%
function M = cmean (X,Y)

if any(abs(X)>pi), error('Input values must be in radians, range -PI to PI'); return; end

if size(X,1)==1, X = X'; end

if nargin<2
    Y = ones(size(X));
else
    if size(Y,1)==1, Y = Y'; end
end

M = atan2(sum(sin(X).*Y),sum(cos(X).*Y));  %  or equivalently: angle(sum(exp(1i * x)))
end
%%
function [x, y] = RectCenterCustom(r)
%   [x,y] = RectCenter(rect);
%
%	RectCenter returns the integer x,y point closest to the center of a rect.
%
%	See also PsychRects, CenterRectOnPoint.

%	9/13/99	Allen Ingling wrote it.
if nargout~=2
    error('Usage: [x, y] = RectCenter(rect);');
end

if PsychNumel(r) == 4
    % Single rect:
    x = round(0.5*(r(1)+r(3)));
    y = round(0.5*(r(2)+r(4)));
else
    % Multi-rect array:
    if (size(r, 2)==4) && (size(r,1)>=1)
        % Multi-row array with one 4-comp. rect per row:
        x = round(0.5*(r(:,1)+r(:,3)));
        y = round(0.5*(r(:,2)+r(:,4)));
    else
        % Something weird and unknown:
        error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
    end
end
end
%%