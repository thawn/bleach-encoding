classdef BleachEvaluationClass < BleachEncodingSim
    %BleachEvaluationClass Evaluate bleach encoded microtubules.
    %
    %See also: BleachEncodingSim
    properties
        ResolutionFactor = 3; %Factor by which the resolution of the average linescan is increased
        PlotRows = 5; %Number of rows for the subplots
        PlotCols = 3; %Number of columns for the subplots
        OrigPos = 1:6; %Position of the original microtubule image
        SuperResPos = 7:9; %Position of the straightened superresolution microtubule image
        SingleLineScanPos = []; %Position of the line scan of the first frame
        ManyLineScansPos = []; %Position of the aligned line scans of all frames
        AveragedLineScanPos = [10 13]; %Position of the full line scan
        CroppedLineScanPos = [11 14]; %Position of the cropped line scan
        FourierPos = [12 15]; %Position of the cropped line scan
        CropLineScan = [0 0]; %Crop the linescan [left right]. Helpful in case findBleachedRegion finds the ends of the microtubule instead of the bleached region.
        Signal; %Signal level of the raw data
        LinescanYLim; %Y-axis limits for line scans
        FourierYLim; %Y-axis limits for fourier plots
        FrameLineScans; %aligned line scans for individual frames
        LineScanX; %X-axis for line scans
        FoundBleachedRegion = false;
        LeftOffset = 0; %Left offset for the bleached region
        RightOffset = 0; %Right offset for the bleached region
        FindRegion = true; %whether or not we want to find a bleached region
        CandidateP = 4:10; %incdices of the wavelenghts in the fft plot that might have been bleached
        Titles = true %turn titles on or off
        FourierRes = struct(); %store results of fourier analysis in case we want to save them to a file for external analysis
    end
    methods
        function BEC = BleachEvaluationClass(varargin)
            BEC@BleachEncodingSim(varargin{:});
        end
        
        
        function processDir(BEC, Dir, varargin)
            %processDir Evaluate all files in a given directory.
            
            p = inputParser;
            p.addParameter('Extension', '.nd2', @ischar);
            p.addParameter('LeftOffset', 0, @isnumeric);
            p.addParameter('RightOffset', 0, @isnumeric);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            Tmp = [fieldnames(p.Unmatched),struct2cell(p.Unmatched)];
            PassThrough = reshape(Tmp',[],1)';
            FilePattern = fullfile(Dir, ['*' p.Results.Extension]);
            ImportedFiles = dir(FilePattern);
            for k = 1 : length(ImportedFiles)
                BEC.FoundBleachedRegion = false;
                BEC.LeftOffset = p.Results.LeftOffset;
                BEC.RightOffset = p.Results.RightOffset;
                try
                    [BEC, SuperResImage, OrigImage, Line] = BEC.importLineScanFile(fullfile(ImportedFiles(k).folder, ...
                        ImportedFiles(k).name), PassThrough{:});
                    if isempty(SuperResImage)
                        BEC.formatPlotsOnly();
                    end
                    Format = BEC.createPlots('SuperResImage', SuperResImage, 'OrigImage', OrigImage, 'Line', Line, ...
                        'name',sprintf('%s Analysis', BEC.DataFileName),'NumberTitle','off');
                    FileName = fullfile(Dir, sprintf('%s_l%d_r%d_%s', ImportedFiles(k).name, BEC.LeftOffset, BEC.RightOffset, Format));
                    saveas(BEC.Fig, [FileName '.pdf'], 'pdf');
                    BEC.saveData(FileName, 'SuperResImage', SuperResImage, 'OrigImage', OrigImage, 'Line', Line);
                    close(BEC.Fig);
                catch e
                    fprintf('Error in %s:\n',BEC.DataFileName);
                    fprintf('%s\n', e.getReport)
                end
            end
        end
        
        
        function Format = createPlots(BEC, varargin)
            p = inputParser;
            p.addParameter('SuperResImage', [], @isnumeric);
            p.addParameter('OrigImage', [], @isnumeric);
            p.addParameter('Line', [], @isnumeric);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            Tmp = [fieldnames(p.Unmatched),struct2cell(p.Unmatched)];
            PassThrough = reshape(Tmp',[],1)';
            BEC.createFigure(PassThrough{:});
            Format = '';
            if ~isempty(BEC.FrameLineScans)
                if ~isempty(BEC.SingleLineScanPos)
                    BEC.createSingleLineScanSubPlot();
                    Format = [Format 'S'];
                end
                if ~isempty(BEC.ManyLineScansPos)
                    BEC.createManyLineScansSubPlot();
                    Format = [Format 'M'];
                end
            end
            if ~isempty(BEC.AveragedLineScanPos)
                BEC.createFullLinescanSubPlot();
                Format = [Format 'A'];
            end
            if ~isempty(p.Results.OrigImage) && ~isempty(BEC.OrigPos)
                BEC.createOriginalImageSubPlot(p.Results.OrigImage, p.Results.Line);
                Format = [Format 'O'];
            end
            if ~isempty(BEC.SuperResPos) && ~isempty(p.Results.SuperResImage)
                BEC.createSuperresImageSubPlot(p.Results.SuperResImage);
                Format = [Format 'R'];
            end
            if ~isempty(BEC.CroppedLineScanPos)
                BEC.createCroppedLinescanSubPlot();
                Format = [Format 'C'];
            end
            if ~isempty(BEC.FourierPos)
                BEC.createFFTSubPlot();
                Format = [Format 'F'];
            end
        end
        
        
        function BEC = formatPlotsOnly(BEC, Direction)
            if nargin <2
                Direction = 'horizontal';
            end
            if strcmpi(Direction, 'horizontal')
                BEC.Aspect = 5;
                BEC.PlotRows = 1;
                BEC.PlotCols = 3;
            else
                BEC.Aspect = 0.4;
                BEC.PlotRows = 3;
                BEC.PlotCols = 1;
            end
            BEC.AveragedLineScanPos = 1;
            BEC.CroppedLineScanPos = 2;
            BEC.FourierPos = 3;
        end
        
        
        function BEC = formatThreeColumn(BEC)
            BEC.Aspect = 2;
            BEC.PlotRows = 5;
            BEC.PlotCols = 3;
            BEC.OrigPos = 1:6;
            BEC.SingleLineScanPos = 7:9;
            BEC.ManyLineScansPos = [];
            BEC.SuperResPos = [];
            BEC.AveragedLineScanPos = [10 13];
            BEC.CroppedLineScanPos = [11 14];
            BEC.FourierPos = [12 15];
        end
        
        
        function BEC = formatWorkflow(BEC)
            BEC.Aspect = 1/3;
            BEC.PlotRows = 6;
            BEC.PlotCols = 1;
            BEC.OrigPos = 1;
            BEC.SingleLineScanPos = 2;
            BEC.ManyLineScansPos = 3;
            BEC.SuperResPos = [];
            BEC.AveragedLineScanPos = 4;
            BEC.CroppedLineScanPos = 5;
            BEC.FourierPos = 6;
        end
        
        
        function BEC = formatMTOnly(BEC)
            BEC.Aspect = 3;
            BEC.PlotRows = 2;
            BEC.PlotCols = 1;
            BEC.OrigPos = 1;
            BEC.SingleLineScanPos = [];
            BEC.ManyLineScansPos = [];
            BEC.SuperResPos = 2;
            BEC.AveragedLineScanPos = [];
            BEC.CroppedLineScanPos = [];
            BEC.FourierPos = [];
        end
        
        
        function BEC = formatMTLinescan(BEC)
            BEC.Aspect = 3/8;
            BEC.PlotRows = 4;
            BEC.PlotCols = 1;
            BEC.OrigPos = 1;
            BEC.SingleLineScanPos = 2;
            BEC.ManyLineScansPos = 3;
            BEC.SuperResPos = [];
            BEC.AveragedLineScanPos = 4;
            BEC.CroppedLineScanPos = [];
            BEC.FourierPos = [];
        end
        
        
        function BEC = formatCompact(BEC)
            BEC.Aspect = 3;
            BEC.PlotRows = 3;
            BEC.PlotCols = 5;
            BEC.OrigPos = 1:3;
            BEC.SingleLineScanPos = [];
            BEC.ManyLineScansPos = [];
            BEC.SuperResPos = [];
            BEC.AveragedLineScanPos = [6, 7, 8];
            BEC.CroppedLineScanPos = [11, 12, 13];
            BEC.FourierPos = [4 5 9 10 14 15];
        end
        
        
        function BEC = formatLinescanFourier(BEC)
            BEC.Aspect = 4/5;
            BEC.PlotRows = 3;
            BEC.PlotCols = 1;
            BEC.OrigPos = [];
            BEC.SingleLineScanPos = [];
            BEC.ManyLineScansPos = [];
            BEC.SuperResPos = [];
            BEC.AveragedLineScanPos = [];
            BEC.CroppedLineScanPos = 1;
            BEC.FourierPos = 2:3;
        end
        
        
        function BEC = formatLinescanFourierSide(BEC)
            BEC.Aspect = 4;
            BEC.PlotRows = 1;
            BEC.PlotCols = 9;
            BEC.OrigPos = [];
            BEC.SingleLineScanPos = [];
            BEC.ManyLineScansPos = [];
            BEC.SuperResPos = [];
            BEC.AveragedLineScanPos = [];
            BEC.CroppedLineScanPos = 1:6;
            BEC.FourierPos = 7:9;
        end
        
        
        function [BEC, SuperResImage, OrigImage, Line] = importLineScanFile(BEC, StackName, varargin)
            %importLineScanFile import a linescan from a file
            %
            %see also: processDir, BleachEncodingSim.importLineScan
            p = inputParser;
            p.addParameter('PixelSize', 160, @isnumeric);
            p.addParameter('CropLineScan', [0 0], @(x) isnumeric(x) && length(x) == 2);
            p.addParameter('SuperResImageHeight', 31, @isscalar);
            p.addParameter('LinescanYLim', [], @(x) isnumeric(x) && (isempty(x) || length(x) == 2));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            Tmp = [fieldnames(p.Unmatched),struct2cell(p.Unmatched)];
            PassThrough = [{'ResolutionFactor', BEC.ResolutionFactor}, reshape(Tmp',[],1)'];
            SuperResImage = [];
            OrigImage = [];
            Line = [];
            BEC.FrameLineScans = [];
            BEC.DataFileName = StackName;
            [Path, Name, Ext] = fileparts(StackName);
            if strcmp(Ext, '.csv')
                PixSize = BEC.PixelSize;
                LineScan = readmatrix(StackName);
                LineScan = LineScan(:, 2);
                BEC.Signal = max(findSignal(LineScan));
            else
                FiestaFiles = dir(fullfile(Path, [Name '*.mat']));
                FiestaFile = fullfile(Path, FiestaFiles(end).name);
                [LineScan, PixSize, BEC.FrameLineScans, SignalLevels] = superResLineScan(FiestaFile, 1, StackName, PassThrough{:});
                BEC.Signal = max(SignalLevels);
                [~, ~, ~, ~, SuperResImage, OrigImage, Line] = superResLineScan(FiestaFile, 1, StackName, PassThrough{:}, 'ScanSize', round((p.Results.SuperResImageHeight - 1) / 2));
                LineScan = LineScan';
            end
            if isempty (p.Results.LinescanYLim)
              BEC.LinescanYLim = [0 BEC.Signal * 1.5];
            else
              BEC.LinescanYLim = p.Results.LinescanYLim;
            end
            BEC.LineScanX = (0:length(LineScan) - 1) .* PixSize / 1000;
            BEC.CropLineScan = p.Results.CropLineScan;
            LineScan = LineScan(1 + BEC.CropLineScan(1):end - BEC.CropLineScan(2));
            BEC.importLineScan(LineScan, PixSize);
            BEC.IntensityProfileX = BEC.IntensityProfileX + (BEC.CropLineScan(1) * BEC.PixelSize / 1000); %shift x axis to compensate for cropping
        end
        
        
        function saveData(BEC, FileName, varargin)
            p = inputParser;
            p.addParameter('SuperResImage', [], @isnumeric);
            p.addParameter('OrigImage', [], @isnumeric);
            p.addParameter('Line', [], @isnumeric);
            p.parse(varargin{:});
            res = struct();
            res.DataFileName = BEC.DataFileName;
            res.OrigImage = p.Results.OrigImage;
            res.SuperResImage = p.Results.SuperResImage;
            res.MTPosLine = p.Results.Line;
            res.PixelSize = BEC.PixelSize;
            res.ResolutionFactor = BEC.ResolutionFactor;
            res.IntensityProfile = BEC.IntensityProfile;
            res.MicrotubuleLength = BEC.MicrotubuleLength;
            res.IntensityProfileX = BEC.IntensityProfileX;
            res.FrameLineScans = BEC.FrameLineScans;
            res.LinescanYLim = BEC.LinescanYLim;
            res.LineScanX = BEC.LineScanX;
            res.FourierRes = BEC.FourierRes;
            save([FileName '.mat'], 'res');
        end
        
        
        function BEC = createFullLinescanSubPlot(BEC)
            %createFullLinescanSubPlot plot the full linescan including the
            %fit used to find the bleached region.
            %
            %see also: processDir, findBleachedRegion, BleachEncodingSim.formatAxes
            
            Ax = subplot(BEC.PlotRows, BEC.PlotCols, BEC.AveragedLineScanPos);
            axes(Ax);
            plot(BEC.IntensityProfileX, BEC.IntensityProfile,'-');
            if BEC.FindRegion
              [BEC, Fitresult] = BEC.findBleachedRegion();
              hold on;
              plot(Fitresult);
              legend off;
            end
            if BEC.Titles
                title (Ax, 'Fitted linescan');
            end
            BEC.formatAxes(Ax);
            ylim(Ax, BEC.LinescanYLim);
            xlim(Ax, [0 max(BEC.LineScanX)]);
        end
        
        
        function createOriginalImageSubPlot(BEC, OrigImage, Line)
            %createOriginalImageSubPlot display the Original image of the
            %microtubule.
            %
            %see also: processDir, plotImage
            
            Ax = subplot(BEC.PlotRows, BEC.PlotCols, BEC.OrigPos, 'Parent', BEC.Fig);
            ImSize = size(OrigImage);
            if ImSize(1) > ImSize(2)
                OrigImage = OrigImage';
                Line = fliplr(Line);
            end
            BEC.plotImage(OrigImage, BEC.PixelSize * BEC.ResolutionFactor, Ax);
            hold on;
            plot(Line(:,1), Line(:,2));
            if BEC.Titles
                title (Ax, 'Fluorescence micrograph of microtubule');
            end
        end
        
        
        function createSuperresImageSubPlot(BEC, SuperResImage)
            %display the superresolution image of the
            %microtubule.
            %
            %see also: processDir, plotImage
            
            Ax = subplot(BEC.PlotRows, BEC.PlotCols, BEC.SuperResPos, 'Parent', BEC.Fig);
            BEC.plotImage(SuperResImage, BEC.PixelSize, Ax, 1/BEC.ResolutionFactor);
            Ax.DataAspectRatio = [BEC.ResolutionFactor, 1, 1];
            if BEC.Titles
                title(Ax, 'Interpolated, Straightened and Averaged Microtubule Image');
            end
        end
        
        
        function createCroppedLinescanSubPlot(BEC)
            %createCroppedLinescanSubPlot plot the cropped linescan of the
            %bleached region.
            %
            %see also: processDir, BleachEncodingSim.formatAxes
            
            if BEC.FindRegion && ~BEC.FoundBleachedRegion
                BEC.findBleachedRegion();
            end
            Ax = subplot(BEC.PlotRows, BEC.PlotCols, BEC.CroppedLineScanPos, 'Parent', BEC.Fig);
            plot(Ax, BEC.IntensityProfileX, BEC.IntensityProfile);
            if BEC.Titles
                title (Ax, 'Intensity profile of bleached region');
            end
            BEC.formatAxes(Ax, 'XLabel', sprintf('distance (%sm)', 181));
            xlim(Ax, [min(BEC.IntensityProfileX), max(BEC.IntensityProfileX)]);
        end
        
        
        function createSingleLineScanSubPlot(BEC)
            %createSingleLineScanSubPlot plot the linescan of the first
            %frame.
            %
            %see also: processDir, BleachEncodingSim.formatAxes
            
            Ax = subplot(BEC.PlotRows, BEC.PlotCols, BEC.SingleLineScanPos, 'Parent', BEC.Fig);
            plot(Ax, BEC.LineScanX, BEC.FrameLineScans(end, :));
            if BEC.Titles
                title (Ax, 'Intensity profile of microtubule in last frame');
            end
            BEC.formatAxes(Ax);
            ylim(Ax, BEC.LinescanYLim);
            xlim(Ax, [0 max(BEC.LineScanX)]);
        end
        
        
        function createManyLineScansSubPlot(BEC)
            %createSingleLineScanSubPlot plot the linescan of the first
            %frame.
            %
            %see also: processDir, BleachEncodingSim.formatAxes
            
            Ax = subplot(BEC.PlotRows, BEC.PlotCols, BEC.ManyLineScansPos, 'Parent', BEC.Fig);
            hold(Ax, 'on')
            for n = 1:size(BEC.FrameLineScans, 1)
                plot(Ax, BEC.LineScanX, BEC.FrameLineScans(n, :));
            end
            if BEC.Titles
                title (Ax, 'Aligned intensity profiles of microtubule in all frames');
            end
            BEC.formatAxes(Ax);
            ylim(Ax, BEC.LinescanYLim);
            xlim(Ax, [0 max(BEC.LineScanX)]);
        end
        
        
        function createFFTSubPlot(BEC)
            %createFFTSubPlot plot the fourier transformation.
            
            if BEC.FindRegion &&~BEC.FoundBleachedRegion
                BEC.findBleachedRegion();
            end
            Ax = subplot(BEC.PlotRows, BEC.PlotCols, BEC.FourierPos, 'Parent', BEC.Fig);
            BEC.plotFFTAnalysisWaveLength('Axes', Ax);
            %if ~isempty(BEC.FourierYLim)
            %    LineY = [0 max(BEC.FourierYLim)];
            %else
            %    LineY = [0 max(ylim(Ax))];
            %end
            %for n=1:length(BEC.ExpectedWavelenghts)
            %    line(Ax, [BEC.ExpectedWavelenghts(n) BEC.ExpectedWavelenghts(n)], LineY, 'Color','r','LineStyle',':');
            %end
            if BEC.Titles
                title (Ax, 'spatial periods detected by FFT');
            end
            if ~isempty(BEC.FourierYLim)
                ylim(Ax, BEC.FourierYLim);
            end
        end
        
        
        function BEC = plotFFTAnalysisWaveLength(BEC, varargin)
            %plotFFTAnalysis plot the fast furier transform analysis of the
            %intensity profile by wavelength.
            
            p = inputParser;
            p.addParameter('Axes', [], @(X) isempty(X) || ishandle(X));
            p.parse(varargin{:});
            
            if isempty(p.Results.Axes)
                [~, Ax] = BEC.createFigure();
            else
                Ax = p.Results.Axes;
            end
            if ~isempty(BEC.FrameLineScans)
                [Fourier, FourierCI, Significant, ControlLimit, WL] = analyseFFTSignificance(BEC);
                errorbar(Ax, WL, Fourier, Fourier - FourierCI(1, :), FourierCI(2, :) - Fourier,  'LineStyle', 'none', 'Marker', '*');
                hold(Ax, 'on');
                plot(Ax, WL(BEC.CandidateP), ControlLimit, 'LineStyle', '-', 'Marker', '.')
                plot(Ax,WL(BEC.CandidateP(Significant)),FourierCI(2,BEC.CandidateP(Significant))+20,'LineStyle', 'none', 'Marker', '*');
                BEC.FourierRes = struct('WaveLengths', WL, 'ControlLimit', ControlLimit, 'Fourier', Fourier, 'FourierCI', FourierCI, 'Significant', Significant, 'CandidateP', BEC.CandidateP);
                hold(Ax, 'off');
            else
                [P1, MTLength] = BEC.fftAnalysis;
                [P1, WL] = BEC.fftWaveLength(P1, MTLength);
                plot(Ax, WL, P1, 'LineStyle', 'none', 'Marker', '*');
                BEC.FourierRes = struct('WaveLengths', WL, 'Fourier', P1);
            end
            BEC.formatAxes(Ax, 'XLabel', 'spatial period (nm)', 'YLabel', 'amplitude (a.u.)');
        end
        
        
        function [Fourier, FourierCI, Significant, ControlLimit, WL, Signal, SignalCI] = analyseFFTSignificance(BEC)
            [Fourier, FourierCI, WL] = BEC.bootciFFT(1000);
            ControlP = 1:length(WL);
            ControlP(BEC.CandidateP) = [];
            LimitFit = fit(WL(ControlP)', FourierCI(2, ControlP)','poly1');
            ControlLimit = LimitFit(WL(BEC.CandidateP)');
            Significant = ControlLimit' < FourierCI(1, BEC.CandidateP);
            Signal = Fourier(1, BEC.CandidateP) - ControlLimit;
            SignalCI = [FourierCI(1, BEC.CandidateP) - ControlLimit; FourierCI(2, BEC.CandidateP) - ControlLimit];
        end
        
        
        function [Fourier, FourierCI, WL] = bootciFFT(BEC, NBoot, Alpha)
            if nargin < 4
                Alpha = 0.05;
            end
            [Fourier, NSamples] = BEC.medianFFT(BEC.FrameLineScans);
            WL = BEC.calculateWaveLengths(size(Fourier, 1), NSamples);
            Discard = BEC.discardWaveLengths(WL);
            Results = zeros(NBoot, length(Fourier));
            NLines = size(BEC.FrameLineScans, 1);
            for n = 1:NBoot
                Rows = randsample(NLines, NLines, true);
                Results(n, :) = BEC.medianFFT(BEC.FrameLineScans(Rows,:));
            end
            FSem = std(Results, 0, 1);
            Bias = mean(bsxfun(@minus, Results, Fourier), 1);
            Za = norminv(Alpha/2);
            Lower = Fourier - Bias + FSem * Za;
            Upper = Fourier - Bias - FSem * Za;
            Fourier(Discard) = [];
            FourierCI = [Lower; Upper];
            FourierCI(:,Discard) = [];
            WL(Discard) = [];
        end
        
        
        function [BEC, Fitresult] = findBleachedRegion(BEC, varargin)
            %findBleachedRegion find the region of a linescan that has been
            %bleached.
            %
            %  See also FIT, CFIT, SFIT.
            
            p = inputParser;
            p.parse(varargin{:});
            % Set up fittype and options.
            Ft = fittype( 'smoothingspline' );
            Opts = fitoptions( 'Method', 'SmoothingSpline' );
            Opts.SmoothingParam = 0.05;
            Fitresult = fit( BEC.IntensityProfileX, BEC.IntensityProfile, Ft, Opts );
            
            %Calculate differentiate
            Fx = differentiate(Fitresult, BEC.IntensityProfileX);
            
            %find min and max
            SortedFx = sort(Fx);
            Found = false;
            LeftBleachborder = find(Fx == SortedFx(1));
            for n=1:length(SortedFx)
                RightBleachborder = find(Fx == SortedFx(end - n + 1));
                if LeftBleachborder < RightBleachborder
                    Found = true;
                    break;
                end
            end
            if ~Found
                RightBleachborder = find(Fx == SortedFx(end));
                for n=1:length(SortedFx)
                    LeftBleachborder = find(Fx == SortedFx(n));
                    if LeftBleachborder < RightBleachborder
                        Found = true;
                        break;
                    end
                end
                if ~Found
                    error('Could not find bleaching borders');
                end
            end
            OrigProfile = BEC.IntensityProfile;
            if BEC.LeftOffset == 0 && BEC.RightOffset == 0
                Cutoffs = combnk(-5:20,2);
                NTries = size(Cutoffs, 1);
                SignalToNoise = zeros(NTries, 1);
                for n = 1:NTries
                    BEC.IntensityProfile = OrigProfile(LeftBleachborder + Cutoffs(n, 1):RightBleachborder - Cutoffs(n, 2));
                    [P1, MTLength] = BEC.fftAnalysis;
                    FFT = BEC.fftWaveLength(P1, MTLength);
                    SignalToNoise(n) = max(FFT)/median(FFT);
                end
                [~, Best] = max(SignalToNoise);
                BEC.LeftOffset = Cutoffs(Best, 1);
                BEC.RightOffset = Cutoffs(Best, 2);
            end
            %cut
            BEC.IntensityProfile = OrigProfile(LeftBleachborder + BEC.LeftOffset:RightBleachborder - BEC.RightOffset);
            BEC.FrameLineScans = BEC.FrameLineScans(:,LeftBleachborder + BEC.LeftOffset + BEC.CropLineScan(1):RightBleachborder - BEC.RightOffset + BEC.CropLineScan(1));
            BEC.calculateIntensityProfileX;
            BEC.MicrotubuleLength = length(BEC.IntensityProfile) * BEC.PixelSize;
            
            %expected wavelenghts
            %if max(BEC.IntensityProfileX) < 19
            %    BEC.calculateExpectedWaveLengths;
            %elseif max(BEC.IntensityProfileX) > 19
            %    BEC.calculateExpectedWaveLengths(2);
            %end
            
            %BEC.ExpectedWavelenghts = round(BEC.ExpectedWavelenghts*1000);
            
            BEC.FoundBleachedRegion = true;
        end
        
        
        function plotSim(BEC, Wavecounts, Intensity, NSim)
            if nargin < 4
                NSim = 20;
            end
            if nargin < 3
                Intensity = 0.1;
            end
            BEC.simulateManyLineScans(NSim, Wavecounts, Intensity)
            BEC.FrameLineScans = zeros(NSim, length(BEC.getLineScan));
            for n = 1:NSim
                BEC.IntensityProfile = BEC.IntensityProfiles(n, :);
                BEC.FrameLineScans(n, :) = BEC.getLineScan;
            end
            BEC.LineScanX = (0:size(BEC.FrameLineScans, 2) - 1) .* BEC.PixelSize / 1000;
            BEC.LinescanYLim = [0 max(BEC.FrameLineScans(:))];
            BEC.IntensityProfile = nanmean(BEC.IntensityProfiles, 1);
            Format = BEC.createPlots;
            Wc = sprintf('%d_', Wavecounts);
            Wc = sprintf('wc[%s]', Wc(1:end-1));
            saveas(BEC.Fig, fullfile('results', filesep, sprintf('simulated_%s_%s.pdf', Wc, Format)), 'pdf');
            close(BEC.Fig);
        end
        
        
        function plotSigSim(BEC)
            BEC.FindRegion = false;
            BEC.CandidateP = 1;
            BEC.MinWaveNumber = 1;
            BEC.MinWaveLengthFactor = 1.5;
            BEC.plotSim(1)
            [~, ~, ~, ~, WL, ~] = BEC.analyseFFTSignificance;
            Signals = zeros(1,length(WL));
            SignalCI = zeros(2, length(WL));
            BEC.formatLinescanFourierSide;
            for n = 1:length(WL)
                BEC.CandidateP = n;
                BEC.plotSim(n - 1 + BEC.MinWaveNumber);
                [~, ~, ~, ~, ~, Signals(n), SignalCI(:, n)] = BEC.analyseFFTSignificance;
            end
            BEC.Aspect = 3/2;
            BEC.createFigure;
            Ax = gca;
            errorbar(Ax, WL, Signals, SignalCI(1,:) - Signals, SignalCI(2,:) - Signals);
            xlabel(Ax, 'spatial period (nm)')
            ylabel(Ax, 'signal above threshold (a.u.)')
            saveas(BEC.Fig, fullfile('results', filesep, sprintf('signal_mtl%d_res%d.pdf', BEC.MicrotubuleLength, BEC.BleachingResolution)), 'pdf');
        end
        
        
        function plotManySim(BEC)
            BEC.FindRegion = false;
            BEC.CandidateP = 1;
            BEC.MinWaveNumber = 7;
            BEC.MinWaveLengthFactor = 1.5;
            BEC.Aspect = 4/3;
            BEC.PlotRows = 1;
            BEC.PlotCols = 1;
            BEC.OrigPos = [];
            BEC.SingleLineScanPos = [];
            BEC.ManyLineScansPos = [];
            BEC.SuperResPos = [];
            BEC.AveragedLineScanPos = [];
            BEC.CroppedLineScanPos = 1;
            BEC.FourierPos = [];
            BEC.plotSim(1)
            [~, ~, ~, ~, WL, ~] = BEC.analyseFFTSignificance;
            for n = 1:length(WL)
                BEC.CandidateP = n;
                BEC.Aspect = 2;
                BEC.FigWidth = 10;
                BEC.CroppedLineScanPos = 1;
                BEC.FourierPos = [];                
                BEC.plotSim(n - 1 + BEC.MinWaveNumber);
                BEC.FigWidth = 5;
                BEC.Aspect = 1;
                BEC.CroppedLineScanPos = [];
                BEC.FourierPos = 1;                
                BEC.plotSim(n - 1 + BEC.MinWaveNumber);
            end
        end
        
        
        function plotBleachInt(BEC, Wavecounts)
            BEC.Aspect = 1/2;
            BEC.createFigure;
            for n=1:length(Wavecounts)
                Ax = subplot(3, 1, n);
                Ax.FontName = 'Arial';
                BEC.plotBleachingIntensityWC(Wavecounts(n), Ax)
            end
            Wc = sprintf('%d_', Wavecounts);
            Wc = sprintf('wc[%s]', Wc(1:end-1));
            saveas(BEC.Fig, fullfile('results', filesep, sprintf('bleach_int_wc%s_res%d.pdf', Wc, BEC.BleachingResolution)), 'pdf');
        end
        
        
        function plotSignals(BEC, Wavecounts)
            BEC.Aspect = 3/2;
            BEC.createFigure;
            for n=1:length(Wavecounts)
                Ax = subplot(1, 3, n);
                BEC.plotBleachingIntensityWC(Wavecounts(n), Ax)
            end
            Wc = sprintf('%d_', Wavecounts);
            Wc = sprintf('wc[%s]', Wc(1:end-1));
            saveas(BEC.Fig, fullfile('results', filesep, sprintf('bleach_int_wc%s_res%d.pdf', Wc, BEC.BleachingResolution)), 'pdf');
        end
        
        
    end
    methods (Static)
        
        
        function plotImage(Image, PixelSize, Parent, ScaleBarHFactor)
            if nargin < 4
                ScaleBarHFactor = 1;
            end
            axes(Parent);
            imshow(Image, [prctile(Image(:), 0.1) prctile(Image(:), 99.8)]);
            BarSize = ceil(size(Image, 2) * PixelSize / 50000) * 5;
            ScaleBarW = BarSize * 1000 / PixelSize;
            ScaleBarH = 200 / PixelSize * ScaleBarHFactor;
            ScaleBarX = size(Image, 2) - ScaleBarW - ScaleBarH;
            ScaleBarY = size(Image, 1) - 2 * ScaleBarH;
            rectangle('Position', [ScaleBarX,ScaleBarY , ScaleBarW, ScaleBarH], 'EdgeColor', 'none', 'FaceColor', 'w');
            text('Position', [ScaleBarX + ScaleBarW / 2, ScaleBarY - 2 * ScaleBarH], 'HorizontalAlignment', 'center', 'Color', 'w',...
                'String',sprintf('%g %cm', BarSize, 181), 'FontUnits', 'points', 'FontSize', 10);
        end
        
        
        function [Result, NSamples] = medianFFT(LineScans)
            Profile = nanmedian(LineScans);
            [Result, NSamples] = BleachEvaluationClass.performFFT(Profile);
        end
        
        
        function [Result, NSamples] = performFFT(Profile)
            NSamples = length(Profile) - 1;
            Fourier = fft(Profile, NSamples);
            Result = BleachEvaluationClass.processFFT(Fourier, NSamples);
        end
        
        
    end
end
