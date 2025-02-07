classdef BleachEncodingSim < handle
    %BleachEncodingSim Simulate bleaching information into microtubules.
    % Example:
    properties
        Microtubule = []; %Vector containing the number of fluorescent dye molecules per SimRes nm.
        MicrotubuleLength = 10000; %Length of the bleached region of microtubule in nm.
        MTOffset = 200; %Length of microtubule to each side of the bleached region in SimRes steps. Set to 0 to disable the offset.
        SimRes = 10; %Resolution of the simulation. The microtubule will be separated into bins of SimRes width.
        LabelingRatio = 0.25 %Ratio of tubulin dimers that are labeled with an average of 2 dye molecules each.
        DyeIntensity = 100; %Average photons emitted per dye molecule per frame.
        DyeLifetime = 0.1; %Lifetime of the dye molecules during bleaching in s.
        ImagingResolution = 400 %Optical resolution of the objective in nm.
        BleachingResolution = 600 %Optical resolution of the bleaching laser in nm.
        PixelSize = 177; %Pixel size used for read out.60x(226 at ov 1, 177 at ov2), 100x (105)
        IntensityProfile = []; %Vector containing the intensity profile.
        IntensityProfileX = []; %X axix for IntensityProfile.
        IntensityProfiles = []; %Matrix containing many intensity profiles
        PSFImage = []; %Vector containing the point spread function of the objective.
        PSFBleach = []; %Vector containing the point spread function of the bleaching laser.
        Seed = []; %Random seed for simulation
        ImportedData = false; %Whether or not the intensity data was imported from a file
        DataFileName = ''; %File name from where data is imported
        Fig; %Figure handle
        FigWidth = 12 %Width of figures
        Aspect = 2; %Aspect ratio (width/height) for the figure
        ExpectedWavelenghts = [];
        ExpectedWaveNumbers = [10 12 14 17];
        MinWaveNumber = 5 % Minimum wave number (or longest wavelength) that is analyzed by FFT.
        MinWaveLengthFactor = 1.8 %Factor by which the imaging resolution is multiplied to determine the minimum wavelength that is analyzed by FFT
    end
    methods
        function BES = BleachEncodingSim(varargin)  %encodes bleaching on MT
            if nargin>0
                BES.Seed = rng('shuffle');
                p=inputParser;
                
                p.KeepUnmatched=true;
                
                p.parse(varargin{:});
                
                Names=fieldnames(p.Unmatched);
                for i=1:length(Names)
                    if isempty(BES.findprop(Names{i}))
                        warning(strcat('BleachEncodingSim does not support the parameter: ',Names{i}));
                    else
                        BES.(Names{i})=p.Unmatched.(Names{i});
                    end
                end
            end
            if isempty(BES.Microtubule)
                BES.resetMicrotubule;
            else
                BES.Microtubule = uint16(round(BES.Microtubule));
                BES.MicrotubuleLength = length(BES.Microtubule) * BES.SimRes;
            end
            if isempty(BES.PSFImage)
                BES.PSFImage = BES.createPSF(BES.ImagingResolution);
            end
            if isempty(BES.PSFBleach)
                BES.PSFBleach = BES.createPSF(BES.BleachingResolution);
            end
            BES.calculateIntensityProfile;
        end
        
        
        function BES = resetMicrotubule(BES)
            %resetMicrotubule rese
            BES.Microtubule = uint16(round(rand(1,round(BES.MicrotubuleLength / BES.SimRes + 2 * BES.MTOffset)) * 13 * BES.LabelingRatio * 2)); % 13 is the number of protofilaments and 2 is the average number of dye molecules per labeled tubulin dimer
        end
        
        
        function BES = setBleachingResolution(BES, Res)
            BES.BleachingResolution = Res;
            BES.PSFBleach = BES.createPSF(BES.BleachingResolution);
        end
        
        
        function BES = setImagingResolution(BES, Res)
            BES.ImagingResolution = Res;
            BES.PSFImage = BES.createPSF(BES.ImagingResolution);
        end
        
        
        function PSF = createPSF(BES, Res)
            PSF = fspecial('gaussian',[1, round(Res * 3 / BES.SimRes)], Res / (BES.SimRes * 2 * sqrt(2 * log(2))));
            PSF = single(PSF / max(PSF));
        end
        
        
        function BES = calculateIntensityProfile(BES)
            %calculateIntensityProfile Calculate the intensity profile of the
            %microtubule.
            %
            % See also: BleachEncodingSim plotLineScan
            
            BES.IntensityProfile = conv(single(BES.Microtubule * BES.DyeIntensity), BES.PSFImage);
            BES.calculateIntensityProfileX;
            ShotNoise = randn(1, length(BES.IntensityProfile));
            BES.IntensityProfile = BES.IntensityProfile + (sqrt(BES.IntensityProfile) .* ShotNoise);
            BES.IntensityProfile = BES.IntensityProfile ./ max(BES.IntensityProfile) * 10;
        end
        
        function Profile = getLineScan(BES)
            Profile = resample(double(BES.removeTails(BES.IntensityProfile, false)), BES.SimRes, BES.PixelSize);
        end
        
        
        function Pattern = getEmptyBleachPattern(BES)
            %getEmptyBleachPattern returns a zero vector with the same size as
            %the microtubule vector
            %
            % See also: BleachEncodingSim getPatternWithSpacing
            
            Pattern = zeros(size(BES.Microtubule));
        end
        
        
        function Pattern = getPatternBarcode(BES, ~)
            %getPatternBarcode returns a vector containing a barcode encoding the
            %number specified by Code.
            %
            % Code: number to be encoded into the barcode.
            
            Pattern = BES.getEmptyBleachPattern;
            
        end
        
        
        function Pattern = getPatternWavelength(BES, Wavelength)
            %getPatternWavelength returns a vector containing a regular
            %bleaching pattern with spacing of Wavelength nm.
            %
            % Wavelength: double pattern wavelength in nanometer.
            %
            % See also: BleachEncodingSim getEmptyBleachPattern
            
            MultiBleachingDistance = 2; %0.75;
            
            BleachesPerMaximum = floor(Wavelength / (2 * MultiBleachingDistance * BES.BleachingResolution));
            if BleachesPerMaximum < 1
                warning('The wavelength is less than the bleaching resolution the bleaching will likely not work very well.');
                BleachesPerMaximum = 1;
            end
            NMaxima = round(BES.MicrotubuleLength / Wavelength);
            SimWL = Wavelength / BES.SimRes;
            Start = BES.MTOffset + 1 - round(length(BES.PSFBleach) / 2); %round(SimWL / 2);
            Maxima = round(Start:SimWL:NMaxima * SimWL + Start);
            if BleachesPerMaximum > 1
                Maxima1 = Maxima;
                for n = 1:BleachesPerMaximum - 1
                    Maxima1 = [Maxima1 Maxima + round( n * MultiBleachingDistance * BES.BleachingResolution / BES.SimRes)]; %#ok<AGROW>
                end
                Maxima = Maxima1;
            end
            Pattern = zeros(1, max(Maxima));
            Pattern(Maxima) = 1;
        end
        
        
        function Pattern = getPatternWaveCount(BES, Wavecount)
            Wavelength = BES.MicrotubuleLength/Wavecount;
            Pattern = BES.getPatternWavelength(Wavelength);
        end
        
        
        function Profile = removeTails(BES,Profile, Bleaching)
            %removeTails removes the tails of intensity profiles such that the
            %intensity profile has the same length as the Microtubule vector.
            
            if nargin < 3
                Bleaching = true;
            end
            MTLength = length(BES.Microtubule) - 2 * BES.MTOffset;
            Overhead = length(Profile) - MTLength;
            if BES.MTOffset == 0
                if Bleaching
                    Start = round(length(BES.PSFBleach) / 2);
                else
                    Start = round(length(BES.PSFImage) / 2) + BES.MTOffset;
                end
            else
                if Bleaching
                    Start = BES.MTOffset;
                else
                    Start = round(length(BES.PSFImage) / 2) + BES.MTOffset;
                end
            end
            if Overhead > Start
                Profile = Profile(Start + 1: Start + MTLength);
            elseif Overhead > 0
                Profile = Profile(floor(Overhead / 2) + 1: end - ceil(Overhead / 2));
            else % Overhead < 0
                Profile = [zeros(1, abs(floor(Overhead / 2))), Profile, zeros(1, abs(ceil(Overhead / 2)))];
            end
        end
        
        
        function BleachIntensityProfile = createBleachIntensityProfile(BES, Pattern)
            %createBleachIntensityProfile creates the intensity profile of the
            %bleaching laser based on the bleaching resolution from an idealized Pattern.
            
            BleachIntensityProfile = conv(single(Pattern), BES.PSFBleach);
        end
        
        
        function BES = bleach(BES, Pattern, Duration)
            %bleach Applies a bleaching pattern defined in Pattern.
            %
            % Calculates how the pattern is affected by the point spread function
            % and applies it to the microtubule.
            %
            % Pattern: Vector containing zero where we don't want to bleach and 1 where we want to bleach
            % Duration: Number in s. How long we want to bleach.
            % Intensity: Number between 0 and 1 what fraction of bleaching
            % intensity we want to apply.
            
            MTLength = length(BES.Microtubule);
            BleachIntensity = BES.createBleachIntensityProfile(Pattern);
            if length(BleachIntensity) > MTLength
                BleachIntensity = BES.removeTails(BleachIntensity);
            end
            DecayProb = exp(-Duration / BES.DyeLifetime);
            NumDyes = max(BES.Microtubule);
            BleachIntensityProfile = zeros(NumDyes, MTLength);
            BleachIntensityProfile(:,1:length(BleachIntensity)) = repmat(BleachIntensity,NumDyes,1);
            Bleach  = rand(NumDyes, MTLength);
            Bleach = Bleach .* BleachIntensityProfile;
            for n=1:MTLength
                BES.Microtubule(n) = BES.Microtubule(n) - sum(Bleach(1:BES.Microtubule(n),n) > DecayProb);
            end
            BES.calculateIntensityProfile;
        end
        
        
        function BES = bleachWavelength(BES, Wavelength, Duration)
            %bleachWavelength create a pattern with a certain wavelength (in nm)
            %and bleach it into the microtubule.
            %
            %  Parameters:
            %    Wavelength: wavelength of the pattern in nm. If this
            %      is a vector, each element of the vector is considered as one
            %      wavelength to be bleached into the microtubule.
            %    Duration: duration of bleaching in seconds.
            
            for n = 1:length(Wavelength)
                if Wavelength(n) > 0
                    Pattern = BES.getPatternWavelength(Wavelength);
                    BES.bleach(Pattern, Duration);
                end
            end
        end
        
        
        function BES = bleachWaveCounts(BES, Wavecount, Duration)
            %bleachWavelength create a pattern with a wavelnength that fits a
            %certain number of maxima into the microtubule length, then bleach
            %this pattern into the microtubule.
            %
            %  Parameters:
            %    Wavecount: number of maxima to fit into the microtubule. If this
            %      is a vector, each element of the vector is considered as one
            %      wavecount to be bleached into the microtubule.
            %    Duration: duration of bleaching in seconds.
            
            for n = 1:length(Wavecount)
                if Wavecount(n) > 0
                    Wavelength = BES.MicrotubuleLength/Wavecount(n);
                    BES.bleachWavelength(Wavelength, Duration);
                end
            end
        end
        
        
        function BES = calculateExpectedWaveLengths(BES, Iterations)
            if nargin < 2
                Iterations = 1;
            end
            BES.ExpectedWavelenghts = max(BES.IntensityProfileX) ./ (BES.ExpectedWaveNumbers .* Iterations);
        end
        
        
        function BES = calculateIntensityProfileX(BES)
            BES.IntensityProfileX = (1:length(BES.IntensityProfile))' .* BES.SimRes ./ 1000;
        end
        
        
        function BES = importLineScan(BES, LineScanData, PixelSize)
            %importLineScan import Linescan Data
            %
            %see also: importLineScanFile
            
            BES.IntensityProfile = LineScanData;
            BES.PixelSize = PixelSize;
            BES.SimRes = PixelSize;
            BES.MicrotubuleLength = length(BES.IntensityProfile) * PixelSize;
            BES.calculateIntensityProfileX;
            BES.ImportedData = true;
        end
        
        
        function BES = simulateManyLineScans(BES, NFrames, Wavecount, Intensity)
            if nargin < 4
                Intensity = 0.2;
            end
            BES.IntensityProfiles = zeros(NFrames, length(BES.IntensityProfile));
            for n = 1:NFrames
                BES.resetMicrotubule;
                BES.bleachWaveCounts(Wavecount, Intensity);
                BES.IntensityProfiles(n,:) = BES.IntensityProfile;
            end
        end
        
        
        function [Result, MTLength] = fftAnalysis(BES)
            %fftAnalysis perform a fast furier transform analysis of the
            %intensity profile.
            
            if ~BES.ImportedData
                Profile = BES.removeTails(BES.IntensityProfile, false);
            else
                Profile = BES.IntensityProfile;
            end
            Profile = resample(double(Profile), round(BES.SimRes), round(BES.PixelSize));
            MTLength = length(Profile) - 1;
            Fourier = fft(Profile, MTLength);
            Result = BES.processFFT(Fourier, MTLength);
        end
        
        
        function [P1, WL] = fftWaveLength(BES, P1, MTLength)
            WL = BES.calculateWaveLengths(size(P1, 1), MTLength);
            Discard = BES.discardWaveLengths(WL);
            P1(Discard) = [];
            WL(Discard) = [];
        end
        
        
        function WL = calculateWaveLengths(BES, NPix, MTLength)
            L = ones(1, NPix);
            WL = ((L .* MTLength) ./ (1:round(MTLength / 2))) .* BES.PixelSize;
        end
        
        
        function Discard = discardWaveLengths(BES, WL)
            MaxWaveLength = BES.MicrotubuleLength / BES.MinWaveNumber;
            Discard = WL > MaxWaveLength | WL < BES.MinWaveLengthFactor * BES.ImagingResolution;
        end
        
        
        function Result = multiBleachWavecounts(BES, Iterations, varargin)
            %multiBleachWavecounts simulate multiple rounds of bleaching. Return
            %the fft analysis results for different wavecounts in the columns of
            %a matrix.
            
            [~, MTLength] = BES.fftAnalysis;
            Result = zeros(Iterations, round(MTLength / 2));
            for n = 1:Iterations
                BES.resetMicrotubule;
                BES.calculateIntensityProfile;
                BES.bleachWaveCounts(varargin{:});
                Result(n,:) = BES.fftAnalysis;
            end
        end
        
        
        function plotFFTAnalysisWaveNumbers(BES, varargin)
            %plotFFTAnalysis plot the fast furier transform analysis of the
            %intensity profile by wave number.
            
            Result = BES.fftAnalysis;
            MaxWaveNumber = round(BES.MicrotubuleLength/BES.BleachingResolution);
            if MaxWaveNumber > length(Result)
                MaxWaveNumber = length(Result);
            end
            Ax = BES.getAxes(varargin{:});
            plot(Ax, (6:MaxWaveNumber), Result(6:MaxWaveNumber), 'LineStyle', 'none', 'Marker', '*');
            xlabel(Ax, sprintf('wave number (1/%d %sm)', BES.MicrotubuleLength/1000, 181));
            ylabel(Ax, 'amplitude (a.u.)');
        end
        
        
        function plotFFTAnalysisWaveLength(BES, varargin)
            %plotFFTAnalysis plot the fast furier transform analysis of the
            %intensity profile by wavelength.
            
            [P1, MTLength] = BES.fftAnalysis;
            [P1, F] = BES.fftWaveLength(P1, MTLength);
            Ax = BES.getAxes(varargin{:});
            plot(Ax, F, P1, 'LineStyle', 'none', 'Marker', '*');
            BES.formatAxes(Ax, 'XLabel', 'spatial period (nm)', 'YLabel', 'amplitude (a.u.)');
        end
        
        
        function plotMultiFFTAnalysisWaveLength(BES, Iterations, Wavecount, Duration, varargin)
            %plotFFTAnalysis plot the fast furier transform analysis of the
            %intensity profile by wavelength.
            
            Results = BES.multiBleachWavecounts(Iterations, Wavecount, Duration);
            [P1, MTLength] = BES.fftAnalysis;
            [P1, F] = BES.fftWaveLength(P1, MTLength);
            Res = zeros(Iterations, length(P1));
            for n=1:Iterations
                Res(n,:) = BES.fftWaveLength(Results(n,:), MTLength);
            end
            S = std(Res);
            M = mean(Res);
            Ax = BES.getAxes(varargin{:});
            errorbar(Ax, F, M, S, 'LineStyle', 'none', 'Marker', '*');
            BES.formatAxes(Ax, 'XLabel', 'spatial period (nm)', 'YLabel', 'amplitude (a.u.)');
            TargetWL = BES.MicrotubuleLength / Wavecount;
            TargetSignal = M(abs(F-TargetWL)==min(abs(F-TargetWL)));
            hold(Ax, 'on');
            plot(Ax, TargetWL, TargetSignal + 0.2,'*', 'Color', [0.9290, 0.6940, 0.1250])
            hold(Ax, 'off');
        end
        
        
        function plotIntensityProfile(BES, varargin)
            %plotLineScan Plot a linescan of the raw intensity profile.
            %
            % See also: BleachEncodingSim

            Ax = BES.getAxes(varargin{:});
            plot(Ax, BES.IntensityProfileX, BES.IntensityProfile);
            BES.formatAxes(Ax)
            xlim(Ax, [0, max(BES.IntensityProfileX)]);
            ylim(Ax, [0 max(BES.IntensityProfile) * 1.1]);
        end
        
        
        function plotLineScan(BES, varargin)
            %plotLineScan Plot a linescan of the intensity profile of the bleached region.
            %
            % See also: BleachEncodingSim
            
            Profile = BES.getLineScan;
            DistanceMicrometer = BES.PixelSize/1000;
            X = 0:DistanceMicrometer:(length(Profile)-1)*DistanceMicrometer;
            Ax = BES.getAxes(varargin{:});
            plot(Ax, X, Profile);
            BES.formatAxes(Ax)
            xlim(Ax, [0, max(X)]);
            ylim(Ax, [0 max(Profile) * 1.1]);
        end
        
        
        function plotBleachingIntensity(BES, Pattern, varargin)
            %plotBleachingIntensity plot the bleaching profile including the
            %gaussian intensity profiles of the individual bleaching positions
            %
            % Pattern: Vector containing zero where we don't want to bleach and 1 where we want to bleach
            
            Maxima = find(Pattern);
            DistanceMicrometer = BES.SimRes/1000;
            OffsetMicrometer = BES.MTOffset * BES.SimRes / 1000;
            MTLengthMicrometer = BES.MicrotubuleLength / 1000;
            Ax = BES.getAxes(varargin{:});
            cla(Ax);
            hold(Ax, 'on');
            HalfLength = (length(BES.PSFBleach) - 1) / 2;
            for n = 1:length(Maxima)
                plot(Ax, (Maxima(n) + HalfLength + (-HalfLength:HalfLength)) .* DistanceMicrometer - OffsetMicrometer, BES.PSFBleach,'r');
            end
            BleachIntensityProfile = BES.createBleachIntensityProfile(Pattern);
            X = (0:DistanceMicrometer:(length(BleachIntensityProfile)-1)*DistanceMicrometer) - OffsetMicrometer;
            plot(Ax, X, BleachIntensityProfile,'b');
            BES.formatAxes(Ax, 'YLabel', 'bleaching intensity (a.u.)');
            xlim(Ax, [0, MTLengthMicrometer]);
            ylim(Ax, [0 1.1]);
        end
        
        
        function plotBleachingIntensityWC(BES, Wavecount, varargin)
            %plotBleachingIntensityWC plot the bleaching intensity for Wavecount.
            %
            % See also: plotBleachingPrediction,getPatternWaveCount

            Pattern = BES.getPatternWaveCount(Wavecount);
            BES.plotBleachingIntensity(Pattern, varargin{:});
        end
        
        
        function plotBleachingPrediction(BES, Wavecount, Duration, varargin)
            p = inputParser;
            p.addParameter('LinescanYlim', [0,11], @isnumeric);
            p.addParameter('FourierYlim', [0,2.8], @isnumeric);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            Tmp = [fieldnames(p.Unmatched),struct2cell(p.Unmatched)];
            PassThrough = reshape(Tmp',[],1)';
            
            BES.createFigure(PassThrough{:})
            Rows = 1;
            Cols = 18;
            BleachAx = subplot(Rows, Cols, 1:6);
            axes(BleachAx)
            BES.plotBleachingIntensityWC(Wavecount, 'Axes', BleachAx);
            BES.resetMicrotubule;
            BES.bleachWaveCounts(Wavecount, Duration)
            LineAx = subplot(Rows, Cols, 8:13);
            BES.plotLineScan('Axes', LineAx);
            ylim(LineAx, p.Results.LinescanYlim);
            FFTAx = subplot(Rows, Cols, 15:18);
            BES.plotMultiFFTAnalysisWaveLength(10, Wavecount, Duration, 'Axes', FFTAx);
            ylim(FFTAx, p.Results.FourierYlim);
        end
                    
        
        function [BES, Ax] = createFigure(BES, varargin)
            %createFigure create a new figure if necessary.
            %
            % See also: getAxes
            
            p = inputParser;
            p.addParameter('Width', BES.FigWidth, @isnumeric);
            p.addParameter('Aspect', BES.Aspect, @isnumeric);
            p.addParameter('XLabel', sprintf('distance along microtubule (%sm)', 181), @ischar);
            p.addParameter('YLabel', 'fluorescence intensity (a.u.)', @ischar);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            Tmp = [fieldnames(p.Unmatched),struct2cell(p.Unmatched)];
            FigArgs = reshape(Tmp',[],1)';
            Width=p.Results.Width;
            Height=Width/p.Results.Aspect;
            if p.Results.Aspect>1
                Orientation='landscape';
            else
                Orientation='portrait';
            end
            Units=get(0,'Units');
            set(0,'Units','centimeters');
            ScreenPos=get(0,'screensize');
            set(0,'Units',Units);
            Screenfill=[Width/ScreenPos(3) Height/ScreenPos(4)];
            if any(Screenfill>1)
                if find(Screenfill==max(Screenfill))==1
                    SWidth=Width/Screenfill(1)*0.95;
                    SHeight=SWidth/p.Results.Aspect;
                else
                    SHeight=Width/Screenfill(2)*0.95;
                    SWidth=SHeight*p.Results.Aspect;
                end
            else
                SWidth=Width;
                SHeight=Height;
            end
            if ishandle (BES.Fig)
                set(BES.Fig, 'PaperOrientation',Orientation);
            else
                BES.Fig = figure('Units','centimeters','PaperType','A4',...
                    'PaperOrientation',Orientation,'PaperUnits','centimeters',...
                    'InvertHardcopy','off','Color','w');
            end
            set(BES.Fig,'Position',[1 1 SWidth SHeight],...
                'PaperPosition',[0 0 Width Height],...
                'PaperSize',[Width Height],...
                'Renderer', 'painters',...
                FigArgs{:});
            if ishandle(BES.Fig.CurrentAxes)
                Ax = BES.Fig.CurrentAxes;
            else
                Ax = axes(BES.Fig);
            end
            Ax.FontName = 'Arial';
        end
        
        
        function Ax = getAxes(BES, varargin)
            %getAxes get axes for plotting.
            %
            %returns Axes object from BEC.Fig (newly creted if necessary).
            %If you pass an Axes object to getAxes, it will use that one
            %instead.
            %
            % See also: createFigure
            
            p = inputParser;
            p.addParameter('Axes', [], @(X) isempty(X) || ishandle(X));
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            Tmp = [fieldnames(p.Unmatched),struct2cell(p.Unmatched)];
            PassThrough = reshape(Tmp',[],1)';
            if isempty(p.Results.Axes)
                [~, Ax] = BES.createFigure(PassThrough{:});
            else
                Ax = p.Results.Axes;
            end
            axes(Ax);
        end
        
        
        function saveFig(BES, Filename)
            saveas(BES.Fig, Filename);
        end
    end
    methods (Static)
        
        
        function Ax = formatAxes(Ax, varargin)
            p = inputParser;
            p.addParameter('XLabel', sprintf('distance along microtubule (%sm)', 181), @ischar);
            p.addParameter('YLabel', 'intensity (a.u.)', @ischar);
            p.parse(varargin{:});
            set(Ax, 'FontName', 'Arial', 'FontSize', 8, 'LabelFontSizeMultiplier', 1.3);
            set(Ax, 'FontName', 'Arial', 'FontSize', 8, 'TitleFontSizeMultiplier', 1.3);
            xlabel(Ax, p.Results.XLabel, 'FontName', 'Arial', 'FontSize', 11);
            ylabel(Ax, p.Results.YLabel, 'FontName', 'Arial', 'FontSize', 11);
            box(Ax, 'on')
        end
        
        
        function Result = processFFT(Fourier, NSamples)
            P2 = abs(Fourier / NSamples);
            Result = P2(1:round(NSamples / 2) + 1);
            Result(2:end - 1) = 2 * Result(2:end - 1);
            Result(1:1) = [];
        end
        
        
    end
end