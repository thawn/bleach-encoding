function [MeanLineScan, PixelSize, LineScans, Signal, SuperResImage, OrigImage, Line] = superResLineScan(FiestaFile, FilamentNo, StackName, varargin)
p = inputParser;
p.addParameter('ScanSize', 0, @isnumeric); %Width of linescan
p.addParameter('ColorIdx', 1, @isnumeric);
p.addParameter('ResolutionFactor', 3, @isscalar);
p.addParameter('Average', true, @islogical);
p.addParameter('OrigImagePadding', 5, @isscalar);
p.addParameter('Debug', false, @islogical)
p.parse(varargin{:});
FiestaData = load(FiestaFile, 'Filament');
Background = double(median(FiestaData.Filament.Data{1}(:,7)));
if p.Results.Average
    Frames = FiestaData.Filament(FilamentNo).Results(:,1)';
else
    Frames = FiestaData.Filament(FilamentNo).Results(1,1)';
end
NFrames = length(Frames);
Stack = fReadND2(StackName);
I = cell(NFrames, 1);
Lengths = zeros(NFrames, 1);
Count = 1;
if p.Results.Debug
    figure();
    Ax1 = subplot(2, 1, 1);
    hold(Ax1, 'on');
    Ax2 = subplot(2, 1, 2);
end
if nargout > 3
    WideLines = cell(NFrames, 1);
end
for Frame = Frames
    nX = (double(FiestaData.Filament(FilamentNo).Data{Count}(:,1)) / FiestaData.Filament(FilamentNo).PixelSize);
    nY = (double(FiestaData.Filament(FilamentNo).Data{Count}(:,2)) / FiestaData.Filament(FilamentNo).PixelSize);
    d = [0; cumsum(sqrt((nX(2:end) - nX(1:end - 1)) .^ 2 + (nY(2:end) - nY(1:end - 1)) .^ 2))];
    dt = max(d) / round(max(d) * p.Results.ResolutionFactor);
    id = (0:round(max(d) * p.Results.ResolutionFactor))' * dt;
    scan_length = length(id);
    idx = nearestpoint(id, d);
    X = zeros(scan_length, 1);
    Y = zeros(scan_length, 1);
    dis = id - d(idx);
    dis(1) = 0;
    dis(end) = 0;
    X(dis == 0) = nX(idx(dis == 0));
    Y(dis == 0) = nY(idx(dis == 0));
    X(dis > 0) = nX(idx(dis > 0)) + (nX(idx(dis > 0) + 1) - nX(idx(dis > 0))) ./ (d(idx(dis > 0) + 1) - d(idx(dis > 0))) .* dis(dis > 0);
    Y(dis > 0) = nY(idx(dis > 0)) + (nY(idx(dis > 0) + 1) - nY(idx(dis > 0))) ./ (d(idx(dis > 0) + 1) - d(idx(dis > 0))) .* dis(dis > 0);
    X(dis < 0) = nX(idx(dis < 0)) + (nX(idx(dis < 0) - 1) - nX(idx(dis < 0))) ./ (d(idx(dis < 0) - 1) - d(idx(dis < 0))) .* dis(dis < 0);
    Y(dis < 0) = nY(idx(dis < 0)) + (nY(idx(dis < 0) - 1) - nY(idx(dis < 0))) ./ (d(idx(dis < 0) - 1) - d(idx(dis < 0))) .* dis(dis < 0);
    iX = zeros(2 * p.Results.ScanSize + 1, scan_length);
    iY = zeros(2 * p.Results.ScanSize + 1, scan_length);
    n = zeros(scan_length, 3);
    for i = 1:length(X)
        if i == 1
            v = [X(i + 1) - X(i), Y(i + 1) - Y(i), 0];
            n(i, :) = [v(2), -v(1), 0] / norm(v);
        elseif i == length(X)
            v = [X(i) - X(i - 1), Y(i) - Y(i - 1), 0];
            n(i, :) = [v(2), -v(1), 0] / norm(v);
        else
            v1 = [X(i + 1) - X(i), Y(i + 1) - Y(i), 0];
            v2 = -[X(i) - X(i - 1), Y(i) - Y(i - 1), 0];
            n(i, :) = v1 / norm(v1) + v2 / norm(v2);
            if norm(n(i, :)) == 0
                n(i, :) = [v1(2), -v1(1), 0]/norm(v1);
            else
                n(i,:) = n(i, :) / norm(n(i, :));
            end
            z=cross(v1, n(i, :));
            if z(3) > 0
                n(i, :) = -n(i, :);
            end
        end
        iX(:, i) = linspace(X(i) + p.Results.ScanSize * n(i, 1), X(i) - p.Results.ScanSize * n(i, 1), 2 * p.Results.ScanSize + 1)';
        iY(:, i) = linspace(Y(i) + p.Results.ScanSize * n(i, 2), Y(i) - p.Results.ScanSize * n(i, 2), 2 * p.Results.ScanSize + 1)';
    end
    d = [0; cumsum(sqrt(((X(2:end) - X(1:end - 1)) .^ 2) + ((Y(2:end) - Y(1:end - 1)) .^ 2)))];
    Z = interp2(double(Stack{p.Results.ColorIdx}(:,:,Frame)),iX,iY);
    I{Count} = mean(Z,1);
    Lengths(Count) = length(I{Count});
    if nargout > 3
        WideLines{Count} = Z;
    end
    if nargout > 3 && Frame == Frames(1)
        StartY = floor(min(nY)) - p.Results.OrigImagePadding;
        StartY(StartY < 1) = 1;
        EndY = ceil(max(nY)) + p.Results.OrigImagePadding;
        EndY(EndY > size(Stack{p.Results.ColorIdx}, 1)) = size(Stack{p.Results.ColorIdx}, 1);
        StartX = floor(min(nX)) - p.Results.OrigImagePadding;
        StartX(StartX < 1) = 1;
        EndX = ceil(max(nX)) + p.Results.OrigImagePadding;
        EndX(EndX > size(Stack{p.Results.ColorIdx}, 2)) = size(Stack{p.Results.ColorIdx}, 2);
        OrigImage = Stack{p.Results.ColorIdx}(StartY:EndY, StartX:EndX, Frame);
        Line = [nX - StartX + 1, nY - StartY + 1];
    end
    if p.Results.Debug
        plot(Ax1, d * FiestaData.Filament(FilamentNo).PixelSize / 1000, I{Count}, 'DisplayName', sprintf('Frame %d', Frame));
    end
    Count = Count + 1;
end
MaxLen = max(Lengths);
LineScans = NaN(NFrames, MaxLen);
if nargout > 3
    Images = zeros(2 * p.Results.ScanSize + 1, MaxLen, NFrames);
end
for k = 1:NFrames
    if k == 1
        Offset = floor((MaxLen - length(I{k})) / 2) + 1;
    else
        [XC, Lag] = xcorr(I{k-1}, I{k});
        [~, MaxCorr] = max(XC);
        Offset = Lag(MaxCorr) + Offset;
        if Offset < 1
            I{k}(1:1-Offset) = [];
            if nargout > 3
                WideLines{k}(:, 1:1-Offset) = [];
            end
            Offset = 1;
        end
    end
    LineScans(k,Offset:Offset + length(I{k}) - 1) = I{k};
    if nargout > 3
        Images(:, Offset:Offset + length(I{k}) - 1, k) = WideLines{k};
    end
end
LineScans = LineScans - Background;
Signal = findSignal(LineScans(:));
MeanLineScan = nanmean(LineScans, 1);
LineScans(:, isnan(MeanLineScan)) = [];
MeanLineScan(isnan(MeanLineScan)) = [];
if nargout > 3
    SuperResImage = mean(Images, 3);
end
PixelSize = FiestaData.Filament(FilamentNo).PixelSize / p.Results.ResolutionFactor;
if p.Results.Debug
    plot(Ax2, (0:length(MeanLineScan) - 1) .* PixelSize / 1000, MeanLineScan);
end
end
