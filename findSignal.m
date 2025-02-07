function Signal = findSignal(Intensities)
NBins=round(length(Intensities)/250)*5; %we want on average 50 intensities per bar
if NBins<7 %make sure we have at least 7 bins so that confint works can calculate the confidence intervals
    NBins=7;
end
BinSize=ceil((max(Intensities)-min(Intensities))/(NBins*10))*10;
MinS=floor(min(Intensities)/BinSize)*BinSize;
MaxS=ceil(max(Intensities)/BinSize)*BinSize;
HistX=MinS:BinSize:MaxS;
[Count,HistX] = hist(Intensities, HistX);
Cfun=fit(double(HistX'),double(Count'),'gauss2');
Signal = [Cfun.b1 Cfun.b2];
if any(isnan(Signal))
    Signal = nanmedian(Intensities);
end
end