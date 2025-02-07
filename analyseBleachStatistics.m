function Significant = analyseBleachStatistics(Data)
ControlCI = bootci(1000, @mean, [Data(:,[10 19]); Data(:,[9 20])]);
DataCI = bootci(1000, @mean, Data(:,11:18));
ControlLimit = interp1([9.5 19.5],ControlCI(2,:),[9.5 11:18 19.5]);
Limits = ControlLimit(2:9);
Significant = DataCI(1,:)>Limits;
end