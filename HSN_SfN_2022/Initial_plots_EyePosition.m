%% Data location
pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-FAGRV5G'
        mainDIR = 'E:\Dropbox\SfN_2022\dataOut_AMC';
        dataDIR = [mainDIR, filesep , 'dataOut_AMC\eyeDATA'];
end


cd(dataDIR)

%% Load the data

load("eyeData_AMC_PY22NO12.mat")

%% Plot types
plottype = 1;

switch plottype
    case 1


SD = chart.ScatterDensity("XData", x, ...
    "YData", y );

SD.SizeData = 400;
grid( SD, "off" )

    case 2



    case 3
        load fisheriris.mat;
x = meas(:,1);
y = meas(:,2);

scatterhist(x,y,'Group',species,'Kernel','on')


    case 4



end