function [] = eyeTrackPlots(cleanedDataLOC, ptID , plotNUM)


cd(cleanedDataLOC)
load(ptID, 'variantS')



switch plotNUM
    case 1 % learning categories

        variNUM = fieldnames(variantS);

        for vi = 1:length(variNUM)

            learnLeft = variantS.(variNUM{vi}).curVariant.Left_L_oT_pupilS_rawCL;
            learnRight = variantS.(variNUM{vi}).curVariant.Right_L_oT_pupilS_rawCL;





        end

    case 2











end

end