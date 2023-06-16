



for si = 1:4
    switch si
        case 1
            stat = statSummaryR1;
        case 2
            stat = statSummaryR2;
        case 3
            stat = statSummaryR3;
        case 4
            stat = statSummaryR4;
    end
    varFields = fieldnames(stat);

    for vi = 1:length(varFields)

        pvalue = stat.(varFields{vi}).Condition.Anova.('Prob>F'){1};

        if pvalue < 0.05
            formatSpec = "Var %s from RecogS %s is significant";
            A1 = num2str(vi);
            A2 = num2str(si);
            str = sprintf(formatSpec,A1,A2);
            disp(str)
        % else
        %     formatSpec = "P value is %d ";
        %     A1 = pvalue;
        %     str = sprintf(formatSpec,A1)
        end
    end
end