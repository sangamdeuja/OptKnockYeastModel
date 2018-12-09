function [carbon7orMore]=hiCarbonMetIndex(model)
    %%sums up the carbon atoms from molecular formula, returns metabolite index
    %%having more carbon number
    carbon7orMore=[];
    for item=1:length(model.mets)
        a=char(model.metFormulas(item));
        carbonSum=carbonCount(a);
        if carbonSum>=7
            carbon7orMore=[carbon7orMore item];
        end
    end
    %higherRxns = model.rxns(find(sum(model.S(notToKeep,:) ~= 0) > 0));
    %modelOut = removeRxns(model,higherRxns);
end