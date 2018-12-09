function modelOut =preProcess(model,fba)
    % modelOut = model without non-relevant processes and reactions
    % that doesn't contain high carbon metabolite
    % fba = flux values of wild model
    %listing for reaction match that involve below organelle and pathways
    rxnWith={'cell envelope','ER memberane','Golgi','Golgi membrane','lipid particle','phospholipase C','mitochondrial membrane','isa complex sphingolipid','IPC','MIPC','M(IP)2C','vacuolar','nucleus'};
    rxnIndex=[];
    for comp=1:length(rxnWith)
        %getRxnIndex function returns the rxn index from rxn name
        rxnId=getRxnIndex(model,rxnWith{comp});
        rxnIndex=[rxnIndex rxnId];
    end
    %finds metabolite index for higher carbon metabolite
    carbon7orMore=hiCarbonMetIndex(model);
    %finds reactions and reaction index that involve higher carbon metabolites
    rxnList=findRxnsFMet(model,model.mets(carbon7orMore));
    carbonIndex=findRxnIDs(model,rxnList);
    rxnIndex=[rxnIndex carbonIndex'];
    rxnIndex=unique(rxnIndex);
    %storing all the high carbon involving reactions that doesn't carry
    %signicicant flux
    nonReq=rxnIndex(fba(rxnIndex)==0);
    modelOut=removeRxns(model,model.rxns(nonReq));
    [modelOut,removedMets,removedRxns] = removeDeadEnds(modelOut);
end