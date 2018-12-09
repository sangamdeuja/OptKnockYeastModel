function [modelOut2,koCandidates]=preprocessModel(model,objRxn,koCandidates)
    %modelOut2=model without non-relevent reactions, without non-flux
    % carrying reactions, without higher carbon involving reactions

    % [m,n]=size(model.S);
    modelOut=model;
    nonReq=find((modelOut.lb==0)&(modelOut.ub==0));%length(nonReq)=31
    display(['Number of reaction without flux in original model is ',num2str(length(nonReq))])
    remRxns=modelOut.rxns(nonReq);
    koCandidates=setdiff(koCandidates,remRxns);
    %display(['Number of KO candidates is ',num2str(length(koCandidates))])
    modelOut=removeRxns(modelOut,remRxns);
    %%Removing the metabolic process(glycerophospholipid metabolic process, glycerophospholipid biosynthetic process, membrane lipid metabolic process, membrane lipid biosynthetic process, tRNA aminoacylation, tRNA aminoacylation for protein translation, tRNA aminoacylation for mitochondrial protein translation, tRNA metabolic process, inorganic ion transmembrane transport, cation transmembrane transport, inorganic cation transmembrane transport, anion transmembrane transport, transmembrane transport and ion transmembrane transport.) from the model
    [num,txt,raw_genes]=xlsread('genes_GO.xlsx');
    optim=optimizeCbModel(model);
    model1=changeObjective(model,objRxn);
    optim1=optimizeCbModel(model1);
    all_rxns={};
    length_rxns=[];
    for count=1:length(raw_genes(:,1))
        genesPerGO= genes_GOterm(raw_genes{count,2});
        [results ListResults] = findRxnsFromGenes(model,genesPerGO,0,1);
        % length_rxns=[length_rxns length(ListResults(:,1)')];
        % all_genes=[all_genes genesPerGO];

        all_rxns=[all_rxns ListResults(:,1)'];
    end
    [hiCarbonRxns,zeroCarbonRxns,nCarbon] = findCarbonRxns(model,6);
    all_rxns=[all_rxns hiCarbonRxns']
    all_rxns=unique(all_rxns);
    z=optim.x;
    non_zeroz=find(z~=0);
    y=optim1.x;
    non_zeroy=find(y~=0);
    biomassProduct=[model.rxns(non_zeroz)' model.rxns(non_zeroy)'];
    toRemove=setdiff(all_rxns,biomassProduct);
    modelOut1=removeRxns(modelOut,toRemove);
    koCandidates=setdiff(koCandidates,toRemove);
    optimOut1=optimizeCbModel(modelOut1);
    [modelOut2,removedMets,removedRxns] = removeDeadEnds(modelOut1);
    koCandidates=setdiff(koCandidates,removedRxns);
    optimOut2=optimizeCbModel(modelOut2);
end