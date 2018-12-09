function[reducedModel,optKnockSolutions]=optKnockImplementation(model,objectiveRxnId,substrateRxnId,substrateUptake,modelReduction,numDel,koListIdentifier)
    % Inputs: model - yeast model
    % objectiveRxnId - reaction id of target production
    % substrateRxnId - substrate rxn id
    % substrateUptake- substrate flux
    % modelReduction - 1 for preprocess (produce reduced model I)
    % 2 for preprocessModel (produce reduced model II)
    % koListIdentifier- 0 for all reactions
    % 1 CC reactions except oxidative phosphorylation
    % [] for all CC reactions (default)

    %%addition of ATP maintenance reaction in our model
    model = addReaction(model,'ATP maintenance',{'s_0434','s_0803','s_0394','s_0794','s_1322'},[-1 -1 1 1 1 ],false);
    %%ATP requirement setting based on yeast research paper
    model = changeRxnBounds(model,'ATP maintenance',1,'b');
    %%Glucose uptake to 10 unit flux
    model = changeRxnBounds(model,substrateRxnId,substrateUptake,'b');
    %%0.01 is the minimal growth condition
    model = changeRxnBounds(model,'r_2111',0.05,'l');
    %CC reactions
    koCandidates={'r_0163','r_2115','r_2116','r_0165','r_0174','r_0175','r_0173','r_0112','r_1838','r_0113','r_0111','r_1054','r_0467','r_0450','r_0486','r_0533','r_0534','r_0886','r_0892','r_0893','r_0962','r_0366','r_0884','r_0303','r_0301','r_1000','r_0455','r_0452','r_0659','r_0661','r_0662','r_0714','r_0715','r_0716','r_0717','r_2305','r_0734','r_1048','r_1049','r_1050','r_0466','r_2126','r_0091','r_0889','r_0887','r_0984','r_0982','r_0990','r_0958','r_0959','r_0718','r_0719','r_0961','r_0280','r_0302','r_0300','r_1022','r_1021','r_0454','r_0451','r_0658','r_0713','r_0831','r_0832','r_2131','r_0449','r_0226','r_0438','r_0439','r_0773','r_4039'};
    if (nargin < 7)
        display('Taking central carbon reactions for knockouts')
    end
    if exist('koListIdentifier','var')
        if koListIdentifier==0
            koCandidates=model.rxns;
            rxns={objectiveRxnId,'r_2111',substrateRxnId,'ATP maintenance'};
            koCandidates=setdiff(koCandidates,rxns);
        end
        if koListIdentifier==1
            %list of oxidative phosphoryaltion reactions
            of={'r_0226','r_0438','r_0439','r_0773','r_4039'};
            koCandidates=setdiff(koCandidates,of);
        end

    end
    growthRxnId='r_2111';
    %%constraints for OptKnock
    constrOpt.rxnList={'r_2111',substrateRxnId,'ATP maintenance'};
    constrOpt.values=[0.05,substrateUptake,1];
    constrOpt.sense=['G','E','E'];
    wildFba=optimizeCbModel(model);
    fba=wildFba.x;
    if modelReduction==1
        modelOut=preProcess(model,fba)
    end
    if modelReduction==2
        [modelOut,koCandidates]=preprocessModel(model,objectiveRxnId,koCandidates);
    end
    reducedModel=modelOut;
    optKnockSolutions=[];
    solution=[];
    title={'wildgrowth' 'wildobj' 'wildoxygen' 'optknockgrowth' 'optknockobj' 'optknockoxygen' 'Rmutgrowth' 'Rmutobj' 'Rmutoxygen' 'Fmgrowth' 'Fmobj' 'Fmoxygen'};
    opt.targetRxn=objectiveRxnId;
    opt.vMax=1000;
    opt.numDelSense='L';
    for deletion=1:numDel
        opt.numDel=deletion;
        optKnockSol=OptKnock(modelOut,koCandidates,opt,constrOpt)
        optKnockSolutions=[optKnockSolutions optKnockSol];
        mutant=modelOut;
        xModel=model;
        for item=1:length(optKnockSol.rxnList)
            rxnInfo(modelOut,strmatch(optKnockSol.rxnList{item},modelOut.rxns))
            mutant=changeRxnBounds(mutant,optKnockSol.rxnList{item},0,'b');
            xModel=changeRxnBounds(xModel,optKnockSol.rxnList{item},0,'b');
        end
        display('wild model')
        display(['wild growth= ',num2str(wildFba.f)])
        display(['wild objectiveRxn= ',num2str(wildFba.x(findRxnIDs(model,objectiveRxnId)))])
        display(['wild oxygen uptake= ',num2str(wildFba.x(findRxnIDs(model,'r_1992')))])
        display('---------------------------------');
        if optKnockSol.stat~=1
            optKnockSol.fluxes(1:length(modelOut.rxns))=0;
        end
        display(['Optknock growth= ',num2str(optKnockSol.fluxes(findRxnIDs(modelOut,'r_2111')))])
        display(['OptKnock objectiveRxn= ',num2str(optKnockSol.obj)])
        display(['Optknock oxygen uptake= ',num2str(optKnockSol.fluxes(findRxnIDs(modelOut,'r_1992')))])
        display('---------------------------------');
        display('validation with reduced model')

        postFba=optimizeCbModel(mutant);
        if postFba.stat~=1
            postFba.x(1:length(modelOut.rxns))=0;
        end
        display(['mutant growth= ',num2str(postFba.f)])
        display(['mutant objectiveRxn= ',num2str(postFba.x(findRxnIDs(modelOut,objectiveRxnId)))])
        display(['mutant oxygen uptake= ',num2str(postFba.x(findRxnIDs(modelOut,'r_1992')))])
        display('---------------------------------');
        display('validation with wild model')
        xModelFba=optimizeCbModel(xModel);
        if xModelFba.stat~=1
            xModelFba.x(1:length(xModel.rxns))=0;
        end
        display(['mutant growth= ',num2str(xModelFba.f)])
        display(['mutant objectiveRxn= ',num2str(xModelFba.x(findRxnIDs(xModel,objectiveRxnId)))])
        display(['mutant oxygen uptake= ',num2str(xModelFba.x(findRxnIDs(xModel,'r_1992')))])
        display('*******************************************************************************');
        toWrite=[wildFba.f wildFba.x(findRxnIDs(model,objectiveRxnId)) wildFba.x(findRxnIDs(model,'r_1992')) optKnockSol.fluxes(findRxnIDs(modelOut,'r_2111')) optKnockSol.obj optKnockSol.fluxes(findRxnIDs(modelOut,'r_1992')) postFba.f postFba.x(findRxnIDs(modelOut,objectiveRxnId)) postFba.x(findRxnIDs(modelOut,'r_1992')) xModelFba.f xModelFba.x(findRxnIDs(xModel,objectiveRxnId)) xModelFba.x(findRxnIDs(xModel,'r_1992'))]
        solution=[solution;toWrite];
    end
    solution=[title;num2cell(solution)]
    xlswrite('optKnockSolution.xls',solution);
end