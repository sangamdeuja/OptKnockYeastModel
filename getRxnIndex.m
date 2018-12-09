function rxnIndex=getRxnIndex(model,rxnName)
    % returns reaction index by reaction name in the model
    reactionNames=[];
    ind=[];
    for item=1:length(model.c)
        reactionNames=[reactionNames,model.rxnNames(item)];
    end
    for item=1:length(model.c)
        rxnCellToChar=char(reactionNames(1,item));
        subStrPos=findstr(rxnName,rxnCellToChar);
        if(subStrPos>=1)
            ind=[ind item];
        end
    end
    rxnIndex=ind
end