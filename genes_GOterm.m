function [geneList]=genes_GOterm(x)
    %generating genelists per GO term from the string
    y=find(x=='/');
    geneList={};
    for item=1:length(y)
        gene='';
        for count=y(item)+1:length(x)
            if x(count)==':' | x(count)==''
                break
            else
                gene=[gene x(count)];
            end
        end
        geneList{item}=gene;
    end
end