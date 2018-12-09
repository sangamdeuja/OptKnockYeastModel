function carbonSum=carbonCount(mystring)
    %finds carbon number from molecular formula
    x=find(mystring=='C');
    noC=length(x);
    carbonSum=0;
    for count=1:noC
        saveNum='';
        for anotherCount=x(count)+1:length(mystring)
            check=str2num(mystring(anotherCount));
            if check | mystring(anotherCount)=='0'
                saveNum=[saveNum mystring(anotherCount)];
            else
                break
            end
        end
        carbonSum=carbonSum+str2num(saveNum);
    end
end