function thestr=cell2str_2(thecell)
if ~iscell(thecell)
    thecell={thecell};
end
    
thestr='';
if ~isempty(thecell)
for p=1:length(thecell)
    if p==1
        thestr=[thestr,thecell{p}];
    else
        thestr=[thestr,',',thecell{p}];
    end
end
end
