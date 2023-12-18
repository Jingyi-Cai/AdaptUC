function [solution,LP]=multi_condition_design0(models,outconstr,delcand,optionstcf)

nc=length(models);
model1=models{1};
thedelCand.rxnass_del1=findRxnIDs(model1,delcand);



delCandrxns=model1.rxns(thedelCand.rxnass_del1);
[~,delCand1]=ismember(delCandrxns,model1.rxns);

if exist('optionstcf')
    if ~isfield(optionstcf,'totalTime')
        optionstcf.totalTime=3600;
    end
    if ~isfield(optionstcf,'totalCasesNeeded')
        optionstcf.totalCasesNeeded=20;
    end
else
    if ~isfield(outconstr,'maxdel')
        outconstr.maxdel=6;
    end
    optionstcf.totalTime=3600;
    optionstcf.totalCasesNeeded=20;
end
if ~exist('outconstr')
    outconstr=[];
    defaultobjflag=1;
elseif ~isfield(outconstr,'objind')|~isfield(outconstr,'sense')
    defaultobjflag=1;
else
    defaultobjflag=0;
end
if ~defaultobjflag&~isfield(outconstr,'objcoef')
         outconstr.objcoef=1;
end

nDelCand=length(delCand1);


LP= Cplex('FBA');

lb{1}=models{1}.lb;
ub{1}=models{1}.ub;
n1=length(models{1}.rxns);m1=length(models{1}.mets);
m1=length(models{1}.mets);    
history_length=length(LP.Model.lb);
LP.addCols(zeros(n1,1), [], lb{1}, ub{1},[],char(strcat('v_',num2str(1),models{1}.rxns)));
history_length_v(1)=history_length;
current_length_v(1)=history_length_v(1)+n1;       
LP.addRows(zeros(m1,1), [sparse(m1,history_length),models{1}.S], zeros(m1,1),char(strcat('sv0_',num2str(1),'_',models{1}.mets))); 

lb{2}=models{2}.lb;
ub{2}=models{2}.ub;
n1=length(models{2}.rxns);m1=length(models{2}.mets);
m1=length(models{2}.mets);    
history_length=length(LP.Model.lb);
LP.addCols(zeros(n1,1), [], lb{2}, ub{2},[],char(strcat('v_',num2str(2),models{2}.rxns)));
history_length_v(2)=history_length;
current_length_v(2)=history_length_v(2)+n1;       
LP.addRows(zeros(m1,1), [sparse(m1,history_length),models{2}.S], zeros(m1,1),char(strcat('sv0_',num2str(2),'_',models{2}.mets)));     

lb{3}=models{3}.lb;
ub{3}=models{3}.ub;
n1=length(models{3}.rxns);m1=length(models{3}.mets);
m1=length(models{3}.mets);    
history_length=length(LP.Model.lb);
LP.addCols(zeros(n1,1), [], lb{3}, ub{3},[],char(strcat('v_',num2str(3),models{3}.rxns)));
history_length_v(3)=history_length;
current_length_v(3)=history_length_v(3)+n1;       
LP.addRows(zeros(m1,1), [sparse(m1,history_length),models{3}.S], zeros(m1,1),char(strcat('sv0_',num2str(3),'_',models{3}.mets))); 
    
lb{4}=models{4}.lb;
ub{4}=models{4}.ub;
n1=length(models{4}.rxns);m1=length(models{4}.mets);
m1=length(models{4}.mets);    
history_length=length(LP.Model.lb);
LP.addCols(zeros(n1,1), [], lb{4}, ub{4},[],char(strcat('v_',num2str(4),models{4}.rxns)));
history_length_v(4)=history_length;
current_length_v(4)=history_length_v(4)+n1;       

LP.addRows(zeros(m1,1), [sparse(m1,history_length),models{4}.S], zeros(m1,1),char(strcat('sv0_',num2str(4),'_',models{4}.mets))); 


if ~isempty(outconstr)
    LP.addRows(outconstr.lb,outconstr.S,outconstr.ub);
end

lb{1}=models{1}.lb;
ub{1}=models{1}.ub;
n1=length(models{1}.rxns);
m1=length(models{1}.mets);


history_length=length(LP.Model.lb);
LP.addCols(zeros(n1,1), [], zeros(n1,1), inf(n1,1),[],char(strcat('dualV',num2str(1),'_lb_',models{1}.rxns)));
history_length_dualvarl(1)=history_length;
current_length_dualvarl(1)=history_length_dualvarl(1)+n1;

history_length=current_length_dualvarl(1);    
LP.addCols(zeros(n1,1), [], zeros(n1,1), inf(n1,1),[],char(strcat('dualV',num2str(1),'_ub_',models{1}.rxns)));
history_length_dualvaru(1)=history_length;
current_length_dualvaru(1)=history_length_dualvaru(1)+n1;

history_length=current_length_dualvaru(1);    
LP.addCols(zeros(m1,1), [], -inf(m1,1), inf(m1,1),[],char(strcat('lambda_',num2str(1),'_',models{1}.mets)));
history_length_dualm(1)=history_length;
current_length_dualm(1)=history_length_dualm(1)+m1;


LP.addRows(models{1}.c(:), [sparse(n1,history_length_dualvarl(1)), -speye(n1), speye(n1),...
    sparse(n1,history_length_dualm(1)-current_length_dualvaru(1)),models{1}.S',...
    ], models{1}.c(:),char(strcat('dual_',num2str(1),models{1}.rxns)));


lb{2}=models{2}.lb;
ub{2}=models{2}.ub;
n1=length(models{2}.rxns);
m1=length(models{2}.mets);


history_length=length(LP.Model.lb);
LP.addCols(zeros(n1,1), [], zeros(n1,1), inf(n1,1),[],char(strcat('dualV',num2str(2),'_lb_',models{2}.rxns)));
history_length_dualvarl(2)=history_length;
current_length_dualvarl(2)=history_length_dualvarl(2)+n1;

history_length=current_length_dualvarl(2);    
LP.addCols(zeros(n1,1), [], zeros(n1,1), inf(n1,1),[],char(strcat('dualV',num2str(2),'_ub_',models{2}.rxns)));
history_length_dualvaru(2)=history_length;
current_length_dualvaru(2)=history_length_dualvaru(2)+n1;

history_length=current_length_dualvaru(2);    
LP.addCols(zeros(m1,1), [], -inf(m1,1), inf(m1,1),[],char(strcat('lambda_',num2str(2),'_',models{2}.mets)));
history_length_dualm(2)=history_length;
current_length_dualm(2)=history_length_dualm(2)+m1;

LP.addRows(models{2}.c(:), [sparse(n1,history_length_dualvarl(2)), -speye(n1), speye(n1),...
    sparse(n1,history_length_dualm(2)-current_length_dualvaru(2)),models{2}.S',...
    ], models{2}.c(:),char(strcat('dual_',num2str(2),models{2}.rxns)));

lb{3}=models{3}.lb;
ub{3}=models{3}.ub;
n1=length(models{3}.rxns);
m1=length(models{3}.mets);

history_length=length(LP.Model.lb);
LP.addCols(zeros(n1,1), [], zeros(n1,1), inf(n1,1),[],char(strcat('dualV',num2str(3),'_lb_',models{3}.rxns)));
history_length_dualvarl(3)=history_length;
current_length_dualvarl(3)=history_length_dualvarl(3)+n1;

history_length=current_length_dualvarl(3);    
LP.addCols(zeros(n1,1), [], zeros(n1,1), inf(n1,1),[],char(strcat('dualV',num2str(3),'_ub_',models{3}.rxns)));
history_length_dualvaru(3)=history_length;
current_length_dualvaru(3)=history_length_dualvaru(3)+n1;

history_length=current_length_dualvaru(3);    
LP.addCols(zeros(m1,1), [], -inf(m1,1), inf(m1,1),[],char(strcat('lambda_',num2str(3),'_',models{3}.mets)));
history_length_dualm(3)=history_length;
current_length_dualm(3)=history_length_dualm(3)+m1;
LP.addRows(models{3}.c(:), [sparse(n1,history_length_dualvarl(3)), -speye(n1), speye(n1),...
    sparse(n1,history_length_dualm(3)-current_length_dualvaru(3)),models{3}.S',...
    ], models{3}.c(:),char(strcat('dual_',num2str(3),models{3}.rxns)));


lb{4}=models{4}.lb;
ub{4}=models{4}.ub;
n1=length(models{4}.rxns);
m1=length(models{4}.mets);
history_length=length(LP.Model.lb);
LP.addCols(zeros(n1,1), [], zeros(n1,1), inf(n1,1),[],char(strcat('dualV',num2str(4),'_lb_',models{4}.rxns)));
history_length_dualvarl(4)=history_length;
current_length_dualvarl(4)=history_length_dualvarl(4)+n1;
history_length=current_length_dualvarl(4);    
LP.addCols(zeros(n1,1), [], zeros(n1,1), inf(n1,1),[],char(strcat('dualV',num2str(4),'_ub_',models{4}.rxns)));
history_length_dualvaru(4)=history_length;
current_length_dualvaru(4)=history_length_dualvaru(4)+n1;
history_length=current_length_dualvaru(4);    
LP.addCols(zeros(m1,1), [], -inf(m1,1), inf(m1,1),[],char(strcat('lambda_',num2str(4),'_',models{4}.mets)));
history_length_dualm(4)=history_length;
current_length_dualm(4)=history_length_dualm(4)+m1;
LP.addRows(models{4}.c(:), [sparse(n1,history_length_dualvarl(4)), -speye(n1), speye(n1),...
    sparse(n1,history_length_dualm(4)-current_length_dualvaru(4)),models{4}.S',...
    ], models{4}.c(:),char(strcat('dual_',num2str(4),models{4}.rxns)));
LP.addCols(zeros(nDelCand,1), [], zeros(nDelCand,1), ones(nDelCand,1),char('B'*ones(1,nDelCand)),char(strcat('del_',models{1}.rxns(delCand1))));
history_length=current_length_dualm(4);
history_length_y=history_length;
current_length_y=history_length_y+nDelCand;
LP.addRows(0,[sparse(1,history_length_y),ones(1,nDelCand)],optionstcf.maxdel);
lbbig{1} = lb{1} <= -999;
ubbig{1} = ub{1} >= 999;

n1=length(models{1}.rxns);
sense = char('E' * ones(1,nDelCand));
rhs = zeros(nDelCand,1);
ind_del = (history_length+1:history_length+nDelCand)';
con = sparse(delCand1,1:nDelCand,ones(nDelCand,1),n1,nDelCand); % c
con1 = [sparse(history_length_v(1),nDelCand);con; sparse(current_length_y-current_length_v(1),nDelCand)]; 
con1 = mat2cell(con1, current_length_y, ones(nDelCand,1));
LP.addIndicators(ind_del, zeros(nDelCand,1), con1, sense, rhs);


delBinf = find(lbbig{1}(delCand1));       
nDelB = numel(delBinf);
ind = history_length+delBinf;
con = sparse(delCand1(delBinf), 1:nDelB , ones(nDelB,1), n1, nDelB);
con1 = [sparse(history_length_dualvarl(1),nDelB); con; sparse(current_length_y-current_length_dualvarl(1),nDelB)];
con1 = mat2cell(con1, current_length_y, ones(nDelB,1));
sense = char('E' * ones(1,nDelB));
rhs = zeros(nDelB,1);
LP.addIndicators(ind, ones(nDelB,1), con1, sense, rhs);

delBinf = find(ubbig{1}(delCand1));
nDelB = numel(delBinf);
ind = history_length+delBinf;
con = sparse(delCand1(delBinf), 1:nDelB , ones(nDelB,1), n1, nDelB);
con1 = [sparse(history_length_dualvaru(1),nDelB); con; sparse(current_length_y-current_length_dualvaru(1),nDelB)];
con1 = mat2cell(con1, current_length_y, ones(nDelB,1));
sense = char('E' * ones(1,nDelB));
rhs = zeros(nDelB,1);
LP.addIndicators(ind, ones(nDelB,1), con1, sense, rhs);


delset0 = setdiff(find(lbbig{1}(1:n1)),delCand1);
LP.Model.ub(history_length_dualvarl(1)+delset0) = 0;
delset0 = setdiff(find(ubbig{1}(1:n1)),delCand1);
LP.Model.ub(history_length_dualvaru(1)+delset0) = 0;


the_ub{1}=-(((~ubbig{1}(:)))' .* ub{1}(:)');
the_lb{1}=((~lbbig{1}(:)))' .* lb{1}(:)';
LP.addRows(0, [sparse(1,history_length_v(1)), models{1}.c',sparse(1,history_length_dualvarl(1)-current_length_v(1)),...
    the_lb{1},the_ub{1},sparse(1,current_length_y-current_length_dualvaru(1))], 0);    

lbbig{2} = lb{2} <= -999;
ubbig{2} = ub{2} >= 999;

n1=length(models{2}.rxns);
sense = char('E' * ones(1,nDelCand));
rhs = zeros(nDelCand,1);
ind_del = (history_length+1:history_length+nDelCand)';
con = sparse(delCand1,1:nDelCand,ones(nDelCand,1),n1,nDelCand); % c
con1 = [sparse(history_length_v(2),nDelCand);con; sparse(current_length_y-current_length_v(2),nDelCand)]; 
con1 = mat2cell(con1, current_length_y, ones(nDelCand,1));
LP.addIndicators(ind_del, zeros(nDelCand,1), con1, sense, rhs);


delBinf = find(lbbig{2}(delCand1));       
nDelB = numel(delBinf);
ind = history_length+delBinf;
con = sparse(delCand1(delBinf), 1:nDelB , ones(nDelB,1), n1, nDelB);
con1 = [sparse(history_length_dualvarl(2),nDelB); con; sparse(current_length_y-current_length_dualvarl(2),nDelB)];
con1 = mat2cell(con1, current_length_y, ones(nDelB,1));
sense = char('E' * ones(1,nDelB));
rhs = zeros(nDelB,1);
LP.addIndicators(ind, ones(nDelB,1), con1, sense, rhs);

delBinf = find(ubbig{2}(delCand1));
nDelB = numel(delBinf);
ind = history_length+delBinf;
con = sparse(delCand1(delBinf), 1:nDelB , ones(nDelB,1), n1, nDelB);
con1 = [sparse(history_length_dualvaru(2),nDelB); con; sparse(current_length_y-current_length_dualvaru(2),nDelB)];
con1 = mat2cell(con1, current_length_y, ones(nDelB,1));
sense = char('E' * ones(1,nDelB));
rhs = zeros(nDelB,1);
LP.addIndicators(ind, ones(nDelB,1), con1, sense, rhs);


delset0 = setdiff(find(lbbig{2}(1:n1)),delCand1);
LP.Model.ub(history_length_dualvarl(2)+delset0) = 0;
delset0 = setdiff(find(ubbig{2}(1:n1)),delCand1);
LP.Model.ub(history_length_dualvaru(2)+delset0) = 0;


the_ub{2}=-(((~ubbig{2}(:)))' .* ub{2}(:)');
the_lb{2}=((~lbbig{2}(:)))' .* lb{2}(:)';
LP.addRows(0, [sparse(1,history_length_v(2)), models{2}.c',sparse(1,history_length_dualvarl(2)-current_length_v(2)),...
    the_lb{2},the_ub{2},sparse(1,current_length_y-current_length_dualvaru(2))], 0);   

    lbbig{3} = lb{3} <= -999;
ubbig{3} = ub{3} >= 999;

n1=length(models{3}.rxns);
sense = char('E' * ones(1,nDelCand));
rhs = zeros(nDelCand,1);
ind_del = (history_length+1:history_length+nDelCand)';
con = sparse(delCand1,1:nDelCand,ones(nDelCand,1),n1,nDelCand); % c
con1 = [sparse(history_length_v(3),nDelCand);con; sparse(current_length_y-current_length_v(3),nDelCand)]; 
con1 = mat2cell(con1, current_length_y, ones(nDelCand,1));
LP.addIndicators(ind_del, zeros(nDelCand,1), con1, sense, rhs);


delBinf = find(lbbig{3}(delCand1));       
nDelB = numel(delBinf);
ind = history_length+delBinf;
con = sparse(delCand1(delBinf), 1:nDelB , ones(nDelB,1), n1, nDelB);
con1 = [sparse(history_length_dualvarl(3),nDelB); con; sparse(current_length_y-current_length_dualvarl(3),nDelB)];
con1 = mat2cell(con1, current_length_y, ones(nDelB,1));
sense = char('E' * ones(1,nDelB));
rhs = zeros(nDelB,1);
LP.addIndicators(ind, ones(nDelB,1), con1, sense, rhs);

delBinf = find(ubbig{3}(delCand1));
nDelB = numel(delBinf);
ind = history_length+delBinf;
con = sparse(delCand1(delBinf), 1:nDelB , ones(nDelB,1), n1, nDelB);
con1 = [sparse(history_length_dualvaru(3),nDelB); con; sparse(current_length_y-current_length_dualvaru(3),nDelB)];
con1 = mat2cell(con1, current_length_y, ones(nDelB,1));
sense = char('E' * ones(1,nDelB));
rhs = zeros(nDelB,1);
LP.addIndicators(ind, ones(nDelB,1), con1, sense, rhs);


delset0 = setdiff(find(lbbig{3}(1:n1)),delCand1);
LP.Model.ub(history_length_dualvarl(3)+delset0) = 0;
delset0 = setdiff(find(ubbig{3}(1:n1)),delCand1);
LP.Model.ub(history_length_dualvaru(3)+delset0) = 0;


the_ub{3}=-(((~ubbig{3}(:)))' .* ub{3}(:)');
the_lb{3}=((~lbbig{3}(:)))' .* lb{3}(:)';
LP.addRows(0, [sparse(1,history_length_v(3)), models{3}.c',sparse(1,history_length_dualvarl(3)-current_length_v(3)),...
    the_lb{3},the_ub{3},sparse(1,current_length_y-current_length_dualvaru(3))], 0);    

lbbig{4} = lb{4} <= -999;
ubbig{4} = ub{4} >= 999;

n1=length(models{4}.rxns);
sense = char('E' * ones(1,nDelCand));
rhs = zeros(nDelCand,1);
ind_del = (history_length+1:history_length+nDelCand)';
con = sparse(delCand1,1:nDelCand,ones(nDelCand,1),n1,nDelCand); % c
con1 = [sparse(history_length_v(4),nDelCand);con; sparse(current_length_y-current_length_v(4),nDelCand)]; 
con1 = mat2cell(con1, current_length_y, ones(nDelCand,1));
LP.addIndicators(ind_del, zeros(nDelCand,1), con1, sense, rhs);


delBinf = find(lbbig{4}(delCand1));       
nDelB = numel(delBinf);
ind = history_length+delBinf;
con = sparse(delCand1(delBinf), 1:nDelB , ones(nDelB,1), n1, nDelB);
con1 = [sparse(history_length_dualvarl(4),nDelB); con; sparse(current_length_y-current_length_dualvarl(4),nDelB)];
con1 = mat2cell(con1, current_length_y, ones(nDelB,1));
sense = char('E' * ones(1,nDelB));
rhs = zeros(nDelB,1);
LP.addIndicators(ind, ones(nDelB,1), con1, sense, rhs);

delBinf = find(ubbig{4}(delCand1));
nDelB = numel(delBinf);
ind = history_length+delBinf;
con = sparse(delCand1(delBinf), 1:nDelB , ones(nDelB,1), n1, nDelB);
con1 = [sparse(history_length_dualvaru(4),nDelB); con; sparse(current_length_y-current_length_dualvaru(4),nDelB)];
con1 = mat2cell(con1, current_length_y, ones(nDelB,1));
sense = char('E' * ones(1,nDelB));
rhs = zeros(nDelB,1);
LP.addIndicators(ind, ones(nDelB,1), con1, sense, rhs);


delset0 = setdiff(find(lbbig{4}(1:n1)),delCand1);
LP.Model.ub(history_length_dualvarl(4)+delset0) = 0;
delset0 = setdiff(find(ubbig{4}(1:n1)),delCand1);
LP.Model.ub(history_length_dualvaru(4)+delset0) = 0;


the_ub{4}=-(((~ubbig{4}(:)))' .* ub{4}(:)');
the_lb{4}=((~lbbig{4}(:)))' .* lb{4}(:)';
LP.addRows(0, [sparse(1,history_length_v(4)), models{4}.c',sparse(1,history_length_dualvarl(4)-current_length_v(4)),...
    the_lb{4},the_ub{4},sparse(1,current_length_y-current_length_dualvaru(4))], 0);  


LP.Model.obj=zeros(length(LP.Model.lb),1);
if defaultobjflag
    LP.Model.obj(history_length_y+1:history_length_y+nDelCand)=1;
    LP.Model.sense='minimize';
else
    LP.Model.obj(outconstr.objind)=outconstr.objcoef;
    LP.Model.sense=outconstr.sense;  
end

ind = (history_length+1:history_length+nDelCand)';
con = sparse(delCand1,1:nDelCand,ones(nDelCand,1),n1,nDelCand);
con1 = [con; sparse(history_length-n1+nDelCand,nDelCand)];
con1 = mat2cell(con1, history_length+nDelCand, ones(nDelCand,1));
sense = char('E' * ones(1,nDelCand));
rhs = zeros(nDelCand,1);
LP.addIndicators(ind, zeros(nDelCand,1), con1, sense, rhs);

rpnt1=length(LP.Model.rhs);
theModel0=LP.Model;


the_sol.x=[];
the_sol.status=1;
hassolution=false;
numsolind=0;
numsol2=0;
solution=struct('fluxes',[],'growth',[],'delCand',[],...
    'rxn_del',[],'delinds',[],'target',[],'status',[]','pool',struct);

solution.delCand1=delCand1;

t1=tic;
solution.timeexceed=0;solution.caseexceed=0;
solution.delCand=delcand;


while the_sol.status~=103&(~solution.caseexceed)&(~solution.timeexceed)
    ts1=tic;
    LP.writeModel('test.lp');
 %  LP.Param.mip.tolerances.mipgap.Cur=optionstcf.relgap;
    the_sol=LP.solve();
  
    ts2=toc(ts1); 
         
    if isfield(the_sol,'x')     
        numsol2=numsol2+1; 
        solution.status{numsol2,1}=the_sol.status;  
        solution.X{numsol2}=the_sol.x;
        delindex=find(the_sol.x(history_length_y+1:history_length_y+nDelCand));
        delinds=delCand1(delindex);
        solution.rxn_del{numsol2,1}=models{1}.rxns(delinds);
        solution.delind{numsol2,1}=delinds;
        solution.status{numsol2,1}=the_sol.status;    
        for i = 1: nc    
                hassolution=true;
                solution.fluxes{numsol2,1}(:,i)=the_sol.x(history_length_v(i)+1:current_length_v(i)); 
                solution.growth(numsol2,i)=solution.fluxes{numsol2}(findRxnIDs(models{i},optionstcf.biomassrxn),i);
%                save('solution_tmp.mat', 'solution','LP')
        end
        LP2=Cplex();
        LP2.Model=theModel0;
        
        they=history_length_y+1:history_length_y+nDelCand;
        they2=they(setdiff(1:nDelCand,delindex));
        LP2.Model.ub(they2)=0;
        the_sol2=LP2.solve();
        if isfield(the_sol2,'x')
            delindex2=find(the_sol2.x(history_length_y+1:history_length_y+nDelCand));
            delinds2=delCand1(delindex2);
            solution.rxn_del2{numsol2,1}=models{1}.rxns(delinds2);        
        end
        totalmod_vec2=sparse(1,length(LP.Model.lb));
        totalmod_vec2(history_length_y+find(the_sol.x(history_length_y+1:history_length_y+nDelCand)))=1;
        rhs=sum(the_sol.x(history_length_y+1:history_length_y+nDelCand))-1;
        LP.addRows(0,totalmod_vec2,rhs);    

    end
    t2=toc(t1);
    solution.timeexceed=(t2>=optionstcf.totalTime);
    solution.caseexceed=(numsol2>=optionstcf.totalCasesNeeded);
    if ~isfield(the_sol,'x')

    end          
        
end








