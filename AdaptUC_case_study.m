
clc;clear;
initCobraToolbox(false)
%setparm
% STEP0: load the model

load iML1515.mat
model=iML1515;

model=changeObjective(model,'BIOMASS_Ec_iML1515_core_75p37M');
%model.csense=model.csense';
% adding methanol gene 

% adding metabolits 
%newmodel =addMetabolite(model,'meoh_c','methanol','CH4O');
%newmodel =addMetabolite(model,'fald_c','formaldehyde','CH2O');
newmodel =addMetabolite(model,'hexulose6p_c','hexulose6p','C6H11O9P');
%newmodel =addMetabolite(newmodel,'meoh_e','methanol','CH4O');
% adding reactions
newmodel =addReaction(newmodel,'FDH','for_c + nad_c --> co2_c + nadh_c');
newmodel =addReaction(newmodel,'H6PI','hexulose6p_c <=> f6p_c');
newmodel =addReaction(newmodel,'H6PS','fald_c + ru5p__D_c  <=> hexulose6p_c');
newmodel =addReaction(newmodel,'MEDH','nad_c + meoh_c  <=> fald_c + nadh_c');
newmodel =addReaction(newmodel,'FRMA','fald_c + nad_c --> for_c + nadh_c');
newmodel=changeRxnBounds(newmodel,{'PFL','OBTFL','DRPA','PAI2T'},0,'b');
%newmodel =addReaction(newmodel,'EX_meoh_e','meoh_e <=>');
% turn off ATPM reaction
newmodel =changeRxnBounds(newmodel,'ATPM',0,'l');
newmodel.csense=ones(1,length(newmodel.mets))*'E';
newmodel=changeRxnBounds(newmodel,'EX_glc__D_e',0,'b');


% STEP1: define task
c1_exchange='EX_meoh_e';
cosub_exchange='EX_xyl__D_e';
load delCand.mat
optionstcf.totalCasesNeeded=1;
optionstcf.totalTime=200;
optionstcf.biomassrxn='BIOMASS_Ec_iML1515_core_75p37M';
optionstcf.bioind0=2669;
optionstcf.c1lb=-20;
optionstcf.maxdel=5;
optionstcf.relgap=0.01;
optionstcf.co_sub_growthRatio=0.2;
optionstcf.c1_growthRatio=0.1;
optionstcf.co_sub_improve=2;

%STEP2: define models

model0=changeRxnBounds(newmodel,cosub_exchange,-10,'b');
model0=changeRxnBounds(model0,c1_exchange,0,'b');
model0=changeObjective(model0,optionstcf.biomassrxn);
sol0=optimizeCbModel(model0);

bioind0=findRxnIDs(newmodel,optionstcf.biomassrxn);
ind_c11=findRxnIDs(newmodel,c1_exchange);
ind_cosub=findRxnIDs(newmodel,cosub_exchange);
% condition 1: the strain can grow in media with methane as the solo carbon
% get the default bounds 
model1=newmodel;
model1=changeRxnBounds(model1,c1_exchange,-10,'l');
model1=changeRxnBounds(model1,cosub_exchange,0,'l');
model1=changeObjective(model1,optionstcf.biomassrxn);
sol1=optimizeCbModel(model1);%check the growth


% condition 2: the strain can grow in media with co-substrate as the solo carbon 
model2=newmodel;
model2=changeRxnBounds(model2,c1_exchange,0,'l');
model2=changeRxnBounds(model2,cosub_exchange,-10,'l');
model2=changeObjective(model2,optionstcf.biomassrxn);


% condition 3 
model3=newmodel;
%ind_bio3=find(strcmp(model3.rxns,optionstcf.biomassrxn));
model3=changeRxnBounds(model3,c1_exchange,-1000,'l');
model3=changeRxnBounds(model3,cosub_exchange,-10,'l');
model3=changeRxnBounds(model3,optionstcf.biomassrxn,0.95*sol0.f,'l');
model3=changeObjective(model3,c1_exchange);


% condition 3: the strain can grow in media with both co-substrat and 

model4=newmodel;
model4=changeRxnBounds(model4,c1_exchange,-10,'l');
model4=changeRxnBounds(model4,cosub_exchange,-10,'l');
model4=changeObjective(model4,optionstcf.biomassrxn);


% models
models=[];
models{1}=model1;models{2}=model2;models{3}=model3;
models{4}=model4;

outconstr=[];
% outer problem 
outconstr.lb=-inf(4,1);outconstr.ub=inf(4,1);
outconstr.S=sparse(4,length(newmodel.rxns)*4);

ns=0;
% for sub problem 1 
ns=ns+1;
outconstr.S(ns,bioind0)=1;
outconstr.lb(ns)=optionstcf.c1_growthRatio*sol1.f;  % obj_model1 >= 0.1*max  for condition 1
% for sub problem 2 growth after deletions 
ns=ns+1;
outconstr.S(ns,length(newmodel.rxns)+bioind0)=1;
outconstr.ub(ns)=optionstcf.co_sub_growthRatio*sol0.f; 
%outconstr.ub(2)=0.05; 
% for sub problem 3
ns=ns+1;
outconstr.S(ns,length(newmodel.rxns)*2+ind_c11)=1;
outconstr.lb(ns)=optionstcf.c1lb;%lb<c1_model2 <= ub  for condition 2 

% for sub problem 4
ns=ns+1;
outconstr.S(ns,bioind0)=-optionstcf.co_sub_improve;
outconstr.S(ns,length(newmodel.rxns)*3+bioind0)=1;  % with co-substrate, its growth is more than 2 times of C1 alone
outconstr.lb(ns)=0; 
outconstr.sense='minimize';
outconstr.objcoef=1;   % for condition 3

%STEP3: solve AdaptUC problem
[solution,LP]=multi_condition_design0(models,outconstr,delCand,optionstcf);


% print the results

for i=1:length(solution.rxn_del)
    fprintf('the reactions to be deleted for solution %d: %s \n',i,cell2str_2(solution.rxn_del{i}));
end






















