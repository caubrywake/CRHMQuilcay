function [kge,r, relvar,bias] = klinggupta(modelled,observed)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% modelled = mod
% observed = obs
modelled(isnan(observed))=NaN;

cflow=[modelled,observed];

sdmodelled=nanstd(modelled);
sdobserved=nanstd(observed);
 
mmodelled=nanmean(modelled);
mobserved=nanmean(observed);

r=corrcoef(cflow,'rows','pairwise'); r=r(1,2);
relvar=sdmodelled/sdobserved;
bias=mmodelled/mobserved;

%KGE timeseries 
kge=1-  sqrt( ((r-1)^2) + ((relvar-1)^2)  + ((bias-1)^2) );
 
end