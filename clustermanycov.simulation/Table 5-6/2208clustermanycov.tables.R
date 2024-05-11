#########################################################################################
## R CODE FOR Cluster Many Covariates
## Cluster Robust Inference in Linear Regression Models with Many Covariates
## DATE: 14-August-2022
## GENERATE TABLES IN LATEX FOR PAPER
#########################################################################################
# source("tables.R")
#########################################################################################
rm(list=ls(all=TRUE)); require("Hmisc");require("xtable")

fixformat = function(x) format(x,format="f",digits=1,nsmall=3)

#########################################################################################
## LOAD dta & CONSTRUCT dta & TABLES LATEX
#########################################################################################
n.cgroup = c(8,8)
colheads = c("","ratio","LZ","HC1","HC2","HC3","LOO","CR","Neg")

q975 = qnorm(.975)

models = 5:6

for (m in models){
  dta = read.csv(paste0("output/output_m",m,".csv"))
  
  EC.tmp = cbind( dta[,"K"],
                  dta[,"ratio"],
                  -q975 <= dta[,"T_LZ"] & dta[,"T_LZ"] <= q975,
                  -q975 <= dta[,"T_HC1"] & dta[,"T_HC1"] <= q975,
                  -q975 <= dta[,"T_HC2"] & dta[,"T_HC2"] <= q975,
                  -q975 <= dta[,"T_HC3"] & dta[,"T_HC3"] <= q975,
                  -q975 <= dta[,"T_LOO"] & dta[,"T_LOO"] <= q975,
                  -q975 <= dta[,"T_CR"] & dta[,"T_CR"] <= q975,
                  dta[,"Neg"]
                  
  )
  out.EC = apply(EC.tmp,2,function(x) {x[is.infinite(x)]=NA; by(x,list(EC.tmp[,1]),mean,na.rm=TRUE)});
  out.EC = out.EC[order(out.EC[,1]),];
  
  IL.tmp = cbind(dta[,"K"],
                 dta[,"ratio"],
                 2*q975*sqrt(dta[,"V_LZ"]),
                 2*q975*sqrt(dta[,"V_HC1"]),
                 2*q975*sqrt(dta[,"V_HC2"]),
                 2*q975*sqrt(dta[,"V_HC3"]),
                 2*q975*sqrt(dta[,"V_LOO"]),
                 2*q975*sqrt(dta[,"V_CR"]),
                 dta[,"Neg"]
  )
  out.IL = apply(IL.tmp,2,function(x) {x[is.infinite(x)]=NA; by(x,list(EC.tmp[,1]),mean,na.rm=TRUE)})
  out.IL = out.IL[order(out.IL[,1]),];
  
  rowname = paste0("K=",fixformat(out.EC[,1])); n.rgroup = length(rowname)
  
  row.names(out.EC) = rowname
  colnames(out.EC) = colheads
  row.names(out.IL) = rowname
  colnames(out.IL) = colheads
  
  out.EC = apply(out.EC[,c(-1)],2, fixformat)
  out.IL = apply(out.IL[,c(-1)],2, fixformat)
  
  out.EC = xtable(out.EC,digits=2)
  out.IL = xtable(out.IL,digits=2)
  
  name_EC_latex = paste0("output/output_m",m,"_EC.tex")
  print(out.EC, type="latex", file=name_EC_latex)
  name_IL_latex = paste0("output/output_m",m,"_IL.tex")
  print(out.IL, type="latex", file=name_IL_latex)
  
  message("Model ",m," Completed.")
}

