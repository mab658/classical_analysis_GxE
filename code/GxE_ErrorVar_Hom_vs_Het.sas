
/* A defined macro that fit  Linear Mixed Model to several analytical variables 
This model assumes homogeneity of variance across environment*/;

%macro GxE_ErrorVar_Hom_vs_Het(dsn);
/* Set up working environment */;
dm 'log;clear;output;clear;odsresults;clear;';
options ls=120 ps=1500 nocenter nodate pageno=1;
ods graphics off;
ods noresults;

/* Sort the raw data on multiple variables */;
proc sort data = &dsn;
	by trait ;
run;

ods html body ="&outpath&dsn._test_homo_vs_hetero_ErrorVar_in_GxE_output.html";
title "Output of fitting reduced model assuming homogeneity of error variance across environmets";
/*Fitting the reduced model to each agronomic traits assuming variance homogeneity */;
ods output Fitstatistics = reduced_fitstat(where=(descr="-2 Res Log Likelihood") rename=(value=reduced_loglik));
ods output Dimensions = reduced_param(rename=(Descr=Descr1 value=reduced_param)where=(Descr1="Covariance Parameters"));
proc mixed data = &dsn noitprint noclprint covtest update;
	by trait;
	class gen env rep;
	model y=env rep(env)/ddfm=satterth;
	random int env/subject=gen;
run;
quit;


title "Output of fitting full model assuming  heterogeneity of error variance across environments";
/*Fitting the full model to agronomic traits assumming error variance heterogeneity */
ods output Fitstatistics = full_fitstat(where=(descr="-2 Res Log Likelihood") rename=(value=full_loglik));
ods output Dimensions = full_param(rename=(Descr=Descr1 value=full_param)where=(Descr1="Covariance Parameters"));
proc mixed data = &dsn noitprint noclprint covtest update scoring=15 convh=1e-5;
	by trait;
	class gen env rep;
	model y=env rep(env)/ddfm=satterth;
	random int env/subject=gen ;
	repeated /group=env;	
run;
quit;

title ;

/* Append fitstatistics of other traits and dm as a single file
and extract out log likelihood statistics*/

data reduced_model(keep=trait reduced_loglik reduced_param);
	merge reduced_fitstat reduced_param;
	by trait;
run;

/* Append parameter dimensions of other traits and dm as a single file */

data full_model(keep = trait full_loglik full_param);
	merge full_fitstat full_param;
	by trait;
run;

/*********************Likelihood Ratio Test *************************/

proc sql;
	create table deviance as
	select *,(reduced_loglik - full_loglik) as chi, (full_param - reduced_param) as dof
	from Reduced_model as red_mod, full_model as full_mod
	where red_mod.trait = full_mod.trait;

quit;


data &dsn._GxE_hom_vs_het_LRT;
	retain trait reduced_loglik full_loglik chi full_param reduced_param dof prob LRT;
	length LRT $15.;
	set deviance;
	prob= 1-probchi(chi,dof);
	if prob < 0.05 then LRT="Significant";
	else LRT="Non-significant";
run;


title "Output of likelihood ratio test homogeneity vs heterogeneity of GxE error variances";
proc print data= &dsn._GxE_hom_vs_het_LRT;run;
title;

%exportdat(&dsn._GxE_hom_vs_het_LRT);


title "Output of likelihood ratio test for homogeneity vs heterogeneity of residual variances across environments";
proc print data= lrtoutput;run;
title ;
/* delete intermediate dataset to free memory space */;
proc datasets lib=work nolist;
	delete reduced_fitstat reduced_param full_fitstat full_param deviance;
run;
quit;
%mend GxE_ErrorVar_Hom_vs_Het;
ods html close;
