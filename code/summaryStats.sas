

/* A defined macro that fit  Linear Mixed Model to several analytical variables */;

%macro summaryStats(dsn);
/* Set up working environment */;
dm 'log;clear;output;clear;odsresults;clear;';
options ls=120 ps=1500 nocenter nodate pageno=1;
ods graphics off;
ods noresults;

/* Sort the raw data on multiple variables */;
proc sort data = &dsn;
	by env trial trait;
run;

/* Compute the summary statistics of the the data */;
proc means data=&dsn nway noprint ;
	class env trait;
	var y;
	output out=&dsn._summ_stats(drop=_type_ _freq_)min=min max=max mean=mean std=std;
	id trial;
run;


title "Model fitting output of analytical variable for individual envs";
/*Fitting the model to each agronomic trait except dry matter (DM) for all envs*/;
/* outp returns conditional internally studentized residual*/;
ods html body ="&outpath&dsn._env_specific_analysis_output.html";
ods output covparms=&dsn._varcomp(drop = StdErr ZValue ProbZ rename=(CovParm = effect Estimate = variance));
ods output solutionR=&dsn._blup(keep = env trial trait gen estimate StdErrPred);
proc mixed data = &dsn noitprint noinfo noclprint covtest update;
     by env trial trait;
     class rep gen;
     model y = rep prop_hav/ddfm=satterth outp = &dsn._studresid(keep=env--Pred StudentResid) residual;
     random gen/solution;
	 ods listing exclude SolutionR;
	 where  trait not in("dm");	 
run;
quit;


/*Fitting the model to dry matter trait (DM) */;
/* outp returns conditional internally studentized residual*/
ods output covparms = &dsn._dm_varcomp(drop = StdErr ZValue ProbZ rename=(CovParm = effect Estimate = variance));
ods output solutionR=&dsn._dm_blup(keep = env trial trait gen estimate StdErrPred);
proc mixed data = &dsn noitprint noinfo noclprint covtest update;
     by env trial trait;
     class rep gen ;
     model y=rep/ddfm=satterth outp = &dsn._dm_studresid(keep=env--Pred StudentResid)residual;
     random gen/solution;
	 ods listing exclude SolutionR;
	  where trait in("dm");
run;
quit;


/* Append variance component estimates  of other traits and dm as a single file */

proc append base=&dsn._varcomp data=&dsn._dm_varcomp force;
run;


/* Append blup estimates  of other traits and dm as a single file */;
proc append base=&dsn._blup data=&dsn._dm_blup force;
run;
/* calculate the prediction error variance */;
data &dsn._blup;
	set &dsn._blup;
	pev = StdErrPred**2;
run;

/* calculate the average prediction error variance by env and trait */;
proc means data=&dsn._blup nway noprint ;
	class env trait;
	var pev;
	output out=&dsn._ave_pev(drop=_type_ _freq_)mean=;
	id trial;
run;


/* Append studentized residual of other traits and dm as a single file */

proc append base = &dsn._studresid data = &dsn._dm_studresid force;
run;


/* remove intermediate dataset to free memory space */;
proc datasets lib=work nolist nodetails;
	delete  &dsn._dm_varcomp &dsn._dm_blup &dsn._dm_studresid;
run;
quit;


/* Sort combined variance component estimate */;

proc sort data = &dsn._varcomp;
	by env trial trait;
run;

Proc transpose data = &dsn._varcomp Out = &dsn._varcomp
      (rename = (gen = gen_var residual = res_var) drop=_name_);
     by env trial trait;
	 Id effect ;
	 var variance;
run;


/* Merge the summary statistics, variance estimates and average PEV by env and trait*/;

data &dsn._varcomp;
	merge &dsn._summ_stats &dsn._varcomp &dsn._ave_pev;
	by env trait;
run;

/* Calculate the total phenotypic variance and broad-sense heritability on plot basis */;
data &dsn._trial_summary;
	set &dsn._varcomp;
	pheno_var = gen_var + res_var;
	if pheno_var in(0,.)then H2 = .;
	else H2 = gen_var/Pheno_var;
	cv = (sqrt(res_var)/mean)*100;
	*GH2 = 1 - (pev/(2 * gen_var));
	Accr = sqrt(1 - (pev/gen_var));
run;

/* Summary statistics across all the envs */;

proc sql;
   create table &dsn._overall_summary as 
   select trait, min(min) as min, max(max) as max, mean(mean) as mean, min(H2) as min_H2, max(H2) as max_H2, 
	min(cv) as min_CV, max(cv) as max_CV,
	min(Accr) as min_Accr, max(Accr) as max_Accr
   from &dsn._trial_summary
   group by trait;
quit;


title "Summary output of analytical variables for individual trials";
proc print data=&dsn._trial_summary;run;
title ;

title "Summary output of analytical variables across all envs";
proc print data=&dsn._overall_summary; run;

/* Export the summary output to files */;
%exportdat(&dsn._trial_summary);
%exportdat(&dsn._overall_summary);
%mend summaryStats;

ods html close;
