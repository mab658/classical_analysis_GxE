/* Set up working environment */;
dm 'log;clear;output;clear;odsresults;clear;';
options ls=120 ps=1500 nocenter nodate pageno=1;
ods graphics off;
ods noresults;

/* A defined macro that fit  Linear Mixed Model to several analytical variables */;

%macro two_way_GE_data(dsn);

/* Sort the raw data on multiple variables */;
proc sort data = &dsn;
	by  env trial trait;
run;

/*Fitting the model to each agronomic traits except dry matter*/;
title "First stage analysis to get adjusted mean from individual trials";
ods html body ="&outpath&dsn._adjusted_mean_output.html";
ods output lsmeans=&dsn._lsmean(keep = loc year env trial gen trait estimate Stderr);
proc mixed data = &dsn noitprint noinfo noclprint covtest update;
     by env trial trait;
     class rep gen ;
     model y = gen prop_hav/ddfm=satterth outp = &dsn._studresid(keep = loc--Pred StudentResid) residual; /* outp returns conditional internally studentized residual*/
     random rep;
	 lsmeans gen;
	 ods exclude lsmeans;
	 where trait ne "dm";
run;
quit;
ods output close;


/*Fitting the model to dry matter*/;
ods output lsmeans=&dsn._dm_lsmean(keep=loc year env trial gen trait estimate Stderr);
proc mixed data = &dsn noitprint noinfo noclprint covtest update;
     by env trial trait;
     class rep gen;
     model y=gen/ddfm=satterth outp=&dsn._dm_studresid(keep=loc--Pred StudentResid) residual; /* outp returns conditional internally studentized residual*/
     random rep;
	 lsmeans gen;
	 ods exclude lsmeans;
	 where trait in("dm");
run; 
quit;
ods output close;

/* Append studentized residual of other traits and dm as a single file */

proc append base = &dsn._studresid data = &dsn._dm_studresid force;
run;


/* Append lsmeans estimates  of other traits and dm as a single file */

proc append base = &dsn._lsmean data = &dsn._dm_lsmean force;
run;

data &dsn._lsmean(rename=(estimate=y));
	set &dsn._lsmean;
	var_mean=Stderr**2;
	w=1/var_mean;
run;

title "One step analysis showing adjusted mean of genotypes and weight for second stage analysis";
proc sort data= &dsn._lsmean;
	by gen env;
run;

/* create two-way table of gen x env for the traits */;
Proc transpose data = &dsn._lsmean Out = gen_env_trait(drop=_name_);
     by gen env;
	 id trait;
	 var y;
run;
proc print data=gen_env_trait;run;


/* Delete intermediate dataset to free memory space raw_pheno_lsmean  */
proc datasets lib=work nolist;
	delete &dsn._dm_studresid &dsn._dm_lsmean;
run;
quit;
/*%exportdat(&dsn._lsmean); */;
%exportdat(gen_env_trait);

%mend two_way_GE_data;
ods html close;
