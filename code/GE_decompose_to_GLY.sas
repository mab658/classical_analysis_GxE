
/* Set up working environment */;
dm 'log;clear;output;clear;odsresults;clear;';
options ls=120 ps=1500 nocenter nodate pageno=1;
ods graphics off;
ods noresults;


/* A defined macro that fit  Linear Mixed Model to several analytical variables 
This model decomposes GE variance into GLY component*/;

%macro GE_decompose_to_GLY(dsn);
/* Sort the raw data on multiple variables */;
proc sort data = &dsn;
	by trait ;
run;

title "Decomposition of GE variance to GL, GY, and GLY Output";
/*Fitting the reduced model to each agronomic traits assuming variance homogeneity */;
/*Fitting mixed model decomposing GE variance for other traits other than dry matter */

ods html body ="&outpath&dsn._Decompose_GE_variance_into_GL_GY_GLY_component.html";
proc mixed data = &dsn noitprint noclprint covtest update;
	by trait;
	class loc year rep gen;
	model y = loc year rep(loc*year) prop_hav /ddfm=satterth;
	random int loc year loc*year/subject=gen;
	repeated /group = loc*year;
	where  trait not in("dm");
run;
quit;


/*Fitting mixed model decomposing GE variance for the dry matter */

proc mixed data = &dsn noitprint noclprint covtest update scoring=18 maxiter=100 convh=1e-4;
	by trait;
	class loc year rep gen;
	model y = loc year rep(loc*year) /ddfm=satterth;
	random int loc year loc*year/subject=gen;
	repeated /group = loc*year;
	where trait in("dm");
run;
quit;
run;
%mend GE_decompose_to_GLY;
ods html close;
