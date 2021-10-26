/* Set up working environment */;
dm 'log;clear;output;clear;odsresults;clear;';
options ls=120 ps=1500 nocenter nodate pageno=1;
ods graphics off;
ods noresults;

/* Define macro variable to set the path for import and export data */
%let dirpath = C:\Users\mab658\Documents\research_project\classical_analysis_GxE\output\;

/* Define macro variable to set the path to write html output to docs directory */
%let outpath = C:\Users\mab658\Documents\research_project\classical_analysis_GxE\docs\;

/* Define a macro to read spreadsheet file */;

%macro importdat(dsn);
proc import out=&dsn
	datafile="&dirpath&dsn..csv"
	dbms=csv replace;
	getnames=yes;
	guessingrows=13042;
run;
%mend importdat;


/* Macro to export some output dataset to an output directory */;
%macro exportdat(dataset);
proc export data= &dataset
	outfile="&dirpath&dataset..csv"
	dbms=csv replace;
run;
%mend exportdat;

/* invoke the macro to import raw data in a narrow fromat */;
%importdat(narrowPhenoDat);

/********Step 1: Compute Summmary statistics for key traits of interest subset from the raw data *******/;
data raw_pheno;
	set narrowPhenoDat;
	where trait in("fyld","dm","dyld","hi","tyld");
run;

%include "C:\Users\mab658\Documents\research_project\classical_analysis_GxE\code\summaryStats.sas";
%summaryStats(raw_pheno); /* invoke macro to carry out summary statistics */;



/* From this point the analysis is based on 17 trials after excluding 3 trials based on data quality control
/***********Step 2: Invoke the macro to carry out hypothesis test for the presence of GxE interaction ***************/;

data raw_pheno17;
	set raw_pheno;
	where trial not in ("18UYT36setAKN","19UYT36setAZA","19UYT36setAMK");
run;

%include "C:\Users\mab658\Documents\research_project\classical_analysis_GxE\code\GxE_YesNo.sas";
%GxE_YesNo(raw_pheno17);


/****** Step 3: Invoke macro to carry out hypothesis test on homogeneity vs heterogeneity of error variance in GxE *****/;

%include "C:\Users\mab658\Documents\research_project\classical_analysis_GxE\code\GxE_ErrorVar_Hom_vs_Het.sas";
%GxE_ErrorVar_Hom_vs_Het(raw_pheno17);


/****** Step 4: Invoke macro to fit model and obtain two-way data of adjusted  means 
to be used for fitting linear bilinear  models: Finlay Wilkinson, AMMI, and GGE model in R *****/;

%include "C:\Users\mab658\Documents\research_project\classical_analysis_GxE\code\two_way_GE_data.sas";
%two_way_GE_data(raw_pheno17);


/***********Step 5:  Decompose GE variance into GLY */;
%include "C:\Users\mab658\Documents\research_project\classical_analysis_GxE\code\GE_decompose_to_GLY.sas";
%GE_decompose_to_GLY(raw_pheno17);
