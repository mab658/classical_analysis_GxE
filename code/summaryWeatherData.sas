/* Set up working environment */;
dm 'log;clear;output;clear;odsresults;clear;';
options ls=120 ps=1500 nocenter nodate pageno=1;
ods graphics off;
ods noresults;

/* Define macro variable to set the path for the output of analysis */
%let dirpath = C:\Users\mab658\Documents\Moshood_PhD_Research_Work\weather_data\;

/* Define a macro to read spreadsheet file */;

%macro importdat(dsn);
proc import out=&dsn
	datafile="&dirpath&dsn..csv"
	dbms=csv replace;
	getnames=yes;
	guessingrows=35269;
run;
%mend importdat;

/* Macro to export some output dataset to an output directory */;
%macro exportdat(dsn);
proc export data= &dsn
	outfile="&dirpath&dsn..csv"
	dbms=csv replace;
run;
%mend exportdat;

/* Import weather data */;
%importdat(weatherData);


/* get the summary across all trials */;
proc means data=weatherData nway noprint;
	class trial_name;
	var crop_growth_cycle Srad--tmean wind--rh_diff;
	output out=weatherSumm_env(drop=_type_ _freq_ rename=(wind=wind_speed)) mean=;
	id env;
run;


/* get the sum of the rain across all the environment */;
proc means data=weatherData nway noprint;
	class trial_name;
	var rain;
	output out=rainSumm_env(drop=_type_ _freq_ rename=(rain=precipitation))sum=;
	id env;
run;

data weatherSumm_env;
 	merge weatherSumm_env rainSumm_env;
 	by trial_name;
run;

/* export the weather dataset */;
%exportdat(weatherSumm_env);


proc print data=weatherSumm_env;run;



