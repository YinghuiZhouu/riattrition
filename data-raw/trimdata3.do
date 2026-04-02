# delimit ;

set more 1;

set mem 400000;

/* this program will keep the final estimation sample */;
/* keep only people who have non-missing outcomes continuously for 208 weeks */;




capture program drop impute; /* delete the program if in memory */;
program define impute; /* write a program that will impute missing values with the average */;
	quietly su `1';
	if $NOBS~=r(N) {; /* if number of obs is different from total number in dataset */; 
		egen `1'_i=mean(`1'); /* create new variable that is the mean */;
		replace `1'_i=`1' if `1'~=.; /* if the orifinal variable is not missing, use that instead */;
		compare `1' `1'_i; /* compare the two variables */;
	}; /* end if */;
end; /* end program  */;

use trimdata2, clear; /* use the data from the previous program */;

count; /* how many observations */;
sort mprid week; /* sort the data */;
keep if hwh~=. & earnh~=.; /* every person-week observation must have non-missing outcoems */;
quietly by mprid: keep if _N==208;
count; /* how many observations now? */;

quietly su treatmnt; /* how many non-missing obs? */;
global NOBS=r(N); /* create global macro */;

tab race_eth, gen(rdum); /* tabulate race, create dummies */;
tab marriage, gen(mardum); /* tabulate marriage, create dummies */;
tab hh_inc, gen(hincdum); /* tabulate household income cat, create dummies */;
tab pers_inc, gen(pincdum); /* tabulate personal income categories, create dummies */;

replace nchld=0 if haschld==0; /* if no children, then number of child=0 */;
replace mosinjob=0 if mosinjob==.; /* if missing, set equal to 0 */;
replace hrswk_jr=0 if hrswk_jr==.; /* if missing, set equal to 0 */;
replace hrwager=0 if hrwager==.; /* if missing, set equal to 0 */;
replace wkearnr=0 if wkearnr==.; /* if missing, set equal to 0 */;


/* impute missing values */;

impute hwh;
impute earnh;
impute treatmnt;
impute female;
impute age_cat;
impute rdum1;
impute rdum2;
impute rdum3;
impute rdum4;
impute haschld;
impute hgc_moth;
impute hgc_fath;
impute mardum1;
impute mardum2;
impute mardum3;
impute mardum4;
impute nchld;
impute hgc;
impute yr_work1;
impute earn_yr;
impute mosinjob;
impute hrswk_jr;
impute hrwager;
impute wkearnr;
impute hincdum1;
impute hincdum2;
impute hincdum3;
impute hincdum4;
impute hincdum5;
impute pincdum1;
impute pincdum2;
impute pincdum3;
impute pincdum4;
impute evarrst1;
impute currjob;

save trimdata3, replace;


