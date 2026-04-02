# delimit ;

set more 1;

set mem 400000;

/* this program will merge on variables from other datasets */;

/* first, merge on treatment/control status */;

tempfile temp1; /* create a temporary dataset name */;

use key_vars, clear; /* use the mathematica data */;
keep mprid treatmnt female age_cat race_eth haschld dsgn_wgt; /* keep only certain variables */;
sort mprid; /* sort the data */;
count if mprid==mprid[_n-1]; /* are there any duplicates? there shouldn't be */;
save `temp1', replace; /* save the data */;

use trimdata1, clear; /* use the data from the previous program */;
sort mprid; /* sort the data */;
merge mprid using `temp1'; /* merge the data by the mprid variable */;
tab _merge ; /* what does the _merge variable look like? */;
keep if _merge==3; /* keep the record only if both datasets contribute to the record */;
drop _merge; /* drop the merge variable */;

save trimdata2, replace;

use baseline, clear; /* use the mathematica data, the baseline data  */;
keep mprid hgc_moth hgc_fath marriage nchld hgc yr_work1 earn_yr mosinjob
hrswk_jr hrwager wkearnr hh_inc pers_inc evarrst1 currjob; /* keep selected variables */;
sort mprid; /* sort the data */;
count if mprid==mprid[_n-1]; /* are there any duplicates? there shouldn't be */;
save `temp1', replace; /* save the data */;

use trimdata2, clear; /* use the data saved above */;
sort mprid; /* sort the data */;
merge mprid using `temp1'; /* merge the data by the mprid variable */;
tab _merge ; /* what does the _merge variable look like? */;
keep if _merge==1 | _merge==3; /* keep the record if trimdata2 contributes, which some will have no baseline info */;
drop _merge; /* drop the merge variable */;

save trimdata2, replace;


