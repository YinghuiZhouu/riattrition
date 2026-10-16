*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	Prepared by: R. Han & R. Wu
*	5/6/2020
*	-------
*	Arab Barometer Data Cleaning
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 
	
	use "${datadir}/raw_data/06_arab_barometer_raw.dta", clear 
	
	//Saudi Arabia
	keep if country == 17
	//males
	keep if q1002 == 1
	
	//A married woman can work outside the home.
	
	//Drop if non-response, or if age=998 ("I Don't Know") or is missing
	//0: Missing; 8: I don't know (do not read); 9: Declined to answer (do not read)
	drop if (inlist(q60102, 0, 8, 9)) | mi(q60102) | q1001==998 | mi(q1001)
	
	gen work_outside_ab = (inlist(q60102,1,2))
	
	keep work_outside_ab wt q1001
	
	label var work_outside_ab 	"Men's opinion towards women working outside the home in Saudi Arabia"
	label var wt 				"Weight"
	label var q1001 			"Age"
	
	save "${datadir}/clean_data/06_arab_barometer_clean.dta", replace
		
	
