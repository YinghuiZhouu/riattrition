*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	Prepared by: R. Han & R. Wu
*	5/6/2020
*	-------
*	Recruitment Experiment Data Cleaning
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	clear all
	set more off 

	use "${datadir}/raw_data/05_recruitment_exp_raw.dta", clear 

* 	First call 

	* Clean arabic variables

	* Q5 Caller's gender 
	gen caller_gender = "female" if Q5_1 == "أنثى"
	replace caller_gender = "male" if caller_gender == ""
	
	* Q7 Respondent answered the phone
	gen valid_response = Q7_1 == "نعم"

	* age_clean
	destring age_clean_1, gen(age)
	
	* Q10 Education
	gen edu = ""
	replace edu = "middle school" if Q10_1 == "المدرسة الإعدادية"
	replace edu = "primary school" if Q10_1 == "المدرسة الابتدائية"
	replace edu = "high school" if Q10_1 == "المدرسة الثانوية"
	replace edu = "vocational training" if Q10_1 == "تدريب مهني"
	replace edu = "professional degree" if Q10_1 == "شهادة اختصاصية"
	replace edu = "4-year bachelor" if Q10_1 == "شهادة البكالوريوس لمدة 4 سنوات"
	replace edu = "6-year bachelor" if Q10_1 == "شهادة بكالوريوس لمدة 6 سنوات"
	replace edu = "master degree" if Q10_1 == "ماجيستير"
	
	gen higher_edu = 1 if inlist(edu, "professional degree", "4-year bachelor", "6-year bachelor", "master degree")
		replace higher_edu = 0 if higher_edu == . 
		
	encode edu, gen(education)
	
	* Q11 Marital Status
	gen married = Q11_1 == "متزوجة"
		replace married = . if Q11_1 == ""
	
	* Q12 Ever worked for FP
	gen previous_work_fp = Q12_1 == "نعم"
		replace previous_work_fp = . if Q12_1 == "" 
		
	* Q13 Ever worked for other svy company 
	gen previous_work_other = Q13_1 == "نعم"
		replace previous_work_other = . if Q13_1 == "" 
		
	* Q14 Currently employed 
	gen employed_now = Q14_1 == "نعم"
		replace employed_now = . if Q14_1 == "" 
		
	* Q16 Currently employed: work outside home 
	gen outside_work_now = Q16_1 == "خارج المنزل"
		replace outside_work_now = . if Q16_1 == "" 
		
	* Q18 Ever worked outside home
	gen outside_work_ever = Q18_1 == "نعم"
		replace outside_work_ever = . if Q18_1 == "" 
		
	* Q19 Number of children
	gen num_child = Q19_1
		replace num_child = "8" if Q19_1 == "8 أطفال أو أكثر"
	destring num_child, replace 
		
	* Q21 Job Takeup 
	gen takeup = Q21_1 == "نعم"
		replace takeup = . if Q21_1 == "" 
	
* 	Second call 
	
	* Clean arabic variables

	* q13 Job offer choice 
	gen outside_job_takeup = regexm(q13, "240") == 1 if !missing(q13)
	
	* fl_8_do treatment status
	gen treated = fl_8_do == "FL_30" if !missing(fl_8_do)

* 	Third call  
		
	gen outside_job_takeupXshow_up = outside_job_takeup * show_up
	
* Sample restriction: women aged 18-45 
	
	keep if age > 17 & age < 46
	
* Keep relevant variables 
	
	keep age higher_edu education married employed_now num_child outside_job_takeup ///
		 treated show_up outside_job_takeupXshow_up
		 
	label var age 							"Age"
	label var higher_edu 					"College education or higher"
	label var education 					"Highest education"
	label var married 						"Married"
	label var employed_now 					"Currently employed"
	label var num_child 					"Number of children in the household"
	label var outside_job_takeup 			"Take up the outside job"
	label var treated 						"Treatment"
	label var show_up 						"Show up for the outside job"
	label var outside_job_takeupXshow_up 	"Take up and show up for the outside job"	

	save "${datadir}/clean_data/05_recruitment_exp_clean.dta" , replace
	
