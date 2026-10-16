*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*
*	Replication for: Misperceived Social Norms
*	L. Bursztyn, A. Gonzalez, D. Yanagizawa-Drott
*	-------
*	Setup Do File: to be run once at the beginning 
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*

	version 15.1
	clear all
	set more off
	
* ------------------------------------------------------------------------------
* Change main directory
* ------------------------------------------------------------------------------
	
	global maindirectory = "C:/Users/arthu/Desktop/Replication_Package_LiangyuanZeng" 

	cd "${maindirectory}"

* ------------------------------------------------------------------------------
* Set up folder structure
* ------------------------------------------------------------------------------
	cap program drop confirmdir
	program define confirmdir, rclass
		local cwd `"`c(pwd)'"'
		qui cap cd `"`1'"'
		local confirmdir=_rc 
		qui cd `"`cwd'"'
		return local confirmdir `"`confirmdir'"'
	end 


	* Define directory to be checked or created
	global codedir  = "${maindirectory}/code"
	global datadir  = "${maindirectory}/data"
	global outdir   = "${maindirectory}/output"
	global adobase  = "${maindirectory}"

	* Create folder if it doesn't exits
	confirmdir "${datadir}"
	if `r(confirmdir)'~=0{
		cap noisily mkdir "${datadir}"
	}
	confirmdir "${datadir}/raw_data"
	if `r(confirmdir)'~=0{
		cap noisily mkdir "${datadir}/raw_data"
	}
	confirmdir "${datadir}/clean_data"
	if `r(confirmdir)'~=0{
		cap noisily mkdir "${datadir}/clean_data"
	}
	
	
	confirmdir "${codedir}/01_cleaning"
	if `r(confirmdir)'~=0{
		cap noisily mkdir "${codedir}/01_cleaning"
	}
	confirmdir "${codedir}/02_analysis"
	if `r(confirmdir)'~=0{
		cap noisily mkdir "${codedir}/02_analysis"
	}

	confirmdir "${outdir}"
	if `r(confirmdir)'~=0{
		cap noisily mkdir "${outdir}"
	}
	
	
	* All required Stata packages are available in the /ado/plus folder
	
	adopath ++ "$adobase/ado/plus"

	mata: mata mlib index

