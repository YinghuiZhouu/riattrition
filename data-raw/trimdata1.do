# delimit ;

cd "/Users/peizan/github/riattrition/data-raw";

set more 1; 

set mem 100000;

/* this program will reshape the data from a wide to a long dataset */;

use empl_tl, clear; /* this is the data from mathematica */;

drop workh*; /* don't need this variable */;

reshape long hwh earnh, i(mprid) j(week);

save trimdata1, replace;

