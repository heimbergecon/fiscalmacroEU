********************************************************************************
* Replication Code: 
* Macroeconomic and distributional effects of fiscal consolidation measures in EU countries
* Philipp Heimberger & Anna Matzner
* February 2026
* Stata v19
********************************************************************************

* Settings
clear 

*******************************************************************************
********* INSERT PATH HERE: location of 00_MasterDoFile ***********************
*******************************************************************************
global PATH ""

cd "$PATH"
global FIGUREDIR "$PATH\figures"
global DATA "$PATH\data"
global TABLEDIR "$PATH\tables"

// Graphstyle
grstyle clear
grstyle init
grstyle set grid
global grfont "P052"
graph set window fontface $grfont 
global mblue "0 114 189"
global mred "217 83 25"
global morange "237 177 32"
global mdred "162 20 47"
global mpurple "126 47 142"
global mgreen "119 172 48"

// Data 

use "$DATA\data_select.dta", clear

gen temp = 100*log_REALGDP
drop log_REALGDP
rename temp log_REALGDP

gen temp = 100*log_PRICES
drop log_PRICES
rename temp log_PRICES

encode ccode, gen(country_id)

// Figure 1
do "$CODEDIR\01_Baseline.do"

// Figure 2, A1, A2, A3
do "$CODEDIR\02_Robustness_AIPW.do"
do "$CODEDIR\03_Robustness_All.do"

// Figure 3, A4
do "$CODEDIR\04_StateDependentLP.do"
