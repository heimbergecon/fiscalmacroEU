********************************************************************************
* Replication Code: 
* Macroeconomic and distributional effects of fiscal consolidation measures in EU countries
* Philipp Heimberger & Anna Matzner
* July 2026
* Stata v19
********************************************************************************

* Settings
clear 

*******************************************************************************
********* INSERT PATH HERE: location of 00_MasterDoFile ***********************
*******************************************************************************
global PATH "C:\Users\matzner\OneDrive - Wiener Institut für internationale Wirtschaftsvergleiche\Projekte\Dezernat Zukunft Projekte\WP1_Fiscal Consolidation\EMPN I Replication Package Submission July 2026"

cd "$PATH"
global FIGUREDIR "$PATH\figures"
global FIGUREDIRC "$PATH\figures_combined"
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
do "$PATH\01_Baseline.do"

// Figure 2
do "$PATH\02_StateDependentLP.do"

// Robustness
do "$PATH\03_Robustness_AIPW.do"
do "$PATH\04_Robustness_All.do"

