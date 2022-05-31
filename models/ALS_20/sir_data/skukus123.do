clear all
* this is just for the model. keep the .csv file since it's quite interesting!

import delimited "owid-covid-data.csv", clear

* keep 3 countries
keep if iso=="GBR" | iso=="KOR" | iso=="USA"

* keep only covid and demographic data
keep iso date *cases *deaths total_tests population aged*

g byte country=1 if iso=="KOR"
replace country=2 if iso=="GBR"
replace country=3 if iso=="USA"
drop iso

rename date sdate
g date = date(sdate,"YMD")
format date %td
drop sdate

order date country *cases *deaths
sort country date
export delimited using "skukus123.csv", replace
