use "C:\STUDIES\allanK\telomere length\FINALSTATS.dta", clear

****PTSD
gen e57=upsetting1
gen e58=dreamss1
recode e57 0/1=0 2/3=1 9=.
recode e58 0/1=0 2/3=1 9=.
tab1 e57 e58
tab e57 mentaldisorder
egen ptsd = rowtotal(e57 e58)
recode ptsd 2=1

ta ptsd mentaldisorder, col chi
*******Telomere length dataset

*******generating chronic stress
*****<1.375
*****1.375,2.375
*****>2.375

gen chronic2=cstr_index

recode chronic2 min/1.375=1 1.3751/2.375=2 2.3751/max=3

ta chronic2



******recent stress reduced to 3 categories

gen recent3 =.

replace recent3=1 if rstr_index<0.3616

replace recent3=2 if rstr_index>=0.3616 & rstr_index<0.6224

replace recent3=3 if rstr_index>0.62241

tab control recent3, col chi

*****
gen mentaldisorder=.

replace mentaldisorder=1 if control==0

replace mentaldisorder=0 if control==1

**********getting odds ratios

logistic mentaldisorder i.recent3

testparm i.recent3


*********merging the telomere length data

sort patid

merge 1:1 patid using "C:\STUDIES\allanK\telomere length\telomere length data.dta"

*********looking at the baseline telomere length data

hist baselinet*, norm

destring baselinet, gen(tlbase)

destring monthstl, gen(tl12)


********removing the outliers*******

********keeping values between 0 andn 2.2
summ tlbase, d


foreach var in tlbase tl12{

replace `var' =. if `var'>2.2

}

***(15 real changes made, 15 to missing)

***(21 real changes made, 21 to missing)

pnorm tlbase
hist tlbase, norm

pnorm tl12
hist tl12, norm


********generating the difference in the telomere lengths

gen tldiff = (tlbase - tl12)


*******testing for associations btn 
*****IMD and chronic stress
logistic mentaldisorder i.chronic2
testparm i.chronic2

logistic mentaldisorder i.chronic2 i.sex1 i.agecat
testparm i.chronic2

*****chronic stress and telomere length
regress tlbase i.chronic2
testparm i.chronic2

regress tlbase i.chronic2 i.sex1 i.agecat
testparm i.chronic2

*******************12 months TL with chronic stress
regress tl12 i.chronic2
testparm i.chronic2

regress tl12 i.chronic2 i.sex1 i.agecat
testparm i.chronic2


*****recent stress and telomere length
regress tlbase i.recent

regress tlbase i.recent i.sex1

oneway tlbase recent

oneway tl12 recent


**************
regress tldiff i.chronic

regress tldiff i.chronic i.sex1


*******stress and the internalising mental disorders

lab def mentaldisorder 1"cases" 0"controls"

lab val mentaldisorder mentaldisorder

ta mentaldisorder chronic, col chi

****chronic is having a category with 1 person (trying a combination)

gen chronic1 = chronic

recode chronic1 3=2

ta mentaldisorder chronic1, col chi


****recent stress

lab def recent3 1"mild stress" 2"moderate stress" 3"severe stress"

lab val recent3 recent3

ta mentaldisorder recent, col chi

ta mentaldisorder recent3, col chi

logistic mentaldisorder i.recent3 i.sex1

testparm i.recent3


logistic mentaldisorder i.recent3 i.agecat

testparm i.recent3


**** looking at the telomere length differences

regress tldiff i.chronic2
testparm i.chronic2

regress tldiff i.chronic2 i.sex1 i.agecat
testparm i.chronic2


******is the difference between baseline and 12 months statistically significant





ttest tlbase = tl12


*******************ANALYSIS ADVISED BY JACQUI WOMMERSLEY*************************

*****correlation between TL and age

*****baseline TL and age

*****getting a rounded off age

gen age1=round(age, 0.1)

corr tlbase age1

corr tlbase age1 if age<=11

corr tlbase age1 if age>=11.1

***** change in telomere length and chronic stress

oneway tldiff chronic

***** change in telomere length and internalising mental disorder

ttest tldiff, by(mentaldisorder)


*****stepwise regression*******

stepwise, pr(0.1): logistic mentaldisorder recent3 chronic2 rs34517220 tldiff httlpr1 rs35531 httlprrs25531_num stin2vntr_num tlbase tl12
**********************************************************************
**********************************************************************
*******************LGC GENOS* from sheet 2****************************

*****Merging with the TenSNP data

duplicates list patid

merge 1:1 patid using "C:\STUDIES\allanK\telomere length\TenSNP.dta" 

****Merging with the OneSNP data

drop _merge

duplicates list patid

merge 1:1 patid using "C:\STUDIES\allanK\telomere length\OneSNP.dta" 


******testing the association with mentaldisorder

sort sex1

foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

ta mentaldisorder `var', col chi

logistic mentaldisorder i.`var' 

testparm i.`var'

logistic mentaldisorder i.`var' i.sex1 i.agecat

testparm i.`var'

}


*******testing the assocition with chronic stress

foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

ta chronic1 `var' , col chi

}


*******testing association with change in telomere length

foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

oneway tldiff `var' 

}


************ODDS RATIOS *****************

foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

logistic mentaldisorder i.`var' 

testparm i.`var'


**********adjusting for sex****************

logistic mentaldisorder i.`var' i.sex1

testparm i.`var'

}


***************chronic stress



***********telomere length change

foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

logistic chronic1 i.`var' 

testparm i.`var'


**********adjusting for sex****************

logistic chronic1 i.`var' i.sex1

testparm i.`var'

}



********testing for HWE on the eleven SNP****************************************

********testing for HWE FOR HTTLPR & STiN2*********
*****httlpr
*****COMBINED CASES & CONTROLS
genhwi 454 209 35, binvar label(LL LS SS)
*****INDEPENDENTLY CASES & CONTROLS
genhwcci 223 115 20 231 94 15, binvar label( LL LS SS)


****STin2
****T=10 W=12********
*****COMBINED CASES & CONTROLS
genhwi 54 238 366, binvar label(TT TW WW)
*****INDEPENDENTLY CASES & CONTROLS
genhwcci 26 126 189 28 112 177, binvar label(TT TW WW)

************************************************************************
ta mental* rs1386494


genhwcci 18 124 210 6 125 202, binvar label(A:A G:A G:G)

ta mental* rs1843809
genhwcci 94 150 105 55 169 113, binvar label(G:G T:G T:T)

ta mental* rs12696304
genhwcci 46 165 133 41 187 104, binvar label(C:C G:C G:G)

ta mental* rs2736100
genhwcci 82 149 106 57 181 95, binvar label(G:G T:G T:T)

ta mental* rs10936599
genhwcci 334 13 0 316 19 0, binvar label(C:C T:C T:T)

ta mental* rs2853669
genhwcci 5 65 265 3 60 257, binvar label(C:C C:T T:T)

ta mental* rs7726159
genhwcci 19 118 219 12 127 198, binvar label(A:A C:A C:C)

ta mental* rs10069690
genhwcci 42 165 138 31 152 148, binvar label(C:C T:C T:T)

ta mental* rs34517220
genhwcci 90 170 91 84 166 84, binvar label(A:A G:A G:G)

ta mental* rs16847897
genhwcci 27 137 186 17 128 187, binvar label(C:C G:C G:G)

ta mental* genotypers4570625
genhwcci 88 177 84 90 156 88, binvar label(G:G T:G T:T)







**************************************************************


*****modelling IDM,TLDIFF, chronic stress and the SNPs

logistic mentaldisorder i.chronic i.sex1 i.agecat i.rs1386494 i.rs1843809 i.rs12696304 i.rs2736100 i.rs10936599 i.rs2853669 i.rs7726159 i.rs10069690 i.rs34517220 i.rs16847897 i.genotypers4570625
testparm i.chronic2

regress tldiff i.chronic i.sex1 i.agecat i.rs1386494 i.rs1843809 i.rs12696304 i.rs2736100 i.rs10936599 i.rs2853669 i.rs7726159 i.rs10069690 i.rs34517220 i.rs16847897 i.genotypers4570625
testparm i.chronic2


***************************************************************************************

******looking at the chronic stress on TL with the genes

foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

regress tldiff i.chronic2
testparm i.chronic2

regress tldiff i.chronic2 i.`var' 
testparm i.chronic2
}

**********looking at TL the disorder with the effect of the gene************************
foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

logistic mentaldisorder i.chronic2
testparm i.chronic2

logistic mentaldisorder i.chronic2 i.`var' 
testparm i.chronic2
}


foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

logistic mentaldisorder i.`var'

}

foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625{

sort mentaldi*

by mentaldisorder: tab chronic2 `var', col chi

}


***here you will need to ass the change in TL over 12 months in youth with IMD compared to those without IMD?

****variables to use: tldiff, with youth in IMD

summ tldiff if mentaldisorder==1

summ tldiff if mentaldisorder==0

ttest tldiff, by(mentaldisorder)

gen tldiff2= tl12 - tlbase

sort patid

merge 1:1 patid using "C:\STUDIES\allanK\telomere length\trial.dta"




*******looking at telomere length at baseline among controls alone and cases alone

sort mentaldisorder

by mentaldisorder: regress tlbase i.chronic2




**********looking at recent stress and TL, TL change

codebook recent3 tlbase tl12 tldiff

foreach var in tlbase tl12 tldiff{
regress `var' i.recent3
testparm i.recent3
}

***For rs2736100 polymorphism, it seems this may be driving the longer TL among cases. Could you kindly check for me the mean TL among GG vs the rest and check whether the difference is statistically significant.
ta rs2736100

********

ttest tldiff, by(mentaldisorder)


*****************

ttest tlbase, by(mentaldisorder)



ttest tlbase=tl12 if mentaldisorder==1

ttest tlbase=tl12 if mentaldisorder==0

ta chronic2

xi: stepwise, pr(.2): regress tldiff i.chronic2 i.studsite1 i.sex1 i.ses_cat i.childeduc1

xi: stepwise, pr(.2): regress chronic2 i.studsite1 i.sex1 i.ses_cat i.childeduc1


*********************looking at the association of socio demographics with change in TL**********

foreach var in sex1 studsite1   ses_cat {
ttest tldiff, by(`var')
}


foreach var in carededuc3 agecatak childeduc1{

oneway tldiff `var', tab

}

******trying out a the stepwise as advised by Soraya Seedat*****

xi: stepwise, pr(.5): regress tldiff i.sex1 i.studsite1 /*i.carededuc3*/ i.agecatak i.ses_cat i.childeduc1 i.chronic2



/*1) Association between baseline TL and Polymorphisms in each of the TERT and TERC genes
2) Association between TL change  and Polymorphisms in each of the TERT and TERC genes
3) Association between Polymorphisms in each of the TERT and TERC genes and IMDS
4) Association between IMDs and TL (We have determined this before) and how this is influenced by Polymorphisms in each of the TERT and TERC genes.
*/

/*
TERC: rs12696304, rs16847897, rs58465486, rs10936599

TERT: rs2736100, rs7726159, rs2853669, rs10069690
*/

tab1 rs12696304 rs16847897 /*rs58465486*/ rs10936599

*********************************************************************************

***Association between baseline TL and Polymorphisms in each of the 
***TERT and TERC genes

foreach var in rs12696304 rs16847897 /*rs58465486*/ rs10936599 {

oneway tlbase `var', tab

}

foreach var in rs2736100 rs7726159 rs2853669 rs10069690 {

oneway tlbase `var', tab

}

*** Association between TL change  and Polymorphisms in each of the
*** TERT and TERC genes

foreach var in rs12696304 rs16847897 /*rs58465486*/ rs10936599 {

oneway tldiff `var', tab

}

foreach var in rs2736100 rs7726159 rs2853669 rs10069690 {

oneway tldiff `var', tab

}

***3) Association between Polymorphisms in each of the 
***TERT and TERC genes and IMDS

foreach var in rs12696304 rs16847897 /*rs58465486*/ rs10936599 {

tab mentaldisorder `var', col chi

}

foreach var in rs2736100 rs7726159 rs2853669 rs10069690 {

tab mentaldisorder `var', col chi

}



***4) Association between IMDs and TL (We have determined 
***this before) and how this is influenced by Polymorphisms
*** in each of the TERT and TERC genes.

regress  tldiff i.mentaldisorder

foreach var in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 {


/*regress  tldiff i.mentaldisorder i.`var'*/

oneway tlbase `var', tab
}


*****I would like to investigate the moderation effect of TERC and TERT genetic variants on the relationship between IMDs and TL change.
*****telomere length at baseline

regress tlbase i.mentaldisorder 
testparm i.mentaldisorder


regress tlbase i.mentaldisorder i.rs34517220 
testparm i.mentaldisorder


regress tlbase i.mentaldisorder i.rs34517220 i.rs1386494 i.rs1843809 i.rs12696304 i.rs2736100 i.rs10936599 i.rs2853669 i.rs7726159 i.rs10069690 i.rs16847897
testparm i.mentaldisorder


*****telomere length at 12 months

regress tl12 i.mentaldisorder 
testparm i.mentaldisorder


regress tl12 i.mentaldisorder i.rs34517220 
testparm i.mentaldisorder

regress tl12 i.mentaldisorder i.rs34517220 i.rs1386494 i.rs1843809 i.rs12696304 i.rs2736100 i.rs10936599 i.rs2853669 i.rs7726159 i.rs10069690 i.rs16847897
testparm i.mentaldisorder


***********************generating continous scores for the disorders**************************
sort mentaldisorder
*****PTSD
foreach var in upsetting1 dreamss1 {

recode `var' 9=.

}
egen ptsd = rowtotal(upsetting1 dreamss1),m
tab ptsd, m
by mentaldisorder: tab ptsd
*****GAD
foreach var in  overconcernn1 worriess1 edgyy1 irritablee1 tensee1 sleepingp1 {

recode `var' 9=.

}

egen gad=rowtotal(overconcernn1 worriess1 edgyy1 irritablee1 tensee1 sleepingp1), m
ta gad
by mentaldisorder: tab gad
*****MDD
foreach var in depressedm1 pleasurabl1 deathsuici1 worthless1 lowenergy1 cofidence1 workright1 {

recode `var' 9=.

}
egen mdd=rowtotal(depressedm1 pleasurabl1 deathsuici1 worthless1 lowenergy1 cofidence1 workright1), m
tab mdd
by mentaldisorder: tab mdd
************************************************************************
foreach var in ptsd gad mdd {

regress tlbase `var' 

}

foreach var in ptsd gad mdd {

regress tl12 `var' 

}


**********combining the 3 continous scores
egen intdis=rowtotal(ptsd gad mdd), m

ta intdis


regress tlbase intdis/*intdis is a continous score of the 3 disorders  - GAD MDD PTSD*/



regress tl12 intdis/*intdis is a continous score of the 3 disorders  - GAD MDD PTSD*/


*******INTERACTION ANALYSIS - 5httlpr

*logistic mentaldisorder i.recent3
*testparm i.recent3
*logistic mentaldisorder i.httlpr1
logistic mentaldisorder i.recent3##i.chronic2
estimates store A
logistic mentaldisorder i.recent3 i.chronic2
estimates store B

lrtest A B, di st
***********************************************
*logistic mentaldisorder i.httlpr1
logistic mentaldisorder i.recent3##i.genotypers4570625
estimates store A
logistic mentaldisorder i.recent3 i.genotypers4570625
estimates store B

lrtest A B, di st
*******************************************************
logistic mentaldisorder i.recent3##i.HTTLPRrs35531_1
estimates store A
logistic mentaldisorder i.recent3 i.HTTLPRrs35531_1
estimates store B

lrtest A B, di st

*************************************************
logistic mentaldisorder i.chronic2##i.HTTLPRrs35531_1
estimates store A
logistic mentaldisorder i.chronic2 i.HTTLPRrs35531_1
estimates store B

lrtest A B, di st


logistic mentaldisorder i.recent3##i.chronic2##i.HTTLPRrs35531_1
estat ic

*************************
logistic mentaldisorder i.recent3##i.chronic2 i.recent3##i.HTTLPRrs35531_1 i.chronic2##i.HTTLPRrs35531_1  if recent3!=. & chronic2!=. & HTTLPRrs35531_1!=. 

logistic mentaldisorder i.recent3##i.chronic2##i.HTTLPRrs35531_1 i.recent3##i.chronic2##i.HTTLPRrs35531_1 if recent3!=. & chronic2!=. & HTTLPRrs35531_1!=.

estimates store A
logistic mentaldisorder i.recent3 i.chronic2 i.HTTLPRrs35531_1 if recent3!=. & chronic2!=. & HTTLPRrs35531_1!=.
estimates store B

lrtest A B 

******
logistic mentaldisorder i.recent3
testparm i.recent3
logistic mentaldisorder i.httlpr1
logistic mentaldisorder i.recent3##i.httlpr1
estimates store A
logistic mentaldisorder i.recent3 i.httlpr1
estimates store B

lrtest A B
**************************************
logistic mentaldisorder i.recent3##i.genotypers4570625
estimates store A
logistic mentaldisorder i.recent3 i.genotypers4570625
estimates store B

lrtest A B

**************************************

*******TLBASE
/*logistic mentaldisorder tlbase
testparm c.tlbase
logistic mentaldisorder i.httlpr1*/
logistic mentaldisorder i.recent3##c.tlbase
estimates store A
logistic mentaldisorder c.tlbase i.recent3
estimates store B

lrtest A B

***********************
logistic mentaldisorder i.chronic2##c.tl12
estimates store A
logistic mentaldisorder c.tl12 i.chronic2
estimates store B

lrtest A B



********

logistic mentaldisorder i.recent3##c.tlbase
estimates store A
logistic mentaldisorder c.tlbase i.recent3
estimates store B

lrtest A B



logistic mentaldisorder i.recent3
logistic mentaldisorder tlbase
logistic mentaldisorder i.recent3 tlbase



****************************

regress tlbase i.mentaldisorder if rs34517220!=. , allbaselevels 


regress tl12 i.mentaldisorder if rs34517220!=. , allbaselevels 
estimates store A

regress tl12 i.mentaldisorder i.rs34517220, allbaselevels
estimates store B

regress tl12  i.rs34517220, allbaselevels
estimates store C

lrtest A B
lrtest B C



********************************
stepwise, pr(0.1): logistic mentaldisorder recent3 chronic2 rs34517220 tldiff httlpr1 rs35531 httlprrs25531_num stin2vntr_num tlbase tl12


************Assocition btn SNPs and tl12, tldiff
foreach v in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 {

oneway tlbase `v', tab

}


foreach v in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 {

oneway tl12 `v', tab

}

foreach v in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 {

oneway tldiff `v', tab

}

*****interaction between each of the 10 SNPs and TLBASE, TL12 & TLDIFF

****************rs16847897**********************************
regress tlbase i.mentaldisorder##i.rs16847897
estimates store A
regress tlbase i.mentaldisorder i.rs16847897

estimates store B

lrtest A B

****************************rs2736100************************
regress tlbase i.mentaldisorder##i.rs2736100
estimates store A
regress tlbase i.mentaldisorder i.rs2736100
estimates store B

lrtest A B


*****TL BASELINE
foreach v in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 {

regress tlbase i.mentaldisorder##i.`v'
estimates store A
regress tlbase i.mentaldisorder i.`v'
estimates store B

lrtest A B

}

****TLDIFF MONTHS
foreach v in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 {

regress tldiff i.mentaldisorder##i.`v'
estimates store A
regress tldiff i.mentaldisorder i.`v'
estimates store B

lrtest A B

}

*****TL 12 MONTHS
foreach v in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 {

regress tl12 i.mentaldisorder##i.`v'
estimates store A
regress tl12 i.mentaldisorder i.`v'
estimates store B

lrtest A B

}

*******focusing on the rs16847897
regress tl12 i.mentaldisorder##i.rs16847897
estimates store A
regress tl12 i.mentaldisorder i.rs16847897
estimates store B

lrtest A B

*****stratum specific coeffs
regress tl12 i.mentaldisorder##i.rs16847897

*******G:C
lincom 1.mentaldisorder + 1.mentaldisorder#2.rs16847897

*******G:G
lincom 1.mentaldisorder + 1.mentaldisorder#3.rs16847897

****************focusing on the rs2736100
regress tl12 i.mentaldisorder##i.recent3
estimates store A
regress tl12 i.mentaldisorder i.recent3
estimates store B

lrtest A B

regress tlbase i.mentaldisorder##i.recent3
estimates store A
regress tlbase i.mentaldisorder i.recent3
estimates store B

lrtest A B
************************
regress tl12 i.mentaldisorder##i.rs2736100
estimates store A
regress tl12 i.mentaldisorder i.rs2736100
estimates store B

lrtest A B

*****stratum specific coeffs
regress tl12 i.mentaldisorder##i.rs2736100

*******T:G
lincom 1.mentaldisorder + 1.mentaldisorder#2.rs2736100

*******T:T
lincom 1.mentaldisorder + 1.mentaldisorder#3.rs2736100

gen 5HTTLPRrs35531=.
replace 5HTTLPRrs35531=1 if httlprrs25531=="LA/LA"
replace 5HTTLPRrs35531=2 if httlprrs25531=="LA/LG"
replace 5HTTLPRrs35531=3 if httlprrs25531=="LA/SA"
replace 5HTTLPRrs35531=4 if httlprrs25531=="LA/SG"
replace 5HTTLPRrs35531=4 if httlprrs25531=="LASG"
replace 5HTTLPRrs35531=5 if httlprrs25531=="LG/LG"
replace 5HTTLPRrs35531=6 if httlprrs25531=="SA/LG"
replace 5HTTLPRrs35531=6 if httlprrs25531=="SALG"
replace 5HTTLPRrs35531=7 if httlprrs25531=="SA/SA"

lab def 5HTTLPRrs35531 1"LA/LA" 2"LA/LG" 3"LA/SA" 4"LA/SG" 5"LG/LG" 6"SA/LG" 7"SA/SA", replace
lab val 5HTTLPRrs35531 5HTTLPRrs35531




5-HTTLPR, rs35531, 5-HTTLPR/rs35531 and STin2.VNTR

ren mental_disorder mentaldisorder
desc mentaldisorder
lab var mentaldisorder "mental disorder"

lab var rs2736100 "rs2736100 genotype"
lab var rs16847897 "rs16847897 genotype"

lab var tl12 "mean relative telomere length at 12 months"

*******


******TWO-WAY ANOVA
foreach v in rs1386494 rs1843809 rs12696304 rs2736100 rs10936599 rs2853669 rs7726159 rs10069690 rs34517220 rs16847897 genotypers4570625 {

anova tlbase mentaldisorder `v' mentaldisorder##`v'

anova tl12 mentaldisorder `v' mentaldisorder##`v'

anova tldiff mentaldisorder `v' mentaldisorder##`v'

}
********chronic stress
anova tlbase mentaldisorder chronicstress mentaldisorder##chronicstress

anovaplot  chronicstress mentaldisorder , scatter(msym(none)) ylabel(0(0.3)1.5) /*title(Means plot - rs2736100genotype mentaldisorder interaction)*/ yscale(r(0.45 1.2))

anova tl12 mentaldisorder chronicstress mentaldisorder##chronicstress

anovaplot  chronicstress mentaldisorder , scatter(msym(none)) ylabel(0(0.3)1.5) /*title(Means plot - rs2736100genotype mentaldisorder interaction)*/ yscale(r(0.45 1.2))


******rs2736100
anova tl12 mentaldisorder rs2736100 mentaldisorder##rs2736100

anovaplot  rs2736100 mentaldisorder , scatter(msym(none)) ylabel(0(0.3)1.5) /*title(Means plot - rs2736100genotype mentaldisorder interaction)*/ yscale(r(0.45 1.2))

anova tl12 mentaldisorder rs16847897 mentaldisorder##rs16847897

anovaplot  rs16847897 mentaldisorder , scatter(msym(none)) ylabel(0(0.3)1.5) /*title(Means plot for rs16847897genotype#mentaldisorder)*/
*******TLBASE
******************rs16847897
anova tlbase mentaldisorder rs16847897 mentaldisorder##rs16847897

**********************************
anova tlbase mentaldisorder rs2736100 mentaldisorder##rs2736100

anovaplot  rs2736100 mentaldisorder , scatter(msym(none)) ylabel(0(0.3)1.5) /*title(Means plot - rs2736100genotype mentaldisorder interaction)*/ yscale(r(0.45 1.2))

***********************************
anova tlbase mentaldisorder rs16847897 mentaldisorder##rs16847897

anovaplot  rs16847897 mentaldisorder , scatter(msym(none)) ylabel(0(0.3)1.5) /*title(Means plot - rs2736100genotype mentaldisorder interaction)*/ yscale(r(0.45 1.2))


**rs1386494 
anova tlbase mentaldisorder rs1386494 mentaldisorder##rs1386494

**rs1843809 
anova tlbase mentaldisorder rs1843809 mentaldisorder##rs1843809

**rs12696304  
anova tlbase mentaldisorder rs12696304  mentaldisorder##rs12696304 

**rs10936599 
anova tlbase mentaldisorder rs10936599  mentaldisorder##rs10936599 
**rs2853669 
anova tlbase mentaldisorder rs2853669  mentaldisorder##rs2853669 
**rs7726159 
anova tlbase mentaldisorder rs7726159   mentaldisorder##rs7726159  
**rs10069690 
anova tlbase mentaldisorder rs10069690   mentaldisorder##rs10069690  
**rs34517220 
anova tlbase mentaldisorder rs34517220   mentaldisorder##rs34517220  

*******TL12
******************rs16847897
anova tl12 mentaldisorder rs16847897 mentaldisorder##rs16847897

**rs1386494 
anova tl12 mentaldisorder rs1386494 mentaldisorder##rs1386494

**rs1843809 
anova tl12 mentaldisorder rs1843809 mentaldisorder##rs1843809

**rs12696304  
anova tl12 mentaldisorder rs12696304  mentaldisorder##rs12696304 

**rs10936599 
anova tl12 mentaldisorder rs10936599  mentaldisorder##rs10936599 
**rs2853669 
anova tl12 mentaldisorder rs2853669  mentaldisorder##rs2853669 
**rs7726159 
anova tl12 mentaldisorder rs7726159   mentaldisorder##rs7726159  
**rs10069690 
anova tl12 mentaldisorder rs10069690   mentaldisorder##rs10069690  
**rs34517220 
anova tl12 mentaldisorder rs34517220   mentaldisorder##rs34517220  
/*predict yhat

sort rs2736100 mentaldisorder

graph twoway scatter yhat mentaldisorder, connect(L)*/

anova tl12 mentaldisorder stin2vntr_   mentaldisorder##stin2vntr_  

anovaplot  stin2vntr_ mentaldisorder , scatter(msym(none)) ylabel(0.0(0.2)1.2)

anova tlbase mentaldisorder stin2vntr_   mentaldisorder##stin2vntr_  

anovaplot  stin2vntr_ mentaldisorder , scatter(msym(none)) ylabel(0.0(0.2)1.2)


*****Frequencies by mental disorder

tab rs16847897 mentaldisorder, col  
tab rs2736100 mentaldisorder, col  



foreach var in rs2736100 rs16847897 {

oneway tlbase `var', tab
oneway tl12 `var', tab
oneway tldiff `var', tab

}

pwcorr tlbase age


graph twoway (lfitci tlbase age1) (scatter tlbase age1)
graph twoway (lfitci tl12 age1) (scatter tl12 age1)

********rs2736100 and rs7726159

tab1 rs2736100  rs7726159

/*
We have found out that rs2736100 and rs7726159 are in LD. with the following haplotypes: CT, CG and AG. Could you kindly run for me the following;

1) Association between the haplotypes and
i) baseline TL
ii) 12 months TL
iii) TL change,

2) Interaction between IMDs and the haplotypes on;
i), ii) and iii).

Kind regards,

Allan.

*/

tab rs2736100 rs7726159

gen rs2736100_rs7726159=.
replace rs2736100_rs7726159=1 if rs2736100==1 & rs7726159==1
replace rs2736100_rs7726159=2 if rs2736100==1 & rs7726159==3
replace rs2736100_rs7726159=3 if rs2736100==3 & rs7726159==3

tab rs2736100_rs7726159, m

lab def rs2736100_rs7726159 1"AG" 2"CG" 3"CT", replace
lab val rs2736100_rs7726159 rs2736100_rs7726159

***************repeated mesures analysis******************
summa tlbase-tl12

tabstat tlbase-tl12, by(mentaldisorder) stat(n mean sd var)

profileplot tlbase-tl12, by(mentaldisorder)


*****Reshaping the data*******
use "C:\STUDIES\allanK\telomere length\FINALSTATS - reshapedfor repeated measures analysis.dta", clear
gen id=_n
ren tl12 tl2
ren tlbase tl1
keep patid id tl1 tl2 mentaldisorder
order patid id tl1 tl2 mentaldisorder
br
reshape long tl, i(id) j(time)


****Repeated measures anova
anova tl mentaldisorder / id|mentaldisorder time mentaldisorder#time  if id<100, repeated(time)

br
