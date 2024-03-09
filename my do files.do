* Convert viralload1 to a string variable
gen str20 viralload1_str = viralload1

* Replace "<20" with missing
replace viralload1_str = "<20" if viralload1_str == "< 20"
replace viralload1_str = "." if viralload1_str == "<20"
replace viralload1_str = "10" if viralload1_str == "<20"

* Remove commas from viralload1_str
gen viralload1_str_no_commas = subinstr(viralload1_str, ",", "", .)

* Convert viralload1_str_no_commas to numeric
gen viralload1_numeric = real(viralload1_str_no_commas) if viralload1_str_no_commas != "."

*Categorise Viral laods categories
recode Viralload (0=1 "Undetectable") (1/1000=2 "Suppressed") (1001/max=3 "Unsuppressed"), generate(Viralload_Category)
(611 differences between Viralload and Viralload_Category)

*drop old viral load
drop viralload1

*FOR ORPHANHOOD
* Change 0.5 to another value (e.g.2) as STATA can't encode, rename 0.5.
replace orphanhood = 2 if orphanhood == 0.5

* Define labels for the modified numeric values
label define orphanhood_labels 0 "Both parents alive" 2 "Single parent alive" 1 "Both parents dead"

* Apply the labels to the numeric variable
label values orphanhood orphanhood_labels

*CREATING CHILD VS ADOLESCENT CATEGORIES FROM AGE GROUPS
* Since the data with ga ecategory is named "agecatak"
gen AgeCategory = ""

* Label categories 13-17 and 9-12 as "Adolescent"
replace AgeCategory = "Adolescent" if agecatak == "13-17" | agecatak == "9-12"

* Label category 5-8 as "Children"
replace AgeCategory = "Children" if agecatak == "5-8"

*UPDATE TOBACCO CATEGORIES
* Generate a new variable named tobacco_status with empty strings
gen tobacco_status = ""

* Recode values based on conditions
replace tobacco_status = "I don't know" if tobacco1 == 9
replace tobacco_status = "Yes" if tobacco1 == 1 | tobacco1 == 2 | tobacco1 == 3
label variable tobacco_status "Smokes tobacco cigarretes"
drop tobacco1

*BMI
gen height_m = heightst1/100 //change height to meters
gen BMI = weightst1/(height_m ^2) //Generate BMI inn kg/m2
* Round BMI to one decimal place to fit categories
gen rounded_bmi = round(BMI, 0.1)
* Recode rounded BMI into categories
recode rounded_bmi (min/18.4=1 "Underweight") (18.5/24.9=2 "Normal weight") (25/29.9=3 "Overweight") (30/max=4 "Obesity"), generate(BMI_category)

*Recode CD4 levels
recode cd4takeoff1 (min/350=1 "Normal CD4 range") (351/500=2 "Possibly immune compromised") (501/max=3 "Advanced immune compromised"), generate(CD4_category)




