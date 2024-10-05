* Create a new variable to store the decoded values
gen New_Age_Category = ""
* Replace encoded values with their corresponding original categories
replace New_Age_Category = "Thirteen to Seventeen Years" if Age_Category == 1
replace New_Age_Category = "Five to Eight" if Age_Category == 2
replace New_Age_Category = "None to Twelve" if Age_Category == 3
replace New_Age_Category = "Thirteen to Seventeen Years" if Age_Category == "13-17"
replace New_Age_Category = "Five to Eight" if Age_Category == "5-8"
replace New_Age_Category = "Nine to Twelve" if Age_Category == "9-12"

order New_Age_Category, after (Age_Category)
drop Age_Category
rename New_Age_Category Age_Category

encode Age_Category, gen (Age_G)
order Age_G, after (Age_Category)
drop Age_Category
rename Age_G Age_Category

rename age Age          
rename religion1 Religion     
rename ses_cat Socio_Economic_Status   
rename chilborhiv Born_with_HIV
rename rs35531 rs35531_polymorphism      
rename tlbase Baseline_telomere_length         
rename sex1 Gender_of_Child         
rename childeduc1 Childs_Education_Level  
rename childartk1 Child_Takes_ART    
rename childworst1 Worst_HIV_Stage 
rename tl12 Telomere_length_12Months          
rename httlpr1 httlpr1_polymorphism  
rename childpremt1 Premature_Birth 
rename rs1843809 rs1843809_polymorphism
rename Rs10482605  Rs10482605_polymorphism
rename rs34517220 rs34517220_polymorphism
rename Rs1360780  Rs1360780_polymorphism  
rename rs1386494 rs1386494_polymorphism
rename rs4570625 rs4570625_polymorphism
