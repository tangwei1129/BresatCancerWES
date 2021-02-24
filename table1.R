library("tableone")

#set working directory
setwd("~/Desktop/WES")

#Read in Data
dat<- read.csv("WES samples clinic database01292021.csv", header = T, na.strings = "")
## add TNBC ###
dat$TNBC <- ifelse(dat$ER == 0 & dat$PR == 0 & dat$HER2 == 0, 1, ifelse(dat$ER == 1 | dat$PR == 1 | dat$HER2 == 1, 0, NA))
#str(dat)


data <- dat
#Data Cleaning and Factor Naming
#data = data[, -(1:3)]
data[is.na(data)] <- "Unknown"

### numeric variables AGE, weight heightin BMI SCHOOL..years. ###
data$AGE <- as.numeric(data$AGE)
data$weight <- as.numeric(data$weight)
data$heightin <- as.numeric(data$heightin)
data$BMI <- as.numeric(data$BMI)
data$SCHOOL..years. <- as.numeric(data$SCHOOL..years.)
###
data$RACE = factor(data$RACE, levels = c("African American","Caucasian","Asian"))
## ER PR HER2 Stage recode ###
#data$ER <- ifelse(is.na(data$ER),"Unknown", as.character(data$ER))
data$PR <- factor(data$PR, levels = c("0","1","U","unclassified","Unknown"), labels = c("0","1","Unknown","Unknown","Unknown"))
data$HER2 <- factor(data$HER2, levels = c("0","1","U","unclassified","Unknown"), labels = c("0","1","Unknown","Unknown","Unknown"))
#unique(data$Stage)
data$Stage <- factor(data$Stage, levels = c("I", "IIA","IIIA","IIIa", "IIIB","IV","IIB", "IIIC"," IIIA", " IIB"," IIA"," IIIA "," IA"," IIIC", "II","IA",  "TIS","phyllodes","U", "Unknown"), 
                                 labels = c("I", "II",  "III", "III", "III", "IV","II",  "III",  "III",    "II", "II",   "III",   "I", "III",  "II", "I",   "I",  "I", "Unknown","Unknown"))
 
### MONEY  ####
data$MONEY <- factor(data$MONEY,levels = c("1 Under $15,000", "2 $15,000-60,000", "3 More than $60,000", "$10,000 to 29,999", "$30,000 to 59,999", "$60,000 to $90,000",  
                                           "betw. $15,000 and 60,000", "D Doesn’t know”", "greater than $60,000", "greater than $90,000", "less than $10,000", "less than $15,000","Unknown"), 
                                labels = c("< $15k","$15 to 60k", "greater than $60k","$15 to 60k", "$15 to 60k", "greater than $60k",
                                           "$15 to 60k", "Unknown","greater than $60k", "greater than $60k", "< $15k", "< $15k","Unknown"))  
 

#### SCHOOL_recode  ####
data$SCHOOL_recode <- factor(data$SCHOOL_recode, levels = c("5th-6th","10th or 11th","10th-11th","High School","College","College degree",
                                                            "Graduate school","Graduate School","junior high","Junior high", 
                                                            "less than 5th","Some College","Technical School","Unknown"),
                                                  labels = c("<High school", "<High school", "<High school", "High school","College/Technical School", "College/Technical School",
                                                            "Graduate school", "Graduate school",  "<High school", "<High school",
                                                            "<High school", "College/Technical School", "College/Technical School","Unknown"))

## diabetes ##
data$Diabetes..Y.N..at.time.of.tumor.collection <- gsub(" ", "", data$Diabetes..Y.N..at.time.of.tumor.collection)
data$Diabetes..Y.N..at.time.of.tumor.collection <- factor(data$Diabetes..Y.N..at.time.of.tumor.collection, levels = c("no","yes","Yes","Unknown"), labels = c("no", "yes", "yes","Unknown"))

summary(CreateTableOne(data=data))

vars2use <- c("SEX","AGE", "ER","PR","HER2","TNBC","Stage" ,"BMI" , "BMI_recode","MONEY","SCHOOL_recode", "SCHOOL..years.","Diabetes..Y.N..at.time.of.tumor.collection")
catVars <- c("SEX","ER","PR","HER2","TNBC","Stage" ,"BMI_recode","MONEY","SCHOOL_recode", "Diabetes..Y.N..at.time.of.tumor.collection")
table1 <- CreateTableOne(data = data, vars = vars2use, factorVars = catVars, test = TRUE, testNormal = oneway.test(data$RACE == "Caucasian"~ data$RACE =="African American"), strata = "RACE")

## reconstruct table2 for Unknown as NA, and remove Asican popluation to get correct pvalues ###
data2= subset(data, RACE %in% c("African American","Caucasian"))
data2$RACE = droplevels(data2$RACE)
data2[data2=="Unknown"] = NA
table2 <- CreateTableOne(data = data2, vars = vars2use, factorVars = catVars, test = TRUE, strata = "RACE")
print(table2, showAllLevels = TRUE, formatOptions = list(big.mark = ","))

### replace table1 pvalue with table2 pvalue to reflect the AA vs EA comparisons ###
attr(table1$CatTable,"pValues") <- attr(table2$CatTable, "pValues")
attr(table1$ContTable,"pValues") <- attr(table2$ContTable, "pValues")

#### output table1 #######
table1.mat <- print(table1, nonnormal = c("AGE","BMI","SCHOOL..years."), formatOptions = list(big.mark = ","),test = T,quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(table1.mat, file = "Table1.csv")
