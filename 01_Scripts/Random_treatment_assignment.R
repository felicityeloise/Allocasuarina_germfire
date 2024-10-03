# Written by Felicity Charles
# 03/04/2024

# This script allows the random assignment of seeds lots to one of 8 treatments that will be used in the germination experiment. The seeds have been split by species but seed lots are randomly assigned to a treatment based on the condition that there are more than 60 seeds in a seed lot. To begin random assignment of the seeds, we want to create two new dataframes for torulosa, as littoralis has too few seeds available for high fire, that contain seeds for high fire and seeds for low fire. We will do the same for littoralis but only create a dataframe for low fire. This will help to ensure that seeds from high fire and low fire should be equally assigned to a treatment. 




# Read in the data
torulosa <- read.csv("./00_Data/Full_experiment/Treatment_grouping/Torulosa.csv", header = T)
head(torulosa)
str(torulosa)
unique(torulosa$Fire_freq_category)

littoralis <- read.csv("./00_Data/Full_experiment/Treatment_grouping/Littoralis.csv", header = T)
head(littoralis)

torulosa_lowfi <- subset(torulosa, Fire_freq_category == "low")
head(torulosa_lowfi)

torulosa_highfi <- subset(torulosa, Fire_freq_category == "High")
head(torulosa_highfi)




# Lets also simplify things by removing any seed lots that have less than 60 seeds. This way we can just assign a random number to the column.


torulosa_lowfi <- subset(torulosa_lowfi, Number_seeds >=60)
head(torulosa_lowfi); dim(torulosa_lowfi)

torulosa_highfi <- subset(torulosa_highfi, Number_seeds >= 60)
head(torulosa_highfi); dim(torulosa_highfi)

littoralis <- subset(littoralis, Number_seeds >=60)
head(littoralis); dim(littoralis)



# Now we can begin randomly assigning variables

# We know that we can assign some of the seed lots to all treatments already 
# Torulosa = WPS55, WPS125, WPS126, WPS134, WPS177, WPS73, WPS113, WPS130, WPS147, WPS149, WPS197
#Aim is to have 30 seed lots in each treatment, above gives us a start with 11 seeds lots in each treatment. Only need 19 more across 9 treatments for torulosa

# The treatments are 1. Control, 2. Temperature one, 3. Temperature two, 4. Smoke, 5. Temperature one + smoke, 6. Temperature two + smoke
library(dplyr)


# Torulosa low fire ----

x <- arrangements::permutations(6, 6, nsample = nrow(torulosa_lowfi))
d <- tibble(ID = torulosa_lowfi$Seed_lot, number_treats_possible = torulosa_lowfi$Number_treatments_possible, Group1 = x[,1], Group2 = x[,2], Group3 = x[, 3], Group4 = x[, 4], Group5 = x[,5], Group6 = x[,6])
d
da <- d %>% 
  mutate(is_equal = if_all(Group2:Group6, `==`, Group1))

# This method worked to ensure none of the rows have the same group in any of the columns - our job is then to manually look at this data, make sure we have enough samples across the treatments. 
write.csv(da, "./03_Results/Treatment_grouping/torulosa_lowfi.csv")

da <- read.csv("./03_Results/Treatment_grouping/torulosa_lowfi.csv")
dim(da)

# See how many individuals have been allocated to each treatment
Group1 <- da %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- da %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- da %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- da %>% count(Group4, sort = F, name = "Treatment frequency")
Group5 <- da %>% count(Group5, sort = F, name = "Treatment frequency")
Group6 <- da %>% count(Group6, sort = F, name = "Treatment frequency")

T1 <- Group1[1,2]+Group2[1,2]+Group3[1,2]+Group4[1,2]+Group5[1,2]+Group6[1,2]+Group7[1,2]+Group8[1,2]
T2 <- Group1[2,2]+Group2[2,2]+Group3[2,2]+Group4[2,2]+Group5[2,2]+Group6[2,2]+Group7[2,2]+Group8[2,2]
T3 <- Group1[3,2]+Group2[3,2]+Group3[3,2]+Group4[3,2]+Group5[3,2]+Group6[3,2]+Group7[3,2]+Group8[3,2]
T4 <- Group1[4,2]+Group2[4,2]+Group3[4,2]+Group4[4,2]+Group5[4,2]+Group6[4,2]+Group7[4,2]+Group8[4,2]
T5 <- Group1[5,2]+Group2[5,2]+Group3[5,2]+Group4[5,2]+Group5[5,2]+Group6[5,2]+Group7[5,2]+Group8[5,2]
T6 <- Group1[6,2]+Group2[6,2]+Group3[6,2]+Group4[6,2]+Group5[6,2]+Group6[6,2]+Group7[6,2]+Group8[6,2]

T1;T2;T3;T4;T5;T6

# Now reduce the number of treatments for each individual based on the number of allowed treatments - in this round, always keep the number assigned to the first group - group 1
dim(da)
View(da)
da[1, 7:10] <- 0
da[2, 6:10] <- 0
da[3, 6:10] <- 0
da[4, 6:10] <- 0
da[5, 8:10] <- 0
da[6, 9:10] <- 0
da[7, 6:10] <- 0
da[8, 5:10] <- 0
da[9, 7:10] <- 0
da[10, 9:10] <- 0
da[12, 4:9] <- 0
da[13, 4:9] <- 0
da[14, 5:10] <- 0
da[15, 7:10] <- 0
da[17, 6:10] <- 0
da[18, 5:10] <- 0
da[20, 6:10] <- 0
da[22, 5:10] <- 0
da[26, 6:10] <- 0
da[27, 5:10] <- 0
da[29, 7:10] <- 0
da[31, c(5, 8:10)] <- 0 # Need to do this one a little different as we would like to keep particular columns
da[32, 6:10] <- 0
da[33, 9:10] <- 0
da[34, 6:10] <- 0
da[35, 6:10] <- 0
da[36, 6:10] <- 0
da[37, 6:10] <- 0
da[38, 5:10] <- 0
da[39, 6:10] <- 0
da[40, 8:10] <- 0
da[42, 7:10] <- 0
da[43, 6:10] <- 0



# Check how many individuals have been allocated to each treatment
Group1 <- da %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- da %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- da %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- da %>% count(Group4, sort = F, name = "Treatment frequency")
Group5 <- da %>% count(Group5, sort = F, name = "Treatment frequency")
Group6 <- da %>% count(Group6, sort = F, name = "Treatment frequency")




T1 <- Group1[2,2]+Group2[2,2]+Group3[2,2]+Group4[2,2]+Group5[2,2]
T2 <- Group1[3,2]+Group2[3,2]+Group3[3,2]+Group4[3,2]+Group5[3,2]
T3 <- Group1[4,2]+Group2[4,2]+Group3[4,2]+Group4[4,2]+Group5[4,2]+Group6[2,2]
T4 <- Group1[5,2]+Group2[5,2]+Group3[5,2]+Group4[5,2]+Group5[5,2]+Group6[3,2]
T5 <- Group1[6,2]+Group2[6,2]+Group3[6,2]+Group4[6,2]+Group5[6,2]+Group6[4,2]
T6 <- Group1[7,2]+Group2[7,2]+Group3[7,2]+Group4[7,2]+Group5[7,2]+Group6[5,2]


T1;T2;T3;T4;T5;T6

View(da) # Look at da and determine which treatment types can be assigned 0 for each column
# Start at the end column and work back towards column 1

unique(da$Group6)
da$Group6[da$Group6 == 3 | da$Group6 == 4 | da$Group6 == 6] <- 0 # Do all but not treatment 5
View(da)

unique(da$Group5)
da$Group5[da$Group5 == 1 | da$Group5 == 2 | da$Group5 == 3 | da$Group5 == 4 | da$Group5 == 6] <- 0 # Do all but treatment 5

View(da)
unique(da$Group4)
da$Group4[da$Group4 == 5 | da$Group4 == 1 | da$Group4 == 2| da$Group4 == 3 | da$Group4 == 6] <- 0 # All but treatment 4

View(da)

da$Group2[da$Group2 == 6] <- 0

# Now do manual removal of numbers for the other columns - choosing to remove those already in many treatments first
da[24, 6] <- 0
da[41, 4] <- 0
da[15, 5] <- 0
da[11, 6] <- 0
da[29, 5] <- 0
da[23, 4] <- 0
da[21, 5] <- 0
da[28, 6] <- 0
da[10, 5] <- 0
da[30, 6] <- 0


# Check how many individuals have been allocated to each treatment
Group1 <- da %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- da %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- da %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- da %>% count(Group4, sort = F, name = "Treatment frequency")
Group5 <- da %>% count(Group5, sort = F, name = "Treatment frequency")
Group6 <- da %>% count(Group6, sort = F, name = "Treatment frequency")


T1 <- Group1[2,2]+Group2[2,2]+Group3[2,2]
T2 <- Group1[3,2]+Group2[3,2]+Group3[3,2]
T3 <- Group1[4,2]+Group2[4,2]+Group3[4,2]
T4 <- Group1[5,2]+Group2[5,2]+Group3[5,2]+Group4[2,2]
T5 <- Group1[6,2]+Group2[6,2]+Group3[6,2]+Group5[2,2]+Group6[2,2]
T6 <- Group1[7,2]+Group3[7,2]


T1;T2;T3;T4;T5;T6


# Write the final output
write.csv(da, './03_Results/Treatment_grouping/Treatment_grouping/torulosa_lowfi_treatments_final.csv')

torlow_final <- read.csv('./03_Results/Treatment_grouping/torulosa_lowfi_treatments_final.csv', header = T)


Control <- torlow_final[torlow_final$Group1 == 1 | torlow_final$Group2 == 1 | torlow_final$Group3 == 1,]
dim(Control)
Control

H1 <- torlow_final[torlow_final$Group1 == 2 | torlow_final$Group2 == 2 | torlow_final$Group3 == 2,]
dim(H1)
H1

H2 <- torlow_final[torlow_final$Group1 == 3 | torlow_final$Group2 == 3 | torlow_final$Group3 == 3,]
dim(H2) 
H2

Smoke <- torlow_final[torlow_final$Group1 == 4 | torlow_final$Group2 == 4 | torlow_final$Group3 == 4 | torlow_final$Group4 == 4, ]
dim(Smoke)
Smoke

H1_smoke <- torlow_final[torlow_final$Group1 == 5 | torlow_final$Group2 == 5 | torlow_final$Group3 == 5 | torlow_final$Group5 == 5 | torlow_final$Group6 == 5, ]
dim(H1_smoke)
H1_smoke


H2_smoke <- torlow_final[torlow_final$Group1 == 6 | torlow_final$Group3 == 6, ]
dim(H2_smoke)
H2_smoke


# Torulosa high fire ----

x1 <- arrangements::permutations(6,6, nsample = nrow(torulosa_highfi))
d1 <- tibble(ID = torulosa_highfi$Seed_lot, number_treats_possible = torulosa_highfi$Number_treatments_possible, Group1 = x1[,1], Group2 = x1[,2], Group3 = x1[, 3], Group4 = x1[, 4], Group5 = x1[,5], Group6 = x1[,6])
d1
da1 <- d1 %>% 
  mutate(is_equal = if_all(Group2:Group6, `==`, Group1))

write.csv(da1, "./03_Results/Treatment_grouping/torulosa_highfi.csv")



# Filter out those treatments that exceed the number allowed for each individual
df <- read.csv("./03_Results/Treatment_grouping/torulosa_highfi.csv")


# See how many individuals have been allocated to each treatment
Group1 <- df %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- df %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- df %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- df %>% count(Group4, sort = F, name = "Treatment frequency")
Group5 <- df %>% count(Group5, sort = F, name = "Treatment frequency")
Group6 <- df %>% count(Group6, sort = F, name = "Treatment frequency")


T1 <- Group1[1,2]+Group2[1,2]+Group3[1,2]+Group4[1,2]+Group5[1,2]+Group6[1,2]
T2 <- Group1[2,2]+Group2[2,2]+Group3[2,2]+Group4[2,2]+Group5[2,2]+Group6[2,2]
T3 <- Group1[3,2]+Group2[3,2]+Group3[3,2]+Group4[3,2]+Group5[3,2]+Group6[3,2]
T4 <- Group1[4,2]+Group2[4,2]+Group3[4,2]+Group4[4,2]+Group5[4,2]+Group6[4,2]
T5 <- Group1[5,2]+Group2[5,2]+Group3[5,2]+Group4[5,2]+Group5[5,2]+Group6[5,2]
T6 <- Group1[6,2]+Group2[6,2]+Group3[6,2]+Group4[6,2]+Group5[6,2]+Group6[6,2]


T1;T2;T3;T4;T5;T6
# Double check which ones had treatment 4, need to remove 2 less from this treatment

View(df)
df[1, 7:10] <- 0
df[2, c(5, 7:10)] <- 0
df[3,5:10] <- 0
df[5, 6:10] <- 0
df[6, 6:10] <- 0
df[8, 6:10] <- 0
df[9, 5:10] <- 0
df[10, 5:10] <- 0
df[11, 8:10] <- 0
df[12, 7:10] <- 0
df[13, 6:10] <- 0
df[14, 9:10] <- 0
df[17, 5:10] <- 0
df[19, 8:10] <- 0
df[22, 6:10] <- 0
df[23, 5:10] <- 0
df[24, 8:10] <- 0
df[27, 7:10] <- 0
df[28, c(5, 7:10)] <- 0
df[29, c(4, 7:10)] <- 0
df[30, 7:10] <- 0
 


# See how many individuals have been allocated to each treatment. Checking that we have 15 in each treatment
Group1 <- df %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- df %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- df %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- df %>% count(Group4, sort = F, name = "Treatment frequency")
Group5 <- df %>% count(Group5, sort = F, name = "Treatment frequency")
Group6 <- df %>% count(Group6, sort = F, name = "Treatment frequency")




T1 <- Group1[2,2]+Group2[2,2]+Group3[2,2]+Group4[2,2]+Group5[2,2]+Group6[2,2]
T2 <- Group1[3,2]+Group2[3,2]+Group3[3,2]+Group4[3,2]+Group6[3,2]
T3 <- Group1[4,2]+Group2[4,2]+Group3[4,2]+Group4[4,2]+Group5[3,2]
T4 <- Group1[5,2]+Group2[5,2]+Group3[5,2]+Group4[5,2]+Group5[4,2]+Group6[4,2]
T5 <- Group1[6,2]+Group2[6,2]+Group3[6,2]+Group4[6,2]+Group5[5,2]+Group6[5,2]
T6 <- Group1[7,2]+Group2[7,2]+Group3[7,2]+Group4[7,2]+Group5[6,2]+Group6[6,2]


T1;T2;T3;T4;T5;T6

View(df)

# We do not need to reduce treatments 2 or 3. And only treatment 6 by 1

unique(df$Group6)
df$Group6[df$Group6 == 1 | df$Group6 == 4 | df$Group6 == 5 | df$Group6 == 6] <- 0

# Now further reduce for each treatment
df$Group5[df$Group5 == 1] <- 0
df$Group4[df$Group4 == 1 | df$Group4 == 4] <- 0
df$Group3[df$Group3 == 5] <- 0


# See how many individuals have been allocated to each treatment. Checking that we have 15 in each treatment
Group1 <- df %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- df %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- df %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- df %>% count(Group4, sort = F, name = "Treatment frequency")
Group5 <- df %>% count(Group5, sort = F, name = "Treatment frequency")
Group6 <- df %>% count(Group6, sort = F, name = "Treatment frequency")




T1 <- Group1[2,2]+Group2[2,2]+Group3[2,2]
T2 <- Group1[3,2]+Group2[3,2]+Group3[3,2]+Group4[2,2]+Group6[2,2]
T3 <- Group1[4,2]+Group2[4,2]+Group3[4,2]+Group4[3,2]+Group5[2,2]
T4 <- Group1[5,2]+Group2[5,2]+Group3[5,2]+Group5[3,2]
T5 <- Group1[6,2]+Group2[6,2]+Group4[4,2]+Group5[4,2]
T6 <- Group1[7,2]+Group2[7,2]+Group3[6,2]+Group4[5,2]+Group5[5,2]


T1;T2;T3;T4;T5;T6


write.csv(df, './03_Results/Treatment_grouping/Treatment_grouping/torulosa_highfi_final_treatments.csv')



torhi_final <- read.csv('./03_Results/Treatment_grouping/torulosa_highfi_final_treatments.csv', header = T)

Control <- torhi_final[torhi_final$Group1 == 1 | torhi_final$Group2 == 1  | torhi_final$Group3 == 1, ]
dim(Control)
Control

H1 <- torhi_final[torhi_final$Group1 == 2 | torhi_final$Group2 == 2 | torhi_final$Group3 == 2 | torhi_final$Group4 == 2 | torhi_final$Group6 == 2,]
dim(H1)
H1

H2 <- torhi_final[torhi_final$Group1 == 3 | torhi_final$Group2 == 3 | torhi_final$Group3 == 3 | torhi_final$Group4 == 3 | torhi_final$Group5 == 3,]
dim(H2)
H2

Smoke <- torhi_final[torhi_final$Group1 == 4 |torhi_final$Group2 == 4 |torhi_final$Group3 == 4 | torhi_final$Group5 == 4,]
dim(Smoke)
Smoke

H1_smoke <- torhi_final[torhi_final$Group1 == 5| torhi_final$Group2 == 5 | torhi_final$Group4 == 5 | torhi_final$Group5 == 5,]
dim(H1_smoke)
H1_smoke

H2_smoke <- torhi_final[torhi_final$Group1 == 6 | torhi_final$Group2 == 6 | torhi_final$Group3 == 6 | torhi_final$Group4 == 6 | torhi_final$Group5 == 6,]
dim(H2_smoke)
H2_smoke





# Littoralis ----

x3 <- arrangements::permutations(6,6, nsample = nrow(littoralis))
d3 <- tibble(ID = littoralis$Seed_lot, number_treats_possible = littoralis$Number_treatments_possible, Group1 = x3[,1], Group2 = x3[,2], Group3 = x3[, 3], Group4 = x3[, 4], Group5 = x3[,5], Group6 = x3[,6])
d3
d3a <- d3 %>% 
  mutate(is_equal = if_all(Group2:Group6, `==`, Group1))

write.csv(d3a, "./03_Results/Treatment_grouping/littoralis.csv")



# Filter out those treatments that exceed the number allowed for each individual
d3 <- read.csv("./03_Results/Treatment_grouping/littoralis.csv")

# See how many individuals have been allocated to each treatment
Group1 <- d3 %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- d3 %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- d3 %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- d3 %>% count(Group4, sort = F, name = "Treatment frequency")
Group5 <- d3 %>% count(Group5, sort = F, name = "Treatment frequency")
Group6 <- d3 %>% count(Group6, sort = F, name = "Treatment frequency")


T1 <- Group1[1,2]+Group2[1,2]+Group3[1,2]+Group4[1,2]+Group5[1,2]+Group6[1,2]
T2 <- Group1[2,2]+Group2[2,2]+Group3[2,2]+Group4[2,2]+Group5[2,2]+Group6[2,2]
T3 <- Group1[3,2]+Group2[3,2]+Group3[3,2]+Group4[3,2]+Group5[3,2]+Group6[3,2]
T4 <- Group1[4,2]+Group2[4,2]+Group3[4,2]+Group4[4,2]+Group5[4,2]+Group6[4,2]
T5 <- Group1[5,2]+Group2[5,2]+Group3[5,2]+Group4[5,2]+Group5[5,2]+Group6[5,2]
T6 <- Group1[6,2]+Group2[6,2]+Group3[6,2]+Group4[6,2]+Group5[6,2]+Group6[6,2]


T1
T2
T3
T4
T5
T6


  
# We can see we have plenty of individuals/treatment so let's start removing treatments by row. Then we will recaculate the number and remove until we only have 15 individuals for each treatment



View(d3) # Look at d3 and determine which samples need to have 0 allocated to particular columns for each row

d3[1,5:10] <- 0
d3[2,9:10] <- 0
d3[4,7:10] <- 0
d3[7,8:10] <- 0
d3[8, 10] <- 0
d3[10, 8:10] <- 0
d3[12, 8:10] <- 0
d3[15, 10:10] <- 0
d3[16, 7:10] <- 0
d3[18, 6:10] <- 0
d3[30, 5:10] <- 0
d3[32, 10:10] <- 0
d3[34, 7:10] <- 0
d3[35, 8:10] <- 0
d3[37, 10:10] <- 0

# NOTE GROUP 1 hasn't been touched. 
Group1 <- d3 %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- d3 %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- d3 %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- d3 %>% count(Group4, sort = F, name = "Treatment frequency")
Group5 <- d3 %>% count(Group5, sort = F, name = "Treatment frequency")
Group6 <- d3 %>% count(Group6, sort = F, name = "Treatment frequency")

T1 <- Group1[1,2]+Group2[2,2]+Group3[2,2]+Group4[2,2]+Group5[2,2]+Group6[2,2]
T2 <- Group1[2,2]+Group2[3,2]+Group3[3,2]+Group4[3,2]+Group5[3,2]+Group6[3,2]
T3 <- Group1[3,2]+Group2[4,2]+Group3[4,2]+Group4[4,2]+Group5[4,2]+Group6[4,2]
T4 <- Group1[4,2]+Group2[5,2]+Group3[5,2]+Group4[5,2]+Group5[5,2]+Group6[5,2]
T5 <- Group1[5,2]+Group2[6,2]+Group3[6,2]+Group4[6,2]+Group5[6,2]+Group6[6,2]
T6 <- Group1[6,2]+Group2[7,2]+Group3[7,2]+Group4[7,2]+Group5[7,2]+Group6[7,2]

T1;T2;T3;T4;T5;T6


# Can still be cut down further
View(d3)

# Now see which treatment types can be removed from whole columns
d3.1 <- d3[,c(1:7)] # Remove group 5 and 6 completely


unique(d3.1$Group4)
d3.1$Group4[d3.1$Group4 == 1 | d3.1$Group4 == 2 | d3.1$Group4 == 3 | d3.1$Group4 == 5 | d3.1$Group4 == 6] <- 0

d3.1$Group2[d3.1$Group2 == 3] <- 0
# Also adjust for those that needed to be partially deleted

d3.1[19, 6] <- 0
d3.1[37, 6] <- 0
d3.1[17, 4] <- 0
d3.1[22, 6] <- 0
d3.1[21, 5] <- 0
d3.1[27, 4] <- 0
d3.1[29, 5] <- 0
d3.1[23, 6] <- 0
d3.1[2, 4] <- 0
d3.1[23, 4] <- 0
d3.1[6, 6] <- 0
d3.1[20, 6] <- 0
d3.1[4, 4] <- 0
d3.1[10, 5] <- 0
d3.1[32, 6] <- 0
d3.1[18, 4] <- 0
d3.1[5, 5] <- 0
d3.1[35, 6] <- 0


Group1 <- d3.1 %>% count(Group1, sort = F, name = "Treatment frequency")
Group2 <- d3.1 %>% count(Group2, sort = F, name = "Treatment frequency")
Group3 <- d3.1 %>% count(Group3, sort = F, name = "Treatment frequency")
Group4 <- d3.1 %>% count(Group4, sort = F, name = "Treatment frequency")



# Keep in mind that we need check which rows correspond to the treatment because there have been zeros added for all but group 1, and some groups do not have all treatments. We also need to remove any groupings that do not contain the treatment number

T1 <- Group1[2,2]+Group2[2,2]+Group3[2,2]
T2 <- Group1[3,2]+Group2[3,2]+Group3[3,2]
T3 <- Group1[4,2]+Group3[4,2]
T4 <- Group1[5,2]+Group2[4,2]+Group3[5,2]+Group4[2,2]
T5 <- Group1[6,2]+Group2[5,2]+Group3[6,2]
T6 <- Group1[7,2]+Group2[6,2]+Group3[7,2]


T1;T2;T3;T4;T5;T6
View(d3.1)
# These are the final treatment allocations for Allocasuarina littoralis individuals

write.csv(d3.1, './03_Results/Treatment_grouping/littoralis_treatments_final.csv')

lit_final <- read.csv('./03_Results/Treatment_grouping/littoralis_treatments_final.csv', header = T)
head(lit_final)


Control <- lit_final[lit_final$Group1 == 1 | lit_final$Group2 == 1 | lit_final$Group3 == 1,]
dim(Control)
Control

Heat1<- lit_final[lit_final$Group1 == 2 | lit_final$Group2 == 2 | lit_final$Group3 == 2,]
dim(Heat1)
Heat1

Heat2 <- lit_final[lit_final$Group1 == 3 | lit_final$Group3 == 3,]
dim(Heat2)
Heat2

Smoke <- lit_final[lit_final$Group1 == 4 | lit_final$Group2 == 4 | lit_final$Group3 == 4 | lit_final$Group4 == 4,]
dim(Smoke)
Smoke

H1_smoke <- lit_final[lit_final$Group1 == 5 | lit_final$Group2 == 5 | lit_final$Group3 == 5,]
dim(H1_smoke)
H1_smoke

H2_smoke <- lit_final[lit_final$Group1 == 6 | lit_final$Group2 == 6 | lit_final$Group3 == 6,]
dim(H2_smoke)
H2_smoke
