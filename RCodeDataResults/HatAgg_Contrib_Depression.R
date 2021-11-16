install.packages("netmeta") #Version 2.0-0 or later
library(netmeta)
library("writexl")

#read in depression data (Rucker and Schwarzer 2014)
dat <- read.csv2("~...\\Rucker2014.txt", header = TRUE, sep = "", quote="\"", dec=".")

#fit NMA model
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, sm = "OR", data = dat)

#aggregate hat matrix:----------------------------------------------------------
hagg <- hatmatrix(net1, method = "Davies", type = "short")

hagg$fixed #aggregate hat matrix for the fixed effect model
hagg$fixed[1, ] #1st row of hagg$fixed (evidence flow for comparison 1-3)
#Save to Excel file
hagg.df <- as.data.frame(hagg$fixed)
write_xlsx(hagg.df,"~...\\Hat_agg_netmeta.xlsx")


#Shortest path algorithm contributions:-----------------------------------------
short <- netcontrib(net1, method = "shortestpath", hatmatrix.F1000 = FALSE)

short$fixed #fixed effect model contribution matrix
#Save to Excel file
short.df <- as.data.frame(short$fixed)
write_xlsx(short.df,"~...\\Contrib_shortestpath_netmeta.xlsx")


#random walk algorithm contributions:-------------------------------------------
rw <- netcontrib(net1, method = "randomwalk", hatmatrix.F1000 = FALSE)

rw$fixed #fixed effect model contribution matrix
#Save to Excel file
rw.df <- as.data.frame(rw$fixed)
write_xlsx(rw.df,"~...\\Contrib_randomwalk_netmeta.xlsx")




