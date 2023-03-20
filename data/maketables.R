#################  WARNING  ##################
#
# Does NOT work in R V-2.5 or later
#   (it may if we re-read and save the md0* files again in R V-2.5) 
#        No... has trouble with new "merge" command in get.naics
#
###############################################################

#---------------------- reading in the data from tables ------------
# Don't need to do again
# Data for 2000 and before do not appear to show suppressed cells
md01<-read.csv(file="cewenb01.csv", header=T)
md02<-read.csv(file="cewenb02.csv", header=T)
md03<-read.csv(file="cewenb03.csv", header=T)
md04<-read.csv(file="cewenb04.csv", header=T)
md05<-read.csv(file="cewenb05.csv", header=T)
md06<-read.csv(file="cewenb06.csv", header=T)
#----------------------------------------------------------------


#------------------internal for:  get.naics --------------------
get.1y.naics<-function(data, naicscode){

  subset(data, data$naics ==naicscode, 
         select =c( "area", "statq1", "wageq1", "statq2", "wageq2", "statq3",
                    "wageq3", "statq4", "wageq4", "annstat", "annwage"))
}
#---------------------------------------------------------------------

#------------------------- get.naics ------------------------
#
##  returns table of all years for the NAICS code listed in the 
##  given vector
##
 
get.naics<-function(naicscode){
  wtab<-get.1y.naics(md01, naicscode)
  wtab<- merge(wtab, get.1y.naics(md02, naicscode), by="area", all.x=T, all.y=T)
  wtab<- merge(wtab, get.1y.naics(md03, naicscode), by="area", all.x=T, all.y=T)
  wtab<- merge(wtab, get.1y.naics(md04, naicscode), by="area", all.x=T, all.y=T)
  wtab<- merge(wtab, get.1y.naics(md05, naicscode), by="area", all.x=T, all.y=T)
  wtab<- merge(wtab, get.1y.naics(md06, naicscode), by="area", all.x=T, all.y=T)
  wtab<-cbind(naics=rep(naicscode, dim(wtab)[1]),wtab)

  names(wtab)<-c("naics", "area",   "s11", "wage01-1", "s12", "wage01-2", "s13", "wage01-3", 
                 "s14", "wage01-4", "s1a", "wage01-a", "s21", "wage02-1", "s22", "wage02-2", 
                 "s23", "wage02-3", "s24", "wage02-4", "s2a", "wage02-a", "s31", "wage03-1", 
                 "s32", "wage03-2", "s33", "wage03-3", "s34", "wage03-4", "s3a", "wage03-a", 
                 "s41", "wage04-1", "s42", "wage04-2", "s43", "wage04-3", "s44", "wage04-4", 
                 "s4a", "wage04-a", "s51", "wage05-1", "s52", "wage05-2", "s53", "wage05-3", 
                 "s54", "wage05-4", "s5a", "wage05-a", "s61", "wage06-1", "s62", "wage06-2", 
                 "s63", "wage06-3", "s64", "wage06-4", "s6a", "wage06-a")

return(wtab)
}
#------------------------- end get.naics ------------------------


#--------------------------make.wtable -----------------------
##
##

make.wtable<-function(v_naics){
 wtab<-NULL
 for(i in v_naics) if(is.null(wtab)) wtab<-get.naics(i) else wtab<-rbind(wtab,get.naics(i)) 

 return(wtab)
}
#--------------------------make.wtable -----------------------




#--------------------------get.wages -----------------------
##
##  gives naics table with just wages.. removes the indicators
##  for whether the data is suppressed or not.
get.wages<-function(wtab){

 wtab[c("naics", "area", "wage01-1", "wage01-2", "wage01-3", "wage01-4", "wage01-a",
        "wage02-1", "wage02-2", "wage02-3", "wage02-4", "wage02-a", "wage03-1", 
        "wage03-2", "wage03-3", "wage03-4", "wage03-a", "wage04-1",  "wage04-2", 
        "wage04-3", "wage04-4", "wage04-a", "wage05-1", "wage05-2", "wage05-3", 
        "wage05-4", "wage05-a", "wage06-1", "wage06-2", "wage06-3", 
        "wage06-4", "wage06-a")]
}
#--------------------------get.wages -----------------------



#------------------------- get.area------------------------------
#
#  returns the wage table for only the given areas
#

get.area<-function(tab, area_v){
table<-NULL
for(i in area_v) table<-rbind(table, subset(tab, area==i))

return(table)

}
#-------------------------- end get.area -----------------------

#------------------------- comp series --------------------------
# gives the percent missing for each area and each naics given in the vector
# input vector of naics
#
get.perc.miss<-function(naics){
  pmis<-NULL

  tab<-get.naics(naics)
  n<-(ncol(tab)-2)/2
  for(i in 1:nrow(tab))
     pmis<-rbind(pmis, c(tab$area[i], sum(tab[i,]=="N")/n))

  return(pmis)
}

comp.ser<-function(nvec){
  tab<-get.perc.miss(nvec[1]) 
  if(length(nvec)==1)  return(tab)
    else for(i in nvec[-1]) 
           tab<-merge(tab, get.perc.miss(i), by=1) 
  colnames(tab)<-c("area", paste(nvec))
  return(tab)  
}
#------------------------- end comp.ser ------------------------



#--------------------- examples ---------------------------------
test<-make.wtable(c(451, 4511))

#res.constr<-make.wtable(c(23611, 236115, 236116, 236117, 236118))

#nonres.constr<-make.wtable(c(2362, 23621, 23622))

#epower<-make.wtable(c(22111, 221111, 221112, 221113, 221119))

#--------------retail clothes ---------------
clothes<-make.wtable(c(448, 4481, 4482, 4483))


comp.ser(c(448, 4481, 4482, 4483))

series1<-get.wages(get.area(clothes, 24037))
series2<-get.wages(get.area(clothes, 24015))
series3<-get.wages(get.area(clothes, 24013))
series4<-get.wages(get.area(clothes, 24999))

save(series1, file="series1")
save(series2, file="series2")
save(series3, file="series3")
save(series4, file="series4")
#-------------------------------------------------------


#--------------Electronic and Appliance retail ---------------
electron<-make.wtable(c(4431, 44311, 44312, 44313))

#Appliance, Computer, Camera
comp.ser(c(4431, 44311, 44312, 44313))

series11<-get.wages(get.area(electron, 24005))
series12<-get.wages(get.area(electron, 24031))
series13<-get.wages(get.area(electron, 24510))

#save(series11, file="series11")
#save(series12, file="series12")
#save(series13, file="series13")
#------------------------------------------------------


#-------Residential Building Construction---------

#23611	Residential Building Co
#236115	New Single-Family Housi
#236116	New Multifamily Housing
#236117	New Housing Operative B
#236118	Residential Remodelers

res.con<-make.wtable(c(23611, 236115, 236116, 236117, 236118))

# Residential Building Construction
comp.ser(c(23611, 236115, 236116, 236117, 236118))


series21<-get.wages(get.area(res.con, 24005))
series22<-get.wages(get.area(res.con, 24025))
series23<-get.wages(get.area(res.con, 24013))

#save(series21, file="series21")
#save(series22, file="series22")
save(series23, file="series23")
#-------------------------------------------

#-----non-resdential construction--------------------

non.res<-make.wtable(c(2362, 23621, 23622))

#Non Residential Construction
comp.ser(c(2362, 23621, 23622))


series31<-get.wages(get.area(non.res, 24033))
series32<-get.wages(get.area(non.res, 24025))
series33<-get.wages(get.area(non.res, 24037))
series34<-get.wages(get.area(non.res, 24999))
series35<-get.wages(get.area(non.res, 24041)) #might be ridiculus

#save(series31, file="series31")
#save(series32, file="series32")
#save(series33, file="series33")
#save(series34, file="series34")
#save(series35, file="series35")
#--------------------------------------------------------

#-----Specialty Food Stores--------------------

food<-make.wtable(c(4452, 44521, 44522, 44523, 44529))

#Specialty Foods, Meat, Fish, Fruit, Other
comp.ser(c(4452, 44521, 44522, 44523, 44529))


series41<-get.wages(get.area(food, 24003))
save(series41, file="series41")

series42<-get.wages(get.area(food, 24027))
save(series42, file="series42")

series43<-get.wages(get.area(food, 24013))
save(series43, file="series43")

#--------------------------------------------------------

#----- Credit Intermediation  -----------------------
#
#5223      Activities Related to Credit Intermediation  
#52231        Mortgage and Nonmortgage Loan Brokers  
#52232        Financial Transactions Processing, Reserve, and Clearinghouse Activities   
#52239        Other Activities Related to Credit Intermediation  

round(comp.ser(c(5223, 52231, 52232, 52239)),3)

credit<-make.wtable(c(5223, 52231, 52232, 52239))

series51<-get.wages(get.area(credit, 24003))
save(series51, file="series51")

series52<-get.wages(get.area(credit, 24005))
save(series52, file="series52")

series53<-get.wages(get.area(credit, 24033))
save(series53, file="series53")


#--------------------------------------------------------


#----- Hair, Nail, and Skin Care Services  -----------------------
#
#81211        Hair, Nail, and Skin Care Services  
#812111          Barber Shops  
#812112          Beauty Salons  
#812113          Nail Salons  

round(comp.ser(c(81211, 812111, 812112, 812113)),3)

hair<-make.wtable(c(81211, 812111, 812112, 812113))

series61<-get.wages(get.area(hair, 24021))
save(series61, file="series61")

series62<-get.wages(get.area(hair, 24027))
save(series62, file="series62")

series63<-get.wages(get.area(hair, 24510))
save(series63, file="series63")

get.wages(get.area(hair, 24017))

get.1y.naics( md01, 61169)
#--------------------------------------------------------




#----- Motor Vehicle and Parts Dealers   -----------------------
#
#441    Motor Vehicle and Parts Dealers  
#4411      Automobile Dealers  
#4412      Other Motor Vehicle Dealers  
#4413      Automotive Parts, Accessories, and Tire Stores    

round(comp.ser(c(441, 4411, 4412, 4413)),3)

cars<-make.wtable(c(441, 4411, 4412, 4413))

series71<-get.wages(get.area(cars, 24001))
save(series71, file="series71")

series72<-get.wages(get.area(cars, 24009))
save(series72, file="series72")

series73<-get.wages(get.area(cars, 24011))
save(series73, file="series73")

series74<-get.wages(get.area(cars, 24023))
save(series74, file="series74")

#--------------------------------------------------------


#----- Furniture and Home Furnishings Stores   -----------------------
#
#442    Furniture and Home Furnishings Stores  
#4421      Furniture Stores  
#4422      Home Furnishings Stores  

round(comp.ser(c(442, 4421, 4422)),3)

furn<-make.wtable(c(442, 4421, 4422))

series81<-get.wages(get.area(furn, 24009))
save(series81, file="series81")

series82<-get.wages(get.area(furn, 24015))
save(series82, file="series82")

series83<-get.wages(get.area(furn, 24037))
save(series83, file="series83")


#--------------------------------------------------------

#----- rencar   -----------------------
#
#5321      Automotive Equipment Rental and Leasing  
#53211        Passenger Car Rental and Leasing  
#53212        Truck, Utility Trailer, and RV 

round(comp.ser(c(5321, 53211, 53212)),3)

rencar<-make.wtable(c(5321, 53211, 53212))

series91<-get.wages(get.area(rencar, 24003))
save(series91, file="series91")

series92<-get.wages(get.area(rencar, 24031))
save(series92, file="series92")

series93<-get.wages(get.area(rencar, 24510))
save(series93, file="series93")


#--------------------------------------------------------

comp.ser(c(5321, 53211, 53212))


