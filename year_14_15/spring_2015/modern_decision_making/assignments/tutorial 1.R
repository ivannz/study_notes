setwd("~/study_notes/year_14_15/spring_2015/modern_decision_making/assignments/")
rm( list = ls( all.names = TRUE ) ) ; invisible( gc() )

Adv_data <- read.csv( "data/advertisement.csv" )

names( Adv_data )

str( Adv_data )
