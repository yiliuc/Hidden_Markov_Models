library(RCurl)
library(jsonlite)

url <- "https://data.cityofchicago.org/resource/ijzp-q8t2.json?$select=date_trunc_ym(occurred_on_date),count(*)&$group=date_trunc_ym(occurred_on_date)&$limit=1000"
crime_data <- fromJSON(url)


library(datasets)
data <- data("UKLungDeaths")
