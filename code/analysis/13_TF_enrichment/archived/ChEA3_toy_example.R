library(httr)
library(jsonlite)

genes <- c("SMAD9", "FOXO1", "MYC", "STAT1", "STAT3", "SMAD3")

url <- "https://maayanlab.cloud/chea3/api/enrich/"
encode <- "json"
payload <- list(query_name = "myQuery", gene_set = genes)

# POST to ChEA3 server
response <- POST(url = url, body = payload, encode = encode)
json <- content(response, "text")

# results as list of R dataframes
results <- fromJSON(json)
