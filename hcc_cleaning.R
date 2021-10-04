## Dealing with missing data
cat.var = which(!(hcc_info$data_type %in% "continuous"))
hcc_temp_imputed = knn.impute(hcc_temp_b4_impute, k= 10, cat.var = cat.var)
hcc_imputed = data.frame(hcc_temp_imputed)
names(hcc_imputed) <- hcc_info$var_names



## Structure setting
nom_nodes = which(hcc_info$data_type == "nominal")
ord_nodes = which(hcc_info$data_type == "ordinal")

hcc_imputed[nom_nodes] <- as.data.frame(lapply(hcc_imputed[nom_nodes], 
                                               function(x) as.factor(x)))
hcc_imputed[ord_nodes] <- as.data.frame(lapply(hcc_imputed[ord_nodes], 
                                               function(x) as.ordered(x)))


## Excluding "alch.s" and "smok.s" since those are duplicated info
hcc_imputed[c("alch.s", "smok.s")] <- NULL

## Transform some of continuous nodes to follow a normal distribution

# age =====
cur_var = "age"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- (temp)^2

# inr =====
cur_var = "inr"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 1/temp

# total_bil ====
cur_var = "tot.bil"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- log(temp)

## alt ==== 
cur_var = "alt"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- log(temp)

# ast ====
cur_var = "ast"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- log(temp)

# ggt ====
cur_var = "ggt"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- log(temp)

# alp ====
cur_var = "alp"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- (temp)^(1/4)

# major ====
cur_var = "major"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- log(temp)

# ferr ====
cur_var = "ferr"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- log(temp+10)


## Discretizing some of continuous nodes

## alch ====
cur_var = "alch"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 0, 99,119, max(temp)), labels = 1:4, ordered_result = T)

## cig ====
cur_var = "cig"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 0, 29, 39, 100, 520), labels = 1:5, ordered_result = T)

## afp  ====
cur_var = "afp"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 10, 100, 1000, 10000, max(temp)), labels = 1:5, ordered_result = T)

## plat ====
cur_var = "plat"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 37999, 253999, max(temp)), labels = 1:3, ordered_result = T)

## tp ====
cur_var = "tp"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 10, max(temp)), labels = 1:2, ordered_result = T)

## crea ====
cur_var = "crea"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 2.18, max(temp)), labels = 1:2, ordered_result = T)

## dir.bil ====
cur_var = "dir.bil"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 5, max(temp)), labels = 1:2, ordered_result = T)

## iron ====
cur_var = "iron"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 35, 95, 400), labels = 1:3, ordered_result = T)

## nodules ====
cur_var = "nodules"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- factor(temp, ordered = T)
#.. converting as an ordinal-level data

## leucocytes ====
cur_var = "leuc"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 2000, max(temp)), labels = 1:2, ordered_result = T)

## sat ====
cur_var = "sat"
temp = hcc_imputed[[cur_var]]
hcc_imputed[[cur_var]] <- 
  cut(temp, c(-1, 19, 43, 60, max(temp)), labels = 1:4, ordered_result = T)
