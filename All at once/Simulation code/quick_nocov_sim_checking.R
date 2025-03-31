
## possible variations of sep: "Low", "Medium", "High"
## possible variations of n: 50, 100, 500
## possible variations of seed: any number

source("Nocov sim func.R")
tt = Sys.time()
ff = final.fit(seed = 12345, n = 50, sep = "Low")
Sys.time() - tt

