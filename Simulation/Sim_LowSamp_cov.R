
Jobid=Sys.getenv('SLURM_ARRAY_TASK_ID')
# Jobid = 5
source("Sim_func_cov.R")

tt = Sys.time()
n = 50
seed = Jobid
fit.BN = cutoff.func(seed = seed, N = n, data.gen = BN, x = c(0, 1),
                     burnin = 5000, numb.iter = 5000)
fit.Skewed = cutoff.func(seed = seed, N = n, data.gen = Skewed, x = c(0, 1),
                         burnin = 5000, numb.iter = 5000)
fit.Mixed = cutoff.func(seed = seed, N = n, data.gen = Mixed, x = c(0, 1),
                        burnin = 5000, numb.iter = 5000)

out = list(fit.BN = fit.BN, fit.Skewed = fit.Skewed, fit.Mixed = fit.Mixed)
Sys.time() - tt

fname2 <- paste("res_LowSamp_cov_",Jobid,".rds",sep="")
saveRDS(out,file=fname2)

