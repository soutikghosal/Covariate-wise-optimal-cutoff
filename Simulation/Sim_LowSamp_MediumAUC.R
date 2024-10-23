
Jobid=Sys.getenv('SLURM_ARRAY_TASK_ID')
source("Sim_func.R")

n = 50
AUCcat = "Medium"
seed = Jobid
fit.BN.equal = cutoff.func(seed = seed, N = n, data.gen = BN.equal, sep=AUCcat)
fit.BN.unequal = cutoff.func(seed = seed, N = n, data.gen = BN.unequal, sep=AUCcat)
fit.Skewed.I = cutoff.func(seed = seed, N = n, data.gen = Skewed.I, sep=AUCcat)
fit.Skewed.II = cutoff.func(seed = seed, N = n, data.gen = Skewed.II, sep=AUCcat)
fit.Skewed.III = cutoff.func(seed = seed, N = n, data.gen = Skewed.III, sep=AUCcat)
fit.Mixed.I = cutoff.func(seed = seed, N = n, data.gen = Mixed.I, sep=AUCcat)
fit.Mixed.II = cutoff.func(seed = seed, N = n, data.gen = Mixed.II, sep=AUCcat)

out = list(fit.BN.equal = fit.BN.equal, fit.BN.unequal = fit.BN.unequal, fit.Skewed.I = fit.Skewed.I,
           fit.Skewed.II = fit.Skewed.II, fit.Skewed.III = fit.Skewed.III, fit.Mixed.I = fit.Mixed.I,
           fit.Mixed.II = fit.Mixed.II)

fname2 <- paste("res_LowSamp_MediumAUC_",Jobid,".rds",sep="")
saveRDS(out,file=fname2)

