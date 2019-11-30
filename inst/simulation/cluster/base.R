args = commandArgs(TRUE)
source("setup.R")

i <- args[1]
set.seed(i)

foo <- replicate(nn, do(P1, "P2", P3, P4))
write.table(foo, file = paste(c("output", P1, "P2", P3, P4, i), collapse = "-"),
            row.names = FALSE, col.names = FALSE)
proc.time()
