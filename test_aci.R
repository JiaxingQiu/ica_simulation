setwd(dirname(rstudioapi::getSourceEditorContext()$path))


path = paste0("./utils")
flst = list.files(path)
sapply(c(paste(path,flst,sep="/")), source, .GlobalEnv)
source("./sim_conditions.R")


# Define matrix P
v1 <- matrix(c(0.626, 0.536, 0.271, 0.715, 0.074, 0.357, 0.317,
              0.374, 0.464, 0.729, 0.285, 0.926, 0.643, 0.683),
            nrow = 7, byrow = F)

# Define matrix Q
v2 <- matrix(c(0.601, 0.410, 0.571, 0.562, 0.365, 0.347, 0.499,
              0.399, 0.590, 0.429, 0.438, 0.635, 0.653, 0.501),
            nrow = 7, byrow = F)


# Get EP and EQ
e1 <- pw_l1_norm(v1)
e2 <- pw_l1_norm(v2)

# Get NDC 
ndc = NDC(v1, v2)
ndc # correct

endc <- c()
for(i in 1:100){
  endc = c(endc, ENDC(v1, v2))
  
}
mean(endc) # correct ~ 0.8189

ACI(v1, v2) # -0.004829699
#if diag included, wrong! 0.2242986


