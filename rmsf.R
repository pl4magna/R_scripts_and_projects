
'''
This script compares rmsf of two different protein structures, RBP4 wild type
and RBP4 with G75D mutation, in order to make some further considerations about
their stabilities. Rmsf are formerly plotted individually in relation to 
hdrophobicity and secondary structure, lastly are plotted together.

'''

library(Peptides)
library(ggplot2)
library(bio3d)

# READING FILES

#-------G75D_olo_mut-------------

rmsf <- unlist(readXVG("rmsf.xvg"))[-1] # erase '&' character in the last row

# select only rmsf values
selection <- rep(c(TRUE, FALSE), length(rmsf) / 2)
rmsf_df <- data.frame(rmsf[!selection], rmsf[selection])
colnames(rmsf_df) <- c('atom', 'rmsf')

xaver <- read.pdb('xaver.pdb')
df <- as.data.frame(xaver$atom)

rmsf_df$aa <- paste0(df$resid, df$resno)
rmsf_df$type <- df$resid

#-------G75D_olo_wt-------------

rmsf_wt <- unlist(readXVG("rmsf_bb.xvg"))[-1] # erase '&' character in the last row

selection_wt <- rep(c(TRUE, FALSE), length(rmsf_wt) / 2)
rmsf_df_wt <- data.frame(rmsf_wt[!selection_wt], rmsf_wt[selection_wt])
colnames(rmsf_df_wt) <- c('atom', 'rmsf')

xaver_wt <- read.pdb('xaver_bb.pdb')
df_wt <- as.data.frame(xaver_wt$atom)

rmsf_df_wt$aa <- paste0(df_wt$resid, df_wt$resno)
rmsf_df_wt$type <- df_wt$resid


# PLOTTING

pdb_wt <- read.pdb(file.path("5NU7", "5nu7.pdb"))

sh.start <- pdb_wt$sheet$start
sh.end <- pdb_wt$sheet$end
sheet <- c()
for(n in 1:9){sheet <- append(sheet, seq(sh.start[n], sh.end[n]))}

he.start <- pdb_wt$helix$start
he.end <- pdb_wt$helix$end
helix <- c()
for(n in 1:3){helix <- append(helix, seq(he.start[n], he.end[n]))}

polar <- c('SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO')
hydrophobic <- c('ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP')
positive <- c('ARG', 'HIS', 'LYS')
negative <- c('ASP', 'GLU')

rmsf_df$type <- as.character(rmsf_df$type)

for(x in rmsf_df$type){rmsf_df$type[rmsf_df$type %in% polar] <- 'polar'}
for(x in rmsf_df$type){rmsf_df$type[rmsf_df$type %in% hydrophobic] <- 'hydrophobic'}
for(x in rmsf_df$type){rmsf_df$type[rmsf_df$type %in% positive] <- 'positive'}
for(x in rmsf_df$type){rmsf_df$type[rmsf_df$type %in% negative] <- 'negative'}

rmsf_df$ss <- 'coil'
rmsf_df$resid <- df$resno

for(x in rmsf_df$ss){rmsf_df$ss[rmsf_df$resid %in% helix] <- 'helix'}
for(x in rmsf_df$ss){rmsf_df$ss[rmsf_df$resid %in% sheet] <- 'sheet'}

# rmsf and hydrophobicity

ggplot(data=rmsf_df, aes(x=resid, y=as.numeric(rmsf))) +
  geom_col(aes(color=type), width = 0.5)

# rmsf and secondary structure

ggplot(data=rmsf_df, aes(x=resid, y=as.numeric(rmsf))) +
  geom_col(aes(color=ss), width = 0.5)


#-----------------------OVERLAY-----------------------------

rmsf_df_wt$resid <- df_wt$resno

ggplot(rmsf_df, aes(x=resid, y=as.numeric(rmsf))) + 
  geom_line() + 
  geom_line(data=rmsf_df_wt, color='red')


  
  
  
  
  


