
'''
This script has been used to calculate some physical properties of proteins and save a backbone pdb 

'''

library(bio3d)
library(Peptides)
library(magrittr)
library(ggplot2)

pdb <- read.pdb('protein.pdb')
pdb_df <- pdb$atom %>% as.data.frame()

t <- torsion.pdb(pdb)  # angle torsions


# GET AA SEQUENCE

seq <- pdb_df[which(pdb_df$elety == 'N'), 'resid']  # extract 3 letters aa sequence 

conversion <- '
ALA A
ARG R
ASN N
ASP D
CYS C
GLU E
GLN Q
GLY G
HIS H
ILE I
LEU L
LYS K
MET M
PHE F
PRO P
SER S
THR T
TRP W
TYR Y
VAL V'

conversion <- read.delim(text=conversion,
                         sep = ' ',
                         col.names=c('3let','1let'),
                         header = F)

# function to replace 3 letters aa sequence with 1 letter 
replace <- function(seq){
  vec <- c()
  for(x in seq){
    vec <- append(vec, conversion[which(conversion$X3let == x), 'X1let'])}
  return(vec)
}
seqq <- replace(seq)


# HYDROPHOBICITY PLOT

hyd <- data.frame(x=seqq, y=hydrophobicity(seqq))
hyd$x <- paste0(hyd$x, 1:174)

hyd <- hyd %>% 
  dplyr::mutate(mycolor = ifelse(y>0, "type1", "type2"))


ggplot(hyd, aes(x=x, y=y)) +
  geom_segment(aes(x=x, xend=x, y=0, yend=y, color=mycolor), size=1.3, alpha=0.9) +
  theme_light() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
  ) +
  xlab("aa sequence") +
  ylab("hydrophobicity")


# TOTAL CHARGE AND pI

charge(seqq)
sum(charge(seqq))
pI(seqq)
mean(pI(seqq))


# SAVING BACKBONE PDB

b.inds <- atom.select(pdb, elety = c("CA", "C", "N"))
backpdb <- trim.pdb(pdb, b.inds)
write.pdb(backpdb, file="backbone")




