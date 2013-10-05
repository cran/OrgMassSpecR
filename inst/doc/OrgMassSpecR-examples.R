### R code from vignette source 'OrgMassSpecR-examples.Rnw'

###################################################
### code chunk number 1: OrgMassSpecR-examples.Rnw:25-26
###################################################
library(OrgMassSpecR)


###################################################
### code chunk number 2: OrgMassSpecR-examples.Rnw:33-34
###################################################
MonoisotopicMass(formula = list(C=14, H=8, Cl=4))


###################################################
### code chunk number 3: OrgMassSpecR-examples.Rnw:39-42
###################################################
MonoisotopicMass(formula = list(C=14, H=8, Cl=3))
MonoisotopicMass(formula = list(C=14, H=8, Cl=2))
# etc, ...


###################################################
### code chunk number 4: OrgMassSpecR-examples.Rnw:47-49
###################################################
MonoisotopicMass(formula = list(C=2, H=8, Cl=4, x = 12), 
                 isotopes = list(x = 13.0033548378))


###################################################
### code chunk number 5: OrgMassSpecR-examples.Rnw:54-56
###################################################
dde.dist <- IsotopicDistribution(formula = list(C=14, H=8, Cl=4))
dde.dist


###################################################
### code chunk number 6: OrgMassSpecR-examples.Rnw:60-69
###################################################
# plot
library(lattice)
print(xyplot(percent ~ mz,
  data = dde.dist,
  type = "h",
  xlab = "m/z",
  ylab = "intensity (%)",
  main = "Isotopic Distribution, DDE")
)


###################################################
### code chunk number 7: OrgMassSpecR-examples.Rnw:81-83
###################################################
hsa <- Digest(example.sequence)
head(hsa)


###################################################
### code chunk number 8: OrgMassSpecR-examples.Rnw:88-90
###################################################
hsa.sub <- subset(hsa, nchar(hsa$peptide) >= 5 & nchar(hsa$peptide) <= 12)
head(hsa.sub)


###################################################
### code chunk number 9: OrgMassSpecR-examples.Rnw:95-97
###################################################
transitions <- FragmentPeptide(c("YLYEIAR", "AEFAEVSK"))
head(transitions)


###################################################
### code chunk number 10: OrgMassSpecR-examples.Rnw:108-112
###################################################
c13.labeled <- FragmentPeptide("YLYEIAr", custom = list(code = "r", 
  mass = MonoisotopicMass(formula = list(C=6, H=12, N=4, O=1), 
                          isotopes = list(C=13.0033548378))))
head(c13.labeled)


###################################################
### code chunk number 11: OrgMassSpecR-examples.Rnw:117-119
###################################################
n15.labeled <- FragmentPeptide("YLYEIAR", N15 = TRUE)
head(n15.labeled)


###################################################
### code chunk number 12: OrgMassSpecR-examples.Rnw:129-137
###################################################
theoretical.dist <- IsotopicDistributionN("YEVQGEVFTKPQLWP", incorp = 0.99)
print(xyplot(percent ~ mz,
  data = theoretical.dist,
  type = "h",
  xlab = "m/z",
  ylab = "intensity (%)",
  main = "Theoretical Isotopic Distribution, YEVQGEVFTKPQLWP, 99% 15N")
)


###################################################
### code chunk number 13: OrgMassSpecR-examples.Rnw:144-153
###################################################
example.spectrum.labeled$percent <- with(example.spectrum.labeled, 
  intensity / max(intensity) * 100)
print(xyplot(percent ~ mz,
  data = example.spectrum.labeled,
  type = "l",
  xlab = "m/z",
  ylab = "intensity (%)",
  main = "Measured Isotopic Distribution, YEVQGEVFTKPQLWP")
)


