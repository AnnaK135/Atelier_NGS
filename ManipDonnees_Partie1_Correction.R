#####
# Gaëlle LELANDAIS < gaelle.lelandais@u-psud.fr >
# Mélina GALLOPIN  < melina.gallopin@u-psud.fr >
#####


#-------------------------------------------------------
#-------------------------------------------------------
#
# Ce script R est pour les étudiants de Master 2. Il comporte des commandes R utiles pour analyser 
# les données RNAseq Ostreococcus tauri.
#
# Ces données sont présentées en détail dans l'article suivant :
#-------------------------------------------------------
# Ostreococcus tauri is a new model green alga for studying iron metabolism in eukaryotic phytoplankton.
#
# Lelandais G, Scheiber I, Paz-Yepes J, Lozano JC, Botebol H, Pilátová J, Žárský V, Léger T, Blaiseau PL, 
# Bowler C, Bouget FY, Camadro JM, Sutak R, Lesuisse E.
#
# BMC Genomics. 2016 May 3;17:319. doi: 10.1186/s12864-016-2666-6.
#-------------------------------------------------------

#-------------------------------------------------------
# Lecture des fichiers de données
#-------------------------------------------------------

# 1) Données de comptage
allData = read.table("mapping_rawdata_allGenes.txt", header = T, row.names = 1)
# --> les noms des gènes servent d'étiquettes (argument : "row.names")
# --> A noter que la longueur des gènes est indiquée en première colonne (nommée "length")
# --> Les autres colonnes correspondent aux différents échantillons d'ARN qui ont été séquencé.

# 2) Description des échantillons
sampleInfo = read.table("sample_info.txt", header = T, row.names = 1)
# --> Chaque condition expérimentale associe 4 états : 
# -- Fer, 
# -- Luminosité, 
# -- Temps de mesure, 
# -- Type de carence en Fer.
# --> 3 réplicats biologiques sont disponibles par état testé (sauf HCA-24 et HC-25 : S7)

#-------------------------------------------------------
# Manipulation des informations
#-------------------------------------------------------

#######
# Question : Quelles sont les conditions de croissance des cellules en terme de :
#
# 1- luminosité ?
# 2- contenu en Fer du milieu de culture ?
# 3- temps de cinétique ?
# 4- type de réponse étudiée "Short term" ou "Long term" ?
#
# Répondez pour les échantillons nommés :
# HCA.10 :
# HCA.20 :
# HCA.30 :
# HCA.40 :
#######

# REPONSE :
sampleInfo[c("HCA.10", "HCA.20", "HCA.30", "HCA.40"),]

#######
# Question : Quels sont les échantillons pour lesquelles les conditions de croissance
# des cellules sont :
#
# 1- luminosité : LIGHT
# 2- contenu en Fer du milieu de culture : YES
# 3- temps de cinétique : 3H
# 4- type de réponse étudiée "Short term" ou "Long term" : SHORT TERM
#
# Vous associez les noms des échantillons dans un vecteur nommé "vec1"
#######

# REPONSE :
vec1 = row.names(sampleInfo[(sampleInfo[, "Light"] == "LIGHT") & (sampleInfo[, "Iron"] == "YES") & 
                            (sampleInfo[, "Time"] == "3H") & (sampleInfo[, "Condition"] == "ST"),])

#######
# Question : Quels sont les échantillons pour lesquelles les conditions de croissance
# des cellules sont :
#
# 1- luminosité : LIGHT
# 2- contenu en Fer du milieu de culture : NO
# 3- temps de cinétique : 3H
# 4- type de réponse étudiée "Short term" ou "Long term" : SHORT TERM
#
# Vous associez les noms des échantillons dans un vecteur nommé "vec2"
#######

# REPONSE :
vec2 = row.names(sampleInfo[(sampleInfo[, "Light"] == "LIGHT") & (sampleInfo[, "Iron"] == "NO") & 
                            (sampleInfo[, "Time"] == "3H") & (sampleInfo[, "Condition"] == "ST"),])

#######
# Question : Quelles sont les valeurs de comptage, pour tous les gènes, associées aux noms des
# échantillons identifiés aux deux questions précédentes ?
#
# Vous associez ces valeurs dans un tableau nommé "subData".
#######

# REPONSE :
subData = allData[,c(vec1, vec2)]

#######
# Question : Combien de gènes et combien d'échantillons avez vous sélectionné dans votre tableau "subData" ?
#######

# REPONSE : 7699 gènes et 6 échantillons
dim(subData)

#-------------------------------------------------------
# Ecriture des résultats
#-------------------------------------------------------

# Les colonnes du tableau sont renommée, afin de distinguer les réplicats de la condition +Fe des réplicats de 
# la conditions -Fe.

colnames(subData) = c(rep("Fe", 3), rep("noFe", 3))

write.table(subData, file = "mapping_rawData_allGenes_sixSamples.txt",
            quote = F, sep = "\t")
