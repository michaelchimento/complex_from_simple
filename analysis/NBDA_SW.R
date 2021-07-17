
# PHIL TRANS 2021 ------------------------------------

# 1) Install packages and load data  -----------------------------------------------------
# install NBDA package from github via devtools
install.packages("devtools")
# install.packages("rtools40")
library(devtools)
devtools::install_github("whoppitt/NBDA")
library(NBDA)

# load the cleaned up solve data for all three diffusions (dial, complex 1st gen, complex 2nd gen)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("../data/df_solves_cleaned.Rda")
load("../data/df_visitors_all.Rda")

######################################################################### #


# 2) DIAL DIFFUSION -------------------------------------------------------

# 2.1 Read in network data and individual attributes  ------------------------------------------------
load("../data/NBDA_input_data.Rda")

# inspect the gbi (group by individual matrix)
gbi <- NBDA_input_data$data.objects.2016$gbi
gbi[1:10, 1:10]

# load library asnipe for network calculation
# install.packages("asnipe")
library(asnipe)

# 2.2. modify individual attributes  --------------------------------------------------------
# subset to great tits only
ids <- NBDA_input_data$data.objects.2016$ids
ids.mod <- ids
ids.mod <- subset(ids.mod, ids.mod$species=="greti")

# subset the gbi to gretis only
gbi <- gbi[,colnames(gbi)%in% ids.mod$id]
gbi[1:10, 1:10]
# each row represents a group, each column an individual.
# individuals have an entry of 1 if part of the group, a 0 if not

# re-assign codes for sex
ids.mod$sex[ids.mod$sex=="U"] <- 0 # birds of unknown sex get a 0
ids.mod$sex[ids.mod$sex=="M"] <- 0.5 # males are assigned 0.5
ids.mod$sex[ids.mod$sex=="F"] <- -0.5 # females are assigned -0.5

# re-assign codes for age classes
ids.mod$age[ids.mod$age=="ad"] <- 0.5 # adults are assigned 0.5
ids.mod$age[ids.mod$age=="juv"] <- -0.5 # juveniles are assigned -0.5

group_info <- NBDA_input_data$data.objects.2016$group_info

# 2.3. prepare NBDA data objects ------------------------------------------
# create a function that automatically creates the NBDA data objects from the input data
prepare.NBDA.data <- function(input.file, locations, label, start.date, end.date, patch){

  # subset the input file to the correct experiment (dial_diffusion)
  input.file.sub <- subset(input.file, input.file$PATCH==patch & input.file$SPECIES=="GRETI" & input.file$experiment=="dial_diffusion")

  # extract who the demonstrators are (in this case trained in captivity)
  demos <- as.character(unique(subset(input.file.sub$RING, input.file.sub$DEMO==TRUE)))

  # a bird is considered a learner if it has solved the dial task three times
  # the time of its first solve is taken as the time of acquisition

  # prepare empty vectors to store the ID of the learner, the number of solves and the time of the first solve
  learners.list <- NULL
  num_solves_all <- NULL
  time_solves <- NULL

  for (i in as.vector(unique(input.file.sub$RING))){
    sub <- subset(input.file.sub, input.file.sub$RING==i)
    num.solves <- length(sub[,1])
    solve.time <- as.numeric(difftime(min(sub$time_stamp), min(input.file.sub$time_stamp), units="days"))
    learners.list[which(unique(input.file.sub$RING)==i)] <- i
    num_solves_all[which(unique(input.file.sub$RING)==i)] <- num.solves
    time_solves[which(unique(input.file.sub$RING)==i)] <- solve.time
  }


  # reassign column names
  learners.list_df <- na.omit(cbind.data.frame(learners.list, num_solves_all,  time_solves))
  colnames(learners.list_df) <- c("RING", "num_solves", "time_acq")
  learners.list_df <- learners.list_df[order(learners.list_df$time_acq),]

  # remove those with fewer than 3 solves
  learners.list_df <- subset(learners.list_df, learners.list_df$num_solves>=3)

  # restrict the data frame to great tits only that are in the ILV file
  learners.list_df <- subset(learners.list_df, learners.list_df$RING%in%ids.mod$id)

  # remove the demonstrators
  learners.list_df <- subset(learners.list_df , !(learners.list_df$RING%in%demos))

  # create network from group by individual matrix
  # subset the gbi to the correct locations and times for each diffusion
  which.groups <- subset(group_info, group_info$location%in% locations & group_info$start_time>=start.date & group_info$stop_time<=end.date)
  which.groups.num <- rownames(which.groups)

  sub.gbi <- gbi[which.groups.num,] # subset the gbi to the right locations and times
  sub.gbi <- sub.gbi[, colSums(sub.gbi)>9] # remove birds that were seen fewer than 10 times

  ## remove learners.list that were never seen on network feeders (as they can only only by asocial learning)
  learners.list_df <- subset(learners.list_df, learners.list_df$RING%in% colnames(sub.gbi))

  # recalculate the time of acquisition by removing the days where the puzzle boxes weren't out (network days)
  # days 0-4 puzzle box
  # days 5+6 network data
  # days 7-11 puzzle box (so subtract two days from the acquisition of birds that learned between day 7-11)
  # days 12+13 network data
  # day 14-18 (subtract 4 days from the acquisition time for birds that have learned between day 15-19)
  # days 19+20
  # day 21-25 (subtract 6 days)

  learners.list_df$time_acq_adj <- NA

  for (i in 1:length(learners.list_df[,1])){
    time <- learners.list_df[i, "time_acq"]
    if(time > 7 & time < 12){
      time.adj <- time-2
    } else if(time > 14 & time < 19){
      time.adj <- time-4
    } else if(time > 21 & time < 26){
      time.adj <- time-6
    } else {time.adj <- time}
    learners.list_df[i, "time_acq_adj"] <- time.adj
  }


  # generate association network
  net <- get_network(association_data = sub.gbi, data_format="GBI", association_index = "SRI")
  net <- net[order(rownames(net)), order(colnames(net))]

  ## the association matrix needs to be in a matrix array
  assocMatrix <- array(data = net, dim=c(nrow(net), ncol(net), 1))
  class(assocMatrix)

  # get order of acquisition for learners.list
  learners.list_OAc <- NULL

  for (i in learners.list_df$RING){
    pos <- which(colnames(net)==i)
    learners.list_OAc <- c(learners.list_OAc, pos)
  }
  # the values refers to the position of each learner in the network in the order they have learned

  # prepare demonstrator vector
  # each demonstrator gets a 1, all naive individuals get a 0
  demo_vec <- rep(0, length(colnames(net)))
  demo_vec[which(colnames(net)==demos)] <- 1


  ## get ILVs ready
  # subset the ILVs to the those IDs that were present in the network
  ids.sub <- subset(ids.mod, ids.mod$id %in% rownames(net))

  sex <- as.matrix(as.numeric(ids.sub$sex))
  age <- as.matrix(as.numeric(ids.sub$age))

  assign(paste("age", label, sep="_"), age, envir = .GlobalEnv)
  assign(paste("sex", label, sep="_"), sex, envir = .GlobalEnv)

  ILVs <- c(paste("age", label, sep="_"), paste("sex", label, sep="_"))

  # ILVs need to be defined as a character vector
  assign(paste("ILVs", label, sep="_"), ILVs)


  # create NBDA Data Object
  # we are defining unconstrained models (ILVs can influence social and asocial learning rate independently)
  nbdaData <- nbdaData(label=label,
                       assMatrix = assocMatrix,
                       asoc_ilv = get(paste("ILVs", label, sep="_")),
                       int_ilv = get(paste("ILVs", label, sep="_")),
                       multi_ilv = "ILVabsent",
                       orderAcq = learners.list_OAc,
                       timeAcq = learners.list_df$time_acq_adj,
                       demons= demo_vec)


  object <- NULL

  object$nbdadata <- nbdaData
  object$label <- label

  object$learners.list <- learners.list_df
  object$num.birds <- length(rownames(net))
  object$num.learners.list <- length(learners.list_df$RING)
  object$IDs <- rownames(net)

  return(object)

}

# prepare NBDA Data objects for each site
# by specifying the start and end date, and the locations of the network feeder sites


BB_NBDA_DATA <- prepare.NBDA.data(input.file = as.data.frame(df_solves),
                                  patch = "BB",
                                  label = "BB_dial",
                                  start.date = as.POSIXct("26/12/2015", format="%d/%m/%Y"),
                                  end.date = as.POSIXct("24/01/2016", format="%d/%m/%Y"),
                                  locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f")
)

GW_NBDA_DATA <- prepare.NBDA.data(input.file = as.data.frame(df_solves),
                                  patch = "GW",
                                  label = "GW_dial",
                                  start.date = as.POSIXct("09/01/2016", format="%d/%m/%Y"),
                                  end.date = as.POSIXct("07/02/2016", format="%d/%m/%Y"),
                                  locations <- c("2c", "2b", "2d", "2h", "2g")
)


MA_NBDA_DATA <- prepare.NBDA.data(input.file = as.data.frame(df_solves),
                                  patch = "MA",
                                  label = "MA_dial",
                                  start.date = as.POSIXct("02/01/2016", format="%d/%m/%Y"),
                                  end.date = as.POSIXct("31/01/2016", format="%d/%m/%Y"),
                                  locations <- c("6h", "6g", "7a", "6i", "6f", "7b", "7c", "7d")
)

## each of the objects contains a slot with the NBDA data object ($nbdadata)
# one with the label $label (e.g. BB_dial)
# one with a data frame with details on the learners.list (great tits only, without demonstrators) ($learners.list)
# one with the number of birds in total ($num.birds)
# one with the number of learners.list ($num.learners.list)
# and finally, one with the IDs (ring numbers) of the birds (alphabetically)


# 2.4. Summary ------------------------------------------------------------


# we extract numbers for the summary results part:
BB_NBDA_DATA$num.learners.list
# [1] 29 learners
BB_NBDA_DATA$num.birds-1 # minus the demonstrator
# out of [1] 72 naive birds seen at network feeders
BB_NBDA_DATA$num.learners.list/(BB_NBDA_DATA$num.birds-1)
# [1] 0.4027778

GW_NBDA_DATA$num.learners.list
# [1] 12 learners
GW_NBDA_DATA$num.birds-1 # minus the demonstrator
# out of [1] 81 naive birds seen at network feeders
GW_NBDA_DATA$num.learners.list/(GW_NBDA_DATA$num.birds-1)
# [1] 0.1481481

MA_NBDA_DATA$num.learners.list
# [1] 8 learners
MA_NBDA_DATA$num.birds-1
# out of [1] 119 naive birds seen at network feeders
MA_NBDA_DATA$num.learners.list/(MA_NBDA_DATA$num.birds-1)
# [1] 0.06722689
mean(c(0.4027778, 0.1481481, 0.06722689))


# extract the number of dial solutions (right, left, total) together
left <- length(as.data.frame(subset(df_solves, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_l"))[,1])
right <- length(as.data.frame(subset(df_solves, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_r"))[,1])

left+right
# [1] 18421 total number of dial solves

# extract the number of left solves and right solves at each site to test whether there is a significant bias to one side
# CP: right
# GW: left
# BB: left
# MA: right

# CP: note that for CP, we do not have network data available
CP.num.left <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_l" & df_solves$PATCH=="CP"))
CP.num.right <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_r" & df_solves$PATCH=="CP"))

# GW:
GW.num.left <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_l" & df_solves$PATCH=="GW"))
GW.num.right <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_r" & df_solves$PATCH=="GW"))

# BB
BB.num.left <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_l" & df_solves$PATCH=="BB"))
BB.num.right <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_r" & df_solves$PATCH=="BB"))

# MA
MA.num.left <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_l" & df_solves$PATCH=="MA"))
MA.num.right <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_r" & df_solves$PATCH=="MA"))

# BO
BO.num.left <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_l" & df_solves$PATCH=="BO"))
BO.num.right <- length(subset(df_solves$solution, df_solves$experiment=="dial_diffusion" & df_solves$solution=="dial_r" & df_solves$PATCH=="BO"))


dial.max <- c(CP.num.right/(CP.num.right+CP.num.left), GW.num.left/(GW.num.left+GW.num.right), BB.num.right/(BB.num.right+BB.num.left), MA.num.right/(MA.num.left+MA.num.right))
dial.min <- c(CP.num.left/(CP.num.right+CP.num.left), GW.num.right/(GW.num.left+GW.num.right), BB.num.left/(BB.num.right+BB.num.left), MA.num.left/(MA.num.left+MA.num.right))

# perform a t test on the proportion of solves
t.test(dial.min, dial.max)



# number of demonstrators in total
BB.demo <-length(unique(as.data.frame(subset(df_solves$RING, df_solves$PATCH=="BB" & df_solves$experiment=="dial_diffusion" & df_solves$DEMO==TRUE))))
GW.demo <-length(unique(as.data.frame(subset(df_solves$RING, df_solves$PATCH=="GW" & df_solves$experiment=="dial_diffusion" & df_solves$DEMO==TRUE))))
MA.demo <-length(unique(as.data.frame(subset(df_solves$RING, df_solves$PATCH=="MA" & df_solves$experiment=="dial_diffusion" & df_solves$DEMO==TRUE))))

BB.demo+GW.demo+MA.demo
# 3 demonstrator birds

# we now calculate the number of non-demonstrator birds
length(unique(c(BB_NBDA_DATA$IDs, GW_NBDA_DATA$IDs, MA_NBDA_DATA$IDs)))
# [1] 270 - 3 demonstrators = 267

# num learners.list
BB_NBDA_DATA$num.learners.list+GW_NBDA_DATA$num.learners.list+ MA_NBDA_DATA$num.learners.list
# [1] 49

# proportion of learners.list out of non-demonstrators
49/267


# 2.5. Prepare Constraints Vector Matrix ----------------------------------
# this function allows the user to create the matrix specifying which models should be run (specific to an NBDA data object)

create.constraints.Vect.Matrix <- function(NBDA_data_object, num_networks, num_ILVs){
  suppressWarnings(
    if(NBDA_data_object@asoc_ilv=="ILVabsent"){
      num.ILV.asoc <- 0
    } else {num.ILV.asoc <- length(NBDA_data_object@asoc_ilv)})

  suppressWarnings(
    if(NBDA_data_object@int_ilv=="ILVabsent"){
      num.ILV.int<- 0
    } else {num.ILV.int<- length(NBDA_data_object@int_ilv)})

  suppressWarnings(
    if(NBDA_data_object@multi_ilv=="ILVabsent"){
      num.ILV.multi <- 0
    } else {num.ILV.multi <- length(NBDA_data_object@multi_ilv)})

  vector <- seq(1:(num_networks+num.ILV.asoc+num.ILV.int+num.ILV.multi))

  count <- 0 # create an object 'count', which starts on 0

  constraintsVect <- matrix(nrow = 10000000, ncol=(num_networks+num.ILV.asoc+num.ILV.int+num.ILV.multi)) # create a matrix to save the combination of parameters in
  constraintsVect[1,] <- vector # the first row gets filled with a sequence from 1:8 (all parameters will be estimated, none are set to 0)

  for (i in 1:(length(vector)-1)){ # a loop for each number of parameters to be estimated
    array <- combn(vector, i, FUN = NULL, simplify = TRUE) # for each number of paramters to be estiamted (e.g. 2) create all possible combinations of numbers between 1:12 (e.g. 2&8, 1&5 etc)

    for (j in 1:length(array[1,])){ # for each of those combinations
      vector2 <- seq(1:((num_networks+(num.ILV.asoc+num.ILV.int+num.ILV.multi))-i)) # create a second vector with 11-i free spaces
      position <- array[,j] # for each created combination
      count <- count+1 # add +1 to the count

      for (k in position){ # at each possible position
        vector2 <- append(vector2, 0, after=k-1) # add a 0 (e.g. 1 0 2 3 ...; 1 2 0 3 4 5 ...; 1 2 3 0 4 5 ....)
      }
      constraintsVect[count+1,] <- vector2 # and save the resulting order in a matrix
    }
  }


  constraintsVect <- na.omit(constraintsVect) # remove all NAs from the matrix

  # extract which columns are networks
  col.networks <- c(1:num_networks)

  col.names <- NULL

  if(num.ILV.asoc!=0){
    col.names <- rep("asoc", num.ILV.asoc)
  }

  if(num.ILV.int!=0){
    col.names <- c(col.names, rep("int", num.ILV.int))
  }

  if(num.ILV.multi!=0){
    col.names <- c(col.names, rep("multi", num.ILV.multi))
  }

  colnames(constraintsVect) <- c(rep("network", num_networks), col.names)

  constraintsVect <- as.matrix(as.data.frame(constraintsVect))

  # extract the models containing any social network

  social.models <- rep(NA, length(constraintsVect[,1]))

  for (k in 1:length(constraintsVect[,1])){
    sum <- sum(constraintsVect[k,1:num_networks])
    if(sum!=0){
      social.models[k] <- k
    }
  }
  social.models <- as.vector(na.omit(social.models))

  social.models.matrix <- constraintsVect[social.models,]

  # if multiplicative models are fit, we need to adjust the matrix
  # if the multiplicative slots are filled, it automatically fits the parameter for asoc and social (just constrained to be the same)
  # meaning that we can remove it from the asoc and int slot

  if(num.ILV.multi!=0){
    social.models.retain <- rep(NA, length(social.model.matrix[,1]))
    multi.models <- rep(NA, length(social.models.matrix[,1]))
    for (k in 1:length(social.models.matrix[,1])){
      sum <- sum(social.models.matrix[k,which(colnames(social.models.matrix)=="multi")])
      sum2 <- sum(social.models.matrix[k, c(which(colnames(social.models.matrix)=="asoc"),which(colnames(social.models.matrix)=="int"))])
      if(sum!=0 & sum2==0){ # if multi models are fit and int and asoc are set to 0
        multi.models[k] <- k # then retain the model
      } else if (sum==0){
        social.models.retain[k] <- k
      }
    }

    multi.models <- as.vector(na.omit(multi.models))
    social.models.retain <- as.vector(na.omit(social.models.retain))

    models.to.retain <- c(multi.models, social.models.retain)

    # these models are retained
    retain.matrix.soc <- social.models.matrix[models.to.retain,]

    social.models.matrix <- retain.matrix.soc
  }

  # extract the models containing no social network

  asocial.models <- rep(NA, length(constraintsVect[,1]))

  for (k in 1:length(constraintsVect[,1])){
    sum <- sum(constraintsVect[k,1:num_networks])
    if(sum==0){
      asocial.models[k] <- k
    }
  }
  asocial.models <- as.vector(na.omit(asocial.models))

  asocial.models.matrix <- constraintsVect[asocial.models,]

  cols.asoc <- which(colnames(constraintsVect)=="asoc")

  asocial.retain <- rep(NA, length(asocial.models))
  for (k in 1:length(asocial.models)){
    sum <- sum(asocial.models.matrix[k,which(colnames(constraintsVect)!="asoc")])
    if(sum==0){
      asocial.retain[k] <- k
    }
  }


  asocial.retain <- as.vector(na.omit(asocial.retain))

  asocial.models.to.retain <- asocial.models.matrix[asocial.retain, ]
  asocial.models.to.retain.matrix <- as.matrix(asocial.models.to.retain)
  constraintsVectMatrix <- rbind(social.models.matrix,asocial.models.to.retain)

  # add the Null model (without social learning, and no ILVs)
  constraintsVectMatrix <- rbind(constraintsVectMatrix, rep(0, length(constraintsVectMatrix[1,])))

  row.names(constraintsVectMatrix) <- NULL
  return(constraintsVectMatrix)
}


# 2.6. Run cTADA - DIAL ---------------------------------------------------
# we first create the constraints vector matrix for one of the three NBDA objects
# (their structure is the same, 1 network + 2 ILVs)
constraintsVectMatrix <- create.constraints.Vect.Matrix(MA_NBDA_DATA$nbdadata, num_networks=1, num_ILVs=2)
ILVs <- c("age", "sex")
colnames(constraintsVectMatrix) <- c("network", paste("asoc", ILVs , sep="_"), paste("soc", ILVs, sep="_"))
constraintsVectMatrix <- as.matrix(constraintsVectMatrix)
class(constraintsVectMatrix)
head(constraintsVectMatrix)
# Each row corresponds to one model - each column specifies the parameters that are fit in the respective model
# if a parameter is set to 0, it is not estimated in that model. If it is >0, it is included in the model
# If two parameters are set to the same number (e.g 2+2), their effect is constrained to be the same
# we do not have this here - all parameters are independently estimated
# column 1: refers to whether the model is a social model (if 1), or an asocial model (if 0)
# column 2: refers to age (which is the first ILV specified in the nbda data objects) influencing the asocial learning rate (!=0) or not (0)
# column 3: refers to sex (2nd ILV specified) influencing the asocial learning rate (!=0) or not (=0)
# column 4: age influnecing the social learning rate
# column 5: sex influencing the social learning rate

# We can now run cTADA with the function AICtabel
AICtable_dial <- tadaAICtable(nbdadata = list(
  BB_NBDA_DATA$nbdadata,
  GW_NBDA_DATA$nbdadata,
  MA_NBDA_DATA$nbdadata
), constraintsVectMatrix = constraintsVectMatrix)

setwd("C:/Users/swild/Desktop/Konstanz/Collaborations/Lucy, Michael - Phil Trans cumulative culture/Analysis/NBDA output/Results DIAL diffusion")

# save TADA object
# save(AICtable_dial, file="AICTableDial.RData")
# load("AICTableDial.RData")

# printing the resulting AIC table with all the models (N=20)
AICtable_dial@printTable

# extracting the network support
network.support.dial <- networksSupport(AICtable_dial)
network.support.dial
write.table(network.support.dial, "network.support.dial.diffusion.txt")

# we can see that the models containing the social network (1) receive much more support than those without social learning (0)
#        support numberOfModels
# 0 1.703305e-16              4
# 1 1.000000e+00             16

# we could look at the support for the different model types in more detail
# here not necessary, as the social models get overwhelming support
# support by type
typeSupport(AICtable_dial)

# extracting variable support
# an ILV is considered important if it has a summed Akaike weight of >0.5
# (more likely to be part of a model than not)
# which is this case are age and sex influencing the social learning rate, but not the asocial learning rate
variable.support.dial <- variableSupport(AICtable_dial)
variable.support.dial
write.table(variable.support.dial, "variable.support.dial.diffusion.txt")
#         s1 ASOC:age_BB_dial ASOC:sex_BB_dial SOCIAL:age_BB_dial SOCIAL:sex_BB_dial
# support  1        0.3050232         0.256732           0.670655          0.8926658

# extract model averaged estimates
MLE.dial  <- modelAverageEstimates(AICtable_dial,averageType = "median")
MLE.dial
write.table(MLE.dial, "MLE.dial.diffusion.txt")


# 2.7. Extract effect sizes -----------------------------------------------
## here, we extract estimates for the social learning parameter s conditional on the best performing (social) model
# we first constrain the NBDA data objects to contain the parameters from the best performing model (see AIC table model 6)
constrained.dial.BB <- constrainedNBDAdata(BB_NBDA_DATA$nbdadata, constraintsVect = constraintsVectMatrix[6,])
constrained.dial.GW <- constrainedNBDAdata(GW_NBDA_DATA$nbdadata, constraintsVect = constraintsVectMatrix[6,])
constrained.dial.MA <- constrainedNBDAdata(MA_NBDA_DATA$nbdadata, constraintsVect = constraintsVectMatrix[6,])

# run TADA fit on the constrained NBDA data
bestModel.dial <- tadaFit(nbdadata = list(
  constrained.dial.BB,
  constrained.dial.GW,
  constrained.dial.MA
), type="social")

cbind.data.frame(bestModel.dial@outputPar, bestModel.dial@varNames)
# the numbers in front of the parameter names (1, 2, 3) are used below for the profile likelihood intervals
bestModel.dial@optimisation
bestModel.dial@se # standard errors are misleading as the profile likelihood intervals are asymmetric
# we often have more information on the lower bound of the parameters than the upper bounds


# therefore, we plot profile likelihood for parameter s (which=1 refers to parameter s)
plotProfLik(which=1, model=bestModel.dial,range=c(0,300), resolution = 10)

# adjust the ranges for upper and lower interval (where the lines cross over)
CIs <- profLikCI(which=1,model=bestModel.dial,upperRange =c(200, 300), lowerRange=c(10,50))
CIs
# Lower CI  Upper CI
# 15.62725 216.48690

# as s is difficult to interpret (as they are relative to a baseline rate of asocial learning),
# we can instead extract how many birds have learned the dial task through social learning
# first for each acquisition event
prop.solve.social.byevent.dial <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.dial.BB,
  constrained.dial.GW,
  constrained.dial.MA
), model=bestModel.dial)
prop.solve.social.byevent.dial
# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.dial <- oadaPropSolveByST( nbdadata = list(
  constrained.dial.BB,
  constrained.dial.GW,
  constrained.dial.MA
), model=bestModel.dial)
prop.solve.social.dial
# overall, an estimated 85.8% of birds have learned the dial task socially

#To get the estimates for the lower bound we should really find the corresponding value of the other parameters to plug in when s1 is constrained to this value
bestModelDataS1LowerBound.dial.BB <- constrainedNBDAdata(
  nbdadata =
    BB_NBDA_DATA$nbdadata,
  constraintsVect = constraintsVectMatrix[6, ],
  offset = c(CIs[1] , rep(0, 4))
)

bestModelDataS1LowerBound.dial.GW <- constrainedNBDAdata(
  nbdadata =
    GW_NBDA_DATA$nbdadata,
  constraintsVect = constraintsVectMatrix[6, ],
  offset = c(CIs[1] , rep(0, 4))
)

bestModelDataS1LowerBound.dial.MA <- constrainedNBDAdata(
  nbdadata =
    MA_NBDA_DATA$nbdadata,
  constraintsVect = constraintsVectMatrix[6, ],
  offset = c(CIs[1] , rep(0, 4))
)

#Now, when we fit an "asocial" model it constrains the value of s1=0, but then the value of s at the lower bound is added to s as an offset
bestModelS1LowerBound.dial <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.dial.BB,
      bestModelDataS1LowerBound.dial.GW,
      bestModelDataS1LowerBound.dial.MA
    ) ,
    type = "asocial"
  )
bestModelS1LowerBound.dial@outputPar
# [1] 281.7277956   0.8630295   1.0775270
#Now plug into the prop solve function
prop.solve.social.lower.dial <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound.dial,
    nbdadata = list(
      bestModelDataS1LowerBound.dial.BB,
      bestModelDataS1LowerBound.dial.GW,
      bestModelDataS1LowerBound.dial.MA
    )
  )
prop.solve.social.lower.dial
# lower bound for % of birds having learned the dial task through social learning is 77.2%

# We repeat the same procedure for the upper limit
bestModelDataS1UpperBound.dial.BB <- constrainedNBDAdata(
  nbdadata =
    BB_NBDA_DATA$nbdadata,
  constraintsVect = constraintsVectMatrix[6, ],
  offset = c(CIs[2]  , rep(0, 4))
)

bestModelDataS1UpperBound.dial.GW <- constrainedNBDAdata(
  nbdadata =
    GW_NBDA_DATA$nbdadata,
  constraintsVect = constraintsVectMatrix[6, ],
  offset = c(CIs[2]  , rep(0, 4))
)

bestModelDataS1UpperBound.dial.MA <- constrainedNBDAdata(
  nbdadata =
    MA_NBDA_DATA$nbdadata,
  constraintsVect = constraintsVectMatrix[6, ],
  offset = c(CIs[2] , rep(0, 4))
)



# We again fit an "asocial" model with the upper bound as an offset
bestModelS1UpperBound.dial <-
  tadaFit(
    list(
      bestModelDataS1UpperBound.dial.BB,
      bestModelDataS1UpperBound.dial.GW,
      bestModelDataS1UpperBound.dial.MA
    ) ,
    type = "asocial"
  )
bestModelS1UpperBound.dial@outputPar
# [1] 3460.6676872    0.6749891    0.8451802
#Now plug into the prop solve function
prop.solve.social.Upper.dial <-
  oadaPropSolveByST(
    model = bestModelS1UpperBound.dial,
    nbdadata = list(
      bestModelDataS1UpperBound.dial.BB,
      bestModelDataS1UpperBound.dial.GW,
      bestModelDataS1UpperBound.dial.MA
    )
  )
prop.solve.social.Upper.dial
# Upper bound for % of birds having learned the dial task through social learning is 92.8%

# we now extract the effect sizes for age and sex influencing social learning as they had a summed Akaike weight of above 0.5

# first for age influencing social learning
# from
bestModel.dial@varNames
# we can see that we have to set which=2 ("2 Social: age_BB_dial")
plotProfLik(which=2, model=bestModel.dial,range = c(-1, 2), resolution = 10)

# adjust the ranges for upper and lower interval (where the lines cross over)
profLikCI(which=2,model=bestModel.dial,upperRange =c(1.1, 1.7), lowerRange=c(-0.5,0.5))
#     Lower CI   Upper CI
# 0.02347124 1.48145026

# back-transform the estimates (note that for the effect size, we take the MLE (median) and for the CI we take profile likelihood intervals from the best model)
exp(c(MLE.dial[4], 0.02347124, 1.48145026))

# SOCIALage_BB_dial
# 1.965231          1.023749          4.399321

# Adults were 2.0 [1.0-4.4] times faster at learning the dial task socially compared to juveniles

# repeat for sex (which=3)
plotProfLik(which=3, model=bestModel.dial,range = c(0, 2), resolution = 10)
profLikCI(which=3,model=bestModel.dial,upperRange =c(1.5, 2), lowerRange=c(0,0.5))

# Lower CI  Upper CI
# 0.230165 1.669527

# back-transform the estimates
exp(c(MLE.dial[5], 0.230165, 1.669527))
# Males were 2.5 [1.3-5.3] times faster at learning the dial task socially compared to females

######################################################################### #

# 3) COMPLEX 1st GEN ------------------------------------------------------
# In this experiment, both the slide and dial components are rewarded, but solving the complex task (slide and dial in combination) results in a bigger reward


# 3.1. Read in file with previous experience ------------------------------
# this file has extracted the numbers of solves for both dial and slide up until the end of the dial experiment
# (so this includes data from a previous experiment - slide - that is not part of this manuscript)
prev_exp <-  NBDA_input_data$data.objects.2016$prev_exp

dial.demos <- subset(prev_exp$RING, prev_exp$DIAL_PAST>=3)
slide.demos <- subset(prev_exp$RING, prev_exp$DOOR_PAST>=3)

# 3.1.1 Summary -----------------------------------------------------------------
num.dial.1st.gen <- length(as.data.frame(subset(df_solves, df_solves$experiment=="complex_prog" & df_solves$solution_category=="dial"))[,1])
num.slide.1st.gen <- length(as.data.frame(subset(df_solves, df_solves$experiment=="complex_prog" & df_solves$solution_category=="slide"))[,1])
num.complex.1st.gen <- length(as.data.frame(subset(df_solves, df_solves$experiment=="complex_prog" & df_solves$solution_category=="complex"))[,1])

# total number of simple solutions
num.dial.1st.gen+num.slide.1st.gen
# [1] 17488

# extract per subpopulation (and experiment)
# length(subset(df_solves, df_solves$experiment=="nextgen" & df_solves$solution_category %in% c("complex") & df_solves$PATCH=="MA")[,1])

# 3.2. Prepare NBDA data objects ------------------------------------------

prepare.NBDA.data.complex.progressive <- function(input.file, locations, label, start.date, end.date, previous_experience, patch, learned.behav, ILVs.include){
  # subset the input file to the correct experiment ("complex_prog")
  input.file.sub <- subset(input.file, input.file$PATCH==patch & input.file$SPECIES=="GRETI" & input.file$experiment=="complex_prog")

  # a bird is considered a learner if it has solved three times
  # the time of its first solve is taken as the time of acquisition
  # prepare empty vectors to store the ID of the learner, the number of solves and the time of the first solve

  learners.list <- NULL
  num_solves_all <- NULL
  time_solves <- NULL
  which.rows <- NULL

  if(learned.behav=="complex"){ # only considering complex solves
    for (i in as.vector(unique(input.file.sub$RING[input.file.sub$solution_category=="complex"]))){
      sub <- subset(input.file.sub, input.file.sub$RING==i & input.file.sub$solution_category=="complex")
      which.row <- rownames(sub)[1]
      num.solves <- length(sub[,1])
      solve.time <- as.numeric(difftime(min(sub$time_stamp), min(input.file.sub$time_stamp), units="days"))
      learners.list[which(unique(input.file.sub$RING)==i)] <- i
      num_solves_all[which(unique(input.file.sub$RING)==i)] <- num.solves
      time_solves[which(unique(input.file.sub$RING)==i)] <- solve.time
      which.rows[which(unique(input.file.sub$RING)==i)] <- which.row
    }
  } else { # if we are looking at the dial or slide diffusion
    for (i in as.vector(unique(input.file.sub$RING[input.file.sub$solution_category %in% c(learned.behav, "complex")]))){
      sub <- subset(input.file.sub, input.file.sub$RING==i & input.file.sub$solution_category %in% c(learned.behav, "complex"))
      which.row <- rownames(sub)[1]
      num.solves <- length(sub[,1])
      solve.time <- as.numeric(difftime(min(sub$time_stamp), min(input.file.sub$time_stamp), units="days"))
      learners.list[which(unique(input.file.sub$RING)==i)] <- i
      num_solves_all[which(unique(input.file.sub$RING)==i)] <- num.solves
      time_solves[which(unique(input.file.sub$RING)==i)] <- solve.time
      which.rows[which(unique(input.file.sub$RING)==i)] <- which.row

    }
  }



  # reassign column names
  learners.list_df <- na.omit(cbind.data.frame(learners.list, num_solves_all,  time_solves, which.rows))
  colnames(learners.list_df) <- c("RING", "num_solves", "time_acq", "which_rows")
  learners.list_df <- learners.list_df[order(learners.list_df$time_acq),]


  # remove those with fewer than 3 solves
  learners.list_df <- subset(learners.list_df, learners.list_df$num_solves>=3)

  # restrict the data frame to great tits only that are in the ILV file
  learners.list_df <- subset(learners.list_df, learners.list_df$RING%in%ids.mod$id)

  # demonstrators - we only have demos for the dial and slide diffusion
  if(learned.behav=="slide"){
    demos.sub <- subset(prev_exp$RING, prev_exp$DOOR_PAST>=3)
  } else if(learned.behav=="dial"){
    demos.sub <- subset(prev_exp$RING, prev_exp$DIAL_PAST>=3)
  } else {demos.sub <- NA} # if complex, then there are no demonstrators


  # create network from group by individual matrix

  # subset the gbi to the corresponding locations and times for each diffusion
  which.groups <- subset(group_info, group_info$location%in% locations & group_info$start_time>=start.date & group_info$stop_time<=end.date)
  which.groups.num <- rownames(which.groups)

  sub.gbi <- gbi[which.groups.num,] # subset the gbi to the right locations and times
  sub.gbi <- sub.gbi[, colSums(sub.gbi)>9] # remove birds that weren't seen at all during this period

  # remove birds from learners.list_df that are not part of the network
  learners.list_df <- subset(learners.list_df , learners.list_df$RING%in% colnames(sub.gbi))

  # adjust the times of acquisition
  # recalculate the time of acquisition by removing the days where the puzzle boxes weren't out (network days)
  # days 0-4 puzzle box
  # days 5+6 network data
  # days 7-11 puzzle box (so subtract two days from the acquisition of birds that learned between day 8-12)
  # days 12+13 network data
  # day 14-18 (subtract 4 days from the acquisition time for birds that have learned between day 15-19)
  # days 19+20 network data
  # day 21-25 (subtract 6)

  learners.list_df$time_acq_adj <- NA

  for (i in 1:length(learners.list_df[,1])){
    time <- learners.list_df[i, "time_acq"]
    if(time > 7 & time < 12){
      time.adj <- time-2
    } else if(time > 14 & time < 19){
      time.adj <- time-4
    } else if(time > 21 & time < 26){
      time.adj <- time-6
    } else {time.adj <- time}
    learners.list_df[i, "time_acq_adj"] <- time.adj
  }

  # generate association network
  net <- get_network(association_data = sub.gbi, data_format="GBI", association_index = "SRI")
  net <- net[order(rownames(net)), order(colnames(net))]


  ## the association matrix needs to be in an array
  assocMatrix <- array(data = net, dim=c(nrow(net), ncol(net), 1))
  class(assocMatrix)

  rownames(assocMatrix) <- rownames(net)
  colnames(assocMatrix) <- rownames(net)
  suppressWarnings(
    if(is.na(demos.sub)){
      demonstrators <- rep(0, nrow(net))
    } else {
      demonstrators <- rep(0, nrow(net))
      demonstrators[which(rownames(net) %in% demos.sub )] <- 1
    }
  )



  ## get ILVs ready
  ids.sub <- subset(ids.mod, ids.mod$id %in% rownames(net))
  learners.list_df <- subset(learners.list_df, learners.list_df$RING%in%rownames(net))

  # as we are coding the previous experience as a time-varying variable that updates with every acquisition event,
  # we also need to create matrices for sex and age, where each column represents one acquisition event

  sex <- matrix(as.numeric(ids.sub$sex), nrow=length(ids.sub$sex), ncol=length(learners.list_df$which_rows), byrow=F)
  age <- matrix(as.numeric(ids.sub$age), nrow=length(ids.sub$age), ncol=length(learners.list_df$which_rows), byrow=F)

  # we now create two matrices that encode for previous experience:
  # previous_experience_one is a 1/0 vector that contains a 1 if a bird has experience with one of the two components, 0 if no experience
  # previous_experience_both is a 1/0 vector that contains a 1 if a birds has experience with both components (but not in combination), and 0 if no experience or experience with just one of the two components

  prev.experience.matrix.one <- as.data.frame(matrix(NA, ncol=length(learners.list_df$which_rows), nrow=length(rownames(net))))
  prev.experience.matrix.both <- as.data.frame(matrix(NA, ncol=length(learners.list_df$which_rows), nrow=length(rownames(net))))

  # at each acquisition event, extract each birds experience level
  for (k in as.numeric(learners.list_df$which_rows)){
    prev.exp.one <- NULL
    prev.exp.both <- NULL

    for (i in rownames(net)){

      sub.prev.exp <- subset(prev_exp, prev_exp$RING==i & prev_exp$SITE==patch)
      exp.at.start.slide <- sum(sub.prev.exp[,"DOOR_PAST"])
      exp.at.start.dial <- sum(sub.prev.exp[,"DIAL_PAST"])

      # subset the input file to dial and slide only
      input.file.sub.slide.dial <- subset(input.file.sub, input.file.sub$solution_category%in% c("slide", "dial"))

      # exract the closest event to the complex solve
      suppressWarnings(
        index.of.closest <- which(abs(as.numeric(rownames(input.file.sub.slide.dial))-k)==min(abs(as.numeric(rownames(input.file.sub.slide.dial))-k)))
      )


      input.file.sub.slide.dial.sub <- input.file.sub.slide.dial[1:index.of.closest,] # subset the input file up to the point of where a bird has learned
      sub.input.i.dial <- subset(input.file.sub.slide.dial.sub, input.file.sub.slide.dial.sub$RING==i & input.file.sub.slide.dial.sub$solution_category=="dial") # subset the data frame only for individual i
      sub.input.i.slide <- subset(input.file.sub.slide.dial.sub, input.file.sub.slide.dial.sub$RING==i & input.file.sub.slide.dial.sub$solution_category=="slid") # subset the data frame only for individual i

      # extract the number of solves up until the point of the acquisition event
      max.i.dial <- length(sub.input.i.dial[,1])
      max.i.slide <- length(sub.input.i.slide[,1])


      if(max.i.dial>0){ # if the birds has gained any experience during the current experiment in solving dial
        exp.at.acquisition.dial <- exp.at.start.dial + max.i.dial # build the sum of the experience at the start and the experience gathered during the experiment
      } else {
        exp.at.acquisition.dial <- exp.at.start.dial # otherwise we just take the experience the bird had at the start of the experiment
      }

      if(max.i.slide>0){ # if the birds has gained any experience during the current experiment in solving slide
        exp.at.acquisition.slide <- exp.at.start.slide + max.i.slide # build the sum of the experience at the start and the experience gathered during the experiment
      } else {
        exp.at.acquisition.slide <- exp.at.start.slide # otherwise we just take the experience the bird had at the start of the experiment
      }

      # next we build two vectors:
      # one that gets a 1 if birds have experience (>=3 solves) with both components
      # one that gets a 1 if birds have experience (>= solves) with one component but not the other
      if(exp.at.acquisition.slide>=3 & exp.at.acquisition.dial>=3){
        prev.exp.both[which(rownames(net)==i)] <- 1
        prev.exp.one[which(rownames(net)==i)] <- 0
      } else if (exp.at.acquisition.slide>=3 & exp.at.acquisition.dial<3 | exp.at.acquisition.slide<3 & exp.at.acquisition.dial>=3){
        prev.exp.one[which(rownames(net)==i)] <- 1
        prev.exp.both[which(rownames(net)==i)] <- 0
      } else {
        prev.exp.one[which(rownames(net)==i)] <- 0
        prev.exp.both[which(rownames(net)==i)] <- 0
      }

    }
    prev.experience.matrix.one[,which(learners.list_df$which_rows==k)] <- prev.exp.one
    prev.experience.matrix.both[,which(learners.list_df$which_rows==k)] <- prev.exp.both

  }

  # ensure they are in matrix format
  prev.experience.matrix.one <- as.matrix(prev.experience.matrix.one)
  prev.experience.matrix.both <- as.matrix(prev.experience.matrix.both)

  # create the ILV vectors
  assign(paste("age", label, sep="_"), age, envir = .GlobalEnv)
  assign(paste("sex", label, sep="_"), sex, envir = .GlobalEnv)
  assign(paste("pre.exp.one", label, sep="_"), prev.experience.matrix.one, envir = .GlobalEnv)
  assign(paste("pre.exp.both", label, sep="_"), prev.experience.matrix.both, envir = .GlobalEnv)

  ILVs <- paste(ILVs.include, label, sep="_")

  assign(paste("ILVs", label, sep="_"), ILVs)

  # get order of acquisition for learners.list
  learners.list_OAc <- NULL


  for (i in unique(learners.list_df$RING)){
    pos <- which(colnames(net)==i)
    learners.list_OAc <- c(learners.list_OAc, pos)
  }


  # create NBDA Data Object


  nbdaData <- nbdaData(label=label,
                       assMatrix = assocMatrix,
                       asoc_ilv = get(paste("ILVs", label, sep="_")),
                       int_ilv = get(paste("ILVs", label, sep="_")),
                       multi_ilv = "ILVabsent",
                       orderAcq = learners.list_OAc,
                       timeAcq = learners.list_df$time_acq_adj,# get time in days
                       asocialTreatment = "timevarying",
                       demons = demonstrators)

  # note that there are no demonstrators in this experiment, as no birds have ever solved the complex task

  object <- NULL

  object$nbdadata <- nbdaData
  object$label <- label
  object$learners.list <- learners.list_df
  object$num.birds <- length(rownames(net))
  object$num.learners.list <- length(learners.list_df$RING)
  object$IDs <- rownames(net)

  return(object)

}

###
# note that warnings can be ignored when creating the nbda data objects
# the warnings can occur when extracting the experience at each acquisition event


# 3.2.1 Create NBDA Data Objects for the dial diffusion -------------------

# Dial diffusion (1st gen complex experiment)
# birds that have previous experience with dial are demonstrators
BB_complex_progressive_dial <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                patch = "BB",
                                                                locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
                                                                start.date = as.POSIXct("06/02/2016", format="%d/%m/%Y"),
                                                                end.date = as.POSIXct("28/02/2016", format="%d/%m/%Y"),
                                                                label = "BB_complex_progressive",
                                                                learned.behav = "dial",
                                                                ILVs.include = c("sex", "age")
)

CP_complex_progressive_dial <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                patch = "CP",
                                                                locations <- c("1b", "1c", "1d", "2e", "2f", "1e", "1f", "1a"),
                                                                start.date = as.POSIXct("23/01/2016", format="%d/%m/%Y"),
                                                                end.date = as.POSIXct("14/02/2016", format="%d/%m/%Y"),
                                                                label = "CP_complex_progressive",
                                                                learned.behav = "dial",
                                                                ILVs.include = c("sex", "age")
)


GW_complex_progressive_dial <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                patch = "GW",
                                                                locations <- c("2c", "2b", "2d", "2h", "2g"),
                                                                start.date = as.POSIXct("13/02/2016", format="%d/%m/%Y"),
                                                                end.date = as.POSIXct("06/03/2016", format="%d/%m/%Y"),
                                                                label = "GW_complex_progressive",
                                                                learned.behav = "dial",
                                                                ILVs.include = c("sex", "age")
)

MP_complex_progressive_dial <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                     patch = "MP",
                                                                     locations = c("7d", "7e", "7f", "7g", "7h", "8f"),
                                                                     start.date = as.POSIXct("30/01/2016", format="%d/%m/%Y"),
                                                                     end.date = as.POSIXct("21/02/2016", format="%d/%m/%Y"),
                                                                     label = "MP_complex_progressive",
                                                                     learned.behav = "dial",
                                                                     ILVs.include = c("sex", "age")
)


## some subpopulations are close - so we have to check for overlap between BB, GW and CP

# overlap between BB and CP - dial
BB_CP_dial.overlap <- intersect(BB_complex_progressive_dial$learners.list$RING, CP_complex_progressive_dial$learners.list$RING)
# 6 birds
# where have they learned first?
cbind.data.frame(BB_complex_progressive_dial$learners.list[BB_complex_progressive_dial$learners.list$RING%in%BB_CP_dial.overlap,], CP_complex_progressive_dial$learners.list[CP_complex_progressive_dial$learners.list$RING%in%BB_CP_dial.overlap,])
# looking at the times of acquisition, they have all learned at BB first and have then solved at CP - meaning we have to filter them out as learners at CP

# extract the position of these three birds in the network of CP
which(CP_complex_progressive_dial$IDs %in% BB_CP_dial.overlap)
CP_complex_progressive_dial$nbdadata@orderAcq

# filter the NBDA data object
CP_complex_progressive_dial$nbdadata <-
  filteredNBDAdata(
    CP_complex_progressive_dial$nbdadata,
    filter = "id",
    exclude = paste(
      "CP_complex_progressive",
      which(CP_complex_progressive_dial$IDs %in% BB_CP_dial.overlap),
      sep = "_"
    )
  )

# also subtract the three birds from the learner data frame and the number of learners.list
CP_complex_progressive_dial$learners.list <- subset(CP_complex_progressive_dial$learners.list, !(CP_complex_progressive_dial$learners.list$RING %in% BB_CP_dial.overlap))
CP_complex_progressive_dial$num.learners.list <- CP_complex_progressive_dial$num.learners.list - length(BB_CP_dial.overlap)



# overlap between BB and GW - dial
BB_GW_dial.overlap <- intersect(BB_complex_progressive_dial$learners.list$RING, GW_complex_progressive_dial$learners.list$RING)
# 1 bird
# where has it learned first?
cbind.data.frame(BB_complex_progressive_dial$learners.list[BB_complex_progressive_dial$learners.list$RING%in%BB_GW_dial.overlap,], GW_complex_progressive_dial$learners.list[GW_complex_progressive_dial$learners.list$RING%in%BB_GW_dial.overlap,])
# learned first in BB, so filtering out at GW
# extract the position of these three birds in the network of GW
which(GW_complex_progressive_dial$IDs %in% BB_GW_dial.overlap)
GW_complex_progressive_dial$nbdadata@orderAcq

# filter the NBDA data object
GW_complex_progressive_dial$nbdadata <-
  filteredNBDAdata(
    GW_complex_progressive_dial$nbdadata,
    filter = "id",
    exclude = paste(
      "GW_complex_progressive",
      which(GW_complex_progressive_dial$IDs %in% BB_GW_dial.overlap),
      sep = "_"
    )
  )

# also subtract the three birds from the learner data frame and the number of learners.list
GW_complex_progressive_dial$learners.list <- subset(GW_complex_progressive_dial$learners.list, !(GW_complex_progressive_dial$learners.list$RING %in% BB_GW_dial.overlap))
GW_complex_progressive_dial$num.learners.list <- GW_complex_progressive_dial$num.learners.list - length(BB_GW_dial.overlap)


# overlap between GW and CP - dial
GW_CP_dial.overlap <- intersect(GW_complex_progressive_dial$learners.list$RING, CP_complex_progressive_dial$learners.list$RING)
# 3 birds
# where have they learned first?
cbind.data.frame(GW_complex_progressive_dial$learners.list[GW_complex_progressive_dial$learners.list$RING%in%GW_CP_dial.overlap,], CP_complex_progressive_dial$learners.list[CP_complex_progressive_dial$learners.list$RING%in%GW_CP_dial.overlap,])
# looking at the times of acquisition, they have all learned at GW first and have then solved at CP - meaning we have to filter them out as learners at CP

# extract the position of these three birds in the network of CP
which(CP_complex_progressive_dial$IDs %in% GW_CP_dial.overlap)
CP_complex_progressive_dial$nbdadata@orderAcq

# filter the NBDA data object
CP_complex_progressive_dial$nbdadata <-
  filteredNBDAdata(
    CP_complex_progressive_dial$nbdadata,
    filter = "id",
    exclude = paste(
      "CP_complex_progressive",
      which(CP_complex_progressive_dial$IDs %in% GW_CP_dial.overlap),
      sep = "_"
    )
  )

# also subtract the three birds from the learner data frame and the number of learners.list
CP_complex_progressive_dial$learners.list <- subset(CP_complex_progressive_dial$learners.list, !(CP_complex_progressive_dial$learners.list$RING %in% GW_CP_dial.overlap))
CP_complex_progressive_dial$num.learners.list <- CP_complex_progressive_dial$num.learners.list - length(GW_CP_dial.overlap)


# 3.2.2 Create NBDA Data Objects for Slide --------------------------------


# slide diffusion (1st gen complex experiment)
# birds that have previous experience with slide are demonstrators
BB_complex_progressive_slide <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                     patch = "BB",
                                                                     locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
                                                                     start.date = as.POSIXct("06/02/2016", format="%d/%m/%Y"),
                                                                     end.date = as.POSIXct("28/02/2016", format="%d/%m/%Y"),
                                                                     label = "BB_complex_progressive",
                                                                     learned.behav = "slide",
                                                                     ILVs.include = c("sex", "age")
)

CP_complex_progressive_slide <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                     patch = "CP",
                                                                     locations <- c("1b", "1c", "1d", "2e", "2f", "1e", "1f", "1a"),
                                                                     start.date = as.POSIXct("23/01/2016", format="%d/%m/%Y"),
                                                                     end.date = as.POSIXct("14/02/2016", format="%d/%m/%Y"),
                                                                     label = "CP_complex_progressive",
                                                                     learned.behav = "slide",
                                                                     ILVs.include = c("sex", "age")
)


GW_complex_progressive_slide <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                     patch = "GW",
                                                                     locations <- c("2c", "2b", "2d", "2h", "2g"),
                                                                     start.date = as.POSIXct("13/02/2016", format="%d/%m/%Y"),
                                                                     end.date = as.POSIXct("06/03/2016", format="%d/%m/%Y"),
                                                                     label = "GW_complex_progressive",
                                                                     learned.behav = "slide",
                                                                     ILVs.include = c("sex", "age")
)


MP_complex_progressive_slide <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                     patch = "MP",
                                                                     locations = c("7d", "7e", "7f", "7g", "7h", "8f"),
                                                                     start.date = as.POSIXct("30/01/2016", format="%d/%m/%Y"),
                                                                     end.date = as.POSIXct("21/02/2016", format="%d/%m/%Y"),
                                                                     label = "MP_complex_progressive",
                                                                     learned.behav = "slide",
                                                                     ILVs.include = c("sex", "age")
)

# again, we need to check for overlapping learners between subpopulations BB, CP and GW

# overlap between BB and CP - slide
BB_CP_slide.overlap <- intersect(BB_complex_progressive_slide$learners.list$RING, CP_complex_progressive_slide$learners.list$RING)
# 5 birds
# where have they learned first?
cbind.data.frame(BB_complex_progressive_slide$learners.list[BB_complex_progressive_slide$learners.list$RING%in%BB_CP_slide.overlap,], CP_complex_progressive_slide$learners.list[CP_complex_progressive_slide$learners.list$RING%in%BB_CP_slide.overlap,])
# looking at the times of acquisition, the first three birds on the list have learned first at CP, the last two first at BB
CP.first.slide <- BB_CP_slide.overlap[1:3]
BB.first.slide <- BB_CP_slide.overlap[4:5]

# we first filter out the first three birds at BB
BB_complex_progressive_slide$nbdadata <-
  filteredNBDAdata(
    BB_complex_progressive_slide$nbdadata,
    filter = "id",
    exclude = paste(
      "BB_complex_progressive",
      which(BB_complex_progressive_slide$IDs %in% CP.first.slide),
      sep = "_"
    )
  )

# also subtract the three birds from the learner data frame and the number of learners.list
BB_complex_progressive_slide$learners.list <- subset(BB_complex_progressive_slide$learners.list, !(BB_complex_progressive_slide$learners.list$RING %in% CP.first.slide))
BB_complex_progressive_slide$num.learners.list <- BB_complex_progressive_slide$num.learners.list - length(CP.first.slide)

# now we filter the last two birds at CP
CP_complex_progressive_slide$nbdadata <-
  filteredNBDAdata(
    CP_complex_progressive_slide$nbdadata,
    filter = "id",
    exclude = paste(
      "CP_complex_progressive",
      which(CP_complex_progressive_slide$IDs %in% BB.first.slide),
      sep = "_"
    )
  )

# also subtract the birds from the learner data frame and the number of learners.list
CP_complex_progressive_slide$learners.list <- subset(CP_complex_progressive_slide$learners.list, !(CP_complex_progressive_slide$learners.list$RING %in% BB.first.slide))
CP_complex_progressive_slide$num.learners.list <- CP_complex_progressive_slide$num.learners.list - length(BB.first.slide)


# overlap between BB and GW - slide
BB_GW_slide.overlap <- intersect(BB_complex_progressive_slide$learners.list$RING, GW_complex_progressive_slide$learners.list$RING)
# no overlap for slide

# overlap between GW and CP - slide
GW_CP_slide.overlap <- intersect(GW_complex_progressive_slide$learners.list$RING, CP_complex_progressive_slide$learners.list$RING)
# no overlap between GW and CP for slide



# extract how many birds have learned dial out of the total number of different birds

# total number of birds
length(unique(c(BB_complex_progressive_dial$IDs, CP_complex_progressive_dial$IDs, GW_complex_progressive_dial$IDs, MP_complex_progressive_dial$IDs)))
# [1] 344




# 3.2.3 Create NBDA data objects for the complex diffusion ----------------

# complex diffusion (1st gen complex experiment)
# no demonstrators in this one
BB_complex_progressive_complex <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                      patch = "BB",
                                                                      locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
                                                                      start.date = as.POSIXct("06/02/2016", format="%d/%m/%Y"),
                                                                      end.date = as.POSIXct("28/02/2016", format="%d/%m/%Y"),
                                                                      label = "BB_complex_progressive",
                                                                      learned.behav = "complex",
                                                                      ILVs.include = c("sex", "age", "pre.exp.one", "pre.exp.both")
)

CP_complex_progressive_complex <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                      patch = "CP",
                                                                      locations <- c("1b", "1c", "1d", "2e", "2f", "1e", "1f", "1a"),
                                                                      start.date = as.POSIXct("23/01/2016", format="%d/%m/%Y"),
                                                                      end.date = as.POSIXct("14/02/2016", format="%d/%m/%Y"),
                                                                      label = "CP_complex_progressive",
                                                                      learned.behav = "complex",
                                                                      ILVs.include = c("sex", "age", "pre.exp.one", "pre.exp.both")
)


GW_complex_progressive_complex <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                      patch = "GW",
                                                                      locations <- c("2c", "2b", "2d", "2h", "2g"),
                                                                      start.date = as.POSIXct("13/02/2016", format="%d/%m/%Y"),
                                                                      end.date = as.POSIXct("06/03/2016", format="%d/%m/%Y"),
                                                                      label = "GW_complex_progressive",
                                                                      learned.behav = "complex",
                                                                      ILVs.include = c("sex", "age", "pre.exp.one", "pre.exp.both")
)

MP_complex_progressive_complex <- prepare.NBDA.data.complex.progressive(input.file = as.data.frame(df_solves),
                                                                      patch = "MP",
                                                                      locations = c("7d", "7e", "7f", "7g", "7h", "8f"),
                                                                      start.date = as.POSIXct("30/01/2016", format="%d/%m/%Y"),
                                                                      end.date = as.POSIXct("21/02/2016", format="%d/%m/%Y"),
                                                                      label = "MP_complex_progressive",
                                                                      learned.behav = "complex",
                                                                      ILVs.include = c("sex", "age", "pre.exp.one", "pre.exp.both")
)

# note that MA does not have enough learners.list to be included
# but we include MP (the control), as some birds have moved from other sites and started solving - improves model fitting

# check for overlaps between sites (BB, GW and CP)


# overlap between BB and CP - complex
BB_CP_complex.overlap <- intersect(BB_complex_progressive_complex$learners.list$RING, CP_complex_progressive_complex$learners.list$RING)
# 3 birds
# where have they learned first?
cbind.data.frame(BB_complex_progressive_complex$learners.list[BB_complex_progressive_complex$learners.list$RING%in%BB_CP_complex.overlap,], CP_complex_progressive_complex$learners.list[CP_complex_progressive_complex$learners.list$RING%in%BB_CP_complex.overlap,])
# looking at the times of acquisition, the first two have learned first at CP, the third first at BB

BB.first_complex <- BB_CP_complex.overlap[3]
CP.first_complex <- BB_CP_complex.overlap[1:2]

# filter the first two birds out at BB
BB_complex_progressive_complex$nbdadata <-
  filteredNBDAdata(
    BB_complex_progressive_complex$nbdadata,
    filter = "id",
    exclude = paste(
      "BB_complex_progressive",
      which(BB_complex_progressive_complex$IDs %in% CP.first_complex),
      sep = "_"
    )
  )

# also subtract the three birds from the learner data frame and the number of learners.list
BB_complex_progressive_complex$learners.list <- subset(BB_complex_progressive_complex$learners.list, !(BB_complex_progressive_complex$learners.list$RING %in% CP.first_complex))
BB_complex_progressive_complex$num.learners.list <- BB_complex_progressive_complex$num.learners.list - length(CP.first_complex)

# filter the third bird out at CP
CP_complex_progressive_complex$nbdadata <-
  filteredNBDAdata(
    CP_complex_progressive_complex$nbdadata,
    filter = "id",
    exclude = paste(
      "CP_complex_progressive",
      which(CP_complex_progressive_complex$IDs %in% BB.first_complex),
      sep = "_"
    )
  )

# also subtract the  bird from the learner data frame and the number of learners.list
CP_complex_progressive_complex$learners.list <- subset(CP_complex_progressive_complex$learners.list, !(CP_complex_progressive_complex$learners.list$RING %in% BB.first_complex))
CP_complex_progressive_complex$num.learners.list <- CP_complex_progressive_complex$num.learners.list - length(BB.first_complex)

# overlap between BB and GW - complex
BB_GW_complex.overlap <- intersect(BB_complex_progressive_complex$learners.list$RING, GW_complex_progressive_complex$learners.list$RING)

# no overlap between BB and GW for complex

# overlap between GW and CP - complex
GW_CP_complex.overlap <- intersect(GW_complex_progressive_complex$learners.list$RING, CP_complex_progressive_complex$learners.list$RING)
# no overlap between GW and CP for complex


# extract how many are producing simple solutions in each subpopulation
## BB
# extract which have solved dial or slide and are not complex solvers
length(which(!(unique(c(BB_complex_progressive_dial$learners.list$RING, BB_complex_progressive_slide$learners.list$RING)) %in% BB_complex_progressive_complex$learners.list$RING)))/length(BB_complex_progressive_dial$IDs)

## CP
length(which(!(unique(c(CP_complex_progressive_dial$learners.list$RING, CP_complex_progressive_slide$learners.list$RING)) %in% CP_complex_progressive_complex$learners.list$RING)))/length(CP_complex_progressive_dial$IDs)

## GW
length(which(!(unique(c(GW_complex_progressive_dial$learners.list$RING, GW_complex_progressive_slide$learners.list$RING)) %in% GW_complex_progressive_complex$learners.list$RING)))/length(GW_complex_progressive_dial$IDs)

## MP
length(which(!(unique(c(MP_complex_progressive_dial$learners.list$RING, MP_complex_progressive_slide$learners.list$RING)) %in% MP_complex_progressive_complex$learners.list$RING)))/length(MP_complex_progressive_dial$IDs)

# number of demonstrators with knowledge of dial at the start of the experiment
length(which(unique(c(BB_complex_progressive_dial$learners.list$RING, CP_complex_progressive_dial$learners.list$RING, GW_complex_progressive_dial$learners.list$RING, MP_complex_progressive_dial$learners.list$RING)) %in% dial.demos))
# number of demonstrators with knowledge of slide at the start of the experiment
length(which(unique(c(BB_complex_progressive_slide$learners.list$RING, CP_complex_progressive_slide$learners.list$RING, GW_complex_progressive_slide$learners.list$RING, MP_complex_progressive_slide$learners.list$RING)) %in% slide.demos))



################ #
# extract the total number of learners.list and total number of birds
length(unique(c(BB_complex_progressive_dial$IDs, CP_complex_progressive_dial$IDs, GW_complex_progressive_dial$IDs, MP_complex_progressive_dial$IDs)))
# [1] 344
BB_complex_progressive_dial$num.learners.list + CP_complex_progressive_dial$num.learners.list + GW_complex_progressive_dial$num.learners.list+MP_complex_progressive_dial$num.learners.list
# DIAL: 111 birds
BB_complex_progressive_slide$num.learners.list + CP_complex_progressive_slide$num.learners.list + GW_complex_progressive_slide$num.learners.list+MP_complex_progressive_slide$num.learners.list
# SLIDE: 76
BB_complex_progressive_complex$num.learners.list + CP_complex_progressive_complex$num.learners.list + GW_complex_progressive_complex$num.learners.list+ MP_complex_progressive_complex$num.learners.list
# [1] 42

#
dial.learners <- unique(c(BB_complex_progressive_dial$learners.list$RING,CP_complex_progressive_dial$learners.list$RING, GW_complex_progressive_dial$learners.list$RING, MP_complex_progressive_dial$learners.list$RING))
slide.learners <- unique(c(BB_complex_progressive_slide$learners.list$RING,CP_complex_progressive_slide$learners.list$RING, GW_complex_progressive_slide$learners.list$RING, MP_complex_progressive_slide$learners.list$RING))
complex.learners <- unique(c(BB_complex_progressive_complex$learners.list$RING,CP_complex_progressive_complex$learners.list$RING, GW_complex_progressive_complex$learners.list$RING, MP_complex_progressive_complex$learners.list$RING))

# extract those that learn simple but not complex
length(unique(which(!(c(dial.learners, slide.learners) %in% complex.learners))))
103/344

length(complex.learners)
42/344

# 3.3. Create constraintsVectMatrix - simple solutions ---------------------------------------


# First on the simple solutions
constraintsVectMatrix <- create.constraints.Vect.Matrix(CP_complex_progressive_dial$nbdadata, num_networks=1, num_ILVs=2)
ILVs <- c("sex", "age")
colnames(constraintsVectMatrix) <- c("network", paste("asoc", ILVs , sep="_"), paste("soc", ILVs, sep="_"))
constraintsVectMatrix <- as.matrix(constraintsVectMatrix)
class(constraintsVectMatrix)
head(constraintsVectMatrix)



# 3.4. Dial - run cTADA ---------------------------------------------------


# DIAL
AICtable_complex_progressive_dial <- tadaAICtable(nbdadata = list(
  BB_complex_progressive_dial$nbdadata,  # this object is filtered
  CP_complex_progressive_dial$nbdadata,
  GW_complex_progressive_dial$nbdadata,
    MP_complex_progressive_dial$nbdadata # included despite control
), constraintsVectMatrix = constraintsVectMatrix)

# inspect the AIC table
AICtable_complex_progressive_dial@printTable
networksSupport(AICtable_complex_progressive_dial)

# support numberOfModels
# 0 3.14715e-14              4
# 1 1.00000e+00             16

variable.support.complex.progressive_dial <- variableSupport(AICtable_complex_progressive_dial)
variable.support.complex.progressive_dial

# s1 ASOC:sex_BB_complex_progressive ASOC:age_BB_complex_progressive SOCIAL:sex_BB_complex_progressive SOCIAL:age_BB_complex_progressive
# support  1                       0.2418332                        0.623821                          0.996681                         0.2088655

# extract model averaged estimates
MLE.complex.progressive_dial  <- modelAverageEstimates(AICtable_complex_progressive_dial,averageType = "median")
MLE.complex.progressive_dial

# ASOCIALsex_BB_complex_progressive ASOCIALage_BB_complex_progressive  SOCIALsex_BB_complex_progressive  SOCIALage_BB_complex_progressive
# 7595.470447                          0.000000                         13.576694                          1.152288                          0.000000

# 3.5. Slide - run cTADA --------------------------------------------------

# SLIDE
AICtable_complex_progressive_slide <- tadaAICtable(nbdadata = list(
  BB_complex_progressive_slide$nbdadata,
  CP_complex_progressive_slide$nbdadata,
  GW_complex_progressive_slide$nbdadata,
  MP_complex_progressive_slide$nbdadata
), constraintsVectMatrix = constraintsVectMatrix)

# inspect the AIC table
AICtable_complex_progressive_slide@printTable

networksSupport(AICtable_complex_progressive_slide)

# support numberOfModels
# 0 1.259776e-07              4
# 1 9.999999e-01             16


variable.support.complex.progressive_slide <- variableSupport(AICtable_complex_progressive_slide)
variable.support.complex.progressive_slide

# s1 ASOC:sex_BB_complex_progressive ASOC:age_BB_complex_progressive SOCIAL:sex_BB_complex_progressive SOCIAL:age_BB_complex_progressive
# support 0.9999999                       0.4517919                       0.8861038                         0.9129344                         0.2242539


# extract model averaged estimates
MLE.complex.progressive_slide  <- modelAverageEstimates(AICtable_complex_progressive_slide,averageType = "median")
MLE.complex.progressive_slide

# ASOCIALsex_BB_complex_progressive ASOCIALage_BB_complex_progressive  SOCIALsex_BB_complex_progressive  SOCIALage_BB_complex_progressive
# 714.423265                          0.000000                         11.675032                          1.370233                          0.000000

# 3.6. Complex - run cTADA ------------------------------------------------

# the run on the complex solution (more ILVs)

constraintsVectMatrix <- create.constraints.Vect.Matrix(CP_complex_progressive_complex$nbdadata, num_networks=1, num_ILVs=4)
ILVs <- c("sex", "age", "prev_exp_one", "prev_exp_both")
colnames(constraintsVectMatrix) <- c("network", paste("asoc", ILVs , sep="_"), paste("soc", ILVs, sep="_"))
constraintsVectMatrix <- as.matrix(constraintsVectMatrix)
class(constraintsVectMatrix)
head(constraintsVectMatrix)

# run cTADA
AICtable_complex_progressive_complex <- tadaAICtable(nbdadata = list(
  BB_complex_progressive_complex$nbdadata,
  CP_complex_progressive_complex$nbdadata,
  GW_complex_progressive_complex$nbdadata,
  MP_complex_progressive_complex$nbdadata
), constraintsVectMatrix = constraintsVectMatrix)

AICtable_complex_progressive_complex@printTable

networksSupport(AICtable_complex_progressive_complex)

# 0 0.5341161             16
# 1 0.4658839            256

setwd("C:/Users/swild/Desktop/Konstanz/Collaborations/Lucy, Michael - Phil Trans cumulative culture/Analysis/NBDA output/Results COMPLEX PROGRESSIVE diffusion")
save(AICtable_complex_progressive_complex, file="AICTable_complex_progressive.RData")
# load("AICTable_complex_progressive.RData")

# we have more support for asocial (0.53) over social models (0.47)
# note that the number of models differs qwith 256 social over only 16 asocial models
# so there appears to be more support for asocial recombination over social learning
network.support.complex.progressive_complex <- networksSupport(AICtable_complex_progressive_complex)
network.support.complex.progressive_complex
write.table(network.support.complex.progressive_complex, "network.support.complex.progressive.diffusion.txt")

# we can also extract the type support to take a closer look
# asocial models have most support
type.support.complex.progressive_complex <- typeSupport(AICtable_complex_progressive_complex)
type.support.complex.progressive_complex

# extract the variable support
variable.support.complex.progressive_complex <- variableSupport(AICtable_complex_progressive_complex)
variable.support.complex.progressive_complex
write.table(variable.support.complex.progressive_complex, "variable.support.complex.progressive.diffusion.txt")

# s1 ASOC:sex_BB_complex_progressive ASOC:age_BB_complex_progressive ASOC:pre.exp.one_BB_complex_progressive ASOC:pre.exp.both_BB_complex_progressive SOCIAL:sex_BB_complex_progressive
# support 0.4658839                       0.9796956                       0.2485206                               0.9999978                                0.9999999                        0.08429036
# SOCIAL:age_BB_complex_progressive SOCIAL:pre.exp.one_BB_complex_progressive SOCIAL:pre.exp.both_BB_complex_progressive
# support                          0.121515                                0.08013701                                 0.08730101


# extract model averaged estimates
MLE.complex.progressive_complex  <- modelAverageEstimates(AICtable_complex_progressive_complex,averageType = "median")
MLE.complex.progressive_complex
write.table(MLE.complex.progressive_complex, "MLE.complex.progressive.diffusion.txt")

# ASOCIALsex_BB_complex_progressive          ASOCIALage_BB_complex_progressive  ASOCIALpre.exp.one_BB_complex_progressive ASOCIALpre.exp.both_BB_complex_progressive
# 0.000000                                   1.140472                                   0.000000                                   3.181277                                   3.572019
# SOCIALsex_BB_complex_progressive           SOCIALage_BB_complex_progressive   SOCIALpre.exp.one_BB_complex_progressive  SOCIALpre.exp.both_BB_complex_progressive
# 0.000000                                   0.000000                                   0.000000                                   0.000000


# 3.5. Extract effect sizes -----------------------------------------------


# 3.5.1 Effect sizes - Dial -----------------------------------------------

## here, we extract estimates for the social learning parameter s
# we first constrain the NBDA data objects to contain the parameters from the best performing model
AICtable_complex_progressive_dial@printTable
# best mode is number 8

# we quickly rerun the constriantsVectMatrix
constraintsVectMatrix <- create.constraints.Vect.Matrix(CP_complex_progressive_dial$nbdadata, num_networks=1, num_ILVs=2)
constraintsVectMatrix[8,]

# create constrained NBDA data objects
constrained.complex.BB.dial <- constrainedNBDAdata(BB_complex_progressive_dial$nbdadata, constraintsVect = constraintsVectMatrix[8,])
constrained.complex.CP.dial <- constrainedNBDAdata(CP_complex_progressive_dial$nbdadata, constraintsVect = constraintsVectMatrix[8,])
constrained.complex.GW.dial <- constrainedNBDAdata(GW_complex_progressive_dial$nbdadata, constraintsVect = constraintsVectMatrix[8,])
constrained.complex.MP.dial <- constrainedNBDAdata(MP_complex_progressive_dial$nbdadata, constraintsVect = constraintsVectMatrix[8,])

# run TADA on the best model to get an estimate for s
bestModel.complex.dial <- tadaFit(nbdadata = list(
  constrained.complex.BB.dial,
  constrained.complex.CP.dial,
  constrained.complex.GW.dial,
  constrained.complex.MP.dial
), type="social")

# check variable names and output parameters
cbind.data.frame(bestModel.complex.dial@outputPar, bestModel.complex.dial@varNames)

bestModel.complex.dial@optimisation

prop.solve.social.byevent.complex.dial <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.complex.BB.dial,
  constrained.complex.CP.dial,
  constrained.complex.GW.dial,
  constrained.complex.MP.dial
), model=bestModel.complex.dial)
prop.solve.social.byevent.complex.dial
# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.complex.dial <- oadaPropSolveByST( nbdadata = list(
  constrained.complex.BB.dial,
  constrained.complex.CP.dial,
  constrained.complex.GW.dial,
  constrained.complex.MP.dial
), model=bestModel.complex.dial)
prop.solve.social.complex.dial
# overall, an estimated 91.5 % of birds have learned the complex task socially

bestModel.complex.dial@varNames

# plot the profile likelihood interval
# first for s (this is based on the second best model, as the top model is an asocial model)
plotProfLik(which=1, model=bestModel.complex.dial, range=c(0, 100), resolution = 10)
# we can see that the confidence intervals approaches infinity, so we can only estimate the lower bound
# as the upper one is 100%
CIs <- profLikCI(which=1, model=bestModel.complex.dial, lowerRange =c(0, 20))


# We extract the Lower limit for s in %
bestModelDataS1LowerBound.complex.BB.dial <- constrainedNBDAdata(
  nbdadata =
    BB_complex_progressive_dial$nbdadata,
  constraintsVect = constraintsVectMatrix[8, ],
  offset = c(CIs[1]    , rep(0, 4))
)

bestModelDataS1LowerBound.complex.CP.dial <- constrainedNBDAdata(
  nbdadata =
    CP_complex_progressive_dial$nbdadata,
  constraintsVect = constraintsVectMatrix[8, ],
  offset = c(CIs[1]    , rep(0, 4))
)


bestModelDataS1LowerBound.complex.GW.dial <- constrainedNBDAdata(
  nbdadata =
    GW_complex_progressive_dial$nbdadata,
  constraintsVect = constraintsVectMatrix[8, ],
  offset = c(CIs[1]   , rep(0, 4))
)

bestModelDataS1LowerBound.complex.MP.dial <- constrainedNBDAdata(
  nbdadata =
    MP_complex_progressive_dial$nbdadata,
  constraintsVect = constraintsVectMatrix[8, ],
  offset = c(CIs[1]   , rep(0, 4))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1LowerBound.complex_dial <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.complex.BB.dial,
      bestModelDataS1LowerBound.complex.CP.dial,
      bestModelDataS1LowerBound.complex.GW.dial,
      bestModelDataS1LowerBound.complex.MP.dial
    ) ,
    type = "asocial"
  )
bestModelS1LowerBound.complex_dial@outputPar
bestModelS1LowerBound.complex_dial@varNames

# [1] 219.446419   1.184576   1.243763
#Now plug into the prop solve function
prop.solve.social.Lower.complex.dial <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound.complex_dial,
    nbdadata = list(
      bestModelDataS1LowerBound.complex.BB.dial,
      bestModelDataS1LowerBound.complex.CP.dial,
      bestModelDataS1LowerBound.complex.GW.dial,
      bestModelDataS1LowerBound.complex.MP.dial
    )
  )
prop.solve.social.Lower.complex.dial
# Lower bound for % of birds having learned the complex task through social learning is 84.6%

## now let's extract effect sizes for age influencing the asocial learning rate
# plot profile likelihood for the influence of sex on social learning (which=2)
plotProfLik(which=2, model=bestModel.complex.dial, range=c(-2, 9), resolution = 10)
plotProfLik(which=2, model=bestModel.complex.dial, range=c(10, 20), resolution = 10) # upper bound not estimable - does not converge
profLikCI(which=2,model=bestModel.complex.dial, lowerRange=c(-1,1))

#    Lower CI   Upper CI
# -0.1801531         NA

exp(c(MLE.complex.progressive_dial[3], -0.1801531, NA))

# back transform the estimates
# Adults were 7.875597e+05 [0.8-NA] times faster at learning the complex task socially compared to first years
# contradictory effect - they could be faster, but may as well be slower

## now let's extract effect sizes for sex influencing the social learning rate
# plot profile likelihood for the influence of sex on social learning (which=3)
# (here which=1 is no longer referring to s, as we are fitting an asocial model)
plotProfLik(which=3, model=bestModel.complex.dial, range=c(0, 4), resolution = 10) # model does not converge
profLikCI(which=3,model=bestModel.complex.dial, lowerRange=c(0,1), upperRange = c(1.5, 2.5))

# Lower CI  Upper CI
# 0.5367547 1.7481495

exp(c(MLE.complex.progressive_dial[4],0.5367547, 1.7481495 ))

# Males were 3.2 x faster at learning dial compared to females [1.7-5.7]

# 3.5.2 Effect sizes - Slide -----------------------------------------------

## here, we extract estimates for the social learning parameter s
# we first constrain the NBDA data objects to contain the parameters from the best performing model
AICtable_complex_progressive_slide@printTable
# which is Model 8
# create constrained NBDA data objects
constrained.complex.BB.slide <- constrainedNBDAdata(BB_complex_progressive_slide$nbdadata, constraintsVect = constraintsVectMatrix[8,])
constrained.complex.CP.slide <- constrainedNBDAdata(CP_complex_progressive_slide$nbdadata, constraintsVect = constraintsVectMatrix[8,])
constrained.complex.GW.slide <- constrainedNBDAdata(GW_complex_progressive_slide$nbdadata, constraintsVect = constraintsVectMatrix[8,])
constrained.complex.MP.slide <- constrainedNBDAdata(MP_complex_progressive_slide$nbdadata, constraintsVect = constraintsVectMatrix[8,])



# run TADA on the best model to get an estimate for s
bestModel.complex.slide <- tadaFit(nbdadata = list(
  constrained.complex.BB.slide,
  constrained.complex.CP.slide,
  constrained.complex.GW.slide,
  constrained.complex.MP.slide
), type="social")

# check variable names and output parameters
cbind.data.frame(bestModel.complex.slide@outputPar, bestModel.complex.slide@varNames)

# in the best model, sex influences the social and asocial learning rate and age influences the asocial learning rate

bestModel.complex.slide@optimisation

prop.solve.social.byevent.complex.slide <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.complex.BB.slide,
  constrained.complex.CP.slide,
  constrained.complex.GW.slide,
  constrained.complex.MP.slide
), model=bestModel.complex.slide)
prop.solve.social.byevent.complex.slide
# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.complex.slide <- oadaPropSolveByST( nbdadata = list(
  constrained.complex.BB.slide,
  constrained.complex.CP.slide,
  constrained.complex.GW.slide,
  constrained.complex.MP.slide
), model=bestModel.complex.slide)
prop.solve.social.complex.slide
# overall, an estimated 71.0% of birds have learned the complex task socially

bestModel.complex.slide@varNames
# [1] "Scale (1/rate):"                       "1 Social transmission 1"               "2 Asocial: age_BB_complex_progressive" "3 Social: sex_BB_complex_progressive"

# plot the profile likelihood interval
# first for s (this is based on the second best model, as the top model is an asocial model)
plotProfLik(which=1, model=bestModel.complex.slide, range=c(0, 100), resolution = 10)
# again, the upper limit seems to approach infitniy
CIs <- profLikCI(which=1, model=bestModel.complex.slide, lowerRange =c(0, 20))


# We extract the Lower limit for s in %
bestModelDataS1LowerBound.complex.BB.slide <- constrainedNBDAdata(
  nbdadata =
    BB_complex_progressive_slide$nbdadata,
  constraintsVect = constraintsVectMatrix[8, ],
  offset = c(CIs[1]    , rep(0, 4))
)

bestModelDataS1LowerBound.complex.CP.slide <- constrainedNBDAdata(
  nbdadata =
    CP_complex_progressive_slide$nbdadata,
  constraintsVect = constraintsVectMatrix[8, ],
  offset = c(CIs[1]    , rep(0, 4))
)


bestModelDataS1LowerBound.complex.GW.slide <- constrainedNBDAdata(
  nbdadata =
    GW_complex_progressive_slide$nbdadata,
  constraintsVect = constraintsVectMatrix[8, ],
  offset = c(CIs[1]   , rep(0, 4))
)

bestModelDataS1LowerBound.complex.MP.slide <- constrainedNBDAdata(
  nbdadata =
    MP_complex_progressive_slide$nbdadata,
  constraintsVect = constraintsVectMatrix[8, ],
  offset = c(CIs[1]   , rep(0, 4))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1LowerBound.complex_slide <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.complex.BB.slide,
      bestModelDataS1LowerBound.complex.CP.slide,
      bestModelDataS1LowerBound.complex.GW.slide,
      bestModelDataS1LowerBound.complex.MP.slide
    ) ,
    type = "asocial"
  )
bestModelS1LowerBound.complex_slide@outputPar
# [1] 288.063122   3.150089   1.466819
#Now plug into the prop solve function
prop.solve.social.Lower.complex.slide <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound.complex_slide,
    nbdadata = list(
      bestModelDataS1LowerBound.complex.BB.slide,
      bestModelDataS1LowerBound.complex.CP.slide,
      bestModelDataS1LowerBound.complex.GW.slide,
      bestModelDataS1LowerBound.complex.GW.slide
    )
  )
prop.solve.social.Lower.complex.slide
# Lower bound for % of birds having learned the complex task through social learning is 56.2%

bestModel.complex.slide@varNames

## now let's extract effect sizes for age influencing the asocial learning rate
# plot profile likelihood for the influence of sex on social learning (which=2)
# (here which=1 is no longer referring to s, as we are fitting an asocial model)
plotProfLik(which=2, model=bestModel.complex.slide, range=c(0, 20), resolution = 10)
# the upper range cannot be reliably obtained
profLikCI(which=2,model=bestModel.complex.slide, lowerRange=c(1,2))
# Lower CI Upper CI
# 1.23667       NA

exp(c(MLE.complex.progressive_slide[3], 1.23667, NA))
# adults 1.175985e+05 times faster at learning compared to juveniles [3.4-NA]
# back transform the estimates


# extract effect of sex influencing asocial learning (which=3)
plotProfLik(which=3, model=bestModel.complex.slide, range=c(0, 5), resolution = 10)
profLikCI(which=3,model=bestModel.complex.slide, upperRange =c(2,3), lowerRange=c(0,1))

# Lower CI  Upper CI
# 0.5192431 2.6579167

exp(c(MLE.complex.progressive_slide[4], 0.5192431, 2.6579167))
# males were 3.9 [1.8-14.3] times faster at learning compared to females

# 3.5.3 Extract effect sizes - complex -------------------------------------------------------------------

constraintsVectMatrix <- create.constraints.Vect.Matrix(CP_complex_progressive_complex$nbdadata, num_networks=1, num_ILVs=4)
AICtable_complex_progressive_complex@printTable
constraintsVectMatrix[259,] # asocial model
constraintsVectMatrix[213,] # corresponding social model

## here, we extract estimates for the social learning parameter s
# we first constrain the NBDA data objects to contain the parameters from the best performing model (model 259)
# the best model is an asocial model, for which the funciton to build a constrained object does not work
# we take the equivalent social model (213), which is also the third best model, to create the constrained NBDa data ojbects
# create constrained NBDA data objects
constrained.complex.BB.complex <- constrainedNBDAdata(BB_complex_progressive_complex$nbdadata, constraintsVect = constraintsVectMatrix[213,])
constrained.complex.CP.complex <- constrainedNBDAdata(CP_complex_progressive_complex$nbdadata, constraintsVect = constraintsVectMatrix[213,])
constrained.complex.GW.complex <- constrainedNBDAdata(GW_complex_progressive_complex$nbdadata, constraintsVect = constraintsVectMatrix[213,])
constrained.complex.MP.complex <- constrainedNBDAdata(MP_complex_progressive_complex$nbdadata, constraintsVect = constraintsVectMatrix[213,])

# run TADA on the second best model to get an estimate for s
bestModel.complex <- tadaFit(nbdadata = list(
  constrained.complex.BB.complex,
  constrained.complex.CP.complex,
  constrained.complex.GW.complex,
  constrained.complex.MP.complex
), type="social")

# check variable names and output parameters
cbind.data.frame(bestModel.complex@outputPar, bestModel.complex@varNames)

# sex influencing asocial learning
# previous experience of one influencing asocial learning
# previous experience of both influencing asocial learning

bestModel.complex@optimisation

prop.solve.social.byevent.complex <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.complex.BB.complex,
  constrained.complex.CP.complex,
  constrained.complex.GW.complex,
  constrained.complex.MP.complex
), model=bestModel.complex)
prop.solve.social.byevent.complex
# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.complex <- oadaPropSolveByST( nbdadata = list(
  constrained.complex.BB.complex,
  constrained.complex.CP.complex,
  constrained.complex.GW.complex,
  constrained.complex.MP.complex
), model=bestModel.complex)
prop.solve.social.complex
# overall, an estimated 3.2% of birds have learned the complex task socially

bestModel.complex@varNames

# plot the profile likelihood interval
# first for s (this is based on the second best model, as the top model is an asocial model)
plotProfLik(which=1, model=bestModel.complex, range=c(0, 50), resolution = 10)
# lower bound spans 0
CIs <- profLikCI(which=1, model=bestModel.complex,  upperRange =c(30, 40))

# Lower CI Upper CI
# 0.00000 36.22077

# ### lower bound
# # We extract the upper limit for s in %
# bestModelDataS1lowerBound.complex.BB <- constrainedNBDAdata(
#   nbdadata =
#     BB_complex_progressive_complex$nbdadata,
#   constraintsVect = constraintsVectMatrix[213, ],
#   offset = c(CIs[1]    , rep(0, 8))
# )
#
# bestModelDataS1lowerBound.complex.CP <- constrainedNBDAdata(
#   nbdadata =
#     CP_complex_progressive_complex$nbdadata,
#   constraintsVect = constraintsVectMatrix[213, ],
#   offset = c(CIs[1]    , rep(0, 8))
# )
#
#
# bestModelDataS1lowerBound.complex.GW <- constrainedNBDAdata(
#   nbdadata =
#     GW_complex_progressive_complex$nbdadata,
#   constraintsVect = constraintsVectMatrix[213, ],
#   offset = c(CIs[1]   , rep(0, 8))
# )
#
# bestModelDataS1lowerBound.complex.MP <- constrainedNBDAdata(
#   nbdadata =
#     MP_complex_progressive_complex$nbdadata,
#   constraintsVect = constraintsVectMatrix[213, ],
#   offset = c(CIs[1]    , rep(0, 8))
# )
#
#
# #Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
# bestModelS1lowerBound.complex <-
#   tadaFit(
#     list(
#       bestModelDataS1lowerBound.complex.BB,
#       bestModelDataS1lowerBound.complex.CP,
#       bestModelDataS1lowerBound.complex.GW,
#       bestModelDataS1lowerBound.complex.MP
#     ) ,
#     type = "asocial"
#   )
# bestModelS1lowerBound.complex@outputPar
# # [1] 843.618715   1.264200   3.475658   3.875846
# #Now plug into the prop solve function
# prop.solve.social.lower.complex <-
#   oadaPropSolveByST(
#     model = bestModelS1lowerBound.complex,
#     nbdadata = list(
#       bestModelDataS1lowerBound.complex.BB,
#       bestModelDataS1lowerBound.complex.CP,
#       bestModelDataS1lowerBound.complex.GW,
#       bestModelDataS1lowerBound.complex.MP
#     )
#   )
# prop.solve.social.lower.complex
# # lower bound for % of birds having learned the complex task through social learning is 4.4%





# We extract the upper limit for s in %
bestModelDataS1UpperBound.complex.BB <- constrainedNBDAdata(
  nbdadata =
    BB_complex_progressive_complex$nbdadata,
  constraintsVect = constraintsVectMatrix[213, ],
  offset = c(CIs[2]    , rep(0, 8))
)

bestModelDataS1UpperBound.complex.CP <- constrainedNBDAdata(
  nbdadata =
    CP_complex_progressive_complex$nbdadata,
  constraintsVect = constraintsVectMatrix[213, ],
  offset = c(CIs[2]    , rep(0, 8))
)


bestModelDataS1UpperBound.complex.GW <- constrainedNBDAdata(
  nbdadata =
    GW_complex_progressive_complex$nbdadata,
  constraintsVect = constraintsVectMatrix[213, ],
  offset = c(CIs[2]   , rep(0, 8))
)

bestModelDataS1UpperBound.complex.MP <- constrainedNBDAdata(
  nbdadata =
    MP_complex_progressive_complex$nbdadata,
  constraintsVect = constraintsVectMatrix[213, ],
  offset = c(CIs[2]    , rep(0, 8))
)

#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1UpperBound.complex <-
  tadaFit(
    list(
      bestModelDataS1UpperBound.complex.BB,
      bestModelDataS1UpperBound.complex.CP,
      bestModelDataS1UpperBound.complex.GW,
      bestModelDataS1UpperBound.complex.MP
    ) ,
    type = "asocial"
  )
bestModelS1UpperBound.complex@outputPar
# [1] 3749.524953    1.260063    4.929133    5.307033

#Now plug into the prop solve function
prop.solve.social.Upper.complex <-
  oadaPropSolveByST(
    model = bestModelS1UpperBound.complex,
    nbdadata = list(
      bestModelDataS1UpperBound.complex.BB,
      bestModelDataS1UpperBound.complex.CP,
      bestModelDataS1UpperBound.complex.GW,
      bestModelDataS1UpperBound.complex.MP
    )
  )
prop.solve.social.Upper.complex
# Upper bound for % of birds having learned the complex task through social learning is 11.1%
#
## we now leave the second best model and rerun the tadaFit to go to the best (asocial) model
bestModel.complex <- tadaFit(nbdadata = list(
  constrained.complex.BB.complex,
  constrained.complex.CP.complex,
  constrained.complex.GW.complex,
  constrained.complex.MP.complex
), type="asocial")

bestModel.complex@varNames

# plot profile likelihood for the influence of sex on asocial learning (which=1)
# (here which=1 is no longer referring to s, as we are fitting an asocial model)
plotProfLik(which=1, model=bestModel.complex, range=c(-1, 3), resolution = 10)
profLikCI(which=1,model=bestModel.complex, upperRange =c(1.5,2.5), lowerRange=c(0,1))

#   Lower CI  Upper CI
# 0.4236277 1.9653743

 exp(c(MLE.complex.progressive_complex[2],0.4236277, 1.9653743 ))

# back transform the estimates
# Males were 3.1 [1.4-6.2] times faster at learning the complex task asocially compared to females

# we do the same for the previous experience of one component on asocial learning
# plot profile likelihood for the influence of previous experience of one component on asocial learning (which=2)
plotProfLik(which=2, model=bestModel.complex, range=c(0, 5), resolution = 10)
profLikCI(which=2,model=bestModel.complex, lowerRange=c(2,3), upperRange = c(4,5))

# Lower CI Upper CI
# 2.233116 4.404639

exp(c(MLE.complex.progressive_complex[4], 2.23311, 4.404639))

# Birds with previous experience of one component were 24.1 [9.3-81.8] times faster at learning the complex task asocially compared to those without any knowledge

# and finally, we extract the effect of having knowledge of both components over having no or partial knowledge
plotProfLik(which=3, model=bestModel.complex, range=c(1.5, 10), resolution = 10)
profLikCI(which=3,model=bestModel.complex, lowerRange=c(2,4), upperRange = c(4,6))

# Lower CI Upper CI
# 2.531724 4.840074

exp(c(MLE.complex.progressive_complex[5],  2.539212, 4.847395 ))

# birds with knowledge of both components were 35.6 [12.6-126.5] times faster at asocial learning compared to those without any knowledge, or knowledge of one
######################################################################### #
# 4) COMPLEX 2nd gen ----------------------------------

# 4.1. assign the data -------------------------------------------------------

gmm.BB <- NBDA_input_data$data.objects.2017$gmm.BB
gmm.CP <- NBDA_input_data$data.objects.2017$gmm.CP
demos <- NBDA_input_data$data.objects.2017$demos # contains IDs of those that have solved the complex task during the previuos experiment
ILVs.next.gen <- NBDA_input_data$data.objects.2017$ILVs.next.gen
greti.list.next.gen <- unique(subset(ILVs.next.gen$Ring, ILVs.next.gen$Species=="greti"))
prev_exp <- NBDA_input_data$data.objects.2016$prev_exp # contains the number of slide and dial solves up until the start of the 1st generation experiment
# so we have to add the ones that have learned in the 1st gen experiment
dial.demos <- subset(prev_exp$RING, prev_exp$DIAL_PAST>=3)
dial.demos <- unique(c(dial.demos, BB_complex_progressive_dial$learners.list$RING,CP_complex_progressive_dial$learners.list$RING,GW_complex_progressive_dial$learners.list$RING,MP_complex_progressive_dial$learners.list$RING))

slide.demos <- subset(prev_exp$RING, prev_exp$DOOR_PAST>=3)
slide.demos <- unique(c(slide.demos, BB_complex_progressive_slide$learners.list$RING,CP_complex_progressive_slide$learners.list$RING,GW_complex_progressive_slide$learners.list$RING,MP_complex_progressive_slide$learners.list$RING))

# 4.2 prepare NBDA data - learning simple and/or complex ------------------

prepare.NBDA.data.next.gen <- function(input.file, locations, label, ILVs.yes.no, start.date, end.date, min.time, gmm, patch, learned.behav, demos.def, ILVs.include, directed){
  # subset the input file to the correct experiment
  input.file.sub <-
    subset(
      input.file,
      input.file$PATCH == patch &
        input.file$SPECIES == "GRETI" &
        input.file$experiment %in% c("nextgen", "complex_prog")
    )
  # we keep all solves in the input file for now, as we will need to extract the experience across the two experiments

  # a bird is considered a learner if it has solved three times
  # the time of its first solve is taken as the time of acquisition
  # prepare empty vectors to store the ID of the learner, the number of solves and the time of the first solve

  learners.list <- NULL
  num_solves_all <- NULL
  time_solves <- NULL
  which.rows <- NULL

  # create a file with only the 2nd gen data
  input.complex.next.gen <-
    subset(
      input.file.sub,
      #input.file.sub$solution_category == "complex" &
        input.file.sub$experiment == "nextgen"
    )

  if(learned.behav=="complex"){
    for (i in as.vector(unique(input.complex.next.gen$RING[input.complex.next.gen$solution_category==learned.behav]))){
      sub <- subset(input.complex.next.gen, input.complex.next.gen$RING==i & input.complex.next.gen$solution_category==learned.behav)
      which.row <- rownames(sub)[1]
      num.solves <- length(sub[,1])
      solve.time <- as.numeric(difftime(min(sub$time_stamp), min(input.file.sub$time_stamp[input.file.sub$experiment=="nextgen"]), units="days"))
      learners.list[which(unique(input.complex.next.gen$RING)==i)] <- i
      num_solves_all[which(unique(input.complex.next.gen$RING)==i)] <- num.solves
      time_solves[which(unique(input.complex.next.gen$RING)==i)] <- solve.time
      which.rows[which(unique(input.complex.next.gen$RING)==i)] <- which.row
    }


    # reassign column names
    learners.list_df <- na.omit(cbind.data.frame(learners.list, num_solves_all,  time_solves, which.rows))
    colnames(learners.list_df) <- c("RING", "num_solves", "time_acq", "which_rows")
    learners.list_df <- learners.list_df[order(learners.list_df$time_acq),]

    # remove those with fewer than 3 solves
    learners.list_df <- subset(learners.list_df, learners.list_df$num_solves>=3)

    # restrict the data frame to great tits only that are in the ILV file
    learners.list_df <- subset(learners.list_df, learners.list_df$RING%in%ILVs.next.gen$Ring)
  } else { # if it is either 'dial' or 'slide'
    for (i in as.vector(unique(input.complex.next.gen$RING[input.complex.next.gen$solution_category %in% c(learned.behav, "complex")]))){ # as there is a simple solution in each complex solution
      sub <- subset(input.complex.next.gen, input.complex.next.gen$RING==i & input.complex.next.gen$solution_category %in% c(learned.behav, "complex"))
      which.row <- rownames(sub)[1]
      num.solves <- length(sub[,1])
      solve.time <- as.numeric(difftime(min(sub$time_stamp), min(input.file.sub$time_stamp[input.file.sub$experiment=="nextgen"]), units="days"))
      learners.list[which(unique(input.complex.next.gen$RING)==i)] <- i
      num_solves_all[which(unique(input.complex.next.gen$RING)==i)] <- num.solves
      time_solves[which(unique(input.complex.next.gen$RING)==i)] <- solve.time
      which.rows[which(unique(input.complex.next.gen$RING)==i)] <- which.row
    }
  }

    # reassign column names
    learners.list_df <- na.omit(cbind.data.frame(learners.list, num_solves_all,  time_solves, which.rows))
    colnames(learners.list_df) <- c("RING", "num_solves", "time_acq", "which_rows")
    learners.list_df <- learners.list_df[order(learners.list_df$time_acq),]

    # remove those with fewer than 3 solves
    learners.list_df <- subset(learners.list_df, learners.list_df$num_solves>=3)

    # restrict the data frame to great tits only that are in the ILV file
    learners.list_df <- subset(learners.list_df, learners.list_df$RING%in%ILVs.next.gen$Ring)


    if(demos.def == "single" & learned.behav=="dial"){
      demos.sub <- dial.demos # extract those with at least 3 dial solves
    } else if(demos.def =="both" & learned.behav=="dial"){
      demos.sub <- dial.demos # extract those with at least 3 dial solves
      demos.sub <- unique(c(demos.sub, demos))
    } else if(demos.def == "single" & learned.behav=="slide"){
      demos.sub <- slide.demos
    } else if(demos.def =="both" & learned.behav=="slide"){
      demos.sub <- slide.demos # extract those with at least 3 dial solves
      demos.sub <- unique(c(demos.sub, demos))
    } else if(learned.behav == "complex" ){
      demos.sub <- demos
    }


  learners.list_df <- subset(learners.list_df, !(learners.list_df$RING %in% demos.sub))

  # create network from group by individual matrix
  gmm$metadata$location <- tolower(substr(gmm$metadata$Location, 1, 2))

  # the start and end times in the gmm object are stored as seconds after the deployment of the feeders
  # so we add another column which transforms start and end time back into date/time format
  gmm$metadata$Start.Date <- as.POSIXct(min(min.time), format="%Y%m%d%H%M%S")+ lubridate::seconds(gmm$metadata$Start)
  gmm$metadata$End.Date <- as.POSIXct(min(min.time), format="%Y%m%d%H%M%S")+ lubridate::seconds(gmm$metadata$End)
  gmm$metadata <- gmm$metadata[order(gmm$metadata$Start.Date),]

  # subset the gbi to the corresponding locations and times for each diffusion
  which.groups <- subset(gmm$metadata, gmm$metadata$location %in% locations & gmm$metadata$Start.Date >= start.date & gmm$metadata$End.Date<=end.date)
  which.groups.num <- rownames(which.groups)
  gbi <- as.matrix(gmm$gbi)

  sub.gbi <- gbi[as.numeric(which.groups.num),] # subset the gbi to the right locations and times
  sub.gbi <- sub.gbi[, colSums(sub.gbi)>9] # remove birds that weren't seen at all during this period

  # subset to gretis only
  sub.gbi <- sub.gbi[, colnames(sub.gbi) %in% greti.list.next.gen]

  # generate association network
  net <- get_network(association_data = sub.gbi, data_format="GBI", association_index = "SRI")
  net <- net[order(rownames(net)), order(colnames(net))]


  ## the association matrix needs to be in an array
  assocMatrix <- array(data = net, dim=c(nrow(net), ncol(net), 1))
  class(assocMatrix)

  ## get ILVs ready
  ids.sub <- subset(ILVs.next.gen, ILVs.next.gen$Ring %in% rownames(net))

  length(unique(ids.sub)[,1])

  # sort according to ring number
  # (need to correspond to order of rownames in matrix)
  ids.sub <- ids.sub[order(ids.sub$Ring),]

  # ensure that we only keep learners.list that are included in the network
  learners.list_df <- subset(learners.list_df, learners.list_df$RING %in% rownames(net))


  # adjust the times of acquisition
  # recalculate the time of acquisition by removing the days where the puzzle boxes weren't out (network days)
  # days 0-4 puzzle box
  # days 5+6 network days
  # days 7-11 puzzle box (so subtract two days from the acquisition of birds that learned between day 8-12)
  # days 12+13 network days
  # day 14-18 (subtract 4 days from the acquisition time for birds that have learned between day 15-19)
  # days 19+20 network days
  # day 21-25 (subtract 6)

  learners.list_df$time_acq_adj <- NA

  for (i in 1:length(learners.list_df[,1])){
    time <- learners.list_df[i, "time_acq"]
    if(time > 7 & time < 12){
      time.adj <- time-2
    } else if(time > 14 & time < 19){
      time.adj <- time-4
    } else if(time > 21 & time < 26){
      time.adj <- time-6
    } else {time.adj <- time}
    learners.list_df[i, "time_acq_adj"] <- time.adj
  }

  # prepare individual level variables

  ids.sub$Age[ids.sub$Age=="adult"] <- 0.5 # adults are assigned 0.5
  ids.sub$Age[ids.sub$Age=="juvenile"] <- -0.5 # juveniles are assigned -0.5

  # as we are coding the preveious experience as a time-varying variable,
  # we also need to create matrices for sex and age, where each column represents one acquisition event

  sex <- matrix(as.numeric(ids.sub$Sex), nrow=length(ids.sub$Sex), ncol=length(learners.list_df$which_rows), byrow=F)
  age <- matrix(as.numeric(ids.sub$Age), nrow=length(ids.sub$Age), ncol=length(learners.list_df$which_rows), byrow=F)

  # extract previous experience (as in the first generation)

  prev.experience.matrix.one <- as.data.frame(matrix(NA, ncol=length(learners.list_df$which_rows), nrow=length(rownames(net))))
  prev.experience.matrix.both <- as.data.frame(matrix(NA, ncol=length(learners.list_df$which_rows), nrow=length(rownames(net))))


  for (k in as.numeric(learners.list_df$which_rows)){
    prev.exp.one <- NULL
    prev.exp.both <- NULL

    for (i in rownames(net)){

      # the file with the previous experience contains the experience up until after the dial experiment (so excludes the complex progressive/1st gen experience)
      # we therefore need to add the experience that the birds gained both during the first gen diffusion + during the second gen diffusion up to the point of each acquisition event
      sub.prev.exp <- subset(prev_exp, prev_exp$RING==i & prev_exp$SITE==patch)
      exp.at.start.slide <- sum(sub.prev.exp[,"DOOR_PAST"])
      exp.at.start.dial <- sum(sub.prev.exp[,"DIAL_PAST"])

      # subset the input file to dial and slide only
      # (from complex progressive and nextgen)
      input.file.sub.slide.dial <- subset(input.file.sub, input.file.sub$solution_category%in% c("slide", "dial"))

      # extrat the closest event to the complex solve
      suppressWarnings(
        index.of.closest <- which(abs(as.numeric(rownames(input.file.sub.slide.dial))-k)==min(abs(as.numeric(rownames(input.file.sub.slide.dial))-k)))
      )


      input.file.sub.slide.dial.sub <- input.file.sub.slide.dial[1:index.of.closest,] # subset the input file up to the point of where a bird has learned
      sub.input.i.dial <- subset(input.file.sub.slide.dial.sub, input.file.sub.slide.dial.sub$RING==i & input.file.sub.slide.dial.sub$solution_category=="dial") # subset the data frame only for individual i
      sub.input.i.slide <- subset(input.file.sub.slide.dial.sub, input.file.sub.slide.dial.sub$RING==i & input.file.sub.slide.dial.sub$solution_category=="slid") # subset the data frame only for individual i

      max.i.dial <- length(sub.input.i.dial[,1])
      max.i.slide <- length(sub.input.i.slide[,1])


      if(max.i.dial>0){
        exp.at.acquisition.dial <- exp.at.start.dial + max.i.dial # build the sum of the experience at the start and the experience gathered during the experiment
      } else {
        exp.at.acquisition.dial <- exp.at.start.dial
      }

      if(max.i.slide>0){
        exp.at.acquisition.slide <- exp.at.start.slide + max.i.slide # build the sum of the experience at the start and the experience gathered during the experiment
      } else {
        exp.at.acquisition.slide <- exp.at.start.slide
      }

      # next we build two vectors:
      # one that gets a 1 if birds have experience (>=3 solves) with both components
      # one that gets a 1 if birds have experience (>= solves) with one component but not the other
      if(exp.at.acquisition.slide>=3 & exp.at.acquisition.dial>=3){
        prev.exp.both[which(rownames(net)==i)] <- 1
        prev.exp.one[which(rownames(net)==i)] <- 0
      } else if (exp.at.acquisition.slide>=3 & exp.at.acquisition.dial<3 | exp.at.acquisition.slide<3 & exp.at.acquisition.dial>=3){
        prev.exp.one[which(rownames(net)==i)] <- 1
        prev.exp.both[which(rownames(net)==i)] <- 0
      } else {
        prev.exp.one[which(rownames(net)==i)] <- 0
        prev.exp.both[which(rownames(net)==i)] <- 0
      }

    }
    prev.experience.matrix.one[,which(learners.list_df$which_rows==k)] <- prev.exp.one
    prev.experience.matrix.both[,which(learners.list_df$which_rows==k)] <- prev.exp.both

  }

  prev.experience.matrix.one <- as.matrix(prev.experience.matrix.one)
  prev.experience.matrix.both <- as.matrix(prev.experience.matrix.both)


  assign(paste("age", label, sep="_"), age, envir = .GlobalEnv)
  assign(paste("sex", label, sep="_"), sex, envir = .GlobalEnv)
  assign(paste("pre.exp.one", label, sep="_"), prev.experience.matrix.one, envir = .GlobalEnv)
  assign(paste("pre.exp.both", label, sep="_"), prev.experience.matrix.both, envir = .GlobalEnv)


  ILVs <- paste(ILVs.include, label, sep="_")
  # ILVs <- c(paste("pre.exp.one", label, sep="_"), paste("pre.exp.both", label, sep="_"))

  assign(paste("ILVs", label, sep="_"), ILVs)

  # get order of acquisition for learners.list
  learners.list_OAc <- NULL


  for (i in unique(learners.list_df$RING)){
    pos <- which(colnames(net)==i)
    learners.list_OAc <- c(learners.list_OAc, pos)
  }

  # prepare demonstrator vector
  # each demonstrator gets a 1, all naive individuals get a 0
  demo_vec <- rep(0, length(colnames(net)))
  demo_vec[which(colnames(net)%in%demos.sub)] <- 1

  ## if we want to include directed network to answer the question whether birds learn simple tasks from complex solves

  if(directed==TRUE){
    # we first need to extract the acquisition events of the complex solvers

    learners.list <- NULL
    num_solves_all <- NULL
    time_solves <- NULL
    which.rows <- NULL

    # create a file with only the 2nd gen data
    input.complex.next.gen <-
      subset(
        input.file.sub,
        #input.file.sub$solution_category == "complex" &
        input.file.sub$experiment == "nextgen"
      )

      for (i in as.vector(unique(input.complex.next.gen$RING[input.complex.next.gen$solution_category=="complex"]))){
        sub <- subset(input.complex.next.gen, input.complex.next.gen$RING==i & input.complex.next.gen$solution_category=="complex")
        which.row <- rownames(sub)[1]
        num.solves <- length(sub[,1])
        solve.time <- as.numeric(difftime(min(sub$time_stamp), min(input.complex.next.gen$time_stamp[input.complex.next.gen$experiment=="nextgen"]), units="days"))
        learners.list[which(unique(input.complex.next.gen$RING)==i)] <- i
        num_solves_all[which(unique(input.complex.next.gen$RING)==i)] <- num.solves
        time_solves[which(unique(input.complex.next.gen$RING)==i)] <- solve.time
        which.rows[which(unique(input.complex.next.gen$RING)==i)] <- which.row
      }


      # reassign column names
      learners.list_df_complex <- na.omit(cbind.data.frame(learners.list, num_solves_all,  time_solves, which.rows))
      colnames(learners.list_df_complex) <- c("RING", "num_solves", "time_acq", "which_rows")
      learners.list_df_complex <- learners.list_df_complex[order(learners.list_df_complex$time_acq),]

      # remove those with fewer than 3 solves
      learners.list_df_complex <- subset(learners.list_df_complex, learners.list_df_complex$num_solves>=3)

      # restrict the data frame to great tits only that are in the ILV file
      learners.list_df_complex <- subset(learners.list_df_complex, learners.list_df_complex$RING%in%ILVs.next.gen$Ring)

    # remove the demos that knew how to solve complex at the start of the experiment
    learners.list_df_complex <- subset(learners.list_df_complex, !(learners.list_df_complex$RING %in% demos))

    # adjust the acquisition times


    learners.list_df_complex$time_acq_adj <- NA

    for (i in 1:length(learners.list_df_complex[,1])){
      time <- learners.list_df_complex[i, "time_acq"]
      if(time > 7 & time < 12){
        time.adj <- time-2
      } else if(time > 14 & time < 19){
        time.adj <- time-4
      } else if(time > 21 & time < 26){
        time.adj <- time-6
      } else {time.adj <- time}
      learners.list_df_complex[i, "time_acq_adj"] <- time.adj
    }

    # now we have all the acquisition times of the complex solves

    # we create two sets of networks - one for directed complex -> simple, the other simple-> simple

    num.acquisition.events <- length(learners.list_OAc)

    #ASIDE:
    #If we have dynamic (time-varying) networks we need to create a four dimensional array of size
    #no. individuals x no.individuals x no.networks x number of time periods
    #and provide an assMatrixIndex vector as shown in Tutorial 1.


    # the association array contains 2x num.acquisition events matrices

    assocMatrix.full <- array(data = net, dim=c(nrow(net), ncol(net), 2, num.acquisition.events))

    # at each acquisition event, we need to state who is knowledgeable of only simple or knowledgeable of complex
    df.knowledgeable.complex <- vector(mode="list", length = num.acquisition.events)
    for(i in 1:num.acquisition.events){
      # who knows to solve complex:
        # all demonstrators
        # and all who have learned up to the acquistion event

      newly.learned.complex <- subset(learners.list_df_complex$RING, learners.list_df_complex$time_acq_adj<learners.list_df$time_acq_adj[i])
      knowledgeable.complex <- unique(c(demos, newly.learned.complex)) # add demos that knew how to solve complex before the start of the experiment and those that learned in this experiment
    df.knowledgeable.complex[[i]] <- knowledgeable.complex

    }

    # now we have to extract who was knowledgeable in just the simple solution at each acquisition event
    df.knowledgeable.simple <- vector(mode="list", length = num.acquisition.events)
    for(i in 1:num.acquisition.events){
      # who knows to solve simple:
      # all demonstrators as extracted above (demos.sub) minus those that solve complex
      # and all who have learned up to the acquistion event
      newly.learned.simple <- subset(learners.list_df$RING, learners.list_df$time_acq_adj<learners.list_df$time_acq_adj[i])

      # extract who has learned the simple solution thus far
      if(learned.behav == "dial"){
        dems <- dial.demos # extract those with at least 3 dial solves
      } else if(learned.behav=="slide"){
        dems <- slide.demos
      }


      knowledgeable.simple <- unique(c(dems, newly.learned.simple)) # add demos that knew how to solve simple before the start of the experiment and those that learned in this experiment
      knowledgeable.simple <- knowledgeable.simple[!(knowledgeable.simple%in%df.knowledgeable.complex[[i]])] # remove those that already know complex at the specific acquisition event
      df.knowledgeable.simple[[i]] <- knowledgeable.simple

    }
    # these lists are mutually exclusive - a bird who is a complex solver cannot be a simple solver at the same time and vice versa

    # now we create directed networks for each acquisition event
    # the first set of networks allows transmission from complex -> simple, but not from simple -> simple
    for(i in 1:num.acquisition.events){

      # we set everybody apart from knowledgeable simple birds to 0 for transmisson complex -> simple
       soc.sub <- assocMatrix.full[,,1, i] # for the first batch of matrices out of two, and for each acquisition event i
      which.rows <- which(rownames(net) %in% df.knowledgeable.simple[[i]])

      soc.sub[which.rows,] <- rep(0, length(soc.sub[,1])) # associations with everybody that solves simple are set to 0
      soc.sub[,which.rows] <- rep(0, length(soc.sub[,1]))

      # put it back in the correct slot of the ass Matrix
      assocMatrix.full[,,1, i] <- as.matrix(soc.sub)

    }

    # now we repeat it to only allow simple to simple, meaning that we set all associations with complex solvers to 0
    for(i in 1:num.acquisition.events){

      # for the second set of networks, we set everybody but complex knowledgeable birds to 0 for transmission complex -> simple
      soc.sub <- assocMatrix.full[,,2, i] # for the second batch of matrices out of two, and for each acquisition event i
      which.rows <- which(rownames(net) %in% df.knowledgeable.complex[[i]])
      soc.sub[which.rows, ] <- rep(0, length(soc.sub[,1])) #
      soc.sub[,which.rows ] <- rep(0, length(soc.sub[,1]))
      # put it back in the correct slot of the ass Matrix
      assocMatrix.full[,,2, i] <- as.matrix(soc.sub)

    }

    # reassign
    assocMatrix <- assocMatrix.full

    # we specify which matrix belongs to which acquisition evet
    assMatrixIndex <- rep(c(1:num.acquisition.events), times=2)


    nbdaData <- nbdaData(label=label,
                         assMatrix = assocMatrix,
                         asoc_ilv = get(paste("ILVs", label, sep="_")),
                         int_ilv = get(paste("ILVs", label, sep="_")),
                         multi_ilv = "ILVabsent",
                         orderAcq = learners.list_OAc,
                         timeAcq = learners.list_df$time_acq_adj,# get time in days
                         asocialTreatment = "timevarying",
                         demons = demo_vec,
                         assMatrixIndex =assMatrixIndex)

    # for testing
    # tfit <- tadaFit(nbdadata = nbdaData, type="social")
    # cbind.data.frame(tfit@outputPar, tfit@varNames)


  } else { # if we do not want to do directed learning


    # as we are running into convergence issues, we here have the option to create the nbda data objects without ILVs
    if(is.na(ILVs.yes.no)){
      nbdaData <- nbdaData(label=label,
                           assMatrix = assocMatrix,
                           asoc_ilv = "ILVabsent",
                           int_ilv = "ILVabsent",
                           multi_ilv = "ILVabsent",
                           orderAcq = learners.list_OAc,
                           timeAcq = learners.list_df$time_acq_adj,# get time in days
                           #  asocialTreatment = "timevarying",
                           demons = demo_vec)
    } else {
      # create NBDA Data Object

      nbdaData <- nbdaData(label=label,
                           assMatrix = assocMatrix,
                           asoc_ilv = get(paste("ILVs", label, sep="_")),
                           int_ilv = get(paste("ILVs", label, sep="_")),
                           multi_ilv = "ILVabsent",
                           orderAcq = learners.list_OAc,
                           timeAcq = learners.list_df$time_acq_adj,# get time in days
                           asocialTreatment = "timevarying",
                           demons = demo_vec)
    }


  }





  object <- NULL

  object$nbdadata <- nbdaData
  object$label <- label
  object$learners.list <- learners.list_df
  object$num.birds <- length(rownames(net))
  object$num.learners.list <- length(learners.list_df$RING)
  object$IDs <- rownames(net)

  return(object)

}





### Let's build NBDA data objects:




# for BB
# learning slide:
# both solving slide and solving complex are considered solving slide
# we only take individuals as demonstrators that have solved slide but not complex
# BB_complex_next_gen.slide.single.dem <- prepare.NBDA.data.next.gen(input.file = as.data.frame(df_solves),
#                                                           locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
#                                                           patch = "BB",
#                                                           start.date = as.POSIXct("14/01/2017", format="%d/%m/%Y"),
#                                                           end.date = as.POSIXct("12/02/2017", format="%d/%m/%Y"),
#                                                           label = "BB_complex_next_gen",
#                                                           ILVs.yes.no = "yes",
#                                                           min.time = "20170114074344",
#                                                           gmm = gmm.BB,
#                                                           ILVs.include = c("sex", "age"),
#                                                           learned.behav = "slide",
#                                                           demos.def = "single",
#                                                           directed=TRUE
# )

# learning dial:
# both solving dial and solving complex are considered solving dial
# we only take individuals as demonstrators that have solved dial, but not complex
# BB_complex_next_gen.dial.single.dem <- prepare.NBDA.data.next.gen(input.file = as.data.frame(df_solves),
#                                                         locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
#                                                         patch = "BB",
#                                                         start.date = as.POSIXct("14/01/2017", format="%d/%m/%Y"),
#                                                         end.date = as.POSIXct("12/02/2017", format="%d/%m/%Y"),
#                                                         label = "BB_complex_next_gen",
#                                                         ILVs.yes.no = "yes",
#                                                         min.time = "20170114074344",
#                                                         gmm = gmm.BB,
#                                                         ILVs.include = c("sex", "age"),
#                                                         learned.behav = "dial",
#                                                         demos.def = "single",
#                                                         directed=TRUE
# )


# 4.2.1 For dial and slide (undirected) -----------------------------------

# UNDIRECTED
# learning slide:
# both solving slide and solving complex are considered solving slide
# demos are those that have solved slide or complex
BB_complex_next_gen.slide.both.dem.undir <- prepare.NBDA.data.next.gen(input.file = as.data.frame(df_solves),
                                                                 locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
                                                                 patch = "BB",
                                                                 start.date = as.POSIXct("14/01/2017", format="%d/%m/%Y"),
                                                                 end.date = as.POSIXct("12/02/2017", format="%d/%m/%Y"),
                                                                 label = "BB_complex_next_gen",
                                                                 ILVs.yes.no = "yes",
                                                                 min.time = "20170114074344",
                                                                 gmm = gmm.BB,
                                                                 ILVs.include = c("sex", "age"),
                                                                 learned.behav = "slide",
                                                                 demos.def = "both",
                                                                 directed = FALSE
)

# learning dial:
# both solving dial and solving complex are considered solving dial
# demos are those that have solved dial or complex
BB_complex_next_gen.dial.both.dem.undir <- prepare.NBDA.data.next.gen(input.file = as.data.frame(df_solves),
                                                                locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
                                                                patch = "BB",
                                                                start.date = as.POSIXct("14/01/2017", format="%d/%m/%Y"),
                                                                end.date = as.POSIXct("12/02/2017", format="%d/%m/%Y"),
                                                                label = "BB_complex_next_gen",
                                                                ILVs.yes.no = "yes",
                                                                min.time = "20170114074344",
                                                                gmm = gmm.BB,
                                                                ILVs.include = c("sex", "age"),
                                                                learned.behav = "dial",
                                                                demos.def = "both",
                                                                directed = FALSE
)



# 4.2.2. For dial and slide - directed ------------------------------------



# learning slide:
# both solving slide and solving complex are considered solving slide
# demos are those that have solved slide or complex
BB_complex_next_gen.slide.both.dem.dir <- prepare.NBDA.data.next.gen(input.file = as.data.frame(df_solves),
                                                                   locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
                                                                   patch = "BB",
                                                                   start.date = as.POSIXct("14/01/2017", format="%d/%m/%Y"),
                                                                   end.date = as.POSIXct("12/02/2017", format="%d/%m/%Y"),
                                                                   label = "BB_complex_next_gen",
                                                                   ILVs.yes.no = "yes",
                                                                   min.time = "20170114074344",
                                                                   gmm = gmm.BB,
                                                                   ILVs.include = c("sex", "age"),
                                                                   learned.behav = "slide",
                                                                   demos.def = "both",
                                                                 directed = TRUE
)

# learning dial:
# both solving dial and solving complex are considered solving dial
# demos are those that have solved dial or complex
BB_complex_next_gen.dial.both.dem.dir <- prepare.NBDA.data.next.gen(input.file = as.data.frame(df_solves),
                                                                  locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
                                                                  patch = "BB",
                                                                  start.date = as.POSIXct("14/01/2017", format="%d/%m/%Y"),
                                                                  end.date = as.POSIXct("12/02/2017", format="%d/%m/%Y"),
                                                                  label = "BB_complex_next_gen",
                                                                  ILVs.yes.no = "yes",
                                                                  min.time = "20170114074344",
                                                                  gmm = gmm.BB,
                                                                  ILVs.include = c("sex", "age"),
                                                                  learned.behav = "dial",
                                                                  demos.def = "both",
                                                                directed = TRUE
)


# 4.2.3. for complex  -----------------------------------------

BB_complex_next_gen.complex <- prepare.NBDA.data.next.gen(input.file = as.data.frame(df_solves),
                                                                    locations <- c("1g", "1h", "3b", "3c", "3d", "1f", "3e", "3f", "3h", "4d", "4e", "4f"),
                                                                    patch = "BB",
                                                                    start.date = as.POSIXct("14/01/2017", format="%d/%m/%Y"),
                                                                    end.date = as.POSIXct("12/02/2017", format="%d/%m/%Y"),
                                                                    label = "BB_complex_next_gen",
                                                                    ILVs.yes.no = "yes",
                                                                    min.time = "20170114074344",
                                                                    gmm = gmm.BB,
                                                                    ILVs.include = c("pre.exp.one", "pre.exp.both"),
                                                                    learned.behav = "complex",
                                                                    demos.def = "both",
                                                                    directed = FALSE
)





# 4.3 Create summary ------------------------------------------------------



# extract how many birds in total and proporiton of learners.list
# num birds total
length(BB_complex_next_gen.dial.both.dem.dir$IDs)
# [1] 108 total
# num learners.list
length(BB_complex_next_gen.dial.both.dem.dir$learners.list$RING)
# dial = 7
length(BB_complex_next_gen.slide.both.dem.dir$learners.list$RING)
# slide = 9
length(BB_complex_next_gen.complex$learners.list$RING)
# complex =9

# extract the number of demonstrators
length(which(dial.demos %in% unique(c(BB_complex_next_gen.dial.both.dem.dir$IDs))))
length(which(slide.demos %in% unique(c(BB_complex_next_gen.slide.both.dem.dir$IDs))))

# extract how many demonstrator birds (that knew the complex task from the previous year) were part of this experiment
length(which(demos %in% unique(c(BB_complex_next_gen.complex$IDs))))
# [1] 11

# in %
9/(108-11)
# 9.3 % of naive birds have learned the complex task at BB

# 4.4. create constraints vector matrix -----------------------------------

# for the simple solves
constraintsVectMatrix <- create.constraints.Vect.Matrix(BB_complex_next_gen.slide.both.dem.undir$nbdadata, num_networks=1, num_ILVs=2)
ILVs <- c("sex", "age")
colnames(constraintsVectMatrix) <- c("network", paste("asoc", ILVs , sep="_"), paste("soc", ILVs, sep="_"))
constraintsVectMatrix <- as.matrix(constraintsVectMatrix)
class(constraintsVectMatrix)
head(constraintsVectMatrix)


# 4.5. Run cTADA ----------------------------------------------------------



# 4.5.1 Slide - undirected ------------------------------------------------

# both demos
AICtable_complex_next_gen.slide.both.dem.undir <- tadaAICtable(nbdadata = list(
  BB_complex_next_gen.slide.both.dem.undir$nbdadata
  #  CP_complex_next_gen$nbdadata
), constraintsVectMatrix = constraintsVectMatrix)

networksSupport(AICtable_complex_next_gen.slide.both.dem.undir)
# support numberOfModels
# 0 0.1034799              4
# 1 0.8965201             16

variableSupport(AICtable_complex_next_gen.slide.both.dem.undir)
# s1 ASOC:sex_BB_complex_next_gen ASOC:age_BB_complex_next_gen SOCIAL:sex_BB_complex_next_gen SOCIAL:age_BB_complex_next_gen
# support 0.8965201                   0.08471093                   0.08487136                      0.1211635                     0.08180491

MLE.complex.next.gen.slide.undir <- modelAverageEstimates(AICtable_complex_next_gen.slide.both.dem.undir, "median")
MLE.complex.next.gen.slide.undir

# ASOCIALsex_BB_complex_next_gen ASOCIALage_BB_complex_next_gen  SOCIALsex_BB_complex_next_gen  SOCIALage_BB_complex_next_gen
# 19.58441996                     0.03836913                    -0.01484823                     0.00000000                     0.00000000

# 4.5.2. Dial - undirected ------------------------------------------------

# both demos
AICtable_complex_next_gen.dial.both.dem.undir <- tadaAICtable(nbdadata = list(
  BB_complex_next_gen.dial.both.dem.undir$nbdadata
  #  CP_complex_next_gen$nbdadata
), constraintsVectMatrix = constraintsVectMatrix)

networksSupport(AICtable_complex_next_gen.dial.both.dem.undir)
# support numberOfModels
# 0 0.07562518              4
# 1 0.92437482             16

variableSupport(AICtable_complex_next_gen.dial.both.dem.undir)
#
# s1 ASOC:sex_BB_complex_next_gen ASOC:age_BB_complex_next_gen SOCIAL:sex_BB_complex_next_gen SOCIAL:age_BB_complex_next_gen
# support 0.9243748                   0.03237092                   0.03863927                     0.02600695                     0.04683006

MLE.complex.next.gen.dial.undir <- modelAverageEstimates(AICtable_complex_next_gen.dial.both.dem.undir, "median")
MLE.complex.next.gen.dial.undir

# ASOCIALsex_BB_complex_next_gen ASOCIALage_BB_complex_next_gen  SOCIALsex_BB_complex_next_gen  SOCIALage_BB_complex_next_gen
# 7.088402e+08                   7.887458e-03                  -1.989386e-02                   0.000000e+00                   0.000000e+00

# 4.5.3. Slide - directed -------------------------------------------------

# needs new constraintsVectMatrix, as we now have two networks
constraintsVectMatrix <- create.constraints.Vect.Matrix(BB_complex_next_gen.slide.both.dem.dir$nbdadata, num_networks=2, num_ILVs=2)
ILVs <- c("sex", "age")
colnames(constraintsVectMatrix) <- c(c("network_complex_simple", "network_simple_simple"), paste("asoc", ILVs , sep="_"), paste("soc", ILVs, sep="_"))
constraintsVectMatrix <- as.matrix(constraintsVectMatrix)
class(constraintsVectMatrix)
head(constraintsVectMatrix)



AICtable_complex_next_gen.slide.both.dem.dir <- tadaAICtable(nbdadata = list(
  BB_complex_next_gen.slide.both.dem.dir$nbdadata
  #  CP_complex_next_gen$nbdadata
), constraintsVectMatrix = constraintsVectMatrix)

networksSupport(AICtable_complex_next_gen.slide.both.dem.dir)

# support numberOfModels
# 0:0 0.22854508              4
# 0:1 0.15185138             16
# 1:0 0.54263912             16
# 1:2 0.07696441             16

variableSupport(AICtable_complex_next_gen.slide.both.dem.dir)

# s1        s2 ASOC:sex_BB_complex_next_gen ASOC:age_BB_complex_next_gen SOCIAL:sex_BB_complex_next_gen SOCIAL:age_BB_complex_next_gen
# support 0.6196035 0.2288158                    0.1097166                   0.09531642                      0.1272142                     0.05823525

MLE.complex.next.gen.slide.dir <- modelAverageEstimates(AICtable_complex_next_gen.slide.both.dem.dir, "median")
MLE.complex.next.gen.slide.dir

# #                            s1                             s2 ASOCIALsex_BB_complex_next_gen ASOCIALage_BB_complex_next_gen  SOCIALsex_BB_complex_next_gen
# 8.70101363                     2.88816630                     0.08092590                    -0.03131701                     0.00000000
# SOCIALage_BB_complex_next_gen
# 0.00000000 _complex_next_gen


# 4.5.4. Dial - directed --------------------------------------------------

AICtable_complex_next_gen.dial.both.dem.dir <- tadaAICtable(nbdadata = list(
  BB_complex_next_gen.dial.both.dem.dir$nbdadata
  #  CP_complex_next_gen$nbdadata
), constraintsVectMatrix = constraintsVectMatrix)

networksSupport(AICtable_complex_next_gen.dial.both.dem.dir)



#        support numberOfModels
# support numberOfModels
# 0:0 0.14317378              4
# 0:1 0.64428538             16
# 1:0 0.17068303             16
# 1:2 0.04185781             16

variableSupport(AICtable_complex_next_gen.dial.both.dem.dir)

# s1        s2 ASOC:sex_BB_complex_next_gen ASOC:age_BB_complex_next_gen SOCIAL:sex_BB_complex_next_gen SOCIAL:age_BB_complex_next_gen
# support 0.2125408 0.6861432                   0.03828958                   0.06026423                     0.02100754                     0.05884849

MLE.complex.next.gen.dial.dir <- modelAverageEstimates(AICtable_complex_next_gen.dial.both.dem.dir, "median")
MLE.complex.next.gen.dial.dir

# s1                             s2 ASOCIALsex_BB_complex_next_gen ASOCIALage_BB_complex_next_gen  SOCIALsex_BB_complex_next_gen
# 4.118573e+07                   6.695483e+08                   1.514526e-02                  -3.819959e-02                   0.000000e+00
# SOCIALage_BB_complex_next_gen
# 0.000000e+00
# 4.5.5 Complex -----------------------------------------------------------


constraintsVectMatrix <- create.constraints.Vect.Matrix(BB_complex_next_gen.complex$nbdadata, num_networks=1, num_ILVs=2)
ILVs <- c("prev_experience1", "prev_experience2")
colnames(constraintsVectMatrix) <- c("network", paste("asoc", ILVs , sep="_"), paste("soc", ILVs, sep="_"))
constraintsVectMatrix <- as.matrix(constraintsVectMatrix)
class(constraintsVectMatrix)
head(constraintsVectMatrix)

AICtable_complex_next_gen.complex <- tadaAICtable(nbdadata = list(
  BB_complex_next_gen.complex$nbdadata
  #  CP_complex_next_gen$nbdadata
), constraintsVectMatrix = constraintsVectMatrix)

AICtable_complex_next_gen.complex@printTable

# the top model is an asocial model (19), followed by a social model (15)

setwd("C:/Users/swild/Desktop/Konstanz/Collaborations/Lucy, Michael - Phil Trans cumulative culture/Analysis/NBDA output/Results COMPLEX NEXT GEN diffusion")
# save(AICtable_complex_next_gen, file="AICTable_complex_next_gen_no_start_vals.RData")
# load("AICTable_complex_next_gen_no_start_vals.RData")

network.support.complex.next_gen <- networksSupport(AICtable_complex_next_gen.complex)
network.support.complex.next_gen
write.table(network.support.complex.next_gen, "network.support.complex.next_gen.diffusion_no_start_vals.txt")

# support numberOfModels
# 0 0.6388354              4
# 1 0.3611646             16

type.support.complex.next_gen <- typeSupport(AICtable_complex_next_gen.complex)
type.support.complex.next_gen

variable.support.complex.next_gen <- variableSupport(AICtable_complex_next_gen.complex)
variable.support.complex.next_gen
write.table(variable.support.complex.next_gen, "variable.support.complex.next_gen.diffusion_no_start_vals.txt")

# s1 ASOC:pre.exp.one_BB_complex_next_gen ASOC:pre.exp.both_BB_complex_next_gen SOCIAL:pre.exp.one_BB_complex_next_gen SOCIAL:pre.exp.both_BB_complex_next_gen
# support 0.3611646                            0.8858901                            0.06914529                              0.1008408                              0.01183816

# extract model averaged estimates
MLE.complex.next_gen  <- modelAverageEstimates(AICtable_complex_next_gen.complex,averageType = "median")
MLE.complex.next_gen
write.table(MLE.complex.next_gen, "MLE.complex.next_gen.diffusion_no_start_vals.txt")

# ASOCIALpre.exp.one_BB_complex_next_gen ASOCIALpre.exp.both_BB_complex_next_gen   SOCIALpre.exp.one_BB_complex_next_gen  SOCIALpre.exp.both_BB_complex_next_gen
# 0.000000                                2.957669                                0.000000                                0.000000                                0.000000


# 4.6. Extract effect sizes -----------------------------------------------



# 4.6.1. dial - undirected -------------------------------------------------------------


AICtable_complex_next_gen.dial.both.dem.undir@printTable
# top model is 16

# re-run constriantsVectMatrix briefly
constraintsVectMatrix <- create.constraints.Vect.Matrix(BB_complex_next_gen.dial.both.dem.undir$nbdadata, num_networks=1, num_ILVs=2)
constraintsVectMatrix[16,]

# create constrained NBDA data objects
constrained.next_gen.BB.dial.undir <- constrainedNBDAdata(BB_complex_next_gen.dial.both.dem.undir$nbdadata, constraintsVect = constraintsVectMatrix[16,])

# run TADA on the best model
bestModel.next_gen.dial.undir <- tadaFit(nbdadata = list(
  constrained.next_gen.BB.dial.undir
), type="social")

# extract the percentage of birds learned through social learning - by event
prop.solve.social.byevent.next.gen.dial.undir <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.next_gen.BB.dial.undir), model = bestModel.next_gen.dial.undir)
prop.solve.social.byevent.next.gen.dial.undir

# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.next.gen.dial.undir <- oadaPropSolveByST( nbdadata = list(
  constrained.next_gen.BB.dial.undir
), model=bestModel.next_gen.dial.undir)
prop.solve.social.next.gen.dial.undir

# 100% of birds were estimated to have learned the complex task through social learning

plotProfLik(which=1, model=bestModel.next_gen.dial.undir, range=c(0,50), resolution = 10)
CI <- profLikCI(which=1, model=bestModel.next_gen.dial.undir, lowerRange = c(0,10))

# extract percentages of lower bound

bestModelDataS1LowerBound.complex.2nd.dial.undir <- constrainedNBDAdata(
  nbdadata =
    BB_complex_next_gen.dial.both.dem.undir$nbdadata,
  constraintsVect = constraintsVectMatrix[16, ],
  offset = c(CI[1]   , rep(0, 4))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1LowerBound.complex_dial <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.complex.2nd.dial.undir

    ) ,
    type = "asocial"
  )
bestModelS1LowerBound.complex_dial@outputPar
# [1] 322.5608

#Now plug into the prop solve function
prop.solve.social.Lower.complex.dial <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound.complex_dial,
    nbdadata = list(
      bestModelDataS1LowerBound.complex.2nd.dial.undir
    )
  )
prop.solve.social.Lower.complex.dial
# Lower bound for % of birds having learned the complex task through social learning is 37.8%


# 4.6.2. slide - undirected -----------------------------------------------
AICtable_complex_next_gen.slide.both.dem.undir@printTable
# top model is 16

# create constrained NBDA data objects
constrained.next_gen.BB.slide.undir <- constrainedNBDAdata(BB_complex_next_gen.slide.both.dem.undir$nbdadata, constraintsVect = constraintsVectMatrix[16,])

# run TADA on the best model
bestModel.next_gen.slide.undir <- tadaFit(nbdadata = list(
  constrained.next_gen.BB.slide.undir
), type="social")

# extract the percentage of birds learned through social learning - by event
prop.solve.social.byevent.next.gen.slide.undir <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.next_gen.BB.slide.undir), model = bestModel.next_gen.slide.undir)
prop.solve.social.byevent.next.gen.slide.undir

# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.next.gen.slide.undir <- oadaPropSolveByST( nbdadata = list(
  constrained.next_gen.BB.slide.undir
), model=bestModel.next_gen.slide.undir)
prop.solve.social.next.gen.slide.undir

# 78.0% of birds were estimated to have learned the complex task through social learning

plotProfLik(which=1, model=bestModel.next_gen.slide.undir, range=c(0,500), resolution = 10)
CI <- profLikCI(which=1, model=bestModel.next_gen.slide.undir, lowerRange = c(0,10), upperRange = c(450,500))

# extract percentages of lower bound

bestModelDataS1LowerBound.complex.2nd.slide.undir <- constrainedNBDAdata(
  nbdadata =
    BB_complex_next_gen.slide.both.dem.undir$nbdadata,
  constraintsVect = constraintsVectMatrix[16, ],
  offset = c(CI[1]   , rep(0, 4))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1LowerBound.complex_slide <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.complex.2nd.slide.undir

    ) ,
    type = "asocial"
  )
bestModelS1LowerBound.complex_slide@outputPar
# [1] 192.425

#Now plug into the prop solve function
prop.solve.social.Lower.complex.slide <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound.complex_slide,
    nbdadata = list(
      bestModelDataS1LowerBound.complex.2nd.slide.undir
    )
  )
prop.solve.social.Lower.complex.slide
# Lower bound for % of birds having learned the complex task through social learning is 46.9%
# extract percentages of upper bound

bestModelDataS1upperBound.complex.2nd.slide.undir <- constrainedNBDAdata(
  nbdadata =
    BB_complex_next_gen.slide.both.dem.undir$nbdadata,
  constraintsVect = constraintsVectMatrix[16, ],
  offset = c(CI[2]   , rep(0, 4))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1upperBound.complex_slide <-
  tadaFit(
    list(
      bestModelDataS1upperBound.complex.2nd.slide.undir

    ) ,
    type = "asocial"
  )
bestModelS1upperBound.complex_slide@outputPar
# [1] 10888.16

#Now plug into the prop solve function
prop.solve.social.upper.complex.slide <-
  oadaPropSolveByST(
    model = bestModelS1upperBound.complex_slide,
    nbdadata = list(
      bestModelDataS1upperBound.complex.2nd.slide.undir
    )
  )
prop.solve.social.upper.complex.slide
# upper bound for % of birds having learned the complex task through social learning is 88.2%




# 4.6.3. Dial - directed --------------------------------------------------

AICtable_complex_next_gen.dial.both.dem.dir@printTable
# top model is 47

# re-run constriantsVectMatrix briefly
constraintsVectMatrix <- create.constraints.Vect.Matrix(BB_complex_next_gen.dial.both.dem.dir$nbdadata, num_networks=2, num_ILVs=2)
constraintsVectMatrix[47,]

# create constrained NBDA data objects
constrained.next_gen.BB.dial.dir <- constrainedNBDAdata(BB_complex_next_gen.dial.both.dem.dir$nbdadata, constraintsVect = constraintsVectMatrix[47,])

# run TADA on the best model
bestModel.next_gen.dial.dir <- tadaFit(nbdadata = list(
  constrained.next_gen.BB.dial.dir
), type="social")

# extract the percentage of birds learned through social learning - by event
prop.solve.social.byevent.next.gen.dial.dir <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.next_gen.BB.dial.dir), model = bestModel.next_gen.dial.dir)
prop.solve.social.byevent.next.gen.dial.dir

# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.next.gen.dial.dir <- oadaPropSolveByST( nbdadata = list(
  constrained.next_gen.BB.dial.dir
), model=bestModel.next_gen.dial.dir)
prop.solve.social.next.gen.dial.dir

# 100% of birds were estimated to have learned the dial task through social learning from others solving dial

plotProfLik(which=1, model=bestModel.next_gen.dial.dir, range=c(0,50), resolution = 10)
CI <- profLikCI(which=1, model=bestModel.next_gen.dial.dir, lowerRange = c(0,10))

# extract percentages of lower bound

bestModelDataS1LowerBound.complex.2nd.dial.dir <- constrainedNBDAdata(
  nbdadata =
    BB_complex_next_gen.dial.both.dem.dir$nbdadata,
  constraintsVect = constraintsVectMatrix[47, ],
  offset = c(CI[1]   , rep(0, 5))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1LowerBound.complex_dial <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.complex.2nd.dial.dir

    ) ,
    type = "asocial"
  )
bestModelS1LowerBound.complex_dial@outputPar
# [1] 238.6029

#Now plug into the prop solve function
prop.solve.social.Lower.complex.dial <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound.complex_dial,
    nbdadata = list(
      bestModelDataS1LowerBound.complex.2nd.dial.dir
    )
  )
prop.solve.social.Lower.complex.dial
# Lower bound for % of birds having learned the complex task through social learning is 54.9%


# 4.6.4. Slide - directed -------------------------------------------------

AICtable_complex_next_gen.slide.both.dem.dir@printTable
# top model is 48

# create constrained NBDA data objects
constrained.next_gen.BB.slide.dir <- constrainedNBDAdata(BB_complex_next_gen.slide.both.dem.dir$nbdadata, constraintsVect = constraintsVectMatrix[48,])

# run TADA on the best model
bestModel.next_gen.slide.dir <- tadaFit(nbdadata = list(
  constrained.next_gen.BB.slide.dir
), type="social")

# extract the percentage of birds learned through social learning - by event
prop.solve.social.byevent.next.gen.slide.dir <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.next_gen.BB.slide.dir), model = bestModel.next_gen.slide.dir)
prop.solve.social.byevent.next.gen.slide.dir

# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.next.gen.slide.dir <- oadaPropSolveByST( nbdadata = list(
  constrained.next_gen.BB.slide.dir
), model=bestModel.next_gen.slide.dir)
prop.solve.social.next.gen.slide.dir

# 59.8% of birds were estimated to have learned the complex task through social learning

plotProfLik(which=1, model=bestModel.next_gen.slide.dir, range=c(0,400), resolution = 10)
CI <- profLikCI(which=1, model=bestModel.next_gen.slide.dir, lowerRange = c(0,10), upperRange = c(150, 250))

# extract percentages of lower bound

bestModelDataS1LowerBound.complex.2nd.slide.dir <- constrainedNBDAdata(
  nbdadata =
    BB_complex_next_gen.slide.both.dem.dir$nbdadata,
  constraintsVect = constraintsVectMatrix[48, ],
  offset = c(CI[1]   , rep(0, 5))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1LowerBound.complex_slide <-
  tadaFit(
    list(
      bestModelDataS1LowerBound.complex.2nd.slide.dir

    ) ,
    type = "asocial"
  )
bestModelS1LowerBound.complex_slide@outputPar
# [1] 157.705

#Now plug into the prop solve function
prop.solve.social.Lower.complex.slide <-
  oadaPropSolveByST(
    model = bestModelS1LowerBound.complex_slide,
    nbdadata = list(
      bestModelDataS1LowerBound.complex.2nd.slide.dir
    )
  )
prop.solve.social.Lower.complex.slide
# Lower bound for % of birds having learned the complex task through social learning is 16.2.9%

# same for upper bound
bestModelDataS1upperBound.complex.2nd.slide.dir <- constrainedNBDAdata(
  nbdadata =
    BB_complex_next_gen.slide.both.dem.dir$nbdadata,
  constraintsVect = constraintsVectMatrix[48, ],
  offset = c(CI[2]   , rep(0, 5))
)


#Now, when we fit an "asocial" model it constrains the value of s1=0, but then an offset is added to s
bestModelS1upperBound.complex_slide <-
  tadaFit(
    list(
      bestModelDataS1upperBound.complex.2nd.slide.dir

    ) ,
    type = "asocial"
  )
bestModelS1upperBound.complex_slide@outputPar
# [1] 3639.354

#Now plug into the prop solve function
prop.solve.social.upper.complex.slide <-
  oadaPropSolveByST(
    model = bestModelS1upperBound.complex_slide,
    nbdadata = list(
      bestModelDataS1upperBound.complex.2nd.slide.dir
    )
  )
prop.solve.social.upper.complex.slide
# upper bound for % of birds having learned the complex task through social learning is 84.8%

# 4.6.5 Complex -----------------------------------------------------------
AICtable_complex_next_gen.complex@printTable
# top model is 19 (asocial), followed by 15 (social)
constraintsVectMatrix <- create.constraints.Vect.Matrix(BB_complex_next_gen.complex$nbdadata, num_networks=1, num_ILVs=2)

# the top model is model 19 - as we cannot use the constrained NBDA data function on an asocial model
# we use the corresponding social model (15) and then set the 'type' to asocial
constraintsVectMatrix[19,]
constraintsVectMatrix[15,]

# create constrained NBDA data objects
constrained.next_gen.BB <- constrainedNBDAdata(BB_complex_next_gen.complex$nbdadata, constraintsVect = constraintsVectMatrix[15,])
# constrained.next_gen.CP <- constrainedNBDAdata(CP_complex_next_gen$nbdadata, constraintsVect = constraintsVectMatrix[15,])

# run TADA on the best model
bestModel.next_gen <- tadaFit(nbdadata = list(
  constrained.next_gen.BB
), type="social")

bestModel.next_gen@varNames
bestModel.next_gen@outputPar

# then averaged across all acquisition events
prop.solve.social.next.gen.complex <- oadaPropSolveByST( nbdadata = list(
  constrained.next_gen.BB
), model=bestModel.next_gen)
prop.solve.social.next.gen.complex

# 44.0%

# extracting CIs for s (conditional on second best model)
plotProfLik(which=1, model=bestModel.next_gen, range=c(0, 200), resolution = 10)
# upper limit cannot be obtained due to convergence issues - lower limit is 0%
# profLikCI(which=1,model=bestModel.next_gen, lowerRange=c(1,2), upperRange=c(3,5))

# rerun the best model, this time asocial
bestModel.next_gen <- tadaFit(nbdadata = list(
  constrained.next_gen.BB
), type="asocial")

# now we extract effect sizes for previous experience of one component on asocial learning (which=1)
plotProfLik(which=1, model=bestModel.next_gen, range=c(0, 8), resolution = 10)
profLikCI(which=1,model=bestModel.next_gen, lowerRange=c(1,3), upperRange=c(4,6))
# Lower CI Upper CI
# 1.625241 4.513349

# back transform the estimates
exp(c(MLE.complex.next_gen[2], 1.625241, 4.513349 ))
# birds with previous experience of one component were 19.3 [5.1-91.2] times faster at learning the complex solution asocially
# compared to those without any previous experience







## we can now extract the social learning parameter s for the best performing social model (model 15)
# which overall is model 2
constrained.next_gen.BB <- constrainedNBDAdata(BB_complex_next_gen$nbdadata, constraintsVect = constraintsVectMatrix[15,])

# run TADA on the best model
bestModel.next_gen <- tadaFit(nbdadata = list(
  constrained.next_gen.BB
), type="social")

# extract the percentage of birds learned through social learning - by event
prop.solve.social.byevent.next.gen <- oadaPropSolveByST.byevent( nbdadata = list(
  constrained.next_gen.BB), model = bestModel.next_gen)
prop.solve.social.byevent.next.gen

# for each acquisition event, we now have an estimate for the likelihood of social learning  (P network)

# then averaged across all acquisition events
prop.solve.social.next.gen <- oadaPropSolveByST( nbdadata = list(
  constrained.next_gen.BB
), model=bestModel.next_gen)
prop.solve.social.next.gen

# 44.3% of birds were estimated to have learned the complex task through social learning


# plot profile likelihood for s to extract confidence intervals
plotProfLik(which=1, model=bestModel.next_gen, range=c(0, 100), resolution = 30)
# the lower confidence interval clearly spans 0, while on the upper level,
# the zigzag lines indicate convergence issues (which is also apparent from the convergence column in the AIC table)
# which means the upper confidence cannot be reliably obtained


# 5. Plot BB networks -----------------------------------------------------


# 5.1. get data ready for plotting ----------------------------------------

# dial diffusion
dial.net <- as.matrix(BB_NBDA_DATA$nbdadata@assMatrix[,,1,1])
rownames(dial.net) <- BB_NBDA_DATA$IDs
colnames(dial.net) <- BB_NBDA_DATA$IDs

demos.dial <- BB_NBDA_DATA$IDs[which(BB_NBDA_DATA$nbdadata@demons==1)]
learners.list.dial <- BB_NBDA_DATA$learners

# complex 1st gen
complex1st.net <- as.matrix(BB_complex_progressive_complex$nbdadata@assMatrix[,,1,1])
rownames(complex1st.net) <- BB_complex_progressive_complex$IDs
colnames(complex1st.net) <- BB_complex_progressive_complex$IDs

learners.list.complex1st <- BB_complex_progressive_complex$learners
# extract a list of birds that managed to solve either dial or slide, but not complex
count <- 1
learners.complex1st.simple <- NULL
for(i in rownames(complex1st.net)){
  i.data <- as.data.frame(subset(df_solves, df_solves$RING==i & df_solves$PATCH=="BB" & df_solves$experiment=="complex_prog" & df_solves$solution_category %in% c("slide", "dial")))
  i.solves.simple <- length(i.data$solution_category)
  if(i.solves.simple>=3){
    learners.complex1st.simple[count] <- i
    count <- count+1
  }
}
# exclude all those who solve complex
learners.complex1st.simple <- subset(learners.complex1st.simple, !(learners.complex1st.simple %in% learners.list.complex1st$RING))


# complex 2nd gen
complex2nd.net <- as.matrix(BB_complex_next_gen.complex$nbdadata@assMatrix[,,1,1])
rownames(complex2nd.net) <- BB_complex_next_gen.complex$IDs
colnames(complex2nd.net) <- BB_complex_next_gen.complex$IDs

learners.list.complex2nd <- BB_complex_next_gen.complex$learners
demos.complex2nd <- BB_complex_next_gen.complex$IDs[which(BB_complex_next_gen.complex$nbdadata@demons==1)]
demos.complex2nd <- demos.complex2nd[-1]
# extract a list of birds that managed to solve either dial or slide, but not complex
count <- 1
learners.complex2nd.simple <- NULL
for(i in rownames(complex2nd.net)){
  i.data <- as.data.frame(subset(df_solves, df_solves$RING==i & df_solves$PATCH=="BB" & df_solves$experiment=="complex_prog" & df_solves$solution_category %in% c("slide", "dial")))
  i.solves.simple <- length(i.data$solution_category)
  if(i.solves.simple>=3){
    learners.complex2nd.simple[count] <- i
    count <- count+1
  }
}
# exclude all those who solve complex
learners.complex2nd.simple <- subset(learners.complex2nd.simple, !(learners.complex2nd.simple %in% c(learners.list.complex2nd$RING, demos.complex2nd)))

# 5.2. plot networks ------------------------------------------------------

library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(igraph)
#install.packages("sna")
#install.packages("network")


plot.network <- function(net, demos, learners.list, experiment, n, solves, diffusion, position, simple.list){
  # create a colour vector that distinguishes between learners.list, naive birds and demonstrators
  cat.all <- NULL
  for (i in rownames(net)){
    if(i %in% demos){
      cat <- "#F0DC6A"
    } else if(i %in% simple.list){
      cat <- "#33CCCC"  # blue for simple learners
    } else if (i %in% learners.list$RING){
      cat <- "#660099"  # purple for complex learners
    }
    else {cat <- "#FFFFFF"} #white for non-learners
    cat.all[which(rownames(net)==i)] <- cat
  }

  # create graph object
  g.net <- graph_from_adjacency_matrix(net, mode = "undirected",
                                        weighted = TRUE, diag = FALSE)
  V(g.net)$colour <- cat.all

  E(g.net)$width <- E(g.net)$weight
  E(g.net)$edge.color <- "gray50"
  # remove edge weights below 0.005
  g.net <- delete_edges(g.net, E(g.net)[weight<0.05])
  # layout fruchtermann reingold
  l <- layout_with_fr(g.net )
  title <-  paste(experiment, paste("(N=", n, ")", sep =
                                    ""), sep = " ")


  plot( g.net,
        vertex.size = 6,
        edge.curved = 0.2,
        edge.color =  E(g.net)$edge.color,
        vertex.color = V(g.net)$colour,
        vertex.label = NA,
        vertex.frame.colour = "black",
        edge.width = E(g.net)$width*10,
        frame = TRUE,
        layout=l,
     #   main= title,
        adj = c(0,-1),
     margin=c(0.02,0.0,0.0,0.0),
     asp = 1,
     rescale = TRUE
        )

  title(title, line=0.5, adj=0.0, cex.main=1.5, font =1, family = "sans")
if(position!="none"){
  legend(x=-1.15, y=-0.7, c("naive","simple", "demonstrator", "complex"), pch=21,

         col="#777777", pt.bg=c("#FFFFFF", "#33CCCC", "#F0DC6A", "#660099"), pt.cex=1.5, cex=1.2, bty="n", ncol=1,
         y.intersp=0.8)
}


}
set.seed(10)

pdf( "networks.pdf", width=12, height=4)
par(mfrow=c(1,3), mai = c(0.0, 0.0, 0.2, 0.07))


p1 <-
  plot.network(
    net = dial.net,
    demos = demos.dial,
    learners.list = learners.list.dial,
    n = BB_NBDA_DATA$num.birds,
    experiment = "A) Dial diffusion",
    solves = "dial",
    diffusion = "dial_diffusion",
    position = "left",
    simple.list = learners.list.dial$RING
  )


p2 <-
  plot.network(
    net = complex1st.net,
    demos = NA,
    learners.list = learners.list.complex1st,
    n = BB_complex_progressive_complex$num.birds,
    experiment = "B) Complex 1st generation",
    solves = "complex",
    diffusion = "complex_prog",
    position = "none",
    simple.list = learners.complex1st.simple
  )

p3 <-
  plot.network(
    net = complex2nd.net,
    demos = demos.complex2nd,
    learners.list = learners.list.complex2nd,
    n = BB_complex_next_gen.complex$num.birds,
    experiment = "C) Complex 2nd generation",
    solves = "complex",
    diffusion = "nextgen",
    position = "none",
    simple.list = learners.complex2nd.simple
  )


dev.off()

