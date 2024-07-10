library(diagram)
library(lubridate)
library(seraphim)

# get path of the script
if (rstudioapi::isAvailable()){
  if (require('rstudioapi') != TRUE){
    install.packages('rstudioapi', repos = r)
  } else {
    library(rstudioapi)
  }
  wd <- dirname(getActiveDocumentContext()$path)
} else {
  wd <- getwd()
}

setwd(wd)

segments = c("L","M", "nss")
nberOfExtractionFiles = 1000

mostRecentSamplingDates = c(decimal_date(ymd("2022-09-12")),
                            decimal_date(ymd("2022-09-12")),
                            decimal_date(ymd("2022-08-25")))

# 1. Preparing all the different environmental factors (rasters) to test

template = raster("All_rasters_1/Croplands.asc")
template[!is.na(template[])] = 1; template[is.na(template[])] = 0
template = raster::aggregate(template, 10, fun=mean)
template[template[]<0.5] = NA; template[template[]>=0.5] = 0
rasterNames = list.files("All_rasters_1")
for (i in 1:length(rasterNames))
	{
		r1 = raster(paste0("All_rasters_1/",rasterNames[i]))
		r2 = raster::aggregate(r1, 10, fun=mean); r2[is.na(template[])] = NA
		print(r2)
		writeRaster(r2, paste0("All_rasters_2/",rasterNames[i]), overwrite=T)
	}

# 2. Extracting spatio-temporal information embedded in 1000 posterior trees

for (h in 1:length(segments))
	{
		localTreesDirectory = paste0("RVFV-",segments[h],"_exts")
		allTrees = scan(paste0(segments[h],"-lineage-C.trees"), 
		                what="", sep="\n", quiet=TRUE)
		burnIn = 1001 # 10% of sampled trees
		randomSampling = FALSE
		nberOfTreesToSample = nberOfExtractionFiles
		mostRecentSamplingDatum = mostRecentSamplingDates[h]
		coordinateAttributeName = "location"
		nberOfCores = 10
		treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample,
						mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
	}

# 3. Investigating the impact of environmental factors on lineage dispersal velocity

for (h in 1:length(segments))
	{
		localTreesDirectory = paste0("RVFV-",segments[h],"_exts")
		pathModel = 3
		fourCells = FALSE
		nberOfRandomisations = 0
		randomProcedure = 3
		outputName = paste0("RVFV-",segments[h],"_CS_preliminary")
		showingPlots = FALSE
		nberOfCores = 10
		OS = "Unix"
		c = 0
		envVariables = list()
		resistances = list()
		avgResistances = list()
		rasterNames = list.files("All_rasters_2")
		for (k in c(10,100,1000))
			{
				for (i in 1:length(rasterNames))
					{
						c = c+1
						rast = raster(paste("All_rasters_2/",rasterNames[i],sep=""))
						rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(gsub(".asc","",rasterNames[i]),"_k",k,sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(gsub(".asc","",rasterNames[i]),"_k",k,sep="")
						resistances[[c]] = TRUE; avgResistances[[c]] = TRUE
					}
				for (i in 1:length(rasterNames))
					{
						c = c+1
						rast = raster(paste("All_rasters_2/",rasterNames[i],sep=""))
						rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
						names(rast) = paste(gsub(".asc","",rasterNames[i]),"_k",k,sep="")
						envVariables[[c]] = rast; names(envVariables[[c]]) = paste(gsub(".asc","",rasterNames[i]),"_k",k,sep="")
						resistances[[c]] = FALSE; avgResistances[[c]] = FALSE
					}
			}
		spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances,
					  fourCells, nberOfRandomisations, randomProcedure, outputName, showingPlots, nberOfCores, OS)
		tab = read.table(paste0("Seraphim_res/RVFV-",segments[h],"_CS_preliminary_linear_regression_results.txt"), head=T)
		Qs_tab = c(); res = c("R","C"); rasterNames = gsub("-",".",rasterNames)
		for (i in 1:length(rasterNames))
			{
				for (j in 1:length(res))
					{
						for (k in c(10,100,1000))
							{
								Qs = tab[,which(colnames(tab)==paste0("Univariate_LR_delta_R2_",gsub(".asc","",rasterNames[i]),"_k",k,"_",res[j]))]
								Qs_tab = rbind(Qs_tab, cbind(round(sum(Qs>0)/length(Qs),2), paste0("Univariate_LR_delta_R2_",gsub(".asc","",rasterNames[i]),"_k",k,"_",res[j])))
								if ((sum(Qs>0)/length(Qs)) >= 0.9) print(paste0("Univariate_LR_delta_R2_",gsub(".asc","",rasterNames[i]),"_k",k,"_",res[j]))
							}
					}
			}
		selected_variables = Qs_tab[which(as.numeric(Qs_tab[,1])>0.9),2]
		c = 0; envVariables = list(); resistances = list(); avgResistances = list()
		for (i in 1:length(selected_variables))
			{
				rasterName = gsub("Univariate_LR_delta_R2_","",selected_variables[i])
				k = as.numeric(gsub("k","",unlist(strsplit(rasterName,"_"))[length(unlist(strsplit(rasterName,"_")))-1]))
				r = unlist(strsplit(rasterName,"_"))[length(unlist(strsplit(rasterName,"_")))]
				rasterName = unlist(strsplit(rasterName,"_k"))[1]; c = c+1
				rast = raster(paste0("All_rasters_2/",rasterName,".asc"))
				rast[rast[]<0] = 0; M = max(rast[], na.rm=T); rast[] = (rast[]*(k/M))+1
				names(rast) = paste0(rasterName,"_k",k)
				envVariables[[c]] = rast; names(envVariables[[c]]) = paste(rasterName,"_k",k,sep="")
				if (r == "R") { resistances[[c]] = TRUE; avgResistances[[c]] = TRUE }
				if (r == "C") { resistances[[c]] = FALSE; avgResistances[[c]] = FALSE }
			}
		nberOfRandomisations = 1
		outputName = paste0("RVFV-",segments[h],"_CS_randomisation")
		spreadFactors(localTreesDirectory, nberOfExtractionFiles, envVariables, pathModel, resistances, avgResistances,
					  fourCells, nberOfRandomisations, randomProcedure, outputName, showingPlots, nberOfCores, OS)
		seraphim_results_1 = read.table(paste0("Seraphim_res/RVFV-",segments[h],"_CS_randomisation_linear_regression_results.txt"), header=T)
		seraphim_results_2 = read.table(paste0("Seraphim_res/RVFV-",segments[h],"_CS_randomisation_Bayes_factor_supports.txt"), header=T)
		allResults = matrix(nrow=length(selected_variables), ncol=5)
		colnames(allResults) = c("environmental factor","regression coefficient","Q statistic","p(Q) > 0","BF")
		for (i in 1:length(selected_variables))
			{
				rasterName = gsub("Univariate_LR_delta_R2_","",selected_variables[i]); allResults[i,1] = rasterName
				index1 = which(grepl("LR_R2",colnames(seraphim_results_1))&grepl(rasterName,colnames(seraphim_results_1)))
				index2 = which(grepl("delta_R2",colnames(seraphim_results_1))&grepl(rasterName,colnames(seraphim_results_1)))
				R2 = seraphim_results_1[,index1]; Qe = seraphim_results_1[,index2]
				allResults[i,2] = paste0(round(median(R2),3)," [",round(hdi(R2)[1:2],3)[1],", ",round(hdi(R2)[1:2],3)[2],"]")
				allResults[i,3] = paste0(round(median(Qe),3)," [",round(hdi(Qe)[1:2],3)[1],", ",round(hdi(Qe)[1:2],3)[2],"]")
				allResults[i,4] = sum(Qe>0)/length(Qe)
				index3 = which(grepl(rasterName,rownames(seraphim_results_2)))
				allResults[i,5] = round(seraphim_results_2[index3,1],1)
			}
		write.table(allResults, paste0("RVFV-",segments[h],"_1.csv"), row.names=F, quote=F, sep=";")
	}

# 4. Investigating the impact of environmental factors on lineage dispersal location

for (h in 1:length(segments))
	{
		localTreesDirectory = paste0("RVFV-",segments[h],"_exts")
		envVariables = list(raster("All_rasters_2/Croplands.asc")); randomProcedure = 3; nberOfCores = 5
		treesRandomisation(localTreesDirectory, nberOfExtractionFiles, envVariables, randomProcedure, nberOfCores)
		envVariables = list(); rasterNames = list.files("All_rasters_2")
		for (i in 1:length(rasterNames)) envVariables[[i]] = raster(paste0("All_rasters_2/",rasterNames[i]))
		for (i in 1:nberOfExtractionFiles)
			{
				obs = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), header=T)
				ran = read.csv(paste0(localTreesDirectory,"/TreeRandomisation_",i,".csv"), header=T)
				envValues_obs = matrix(nrow=dim(obs)[1], ncol=length(envVariables))
				envValues_ran = matrix(nrow=dim(ran)[1], ncol=length(envVariables))
				colnames(envValues_obs) = gsub(".asc","",rasterNames)
				colnames(envValues_ran) = gsub(".asc","",rasterNames)
				for (j in 1:length(envVariables))
					{
						if (dim(envVariables[[j]])[3] > 1)
							{
								time_intervals = matrix(nrow=length(names(envVariables[[j]])), ncol=2)
								for (k in 1:length(names(envVariables[[j]])))
									{
										time_intervals[k,1] = as.numeric(unlist(strsplit(names(envVariables[[j]])[k],"_"))[2])
										time_intervals[k,2] = as.numeric(unlist(strsplit(names(envVariables[[j]])[k],"_"))[3])
									}
							}
						envValues_obs[,j] = raster::extract(envVariables[[j]], SpatialPoints(obs[,c("endLon","endLat")]))
						envValues_ran[,j] = raster::extract(envVariables[[j]], SpatialPoints(ran[,c("endLon","endLat")]))
					}
				write.csv(envValues_obs, paste0(localTreesDirectory,"/EnvValues_obs_",i,".csv"), row.names=F, quote=F)
				write.csv(envValues_ran, paste0(localTreesDirectory,"/EnvValues_ran_",i,".csv"), row.names=F, quote=F)
			}
		BFs = matrix(nrow=length(rasterNames), ncol=2); onlyConsideringTheTipBranches = TRUE
		row.names(BFs) = gsub(".asc","",rasterNames); colnames(BFs) = c("lower","higher")
		meanEnvValues_obs_list = list(); meanEnvValues_ran_list = list()
		for (i in 1:length(rasterNames))
			{
				lowerEnvValues = 0; meanEnvValues_obs_list = list(); meanEnvValues_ran_list = list()
				meanEnvValues_obs = rep(NA, nberOfExtractionFiles); meanEnvValues_ran = rep(NA, nberOfExtractionFiles)
				for (j in 1:nberOfExtractionFiles)
					{
						obs = read.csv(paste0(localTreesDirectory,"/EnvValues_obs_",j,".csv"))[,gsub(".asc","",gsub("-",".",rasterNames[i]))]
						ran = read.csv(paste0(localTreesDirectory,"/EnvValues_ran_",j,".csv"))[,gsub(".asc","",gsub("-",".",rasterNames[i]))]
						if (onlyConsideringTheTipBranches == TRUE)
							{
								tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
								obs = obs[which(tab[,"node2"]%in%tab[,"node1"])]
								tab = read.csv(paste0(localTreesDirectory,"/TreeRandomisation_",i,".csv"), head=T)
								ran = ran[which(tab[,"node2"]%in%tab[,"node1"])]
							}
						meanEnvValues_obs[j] = mean(obs, na.rm=T); meanEnvValues_ran[j] = mean(ran, na.rm=T)
						if (meanEnvValues_obs[j] < meanEnvValues_ran[j]) lowerEnvValues = lowerEnvValues+1				
					}
				p = lowerEnvValues/nberOfExtractionFiles; BFs[i,"lower"] = round((p/(1-p))/(0.5/(1-0.5)),1)
				p = (1-(lowerEnvValues/nberOfExtractionFiles)); BFs[i,"higher"] = round((p/(1-p))/(0.5/(1-0.5)),1)
				meanEnvValues_obs_list[[i]] = meanEnvValues_obs; meanEnvValues_ran_list[[i]] = meanEnvValues_ran
			}
		write.csv(BFs, paste0("RVFV-",segments[h],"_2.csv"), quote=F)
	}

# 5. Visualising the dispersal history of RVFV lineages

source("mccTree_R.r")
world = shapefile("Countries.shp")
study_area = crop(world, extent(4.75,59,-39,33))
for (h in 1:length(segments))
	{
		# 5.1. Extracting spatio-temporal information embedded in the MCC tree

		mcc_tre = readAnnotatedNexus(paste0(segments[h],"-lineage-C.mcc.tree"))
		mcc_tab = mccTree_R(mcc_tre, mostRecentSamplingDates[h])
		write.csv(mcc_tab, paste0("RVFV-",segments[h],".csv"), row.names=F, quote=F)
		
			# 5.2. Estimating the HPD region for each time slice
		
		prob = 0.80; precision = 2
		minYears = rep(NA, nberOfExtractionFiles)
		for (i in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
				minYears[i] = min(tab[,"startYear"])
			}
		startDatum = HDInterval::hdi(minYears)[1] # to get the lower bound of the 95% HPD interval for the starting year to consider
		polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
		
			# 5.3. Defining the different colour scales to use
		
		cols = gsub("FF","",rev(viridis::magma(151))[1:101])
		minYear = startDatum; maxYear = mostRecentSamplingDates[h]
		mcc = read.csv(paste0("RVFV-",segments[h],".csv"), head=T)
		mcc = mcc[order(mcc[,"endYear"],mcc[,"startYear"],decreasing=F),]
		startYear_index = (((min(mcc[,"startYear"])-minYear)/(maxYear-minYear))*100)+1
		startYear_colour = cols[startYear_index]
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices[which(endYears_indices<1)] = NA
		endYears_colours = cols[endYears_indices]
		polygons_colours = rep(NA, length(polygons))
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				if (polygon_index >= 1)
					{
						polygons_colours[i] = paste0(cols[polygon_index],"20")
					}
			}
		
			# 5.4. Co-plotting the HPD regions and MCC tree
				
		pdf(paste0("RVFV-",segments[h],".pdf"), width=6.5, height=8)
		par(mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o", col="gray30")
		plot(study_area, col="ivory2", border="gray30", lwd=0.5)
		for (i in 1:length(polygons))
			{
				plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
			}
		for (i in 1:dim(mcc)[1])
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
							arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		for (i in dim(mcc)[1]:1)
			{
				if (i == 1) # for the root (only correct if the first line of the CSV file corresponds to one of the two most ancestral branches)
					{
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYear_colour, cex=0.7)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
					}
				if (mcc[i,"node2"]%in%mcc[,"node1"]) # for internal nodes (dots)
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.7)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)	
					}	else	{	# for the tip nodes (squares)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.6)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="gray30", lwd=0.2, cex=0.6)				
					}
			}
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
		plot(rast, legend.only=T, add=T, col=cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.90,0.91,0.04,0.50),
			 legend.args=list(text="", cex=0.5, line=0.3, col="gray30", col.lab="gray30"), horizontal=F,
			 axis.args=list(cex.axis=0.55, lwd=0, lwd.tick=0.2, tck=-0.8, col.tick="gray30", col.axis="gray30", line=0, mgp=c(0,0.40,0)))
		dev.off()
	}

