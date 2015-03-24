#organising c subunit ion table data
subunit_ion_table <- read.table("subunit_ion_table", header=TRUE, row.names=1)

rownames(subunit_ion_table) <- c("PM5","PM159","PM50","PM70","PM30","PM91")

ordered_subunit_ion_table <- subunit_ion_table[c("PM159","PM91","PM70","PM50","PM30","PM5"),]

prop_table = matrix(0, ncol = 6, nrow = 3)

for(i in 4:6)
{
	total_c_subunits = ordered_subunit_ion_table[i,1] + ordered_subunit_ion_table[i,2] + ordered_subunit_ion_table[i,3]
	prop_table[,i] = c(ordered_subunit_ion_table[i,1]/total_c_subunits,ordered_subunit_ion_table[i,2]/total_c_subunits,ordered_subunit_ion_table[i,3]/total_c_subunits)
}


#organising BLAST data
blast_table <- read.table("real_PM_c_subunits_vs_NR.tdt", stringsAsFactors = FALSE)

name_vector = vector()

for(i in 1: dim(blast_table)[1])
{
	name_vector[i] = strsplit(blast_table[i,1],"_")[[1]][1]
}

blast_table[,dim(blast_table)[2]+1] = name_vector

blast_averages = vector()

PM_mean_identities = vector()
PM_mean_bitscores = vector()
PM_sd_identities = vector()
PM_sd_bitscores = vector()

#PM5
PM5_identities = vector()
PM5_bitscores = vector()

for(i in 1:dim(blast_table)[1])
{
	if(blast_table[i,13] == "R1qual32.paired")
	{
		PM5_identities[length(PM5_identities)+1] <- blast_table[i,3]
		PM5_bitscores[length(PM5_bitscores)+1] <- blast_table[i,12]
	}
}

PM_mean_identities[1] <- mean(PM5_identities)
PM_mean_bitscores[1] <- mean(PM5_bitscores)
PM_sd_identities[1] <- sd(PM5_identities)
PM_sd_bitscores[1] <- sd(PM5_bitscores)


#PM30
PM30_identities = vector()
PM30_bitscores = vector()

for(i in 1:dim(blast_table)[1])
{
	if(blast_table[i,13] == "PM30")
	{
		PM30_identities[length(PM30_identities)+1] <- blast_table[i,3]
		PM30_bitscores[length(PM30_bitscores)+1] <- blast_table[i,12]
	}
}

PM_mean_identities[2] <- mean(PM30_identities)
PM_mean_bitscores[2] <- mean(PM30_bitscores)
PM_sd_identities[2] <- sd(PM30_identities)
PM_sd_bitscores[2] <- sd(PM30_bitscores)


#PM50
PM50_identities = vector()
PM50_bitscores = vector()

for(i in 1:dim(blast_table)[1])
{
	if(blast_table[i,13] == "PM50")
	{
		PM50_identities[length(PM50_identities)+1] <- blast_table[i,3]
		PM50_bitscores[length(PM50_bitscores)+1] <- blast_table[i,12]
	}
}

PM_mean_identities[3] <- mean(PM50_identities)
PM_mean_bitscores[3] <- mean(PM50_bitscores)
PM_sd_identities[3] <- sd(PM50_identities)
PM_sd_bitscores[3] <- sd(PM50_bitscores)


#PM70
#PM70_identities = vector()
#PM70_bitscores = vector()

#for(i in 1:dim(blast_table)[1])
#{
#	if(blast_table[i,13] == "PM70")
#	{
#		PM70_identities[length(PM70_identities)+1] <- blast_table[i,3]
#		PM70_bitscores[length(PM70_bitscores)+1] <- blast_table[i,12]
#	}
#}

#PM_mean_identities[4] <- mean(PM70_identities)
#PM_mean_bitscores[4] <- mean(PM70_bitscores)
#PM_sd_identities[4] <- sd(PM70_identities)
#PM_sd_bitscores[4] <- sd(PM70_bitscores)


#PM91
#PM91_identities = vector()
#PM91_bitscores = vector()

#for(i in 1:dim(blast_table)[1])
#{
#	if(blast_table[i,13] == "PM91")
#	{
#		PM91_identities[length(PM91_identities)+1] <- blast_table[i,3]
#		PM91_bitscores[length(PM91_bitscores)+1] <- blast_table[i,12]
#	}
#}

#PM_mean_identities[5] <- mean(PM91_identities)
#PM_mean_bitscores[5] <- mean(PM91_bitscores)
#PM_sd_identities[5] <- sd(PM91_identities)
#PM_sd_bitscores[5] <- sd(PM91_bitscores)


#PM159
#PM159_identities = vector()
#PM159_bitscores = vector()

#for(i in 1:dim(blast_table)[1])
#{
#	if(blast_table[i,13] == "PM159")
#	{
#		PM159_identities[length(PM159_identities)+1] <- blast_table[i,3]
#		PM159_bitscores[length(PM159_bitscores)+1] <- blast_table[i,12]
#	}
#}

#PM_mean_identities[6] <- mean(PM159_identities)
#PM_mean_bitscores[6] <- mean(PM159_bitscores)
#PM_sd_identities[6] <- sd(PM159_identities)
#PM_sd_bitscores[6] <- sd(PM159_bitscores)


#Plotting the data
#par(mfrow=c(1,2))
pdf("~/Documents/Bo_review/sodium_atp_figure/figure_y_sodium_atp_fractions.pdf", width=8, height=4)
barplot(prop_table[,4:6], horiz=TRUE, names.arg = row.names(ordered_subunit_ion_table)[4:6])
dev.off()
#barplot(rev(PM_mean_identities),horiz=TRUE,xlim=c(0,110))
#for(i in 1:length(PM_sd_identities))
#{
#	lines(c(rev(PM_mean_identities[i])-rev(PM_sd_identities[i]),rev(PM_mean_identities[i])+rev(PM_sd_identities[i])),c(7.9-1.2*i,7.9-1.2*i))
#}
#for(i in 1:length(PM_sd_identities))
#{
#	lines(c(rev(PM_mean_identities[i])-rev(PM_sd_identities[i]),rev(PM_mean_identities[i])-rev(PM_sd_identities[i])),c(8-1.2*i,7.8-1.2*i))
#}
#for(i in 1:length(PM_sd_identities))
#{
#	lines(c(rev(PM_mean_identities[i])+rev(PM_sd_identities[i]),rev(PM_mean_identities[i])+rev(PM_sd_identities[i])),c(8-1.2*i,7.8-1.2*i))
#}
