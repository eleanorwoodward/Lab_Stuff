# Schwannoma Plotting Code
# Copyright PK Agarwalla
# 8/24/15

# To make comut plot for schwannoma discovery set
# SIF file will have the LOH, Copy Loss of 22, copy-neural LOH, and number of mutations per sample already encoded into it. It will also have the location of the PON'ed and exac filtered mafs (although less important for NF2).

# For now, just focus on adding the additional annotations.

# Current sif will be master_sif_schwannoma_8_24_15.txt



# Code to create co-mut plot of all data including whole exome and whole genoma samples
# Load SIF
master_sif=read.delim('/xchip/beroukhimlab/agarwalla/meningioma/analysis/validation/sif/meningioma_validation_sif_20150908.txt',stringsAsFactors=FALSE)

#Prepare co-mut rows and type of alteration
sig_gene_list=c('NF2', 'CDKN2C', 'NF2')
alt_type=c('mut', 'mut', 'Loss')

#Reverse order
sig_gene_list=rev(sig_gene_list)
alt_type=rev(alt_type)

# set margins at 0 and save it as a global variable
# The double assignment makes the variable global (can be used outside functions too)

old.par<<-par(mar=c(0,0,0,0))

chromlocs=read.delim('/xchip/beroukhimlab/agarwalla/schwannoma/analysis/chromInfo.txt', stringsAsFactors=FALSE)

silent_mutations=c("Silent", "3'UTR","5'UTR",'Intron',"5'Flank","Non-coding_Transcript","IGR", "lincRNA", "RNA")
# update the silent mutation list

# Function to assign the plot colors with input of the maf and whether there is a mut match

assign_plot_col<-function(maf,mut_matches){
	# Set the plotting defaults for color
	null_col='gainsboro'
	nonsense_col=rgb(205, 39, 48, maxColorValue=255)
	missense_col=rgb(61, 117, 164,maxColorValue=255)
	synonymous_col='gray70'
	frameshift_col=rgb(228, 118, 54, maxColorValue=255)
	compound_col='black'
	splice_site_col=rgb(131, 76, 139, maxColorValue=255)
	hotspot_mutation_col=rgb(77, 165, 76,maxColorValue=255)
	in_frame_indel_col=rgb(244, 232, 83,maxColorValue=255)
	
	#hotspots=c('p.P286R','p.V411L','p.R660*')
	hotspots=c('p.P44L','p.H28R','p.R60Q')
	#PPP2R1A_hotspots=c('p.P179R','p.R183W','p.S256F')
	
	if (length(mut_matches)>0){
		pChanges=maf[mut_matches,'Protein_Change']
		vcs=maf[mut_matches,'Variant_Classification']
		nonsilent_vcs=vcs[which(!(vcs %in% silent_mutations))]
	
		if (any(pChanges %in% hotspots)){
			plot_col=hotspot_mutation_col
			return(plot_col)
		}
		if (length(nonsilent_vcs)>1){
			# Now we know that it is a compound mutation, let's color the box accordingly
			plot_col=compound_col
		}
		if (length(nonsilent_vcs)==1){
			if (nonsilent_vcs=='Missense_Mutation'){plot_col=missense_col}
			if (nonsilent_vcs=='Nonsense_Mutation'){plot_col=nonsense_col}
			if (nonsilent_vcs=='Frame_Shift_Del'){plot_col=frameshift_col}
			if (nonsilent_vcs=='Frame_Shift_Ins'){plot_col=frameshift_col}
			if (nonsilent_vcs%in%c('In_Frame_Del','In_Frame_Ins')){plot_col=in_frame_indel_col}
			if (nonsilent_vcs=='Splice_Site'){plot_col=splice_site_col}
			cat(nonsilent_vcs,' ; ',plot_col,'\n')
		}
		if (length(nonsilent_vcs)==0){plot_col=synonymous_col}

		}else{plot_col=null_col}
return(plot_col)
}


# Variables for workspace environment
# NF2 gene location 29999545 to 30094589
NF2_begin=29999545
NF2_end = 30094589
loss_list = master_sif$loss.22q.tuning0.8
#cnLOH_list = master_sif$cnLOH
#loh_list = master_sif$LOH
n_mut = rep(1, length(loss_list))
NF2_mut_list = c()
master.maf = read.delim('/xchip/beroukhimlab/agarwalla/meningioma/analysis/validation/meningioma.all.germline.filtered.20150816.maf', stringsAsFactors = FALSE, comment.char = "#")
# note that sif file now has rate_per_MB which is mutation rate per MB and can be used to normalize the samples.


for ( i in 1:nrow(master_sif)){
	cat('Calculating maf for sample ', i ,'\n')
	individual.index <- grep(master_sif$case_sample[i], master.maf$Tumor_Sample_Barcode)
	temp.maf <- master.maf[individual.index, ]
	NF2.mut.index <- which(temp.maf$'Hugo_Symbol' == 'NF2' & !(temp.maf$Variant_Classification %in% silent_mutations))
	if(length(NF2.mut.index) >0){NF2_mut_list=c(NF2_mut_list,1)}else{NF2_mut_list=c(NF2_mut_list,0)}
}


new_order=rev(order(NF2_mut_list,loss_list,n_mut))
reordered_master_sif=master_sif[new_order,]
unique_individuals=reordered_master_sif$case_sample
box_width=0.9


output_fn='/xchip/beroukhimlab/agarwalla/meningioma/analysis/validation/meningioma_validation_comut_8_30_15.pdf'


#Plot the coMut output

pdf(output_fn)
plot(-10,-10,xlim=c(0,length(unique_individuals)+1),ylim=c(0,length(sig_gene_list)+3),yaxt='n',xaxt='n',xlab='',ylab='')

y_label_positions=seq(from=1,to=length(sig_gene_list))
y_labels=sig_gene_list
axis(side=2,at=y_label_positions,labels=FALSE,cex.axis=0.75)
text(y = y_label_positions, par("usr")[1], labels = y_labels, srt = 0, pos = 2, xpd = TRUE)


x_label_positions=seq(from=1,to=length(unique_individuals))
x_labels=unique_individuals
axis(side=1,at=x_label_positions,labels=x_labels,las=2,cex.axis=0.75)



for(i in 1:nrow(reordered_master_sif)){
	#temp_pair_id=reordered_master_sif[i,'pair_id']
	temp_individual=unique_individuals[i]
	#sample_name=temp_individual

	individual.index <- grep(temp_individual, master.maf$Tumor_Sample_Barcode)
	sample_maf <- master.maf[individual.index, ]

	for (j in 1:length(sig_gene_list)){
			g=sig_gene_list[j]
			
			if(alt_type[j]=='mut'){
							# Check the maf for a nonsilent mutation in the gene.
			#mut_matches=which(merged_fc_maf[,'Hugo_Symbol']==g & !(merged_fc_maf[,'Variant_Classification']%in%silent_mutations))
			mut_matches=which(sample_maf[,'Hugo_Symbol']==g)
			plot_col=assign_plot_col(sample_maf,mut_matches)
			
			# fill the rectangle
			rect(xleft=i-box_width/2,xright=i+box_width/2,ybottom=j-box_width/2,ytop=j+box_width/2,col=plot_col,border=NA)
			}
			
			if(alt_type[j]=='LOH'){
				if (g=='NF2'){
					if(loh_list[new_order[i]]==1){
						plot_col='blueviolet'
					}else{plot_col='gainsboro'}
				}
				rect(xleft=i-box_width/2,xright=i+box_width/2,ybottom=j-box_width/2,ytop=j+box_width/2,col=plot_col,border=NA)					
			}
			if(alt_type[j]=='Loss'){
				if (g=='NF2'){
					if(loss_list[new_order[i]]==1){
						plot_col='dodgerblue4'
					}else{plot_col='gainsboro'}
				}
				rect(xleft=i-box_width/2,xright=i+box_width/2,ybottom=j-box_width/2,ytop=j+box_width/2,col=plot_col,border=NA)					
			}
			if(alt_type[j]=='cnLOH'){
				if (g=='NF2'){
					if(cnLOH_list[new_order[i]]==1){
						plot_col='darkolivegreen'
					}else{plot_col='gainsboro'}
				}
				rect(xleft=i-box_width/2,xright=i+box_width/2,ybottom=j-box_width/2,ytop=j+box_width/2,col=plot_col,border=NA)					
			}


	}
		
	n_mutations = rep(1, length(loss_list)) #reordered_master_sif[i,'rate_per_MB']

	nMutHeight=n_mutations
	rect(xleft=i-box_width/2,xright=i+box_width/2,ybottom=0.1+j+box_width/2,ytop=0.1+j+box_width/2+nMutHeight,col=rgb(61, 117, 164,maxColorValue=255),border=NA)
	
}


# par(mfrow=c(4,1))

# plot(-10,-10,xlim=c(0,length(unique_individuals)+1),ylim=c(0,100),xaxt='n',xlab='',ylab='Purity')
# for ( i in 1:nrow(reordered_master_sif)){
# 	temp_purity=100*reordered_master_sif[i,'purity']
# 	points(i,temp_purity,pch=16)
# }


# par(old.par)

dev.off()

# # The following code puts all of the copy number heatmaps onto the same plot.

# plot.new()
# frame()
# # par function sets the figure space using normalized device coordinates (x1, x2, y1, y2) and adds it to the current plot using the new=TRUE, the margins of (bottom, left, top, right) are also set
# # the following sets the plotting parameters to be only 1/6 of the device up
# par(fig=c(0,1,0,1/6), new=TRUE, mar=c(0,0,0,0))

# # creates a base plot with 1:number as the y values and the x is the index (here x = y), but type = 'n' suppresses the output on the frame
# plot(1:30000,,ylim=c(0,1),type='n',xlim=c(0,2900000000),frame.plot=FALSE,axes=T,xaxt='n',xlab='')
# chromlocs=read.delim('/xchip/beroukhimlab/agarwalla/schwannoma/analysis/chromInfo.txt',stringsAsFactors=FALSE)
# # this loop just adds the chromosome labels to the bottom of the graph at the centromere position
# for(i in 1:nrow(chromlocs)){
# 	  if(i%%2==1){bar_col='gray90'}
# 	  if(i%%2==0){bar_col='gray100'}
# 	#rect(xleft=0,xright=0.5, ybottom=chromlocs[i,'total_chr_start_pos'], ytop=chromlocs[i,'total_chr_end_pos'],col=bar_col)
# 	#abline(h=chromlocs[i,'total_centromere_pos'],lty=3,col='gray35')
# 	text(y=1,x=chromlocs[i,'total_centromere_pos'],i,cex=0.5, las = 2)
# 	}

# # start next plotting element with full width of matte, but from 1/6 to 4/6 of the 
# # the following plots the gray-scale color bar as the key for the heat map
# # segments command makes a line at the same x value and y at the specified values of ybottom and ytop (-0.1 and .05), this is for the centromere dashed line
# # recall that the ylim is set to 1 here
# # lty = 3 is a line type that is short-dashed

# par(fig=c(0,1,1/6,4/6), new=TRUE, mar=c(0,0,0,0))
# plot(1:30000,,xlim=c(0,2900000000),type='n',ylim=c(-0.1,length(unique_individuals)+.1),frame.plot=FALSE,axes=T,xaxt='n',xlab='',yaxt='n')
# for(i in 1:nrow(chromlocs)){
# 		  if(i%%2==1){bar_col='gray90'}
# 		  if(i%%2==0){bar_col='gray100'}
# 		rect(xleft=chromlocs[i,'total_chr_start_pos'],xright=chromlocs[i,'total_chr_end_pos'], ybottom=-1.5, ytop=.05,col=bar_col)
# 		#abline(h=chromlocs[i,'total_centromere_pos'],lty=3,col='gray35')
# 		segments(x0=chromlocs[i,'total_centromere_pos'],y0=-1,x1=chromlocs[i,'total_centromere_pos'],y1=.05,lty=3,col='gray35')
# 		if(FALSE){text(y=0,x=chromlocs[i,'total_centromere_pos'],i,cex=0.5,pos=1,offset=1)}
# 		}

# #colorRampPalette takes colors and returns a function that can then interpolate a spectrum with an integer value so here cn_cols is a function
# 	cn_cols<<-colorRampPalette(c('blue','white','red')) # select color scheme
# 	wg_col_set=cn_cols(100) # creates a 100 panel color set in hexa color scheme
# 	# new function to take the temporary copy-number total (example is 2) or the observed copy number ratio and find the maximum color out of the 100 that goes up to the temporary copy number level
# 	# if nothing satisfies this requirement, then it sets the position to 1 which is blue (0000FF in the wg_col_set)
# 	# function returns the end_col which is the position in the wg_col_set
# 	wg_cols<-function(y,wg_col_set){
# 	split=seq(from=0,to=4,length.out=100)  #modified from split=seq(from=0,to=4,length.out=100)
# 	temp_pos=which(split<y)
	
# 		if (length(temp_pos)>0){
# 			temp_pos=max(temp_pos)
# 		}else{temp_pos=1}
	
# 	end_col=wg_col_set[temp_pos]
# 	return(end_col)
# 	}

# # temp_col_scheme=c('blue','orange','purple','yellow','red','black','pink')

# # load in all segs file (this file has undergone ISAR correction, removal of events < 1 Megabase, and thresholed for 0.2)

# ss_all = read.delim('/xchip/beroukhimlab/agarwalla/schwannoma/analysis/comut/')
# for (j in 1:length(unique_individuals)){
# 		ss_filename = reordered_master_sif[j,'absolute_seg']

# 	if (is.na(ss_filename)){
# 		ss_recapseg = read.delim(paste(recapseg_dir, recapseg_files[grep(unique_individuals[j], recapseg_files)], sep=""), stringsAsFactors = FALSE)
# 		# Remove artifct segmentations (>10)
# 		bad_segment_means <- which(abs(ss_recapseg$Segment_Mean) > 5)
# 		if (length(bad_segment_means) > 0){ss <- ss_recapseg[-bad_segment_means, ]}else{ss <- ss_recapseg}
# 		rect(xleft=0,xright=2900000000,ybottom=j-1+.1,ytop=j-.1,col='white',border=NA)

		
# 	#ISAR - in silico admixture removal - updated equations from Schumacher
# 	# R(x) - Observed copy ratio at locus x
# 	# alpha = purity
# 	# tau = cancer cell ploidy (the average copy number across the genome)
# 	# q(x) = integer copy number in the cancer cells
# 	# D = average ploidy of all cells (including non-cancer cells)

# 	# q = D*(R/alpha) - 2*((1-alpha)/alpha)
	
# 	#  ReCapseg output is already normalized and relativized, so it is a copy ratio that is uncorrected for purity (e.g. Segment Mean = 0 is normal 2n ploidy)
# 	# Therefore, ReCapSeg Segment mean is R(x) and we would like to calculate Rprime which is the purity corrected form

# 	alpha = reordered_master_sif$purity[j]
# 	tau = 2 # (assume that average copy number across genome is 2)
# 	S = (alpha*tau + (1-alpha)*2) # Sample ploidy (Scott's D value)
# 	b = 2*(1-alpha)/S # copy ratio of homozygous deletion
# 	dt = alpha/S # spacing between discrete levels
# 	thresh = 0.7

# 	# R = 2^(log copy ratio)
# 	# qhat = (R-b)/dt
# 	# isar = qhat / T
	
# 	#D = (alpha*tau) + 2*(1-alpha)

# 		for (i in 1:nrow(ss)){
# 			if (!(is.na(ss[i,'Segment_Mean']))){
# 				ss$length[i] = ss$End[i] - ss$Start[i]
# 				temp_chr=ss[i,'Chromosome']
# 				chr_start=chromlocs[temp_chr,'total_chr_start_pos']
# 				start_seg=chr_start+ss[i,'Start']
# 				end_seg=chr_start+ss[i,'End']
# 				ss$obs_copy_ratio[i] = (2^(ss[i, 'Segment_Mean'])) # Recapseg gives log2 copy ratio, so obs copy ratio is regular copy ratio (e.g. 0.5 is LOH in pure sample)
# 				ss$qhat[i] = (ss$obs_copy_ratio[i] - b)/dt # plot qhat for copy number
# 				ss$isar[i] = ss$qhat[i]/tau # is new the copy ratio
# 				ss$new_log_copy_ratio[i] = log(ss$isar[i], base = 2)
				
# 				#Threshold the heat map data on the new log_copy_ratio
# 				if ((abs((ss$qhat[i])-2) < thresh) | (is.na(ss$qhat[i]) == TRUE)){temp_CN_total = 2}else{temp_CN_total = ss$qhat[i]}
# 				#temp_CN_total = 2*(2^(ss[i, 'Segment_Mean']))
# 				temp_CN_col = wg_cols(temp_CN_total,wg_col_set)
# 				rect(xleft=start_seg,xright=end_seg,ybottom=j-1+.1,ytop=j-.1,col=temp_CN_col,border=NA)
# 			}
# 		}

	
# 	#reordered_master_sif$avg_rel_copy_level[j] <- sum(as.numeric(ss[,'length'])*as.numeric(ss[,'Segment_Mean']))/sum(as.numeric(ss[,'length']))



# 	}else{
	
# 		ss=read.delim(reordered_master_sif[j,'absolute_seg'])
# 		#ss<<-read.delim(master_sif$absolute_seg[which(master_sif$case_sample==colnames(mut_matrix)[col_order[j]])],stringsAsFactors=FALSE)
# 		genome_wide_average_expected_copy_level=sum(as.numeric(ss[,'length'])*as.numeric(ss[,'expected_total_cn']))/sum(as.numeric(ss[,'length']))
# 		genome_wide_average_corrected_copy_level=sum(as.numeric(ss[,'length'])*as.numeric(ss[,'corrected_total_cn']))/sum(as.numeric(ss[,'length']))
# 		genome_wide_average_total_copy_level=sum(as.numeric(ss[,'length'])*as.numeric(ss[,'modal.a1']+ss[,'modal.a2']))/sum(as.numeric(ss[,'length']))
# 		rect(xleft=0,xright=2900000000,ybottom=j-1+.1,ytop=j-.1,col='white',border=NA)

# 		for (i in 1:nrow(ss)){

# 			if (!(is.na(ss[i,'copy.ratio']))){
# 			  temp_chr=ss[i,'Chromosome']
# 			  chr_start=chromlocs[temp_chr,'total_chr_start_pos']
# 			  start_seg=chr_start+ss[i,'Start.bp']
# 			  end_seg=chr_start+ss[i,'End.bp']
# 			  temp_CN_minor=ss[i,'modal.a1']
# 			  temp_CN_major=ss[i,'modal.a2']
# 			  temp_CN_total=ss[i,'modal.a1']+ss[i,'modal.a2']
# 			  temp_CN_total_relative=log(temp_CN_total/genome_wide_average_total_copy_level,base=2)
# 			  temp_CN_col =  wg_cols(temp_CN_total, wg_col_set)
# 			  # if (ss[i,'LOH']==1){
# 			  # 	temp_CN_col='green4'
# 			  # }else{ temp_CN_col=  wg_cols(temp_CN_total,wg_col_set)}
# 			  temp_cn_diff<<-temp_CN_major-temp_CN_minor
# 			  rect(xleft=start_seg,xright=end_seg,ybottom=j-1+.1,ytop=j-.1,col=temp_CN_col,border=NA)
# 			  #segments(x0=start_seg,x1=end_seg,y0=j-1,y1=j-1,col=temp_CN_col,lwd=2,lend=1)  
# 			}

# 		}

# 	}

# }


# dev.off()