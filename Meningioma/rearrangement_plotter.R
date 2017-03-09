# load libraries
#source.all()
source("/xchip/gistic/Jeremiah/fiss_get.R")
library(gUtils)
library(gTrack)
FISS='/xchip/tcga/Tools/gdac/bin/fiss'

## FISS GET
system.time(ind <- as.data.table(fiss_get('total_gistic', wkspace='Men_Pairs', type='pair', fuse=TRUE)))
ind[, Sample := pair_id]
#ind[, Sample := gsub("(.*?)-.*", "\\1", Sample)]
setkey(ind, Sample)

snow <- rbindlist(lapply(ind$snowman_annotate_csv[file.exists(ind$snowman_annotate_csv)], fread))

### get all of the events
#snow <- fread("/xchip/gistic/Jeremiah/Projects/MenPairs/EventPlots/v105_160722.csv")
snow[, Sample := gsub("(.*?).snowman.somatic.sv.vcf", "\\1", sample)]
setkey(snow, Sample)

## 
snow <- ind[, .(Tumor_bam_cov_gr_fragcounter, Tumor_bam_cov_gr_fragcounter_capture, Normal_bam_cov_gr_fragcounter, Normal_bam_cov_gr_fragcounter_capture, Sample)][snow]
setkey(snow, Sample)
snow <- snow[SPAN > 200 | SPAN < 0]

## load all the coverages
covs <- lapply(unique(snow$Sample), function(x) {
    
    if (TRUE || grepl("G", x)) {
        nc <- snow[x]$Normal_bam_cov_gr_fragcounter[1]
        tc <- snow[x]$Tumor_bam_cov_gr_fragcounter[1]
    } else {
        nc <- snow[x]$Normal_bam_cov_gr_fragcounter_capture[1]
        tc <- snow[x]$Tumor_bam_cov_gr_fragcounter_capture[1]
    }
    
    if (!file.exists(tc) || !file.exists(nc))
    {
        print(paste("not exist: ", basename(x)))
        return()
    }
    
    ## load cov
    print(basename(x))
    ncov <- readRDS(nc) 
    tcov <- readRDS(tc) 
    
    tcov$ratio = ifelse(ncov$reads > 0, tcov$reads / ncov$reads, 0)
    
    return(tcov)
    
})


##
#names(covs) <- unique(snow$Sample)
#saveRDS(covs, "wg_covs.rds")
covs <- readRDS("/xchip/gistic/Jeremiah/Projects/MenPairs/EventPlots/wg_covs.rds")

## lapply(unique(snow$Sample), function(x) {

##   if (!length(covs[[x]]))
##     return()

##   lks <- with(snow[x], GRanges(c(chr1,chr2), IRanges(c(pos1, pos2), width=1), c(strand1, strand2), id=rep(paste0(pos1,pos2),2)))
##   lks <- split(lks, lks$id)

##   g <- gTrack(covs[[x]], y.field='ratio')
##   pdf(dd <- paste0(x,".pdf"), width=40)
##   print(dd)
##   plot(g, links=lks, windows=sort(streduce(lks, pad=20000)))
##   dev.off()

## })

lapply(unique(snow$Sample), function(x) {
    
    if (!length(covs[[x]]))
        return()
    
    dir.create("full")
    
    #lks <- with(snow[x], GRanges(c(chr1,chr2), IRanges(c(pos1, pos2), width=1), c(strand1, strand2), id=rep(paste0(pos1,pos2),2)))
    lks <- with(snow[x], GRanges(c(seqnames,altchr), IRanges(c(start, altpos), width=1), c(strand, altstrand), id=rep(uid,2)))  
    lks <- split(lks, lks$id)
    
    g <- gTrack(covs[[x]][seq(1,length(covs[[x]]), by=20)], y.field='ratio')
    pdf(file.path('full',dd <- paste0(x,"_full.pdf")), width=40)
    print(dd)
    plot(g, links=lks, windows=si2gr(si))
    dev.off()
    
})

lapply(unique(snow$Sample), function(x) {
    
    if (length(covs[[x]]) == 0)
        return()
    
    print(x)
    lks <- with(snow[x], GRanges(c(chr1,chr2), IRanges(c(pos1, pos2), width=1), c(strand1, strand2), id=rep(paste0(pos1,pos2),2), span=c(SPAN,SPAN)))
    #lks <- with(snow[x], GRanges(c(chr1,chr2), IRanges(c(, altpos), width=1), c(strand, altstrand), id=rep(uid,2)))
    
    lks <- split(lks, lks$id)
    mcols(lks)$span <- snow$SPAN[match(names(lks),snow[x]$uid)]
    mcols(lks)$span <- snow[x][!duplicated(uid), SPAN]
    
    dir.create(file.path("v117",x))
    
    g <- gTrack(covs[[x]], y.field='ratio')
    
    lapply(seq_along(lks), function(y) {
        
        print('...')
        uid = snow[x][y,uid]
        ##print(paste(dd,"rar", y))
        mcols(lks)$col = "black"
        mcols(lks)$col[y] = "red"
        
        sp = lks[[y]]$span[1]
        if (sp < 0 | sp > 100000)
            suppressWarnings(windowr <- streduce(streduce(lks[y], pad=500000)+500000))
        else
            suppressWarnings(windowr <- lks[[y]][1] + max(10e3,min(sp * 7, 10000)))
        
        if (sp > 500 || sp < 0) {
            pdf(file.path(x, dd <- paste0(x,"_rar_", y, "_UID_", uid, ".pdf")), width=10)
            suppressWarnings(plot(g, links=lks, windows=windowr))
            dev.off()
        }
        
    })
})







### Noah's version
tumor.coverage <- readRDS("C:/Users/Noah/Downloads/men0025.tumor.cov.rds")
normal.coverage <- readRDS("C:/Users/Noah/Downloads/men0025.normal.cov.rds")
tumor.coverage.modified <- tumor.coverage[tumor.coverage@seqinfo@seqnames == "22", ]
tumor.coverage.modified <- tumor.coverage.modified[tumor.coverage.modified@ranges@start > 29000000 & tumor.coverage.modified@ranges@start < 31000000, ]
normal.coverage.modified <- normal.coverage[normal.coverage@seqinfo@seqnames == "22", ]
normal.coverage.modified <- normal.coverage.modified[normal.coverage.modified@ranges@start > 29000000 & normal.coverage.modified@ranges@start < 31000000, ]
ratio <- tumor.coverage.modified@elementMetadata@listData$reads / normal.coverage.modified@elementMetadata@listData$reads
location <- tumor.coverage.modified@ranges@start
plot(location, ratio)


