# tnseq.R -- code to manipulate Tn-seq data from ISB

#   as of 19-Oct-2017

#  R package is available at: https://github.com/robertdouglasmorrison/DuffyTools
library( DuffyTools)

# assume MTb, but the top level function can override
setCurrentSpecies( "MT_H37")

# set some defaults for the top level function, but it can override on the fly
WIG_FOLDER <- "Test1.pH5-pH7"

# set an expected minimum read count, warn if not enough per dataset
MIN_COUNTS_PER_WIG <- 1000000

# some constants the help decide with sites get thrown out for not enough reads
MIN_COUNT_DROP <- 8
MIN_MEDIAN_PCT_DROP <- 0.10
DROP_MODE <- "percent"

# an internal constant, for separating the fixed columns about GeneID, location, etc, from the
# variable columns from each data file
# we have added 5 columns of annotation before the insertion counts data
FIRST_COUNT_COLUMN <- 6    

# how to combine all replicates and site fitness scores for one gene
GENE_AVG_FUN <- median
REPLICATE_AVG_FUN <- median

# some constants for making gene plots
TAIL_WIDTH <- 1200
N_HTML_GENES <- 100


# top level function that completely processes one folder of raw data, and puts all results 
# into that same folder
do.all <- function( path=WIG_FOLDER, reload=FALSE, speciesID=getCurrentSpecies(), offset.readCount=NULL,
			nGenePlots=N_HTML_GENES) {

	if ( speciesID != getCurrentSpecies()) {
		setCurrentSpecies( speciesID)
	}
	gmap <<- getCurrentGeneMap()

	# the MTb genome has a lot of alternate "MTxxx" genes too, drop them?
	isAltGenes <- grep( "^MT[0-9]{4}", gmap$GENE_ID)
	if ( length( isAltGenes)) gmap <<- gmap[ -isAltGenes, ]

	if ( path != WIG_FOLDER) WIG_FOLDER <<- path

	# read in the WIG files and create or read the Sample Key
	rawMatrix <<- loadWigFileData( path, reload=reload)
	sampleKey <<- loadSampleDetails( path)
	if ( all( is.na( sampleKey$ExpansionFactor))) {
		cat( "\nWIG files read in successfully")
		cat( "\nSampleKey file created.   Now manually add 'ExpansionFactor' values.\n")
		return()
	}

	# process that data...
	fitData <<- processData( path, rawMatrix, sampleKey, control.genes=NULL, offset.readCount=offset.readCount)

	# making all the plots will reuse many data.frames from disk, set up quick buffering
	fileBuffering( "setup")

	# generate files & plots of results for each condition
	finalAns <<- calculateGeneFitness( path, fitData, sampleKey, nGenePlots=nGenePlots)

	# and do a comparison between conditions, for genes most different
	diffAns <<- compareConditionFitness( path, sampleKey, min.delta=0.1, min.sites=1, min.fits=3, nGenePlots=nGenePlots)
}


processData <- function( path=WIG_FOLDER, tbl=rawMatrix, samples=sampleKey, control.genes=NULL,
				offset.readCount=0) {

	# new idea:  don't do the 'drop' and normalize on everyone;  instead do it on each pair on the fly

	# before we process, toss out any pairs that don't deserve to be used
	toExclude <- which( samples$Exclude)
	if ( length( toExclude)) {
		columnsToDrop <- c( samples$Time1.Name[toExclude], samples$Time2.Name[toExclude])
		cat( "\nDropping ", length(columnsToDrop), "columns from Counts Table\n")
		where <- match( columnsToDrop, colnames(tbl))
		tbl <- tbl[ , -where]
	}

	fitData <- calculateAllFitness( tbl, samples, control.genes=NULL, offset.readCount=offset.readCount)
	write.table( fitData, file.path( path, "FitScores.csv"), sep=",", quote=T, row.names=F)
	return( fitData)
}


loadWigFileData <- function( path=WIG_FOLDER, reload=FALSE) {

	if ( ! file.exists( path)) stop( paste( "Folder of WIG files not found: ", path))

	# we will read (or create from raw WIG files) one matrix of counts
	countFile <- file.path( path, "RawCounts.csv")
	if ( file.exists( countFile) && ! reload) {
		cat( "\nLoaded existing counts file:  ", countFile, "\n")
		return( read.csv( countFile, as.is=T))
	}

	# read in the raw Wig files to create the count data
	m <- readAllWigFiles( path)

	# augment with annotation details
	pos <- as.integer( rownames(m))
	where <- fastFindInterval( pos, gmap$POSITION)
	gid <- gmap$GENE_ID[where]
	gsym <- gmap$NAME[where]
	gprod <- gene2Product(gid)
	out <- data.frame( "GENE_ID"=gid, "SYMBOL"=gsym, "PRODUCT"=gprod, "POSITION"=pos, "DROP"=FALSE,
				m, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	write.table( out, countFile, sep=",", quote=T, row.names=F)
	cat( "\nWrote new counts data to file:  ", countFile, "\n")
	return( out)
}


readAllWigFiles <- function( path=WIG_FOLDER, pattern=".wig$") {

	# read in all the raw Wig files in one folder, return as one matrix
	fset <- dir( path, pattern=pattern, full.name=T)
	ans <- lapply( fset, readOneWigFile)

	# combine, expecting all positions to match exactly
	Nwig <- length( ans)
	Nloc <- nrow( ans[[1]])
	m <- matrix( NA, nrow=Nloc, ncol=Nwig)
	colnames(m) <- sub( pattern, "", basename(fset))

	cat( "\nCombining..")
	for ( i in 1:Nwig) {
		thisAns <- ans[[i]]
		if ( i == 1) {
			allLocations <- thisAns$location
			rownames(m) <- allLocations
		} else {
			if ( any( thisAns$location != allLocations)) {
				cat( "\nWIG Location mis-matches!   File=", fset[i], "\n")
			}
		}
		m[ , i] <- thisAns$depth
	}
	cat( "  Done.")
	return( m)
}


readOneWigFile <- function( f) {

	# these files are "Variable Step", with "one location, one Depth" layout, space delimited
	cat( "\nReading WIG file: ", basename(f))
	txt <- readLines( f)
	cat( "  N_Lines: ", length(txt))

	# line 1 is comment, line 2 is definition, all usable lines start with digits
	isData <- grep( "^[0-9]", txt)

	# one space separates location, depth
	terms <- strsplit( txt[isData], split=" +")
	loc <- as.numeric( sapply( terms, FUN=`[`, 1))
	depth <- as.numeric( sapply( terms, FUN=`[`, 2))

	out <- data.frame( "location"=loc, "depth"=depth, stringsAsFactors=F)
	ord <- order( out$location)
	out <- out[ ord, ]

	totalCounts <- sum( out$depth, na.mr=T)
	cat( "  N_Counts: ", totalCounts)
	if ( totalCounts < MIN_COUNTS_PER_WIG) {
		cat( "\n  Warning:  Low Counts in WIG file!  Perhaps exclude this experiment.\n")
	}
	return( out)
}


loadSampleDetails <- function( path=WIG_FOLDER) {

	# load or create the table of files, samples, times, replicates, etc
	sampleFile <- file.path( path, "SampleKey.csv")
	if ( file.exists( sampleFile)) {
		tbl <- read.csv( sampleFile, as.is=T)
		# if the file exists, we should force the expansion factor 'd' to be set already
		if ( any( is.na( tbl$ExpansionFactor))) {
			cat( "\nSome Expansion Factors not yet set in Sample Key file.")
			cat( "\nHand edit file:  ", sampleFile)
			cat( "\nOnce all expansion factors are in place, re-run script.")
			stop()
		}
		cat( "\nLoaded existing Sample Key:  ", sampleFile, "\n")
		toExclude <- which( tbl$Exclude)
		if ( N <- length(toExclude)) {
			cat( "\n  Experiments flagged for exclusion:  ", N, "\t",tbl$SampleID[ toExclude], "\n")
		} else {
			cat( "\n  No Experiments excluded.\n")
		}
		return( tbl)
	}

	# when not found, create it from the counts matrix column names, using some hueristics
	# from how Eliza named her samples.    This is a bit subjective...
	cat( "\nCreating Sample Key from Count data file..")
	countFile <- file.path( path, "RawCounts.csv")
	if ( !file.exists( countFile)) {
		cat( "\nError:  Missing existing counts file:  ", countFile)
		stop()
	} 
	tbl <- read.csv( countFile, as.is=T)
	sampleNames <- colnames(tbl)[ FIRST_COUNT_COLUMN:ncol(tbl)]

	timePoint <- sub( "(^T[12])(\\.[0-9]+_)(.+)", "\\1", sampleNames)
	replicate <- sub( "(^T[12]\\.)([0-9]+)(_.+)", "\\2", sampleNames)
	condition <- sub( "(^T[12].[0-9]+_)(.+$)", "\\2", sampleNames)

	# try to catch any unexpected filename formats that mess up our naming convention
	anyBad <- FALSE
	badTimePoint <- which( ! (timePoint %in% c( "T1", "T2")))
	if ( length( badTimePoint)) {
		cat( "\n\nWIG filename has unexpected 'TimePoint' term: ", sampleNames[badTimePoint])
		cat( "\nExpected 'T1.replicate_condition' or 'T2.replicate_condition'")
		anyBad <- TRUE
	}
	badReplicate <- which( is.na( as.numeric( replicate)))
	if ( length( badReplicate)) {
		cat( "\n\nWIG filename has unexpected 'Replicate' term: ", sampleNames[badReplicate])
		cat( "\nExpected 'T1.replicate_condition' or 'T2.replicate_condition'")
		anyBad <- TRUE
	}
	if (anyBad) {
		cat( "\n\nFix naming of WIG files before continuing..")
		# also remove any just created files
		file.delete( sampleFile)
		file.delete( countFile)
		stop()
	}

	sampleKeys <- paste( condition, replicate, sep="-")
	sampleIDs <- sort( unique( sampleKeys))
	NS <- length( sampleIDs)
	whoT1 <- whoT2 <- vector( length=NS)
	for ( i in 1:NS) {
		whoT1[i] <- intersect( which( sampleKeys == sampleIDs[i]), which( timePoint == "T1"))[1]
		whoT2[i] <- intersect( which( sampleKeys == sampleIDs[i]), which( timePoint == "T2"))[1]
	}

	out <- data.frame( "SampleID"=sampleIDs, "Condition"=condition[whoT1], "Replicate"=replicate[whoT1],
			"ExpansionFactor"=NA, "Exclude"=FALSE, "Time1.Name"=sampleNames[whoT1], 
			"Time2.Name"=sampleNames[whoT2], stringsAsFactors=FALSE)
	rownames(out) <- 1:NS
	cat( "\nNumber of Samples:          ", NS)
	cat( "\nExperimental Conditions:    ", sort( unique( condition)))
	write.table( out, sampleFile, sep=",", quote=T, row.names=F)
	cat( "\nWrote new Sample Key file:  ", sampleFile, "\n")
	return( out)
}


dropLowCountInsertions <- function( tbl) {

	# we won't actually drop rows, just flag them
	# now we are always called  with exactly 2 columns of data -- the 2 time points being compared right now
	from <- FIRST_COUNT_COLUMN
	to <- ncol(tbl)
	if ( to != from+1) stop( "Expected exactly 2 columns of count data")
	m <- as.matrix( tbl[ , from:to])

	# let's allow using a percent of the mean cutoff instead of a simple count
	cat( "\nAssessing which insertion sites to drop for low counts:")
	if ( DROP_MODE == "percent") {
		# get the mean count for each column separately
		columnMeans <- apply( m, 2, mean, na.rm=T)
		# our cutoff is to be a percentage of that mean
		columnMinCounts <- round( columnMeans * MIN_MEDIAN_PCT_DROP, digits=2)
	} else {
		columnMinCounts <- rep.int( MIN_COUNT_DROP, ncol(m))
	}
	columnMinCounts <- round( columnMinCounts)

	nHighCounts <- apply( m, MARGIN=1, function(x) sum( x >= columnMinCounts))
	drops <- which( nHighCounts == 0)
	tbl$DROP <- FALSE
	if ( length(drops)) {
		if ( DROP_MODE == "count") {
			cat( "\nSites dropped for COUNT <", MIN_COUNT_DROP, "\tN_Drops: ", length(drops))
		} else {
			cat( "\nSites dropped for PERCENT <", MIN_MEDIAN_PCT_DROP, "\tN_Drops: ", length(drops))
		}
		tbl$DROP[ drops] <- TRUE
		m <- m[ -drops, ]
	}
	cat( "\nN_Kept_Insertions:                  ", sum( ! tbl$DROP))
	nTotal <- apply( m, 2, sum, na.rm=T)
	cat( "\nMean Counts per Time Point:         ", round( mean( nTotal)))
	cat( "\nAverage Counts per Insertion Site:  ", round( mean( nTotal/nrow(m)), digits=2), "\n")
	return( tbl)
}


normalizeInsertionCounts <- function( tbl, max.scaling=2) {

	# normalize the counts in the data matrix
	cat( "\nNormalizing Raw Counts data..")
	from <- FIRST_COUNT_COLUMN
	to <- ncol(tbl)
	if ( to != from+1) stop( "Expected exactly 2 columns of count data")
	m <- as.matrix( tbl[ , from:to])

	# only look at the 'not-dropped' rows
	toUse <- which( tbl$DROP == FALSE)
	totalCounts <- apply( m[toUse, ], MARGIN=2, sum)

	targetCount <- median( totalCounts)
	mOut <- m
	for( i in 1:ncol(m)) {
		thisV <- m[ toUse, i]
		thisTotal <- totalCounts[i]
		scaleFactor <- targetCount / thisTotal
		newV <- thisV * scaleFactor
		mOut[ toUse, i] <- round( newV, digits=3)
		cat( "\n  ScaleFactor applied:  ", colnames(m)[i], "\t", round( scaleFactor, digits=2))
		if ( scaleFactor > max.scaling || scaleFactor < (1/max.scaling)) {
			cat( "\n    Warn:  Excessive Scaling was applied.   Consider excluding: ", colnames(m)[i], "\n")
		}
	}

	normTotal <- apply( mOut[ toUse, ], 2, sum, na.rm=T)
	cat( "\nNormalized Counts per Time Point:   ", round( mean( normTotal)))
	cat( "\nAverage Counts per Insertion Site:  ", round( mean( normTotal/length(toUse)), digits=2), "\n")
	out <- data.frame( tbl[ ,1:(from-1)], mOut, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	return( out)
}


calculateAllFitness <- function( tbl, samples, control.genes=NULL, needDropNormalize=TRUE,
				offset.readCount=NULL) {

	# neww idea:  we start from the raw counts, and do the 'drop' and normalize on just the active pair of samples
	# combine the un-normalized count data and the sample key, and do the fitness calc for all samples
	samples <- subset( samples, ! Exclude)
	NS <- nrow(samples)
	cat( "\nCalculating fitness values..")
	out <- tbl[ , 1:(FIRST_COUNT_COLUMN-1)]

	for ( i in 1:NS) {
		thisID <- samples$SampleID[i]
		thisCond <- samples$Condition[i]
		thisRepl <- samples$Replicate[i]
		cat( "\n\nProcessing:    Sample=", thisID, "\tCondition=", thisCond, "\tReplicate=", thisRepl, "\n")

		fitData <- calculateOneFitness( tbl, col1=samples$Time1.Name[i], col2=samples$Time2.Name[i],
					expansion.factor=samples$ExpansionFactor[i], control.genes=control.genes,
					needDropNormalize=needDropNormalize, offset.readCount=offset.readCount)

		cat( "\n\nFitness Summary: ", i, samples$SampleID[i], "\n")
		print( summary( as.numeric(fitData)))

		out <- cbind( out, "FitData"=fitData, stringsAsFactors=FALSE)
		colnames(out)[ ncol(out)] <- samples$SampleID[i]
	}
	cat( "\nDone.\n")
	return( out)
}


calculateOneFitness <- function( tbl, col1, col2, expansion.factor=2, control.genes=NULL, needDropNormalize=TRUE,
				offset.readCount=0) {

	# implement the fitness calculation from "Methods Online Van Opijnem Tn-seq paper 2009"
	# Nature Methods.   2009 October, 6(10) 767-772

	where.col1 <- match( col1, colnames(tbl), nomatch=0)
	where.col2 <- match( col2, colnames(tbl), nomatch=0)
	if ( any( c( where.col1, where.col2) == 0)) {
		cat( "\nError:  named column(s) not found in table:  ", col1, col2)
		stop()
	}

	# get those 2 columns of data
	x1 <- as.numeric( tbl[[where.col1]])
	x2 <- as.numeric( tbl[[where.col2]])

	# if needed, do the 'drop' and normalize here
	from <- FIRST_COUNT_COLUMN
	if (needDropNormalize) {
		tmpTbl <- tbl[ , 1:(FIRST_COUNT_COLUMN-1)]
		tmpTbl$Time1 <- x1
		tmpTbl$Time2 <- x2
		tmpTbl <- dropLowCountInsertions( tmpTbl)
		tmpTbl <- normalizeInsertionCounts( tmpTbl, max.scaling=2)
		tbl <- tmpTbl
		x1 <- as.numeric( tbl[[ FIRST_COUNT_COLUMN ]])
		x2 <- as.numeric( tbl[[ FIRST_COUNT_COLUMN+1 ]])
	}
	to <- ncol(tbl)
	m <- as.matrix( tbl[ , from:to])
	toUse <- which( tbl$DROP == FALSE)
	x1.original <- x1[ toUse]
	x2.original <- x2[ toUse]

	# allow a fixed offset of reads, added to all these non-dropped sites.
	# allows us to estimate read count sensetivity, and control for effects of very low counts
	if ( is.null( offset.readCount)) {
		AUTO_OFFSET <- TRUE
		offset.readCount <- 0
	} else {
		AUTO_OFFSET <- FALSE
	}

	repeat {
		x1 <- x1.original + offset.readCount
		x2 <- x2.original + offset.readCount
		# turn the counts to frequencies
		freq1 <- x1 / sum(x1, na.rm=T)
		freq2 <- x2 / sum(x2, na.rm=T)
		# prevent divide by zeros
		smallF <- 0
		if ( any( c( freq1, freq2) <= 0)) {
			# find the smallest non-zero value, (presumable from N=1 read), 
			# possibly by some multiplier,  and use that as a linear offset
			smallF <- min( freq1[ freq1 > 0], freq2[ freq2 > 0], na.rm=T) * 1.0
			freq1 <- freq1 + smallF
			freq2 <- freq2 + smallF
		}
		N1 <- freq1
		N2 <- freq2
		# scale and make the '1 minus' terms too
		N1.scaled <- expansion.factor / N1
		m1N1.scaled <- expansion.factor / ( 1 - N1)
		m1N2 <- (1 - N2)
		# make the two final log terms
		top.term <- log( N2 * N1.scaled)
		bottom.term <- log( m1N2 * m1N1.scaled)
		ans <- top.term / bottom.term

		# correct for a control set?
		if ( ! is.null( control.genes) && (nCG <- length( control.genes))) {
			isCTRL <- which( tbl$GENE_ID %in% control.genes)
			ctrlFits <- ans[ isCTRL]
			avgCtrl <- median( ctrlFits)
			cat( "\nRe-scale Fitness against ", length(isCTRL), " Insertions from ", nCG, " Genes.")
			if (smallF > 0) cat( "\nRe-Scaling Factor= ", avgCtrl, "\tSmall Linear Offset= ", smallF)
			if (offset.readCount > 0) cat( "\nRe-Scaling Factor= ", avgCtrl, "\tFixed Read Count Offset= ", offset.readCount)
			ans <- ans + (1 - avgCtrl)
		} else {
			# use a global normalization based on all genes/sites
			avgCtrl <- median( ans)
			cat( "\nRe-scale Fitness against ALL Insertions from ALL Genes.")
			if (smallF > 0) cat( "\nRe-Scaling Factor= ", avgCtrl, "\tSmall Linear Offset= ", smallF)
			if (offset.readCount > 0) cat( "\nRe-Scaling Factor= ", avgCtrl, "\tFixed Read Count Offset= ", offset.readCount)
			cat( "\nRe-Scaling Factor= ", avgCtrl, "\tSmall Linear Offset= ", smallF)
			ans <- ans + (1 - avgCtrl)
		}

		# we may be finding the optimal offset
		if ( ! AUTO_OFFSET) break
		# see if less than K % of the data is extreme.  Fitness is centered around 1.0
		v <- ans[ !is.na( ans)]
		nLow <- sum( v < 0)
		nHigh <- sum( v > 2)
		nExtreme <- nLow + nHigh
		pctExtreme <- nExtreme / length(v)
		if ( pctExtreme < 0.05) break
		# not good enough, go arond again
		cat( "\n\nAuto finding the optimal read count offset.  Current extreme fitness pct: ", pctExtreme)
		offset.readCount <- offset.readCount + 5
		if ( offset.readCount > 1000) break
		cat( "\nRaising offset read count to: ", offset.readCount)
	}

	# we can round a bit, as the precision is excessive
	ans <- round( ans, digits=5)
	out <- rep.int( NA, nrow(tbl))
	out[ toUse] <- ans
	return( out)
}


calculateGeneFitness <- function( path=WIG_FOLDER, fitData=fitData, samples=sampleKey, 
					keepIntergenics=FALSE, makeHTML=TRUE, nGenePlots=N_HTML_GENES) {

	# we will build up results for each gene, by combining replicate, for each condition
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", fitData$GENE_ID, fixed=T)
		if ( length( drops)) {
			fitData <- fitData[ -drops, ]
		}
	}

	geneFac <- factor( fitData$GENE_ID)
	allGenes <- levels( geneFac)
	NG <- length( allGenes)
	geneFitScores <- vector( mode="list", length=NG)
	names(geneFitScores) <- allGenes
	whereG <- match( allGenes, fitData$GENE_ID)

	# get the conditions and replicates
	samples <- subset( samples, !Exclude)
	conditions <- sort( unique( samples$Condition))
	NC <- length(conditions)

	ans <- vector( mode="list", length=NC)
	names(ans) <- conditions

	cat( "\nSummarize across replicates for each Gene..")
	for ( ic in 1:NC) {
		thisCondition <- conditions[ic]
		cat( "\n  ", thisCondition)
		myReplicates <- sort( unique( samples$Replicate[ which( samples$Condition == thisCondition)]))
		NR <- length(myReplicates)

		# rezero any gene data for each condition
		for ( ig in 1:NG) geneFitScores[[ig]] <- vector()

		# now accumulate all those scores
		for ( ir in 1:NR) {
			thisReplicate <- myReplicates[ir]
			mySampleRow <- which( samples$Condition == thisCondition & samples$Replicate == thisReplicate)
			myFitScores <- fitData[[ which( colnames(fitData) == samples$SampleID[ mySampleRow]) ]]
	
			tapply( 1:nrow(fitData), geneFac, function(x) {
					myGene <- fitData$GENE_ID[ x[1]]
					where <- match( myGene, allGenes)
					myScores <- myFitScores[x]
					geneFitScores[[where]] <<- c( geneFitScores[[where]], myScores)
					return(NULL)
			})
		}

		# with all replicates combined, we can assess mean, etc.
		meanFit <- sapply( geneFitScores, GENE_AVG_FUN, na.rm=T)
		sdFit <- sapply( geneFitScores, sd, na.rm=T)
		pvalFit <- sapply( geneFitScores, robust.t.test)
		nFit <- sapply( geneFitScores, function(x) sum( !is.na(x)))
		ttlSites <- round( sapply( geneFitScores, length) / NR)
		avgFitsPerGene <- round( nFit / NR, digits=2)

		meanFit <- round( meanFit, digits=4)
		sdFit <- round( sdFit, digits=4)
		pvalFit <- round( pvalFit, digits=6)

		allSymb <- fitData$SYMBOL[ whereG]
		allProd <- fitData$PRODUCT[ whereG]

		#pvalFit <- p.adjust( pvalFit, method="BH")

		out <- data.frame( "GENE_ID"=allGenes, "SYMBOL"=allSymb, "PRODUCT"=allProd, "Total_Sites"= ttlSites, "Avg_Fits_per_Replicate"=avgFitsPerGene,
				"N_Fit_Values"=nFit, "AVG_Fitness"=meanFit, "SD_Fitness"=sdFit, "PVALUE"=pvalFit, stringsAsFactors=F)
		# revert these genes into chromosomal order
		ord <- order( match( allGenes, gmap$GENE_ID))
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)

		outfile <- file.path( path, paste( thisCondition, "GeneFitness.csv", sep="."))
		write.table( out, outfile, sep=",", quote=T, row.names=F)
		if ( makeHTML) makeGeneFitnessHTML( out, outfile, path=path, condition=thisCondition, N=nGenePlots)

		ans[[ic]] <- out
	}
	return( ans)
}


makeGeneFitnessHTML <- function( tbl, csvfile, path=WIG_PATH, condition="", N=N_HTML_GENES, tail.width=TAIL_WIDTH) {

	# make the HTML file name
	htmlfile <- sub( "csv$", "html", csvfile)
	htmlPath <- file.path( path, localPath <- paste( condition, "GenePlots", sep="."))
	if ( ! file.exists( htmlPath)) dir.create( htmlPath)

	# re-order to show the most interesting genes
	ord <- diffExpressRankOrder( tbl$AVG_Fitness, tbl$PVALUE, notDE.value=1.0)

	# tweak the column names to better view in HTML
	colnames(tbl)[4:ncol(tbl)] <- gsub( "_", " ", colnames(tbl)[4:ncol(tbl)])

	tblUP <- tbl[ ord, ]
	tblUP <- tblUP[ 1:min(N,nrow(tbl)), ]
	rownames(tblUP) <- 1:nrow(tblUP)
	upFile <- sub( "GeneFitness", "UP.Genes", htmlfile)
	table2html( tblUP, upFile, title=paste( "Top Genes with Increased fitness:  ", condition), linkPaths=localPath)

	tblDOWN <- tbl[ rev(ord), ]
	tblDOWN <- tblDOWN[ 1:min(N,nrow(tbl)), ]
	rownames(tblDOWN) <- 1:nrow(tblDOWN)
	downFile <- sub( "GeneFitness", "DOWN.Genes", htmlfile)
	table2html( tblDOWN, downFile, title=paste( "Top Genes with Decreased fitness:  ", condition), linkPaths=localPath)

	genesToPlot <- sort( c( tblUP$GENE_ID, tblDOWN$GENE_ID))
	checkX11( width=10, height=7)
	par( mfrow=c(1,1))

	cat( "\nMaking plots..\n")
	for ( i in 1:length(genesToPlot)) {
		thisGene <- genesToPlot[i]
		visualizeGeneFitness( thisGene, condition=condition, path=path, tail.width=tail.width)
		plotFile <- paste( thisGene, "png", sep=".")
		plotFile <- file.path( htmlPath, file.cleanSpecialCharactersFromFileName(plotFile))
		dev.print( png, plotFile, width=1000)
		cat( "\r", i, thisGene)
	}
}


visualizeGeneFitness <- function( gene, condition, path=WIG_FOLDER, tail.width=TAIL_WIDTH, clipFitness=c( -0.2, 3.2)) {

	whereMap <- match( gene, gmap$GENE_ID, nomatch=0)
	if ( ! whereMap) stop( paste( "Gene not found: ", gene))
	leftEdge <- max( gmap$POSITION[whereMap] - tail.width, 1)
	rightEdge <- min( gmap$END[whereMap] + tail.width, max(gmap$POSITION))
	gmap <- subset( gmap, POSITION <= rightEdge & END >= leftEdge)
	whereMap <- match( gene, gmap$GENE_ID, nomatch=0)

	condFile <- file.path( path, paste( condition, "GeneFitness.csv", sep="."))
	if ( ! file.exists( condFile)) stop( paste( "GeneFitness results file not found: ", condFile))
	#geneTbl <- read.csv( condFile, as.is=T)
	geneTbl <- fileBuffering( "read", condFile)
	geneTbl <- subset( geneTbl, GENE_ID %in% gmap$GENE_ID)
	if ( ! nrow( geneTbl)) stop( "No gene fitness results for this gene region.")

	fitFile <- file.path( path, "FitScores.csv")
	if ( ! file.exists( fitFile)) stop( paste( "Fitness scores file not found: ", fitFile))
	#fitTbl <- read.csv( fitFile, as.is=T)
	fitTbl <- fileBuffering( "read", fitFile)
	fitTbl <- subset( fitTbl, GENE_ID %in% gmap$GENE_ID)
	if ( ! nrow( fitTbl)) stop( "No insertion site fitness results for this gene region.")
	
	sampleFile <- file.path( path, "SampleKey.csv")
	if ( ! file.exists( sampleFile)) stop( paste( "Sample Key file not found: ", sampleFile))
	#samples <- read.csv( sampleFile, as.is=T)
	samples <- fileBuffering( "read", sampleFile)
	samples <- subset( samples, Condition == condition)
	samples <- subset( samples, ! Exclude)
	if ( ! nrow( samples)) stop( paste( "No samples found for this condition: ", condition))

	fitM <- as.matrix( fitTbl[ , match( make.names(samples$SampleID), colnames(fitTbl))])
	fitMeans <- apply( fitM, 1, REPLICATE_AVG_FUN, na.rm=T)
	fitMeans <- fitMeans[ ! is.nan( fitMeans)]

	# clip the displayed values to keep things visually sane
	fitMeans[ fitMeans < clipFitness[1]] <- clipFitness[1]
	fitMeans[ fitMeans > clipFitness[2]] <- clipFitness[2]
	yRange <- range( c( -0.3, 0, 2.0, fitMeans), na.rm=T)
	par( mai=c( 1,1,0.8,0.4))
	gSymbol <- gmap$NAME[whereMap]
	if ( gSymbol == gene) gSymbol <- ""
	if ( gSymbol != "") gSymbol <- paste( "   (", gSymbol, ") ", sep="")
	mainText <- paste( "Condition:  ", condition, "      Gene: ", gene, gSymbol, "\n", gmap$PRODUCT[whereMap])

	plot( 1,1, type="n", main=mainText, xlim=c(leftEdge, rightEdge), ylim=yRange, xlab="Chromosomal Location  (bp)",
			ylab="Fitness Score  (W)", xaxs="i", font.axis=2, font.lab=2, cex.lab=1.1)
	barwidth <- max( 16, min( 100, diff( fitTbl$POSITION)))

	# start with a line at zero and one
	lines( c( leftEdge-10, rightEdge+10), c( 0,0), col='gray30', lwd=2, lty=1)
	lines( c( leftEdge-10, rightEdge+10), c( 1,1), col='gray50', lwd=2, lty=2)

	# mark all the posible sites
	fitSites <- fitTbl$POSITION
	nSites <- length( fitSites)
	text( fitSites, rep.int(0,nSites), rep.int("x",nSites), cex=0.8, col='brown')
	nSitesThisGene <- length( fitTbl$POSITION[ fitTbl$GENE_ID == gene])

	# get the gene averages, and P-values, and make up some colors based on both magnitude & significance
	wh <- match( gmap$GENE_ID, geneTbl$GENE_ID, nomatch=NA)
	gmap$AVG_FIT <- geneTbl$AVG_Fitness[wh]
	gmap$PVALUE <- geneTbl$PVALUE[wh]
	gmap$COLOR <- 'brown'
	# if the P-value is terrible, use green for not different
	gmap$COLOR[ gmap$PVALUE >= 0.25] <- 'dodgerblue'
	gmap$COLOR[ gmap$AVG_FIT >= 0.96 & gmap$AVG_FIT <= 1.04] <- 'dodgerblue'
	gmap$COLOR[ gmap$PVALUE < 0.25 & gmap$AVG_FIT >= 1.04] <- 'cyan'
	gmap$COLOR[ gmap$PVALUE < 0.25 & gmap$AVG_FIT <= 0.96] <- 'purple'
	gmap$COLOR[ gmap$PVALUE < 0.05 & gmap$AVG_FIT >= 1.1] <- 'yellow'
	gmap$COLOR[ gmap$PVALUE < 0.05 & gmap$AVG_FIT <= 0.90] <- 'deeppink'
	fitCategories <<- c( "Significant Advantage", "Mild Advantage", "Neutral", "Mild Disadvantage", "Significant Disadvantage")
	fitColors <<- c( "yellow", "cyan", "dodgerblue", "purple", "deeppink")

	# if no valid fitness data, brown it out
	gmap$COLOR[ is.na( gmap$AVG_FIT)] <- 'brown'

	# draw the averges
	gptrs <- which( gmap$REAL_G)
	for ( i in 1:length(gptrs)) {
		ig <- gptrs[i]
		lines( c(gmap$POSITION[ig],gmap$END[ig]), rep.int(gmap$AVG_FIT[ig],2), col=1, lwd=2, lty=1)
		#if ( ig < max(gptrs)) {
		#	# only draw between when its the expected gene layout...
		#	xleft <- gmap$END[ig]
		#	xright <- gmap$POSITION[gptrs[i+1]]
		#	if ( xright-xleft > -10) lines( c(xleft,xright), gmap$AVG_FIT[ gptrs[i:(i+1)]], col=1, lwd=2, lty=1)
		#}
	}

	# draw the bars; and gather what we need to show the P-value on the fly
	halfbar <- barwidth / 2
	nBarsThisGene <- 0
	allFitsThisGene <- vector()
	for ( i in 1:nrow(fitTbl)) {
		thisPos <- fitTbl$POSITION[i]
		thisFitSet <- fitM[ i, ]
		thisFitSet <- thisFitSet[ !is.na( thisFitSet)]
		if ( length( thisFitSet)) {
			if (fitTbl$GENE_ID[i] == gene) {
				nBarsThisGene <- nBarsThisGene + 1
				allFitsThisGene <- c( allFitsThisGene, thisFitSet)
			}
		} else {
			next
		}
		thisAvg <- REPLICATE_AVG_FUN( thisFitSet)
		thisGptr <- match( fitTbl$GENE_ID[i], gmap$GENE_ID)
		rect( thisPos-halfbar, 0, thisPos+halfbar, thisAvg, col=gmap$COLOR[thisGptr], border=gmap$COLOR[thisGptr])
		errorBar( thisFitSet, mode="se", average.FUN=REPLICATE_AVG_FUN, at=thisPos, whisker=halfbar/2)
	}
	avgFitThisGene <- REPLICATE_AVG_FUN( allFitsThisGene)
	pvalThisGene <- robust.t.test( allFitsThisGene)

	# draw the gene arrows, they are fat so trim the end points a bit
	arrowheadlength <- 0.28 / sqrt(length(gptrs))
	arrowwidth <- round( 16 / sqrt(length(gptrs)))
	arrowtextcex <- 2 / sqrt(length(gptrs))
	for (ig in (gptrs <- which( gmap$REAL_G))) {
		if ( gmap$STRAND[ig] == "+") {
			from <- gmap$POSITION[ig] + 10
			to <- gmap$END[ig] - 10
		} else {
			to <- gmap$POSITION[ig] + 10
			from <- gmap$END[ig] - 10
		}
		arrows( from, -0.15, to, -0.15, col=1, lwd=arrowwidth, lty=1, length=arrowheadlength)
		arrows( from, -0.15, to, -0.15, col=gmap$COLOR[ig], lwd=arrowwidth*0.7, lty=1, length=arrowheadlength)
		textAt <- (from+to)/2
		if ( textAt < leftEdge) textAt <- leftEdge
		if ( textAt > rightEdge) textAt <- rightEdge
		text( textAt, -0.2, gmap$GENE_ID[ig], pos=1, cex=arrowtextcex)
	}

	# some info legends
	legend( "top", paste( "Site Usage: ", nBarsThisGene, "of", nSitesThisGene), pch='x', col='brown', bg='white')
	legend( "topright", c( paste( "Gene Fitness: ", round(avgFitThisGene,digits=2)), 
					paste( "P-value: ", formatC( pvalThisGene, format="e", digits=3))), 
					bg='white')
	legend( "topleft", fitCategories, col=fitColors, lwd=3, bg='white', cex=0.8) 

	dev.flush()
}


compareConditionFitness <- function( path=WIG_FOLDER, samples=sampleKey, min.delta=0.1, min.sites=1, min.fits=2,
					makeHTML=TRUE, nGenePlots=N_HTML_GENES, keepIntergenics=FALSE) {

	cat( "\n\nComparing Fitness between conditions..")

	# reload the fitness data
	fitData <- read.csv( file.path( path, "FitScores.csv"), as.is=T)
	# reset the column names to look like Sample IDs
	colnames(fitData) <- gsub( ".", "-", colnames(fitData), fixed=T)

	# do all pairs of conditions
	allConditions <- sort( unique( samples$Condition))
	NC <- length( allConditions)
	for ( i in 1:NC) for ( j in 1:NC) {
		if ( i == j) next
		thisCond <- allConditions[i]
		thatCond <- allConditions[j]
		cat( "\nCond_1=", thisCond, "\tCond_2=", thatCond)

		ans <- compareConditions( fitData, thisCond, thatCond, path1=path, samples=samples, 
					min.delta=min.delta, min.sites=min.sites, min.fits=min.fits, 
					keepIntergenics=keepIntergenics)
		f <- paste( "DeltaFit", thisCond, "v", thatCond, "Genes.csv", sep=".")
		f <- file.path( path, f)
		write.table( ans, f, sep=",", quote=T, row.names=F)
		cat( "\nWrote file:  ", f, "\tN_Genes: ", nrow(ans))

		if ( makeHTML && nrow(ans)) {
			makeCompareFitnessHTML( ans, csvfile=f, path=path, cond1=thisCond, cond2=thatCond,
						N=nGenePlots, tail.width=TAIL_WIDTH)
		}
	}
	# put the display back to normal
	dev.off()
	checkX11( width=10, height=7)
	par( mfrow=c(1,1))
}


# this comparison just used Gene Results...
compareConditions.oldMethod <- function( condition1, condition2, path1=WIG_FOLDER, path2=path1, min.delta=0.1,
				min.sites=1, min.fits=2, keepIntergenics=FALSE) { 

	cond1File <- file.path( path1, paste( condition1, "GeneFitness.csv", sep="."))
	if ( ! file.exists( cond1File)) stop( paste( "GeneFitness results file not found: ", cond1File))
	geneTbl1 <- read.csv( cond1File, as.is=T)
	cond2File <- file.path( path2, paste( condition2, "GeneFitness.csv", sep="."))
	if ( ! file.exists( cond2File)) stop( paste( "GeneFitness results file not found: ", cond2File))
	geneTbl2 <- read.csv( cond2File, as.is=T)

	# the 2 files should be identical layout, but force it
	genes1 <- geneTbl1$GENE_ID
	genes2 <- geneTbl2$GENE_ID
	if ( ! all( genes1 == genes2)) {
		both <- intersect( genes1, genes2)
		geneTbl1 <- geneTbl1[ match( both, genes1), ]
		geneTbl2 <- geneTbl2[ match( both, genes2), ]
		genes1 <- geneTbl1$GENE_ID
		genes2 <- geneTbl2$GENE_ID
	}

	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes1, fixed=TRUE)
		if ( length( drops)) {
			geneTbl1 <- geneTbl1[ -drops, ]
			geneTbl2 <- geneTbl2[ -drops, ]
			genes1 <- geneTbl1$GENE_ID
			genes2 <- geneTbl2$GENE_ID
		}
	}

	deltaFit <- geneTbl1$AVG_Fitness - geneTbl2$AVG_Fitness
	bestPval <- apply( matrix( c(geneTbl1$PVALUE, geneTbl2$PVALUE), length(genes1), 2), MARGIN=1, min)
	ord <- diffExpressRankOrder( deltaFit, bestPval, notDE.value=0, wt.pvalue=2)

	# grab the ID, symbol, product, site count, fit count, and then the mean & p-value from both
	# there are now 6 fixed columns to carry forward
	out <- data.frame( geneTbl1[ ,1:6], "Delta_Fitness"=deltaFit, "Best_PVALUE"=bestPval, geneTbl1[,c("AVG_Fitness","PVALUE")], 
			geneTbl2[,c("AVG_Fitness","PVALUE")], stringsAsFactors=F)
	colnames(out)[9:12] <- paste( rep( c( condition1, condition2), each=2), rep( c("AVG_Fitness","PVALUE"), times=2), sep="_")
	out <- out[ ord, ]
	tmpSAVE <<- out

	# only keep those with some difference
	keep <- which( abs( out$Delta_Fitness) >= min.delta)
	out <- out[ keep, ]

	# this is very sensitive to genes with almost no data, also allow filter by sites and fits
	keepS <- which( out$Avg_Fits_per_Replicate >= min.sites)
	keepF <- which( out$N_Fit_Values >= min.fits)
	keep <- intersect( keepS, keepF)
	out <- out[ keep, ]

	if (nrow(out)) rownames(out) <- 1:nrow(out)
	cat( "\nCompare Conditions:  ", condition1, "vs", condition2)
	cat( "\nN_Genes Kept: ", nrow(out))
	return( out)
}

compareConditions <- function( fitData, condition1, condition2, path1=WIG_FOLDER, path2=path1, samples=sampleKey, 
				min.delta=0.1, min.sites=1, min.fits=2, keepIntergenics=FALSE) { 
	
	# we will build up results for each gene, by combining replicate, for each condition
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", fitData$GENE_ID, fixed=T)
		if ( length( drops)) {
			fitData <- fitData[ -drops, ]
		}
	}

	geneFac <- factor( fitData$GENE_ID)
	allGenes <- levels( geneFac)
	NG <- length( allGenes)
	cond1FitScores <- cond2FitScores <- vector( mode="list", length=NG)
	names(cond1FitScores) <- names(cond2FitScores) <- allGenes
	whereG <- match( allGenes, fitData$GENE_ID)

	# get the conditions and replicates
	samples <- subset( samples, !Exclude)
	isCond1 <- which( samples$Condition == condition1)
	isCond2 <- which( samples$Condition == condition2)
	NR1 <- length( isCond1)
	NR2 <- length( isCond2)

	# now accumulate all those scores in both conditions
	for ( ir in 1:NR1) {
		mySampleRow <- isCond1[ir]
		myFitScores <- fitData[[ which( colnames(fitData) == samples$SampleID[ mySampleRow]) ]]
		tapply( 1:nrow(fitData), geneFac, function(x) {
				myGene <- fitData$GENE_ID[ x[1]]
				where <- match( myGene, allGenes)
				myScores <- myFitScores[x]
				cond1FitScores[[where]] <<- c( cond1FitScores[[where]], myScores)
				return(NULL)
		})
	}
	for ( ir in 1:NR2) {
		mySampleRow <- isCond2[ir]
		myFitScores <- fitData[[ which( colnames(fitData) == samples$SampleID[ mySampleRow]) ]]
		tapply( 1:nrow(fitData), geneFac, function(x) {
				myGene <- fitData$GENE_ID[ x[1]]
				where <- match( myGene, allGenes)
				myScores <- myFitScores[x]
				cond2FitScores[[where]] <<- c( cond2FitScores[[where]], myScores)
				return(NULL)
		})
	}

	# with all replicates combined, we can assess mean, etc.
	meanFit1 <- sapply( cond1FitScores, GENE_AVG_FUN, na.rm=T)
	meanFit2 <- sapply( cond2FitScores, GENE_AVG_FUN, na.rm=T)
	nFit1 <- sapply( cond1FitScores, function(x) sum( !is.na(x)))
	nFit2 <- sapply( cond2FitScores, function(x) sum( !is.na(x)))
	sdFit1 <- sapply( cond1FitScores, sd, na.rm=T)
	sdFit2 <- sapply( cond2FitScores, sd, na.rm=T)
	pvalFit <- mapply( FUN=robust.2.t.test, cond1FitScores, cond2FitScores)
	ttlSites <- round( sapply( cond1FitScores, length) / NR1)
	avgFitsPerGene1 <- round( nFit1 / NR1, digits=2)
	avgFitsPerGene2 <- round( nFit2 / NR2, digits=2)
	deltaFit <- meanFit1 - meanFit2

	deltaFit <- round( deltaFit, digits=4)
	meanFit1 <- round( meanFit1, digits=4)
	meanFit2 <- round( meanFit2, digits=4)
	sdFit1 <- round( sdFit1, digits=4)
	sdFit2 <- round( sdFit2, digits=4)
	pvalFit <- round( pvalFit, digits=6)

	allSymb <- fitData$SYMBOL[ whereG]
	allProd <- fitData$PRODUCT[ whereG]

	#pvalFit <- p.adjust( pvalFit, method="BH")

	out <- data.frame( "GENE_ID"=allGenes, "SYMBOL"=allSymb, "PRODUCT"=allProd, "Total_Sites"= ttlSites, 
			"Delta_Fitness"=deltaFit, "PVALUE"=pvalFit, 
			"Avg_Fits1"=avgFitsPerGene1, "N_Fit_Values1"=nFit1, "AVG_Fitness1"=meanFit1, "SD_Fitness1"=sdFit1, 
			"Avg_Fits2"=avgFitsPerGene2, "N_Fit_Values2"=nFit2, "AVG_Fitness2"=meanFit2, "SD_Fitness2"=sdFit2, 
			stringsAsFactors=F)

	# revert these genes into chromosomal order
	ord <- order( match( allGenes, gmap$GENE_ID))
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	# only keep those with some difference
	if ( ! is.null( min.delta)) {
		keep <- which( abs( out$Delta_Fitness) >= min.delta)
		cat( "\nTest for 'min.delta': ", length(keep), "of", nrow(out))
		out <- out[ keep, ]
	}

	# this is very sensitive to genes with almost no data, also allow filter by sites and fits
	if ( ! is.null( min.sites)) {
		keepS <- which( out$Total_Sites >= min.sites)
		cat( "\nTest for 'min.sites': ", length(keepS), "of", nrow(out))
		out <- out[ keepS, ]
	}
	
	if ( ! is.null( min.fits)) {
		nFits <- pmin( out$N_Fit_Values1, out$N_Fit_Values2)
		keepF <- which( nFits >= min.fits)
		cat( "\nTest for 'min.fits': ", length(keepF), "of", nrow(out))
		out <- out[ keepF, ]
	}
	if (nrow(out)) rownames(out) <- 1:nrow(out)
	cat( "\nCompare Conditions:  ", condition1, "vs", condition2)
	cat( "\nN_Genes Kept: ", nrow(out))
	return( out)
}


makeCompareFitnessHTML <- function( tbl, csvfile, path=WIG_PATH, cond1="", cond2="", N=N_HTML_GENES, 
					tail.width=TAIL_WIDTH) {

	# make the HTML file name
	htmlfile <- sub( "csv$", "html", csvfile)
	htmlPath <- file.path( path, localPath <- paste( cond1, "v", cond2, "GenePlots", sep="."))
	if ( ! file.exists( htmlPath)) dir.create( htmlPath)

	# re-order to show the most interesting genes
	ord <- diffExpressRankOrder( tbl$Delta_Fitness, tbl$PVALUE, notDE.value=0, wt.pvalue=2)

	tblUP <- tbl[ ord, ]
	tblUP <- tblUP[ 1:min(N,nrow(tbl)), ]
	tblUP <- subset( tblUP, Delta_Fitness > 0)
	rownames(tblUP) <- 1:nrow(tblUP)
	colnames(tblUP)[4:ncol(tblUP)] <- gsub( "_", " ", colnames(tblUP)[4:ncol(tblUP)])
	upFile <- sub( "Genes", "UP.Genes", htmlfile)
	table2html( tblUP, upFile, title=paste( "Top Genes with Increased fitness:  ", cond1, " vs ", cond2), linkPaths=localPath)

	genesToPlot <- sort( tblUP$GENE_ID)
	dev.off()
	checkX11( width=10, height=10)

	cat( "\nMaking plots..\n")
	for ( i in 1:length(genesToPlot)) {
		par( mfrow=c(2,1))
		thisGene <- genesToPlot[i]
		visualizeGeneFitness( thisGene, condition=cond1, path=path, tail.width=tail.width)
		visualizeGeneFitness( thisGene, condition=cond2, path=path, tail.width=tail.width)
		plotFile <- paste( thisGene, "png", sep=".")
		plotFile <- file.path( htmlPath, file.cleanSpecialCharactersFromFileName(plotFile))
		dev.print( png, plotFile, width=1000)
		cat( "\r", i, thisGene)
	}
}


robust.t.test <- function( x) {

	#gracefully handle anything that might break t.test
	x <- x[ !is.na(x)]
	if ( length(x) < 2) return(1)
	if ( diff( range( x)) < 0.001) x <- jitter(x)
	return( t.test( x, mu=1.0)$p.value)
}


robust.2.t.test <- function( x, y) {

	#gracefully handle anything that might break t.test
	x <- x[ !is.na(x)]
	y <- y[ !is.na(y)]
	if ( length(x) < 2) return(1)
	if ( length(y) < 2) return(1)
	if ( diff( range( x)) < 0.001) x <- jitter(x)
	if ( diff( range( y)) < 0.001) y <- jitter(y)
	return( t.test( x, y)$p.value)
}


fileBuffering <- function( mode=c("read", "setup"), fname, sep=",") {

	mode <- match.arg( mode)

	if ( mode == "setup") {
		fileBufferList <<- vector( mode="list")
		return()
	}

	if ( mode == "read") {
		if ( ! exists( "fileBufferList")) fileBuffering( "setup")
		curLength <- length( fileBufferList)
		curNames <- names( fileBufferList)
		where <- match( fname, curNames)
		if ( is.na( where)) {
			tmp <- read.delim( fname, as.is=TRUE, sep=sep)
			curLength <- curLength + 1
			fileBufferList[[ curLength]] <<- tmp
			names( fileBufferList)[ curLength] <<- fname
			return(tmp)
		}
		return( fileBufferList[[ where]])
	}
}
