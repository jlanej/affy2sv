#' Function to create compatible txt files for \code{MAD} or for \code{PennCNV}
#' 
#' This function uses the files created with Ax2APT, metrix, summary and genotype, to create a list containing a map file and a SnpMatrix with the genotype of each probe in the array. It can also be used to create a .tped file, compatible with PLINK (is in charge of the user to create the partner .tfam file).
#' 
#' @param calls.files Path until genotype file.
#' @param summary.file Path until summary file.
#' @param annotation.file NetAffx Annotation database file in CSV format.
#' @param output.name Path where the output files will be stored (for 
#' \code{"mad"} as \code{output.type}) or the name of the generated file (for 
#' \code{"penncnv"} as \code{output.type}).
#' @param output.type By default \code{"mad"}. It can be set to 
#' \code{"penncnv"}. Used to determiner the type of file(s) the function will 
#' generate.
#' @param metrix.file If generated, path until the metrix file.
#' @param metrix.tm By default \code{97.0}. If specified metrix file, 
#' threshold to consider a probe well genotyped.
#' @param verbose By default FALSE. If TRUE the function will shows messages 
#' indicating the process.
#' @examples 
#' \dontrun{
#' # To create files compatible with MAD/R-GADA
#' Ax2Mad(
#' calls.file=paste0(output.path, .Platform$file.sep, "AxiomGT1.calls.txt"),
#' summary.file=paste0(output.path, .Platform$file.sep, "AxiomGT1.summary.txt"), 
#' annotation.file=paste0(path.lib, .Platform$file.sep, "Axiom_GW_Hu_SNP.na34.annot.csv"),
#' output.name="ax.mad",
#' output.type="mad"
#' )
#' }
#' \dontrun{
#' # To create a file compatible with PennCNV
#' Ax2Mad(
#' calls.file=paste0(output.path, .Platform$file.sep, "AxiomGT1.calls.txt"),
#' summary.file=paste0(output.path, .Platform$file.sep, "AxiomGT1.summary.txt"), 
#' annotation.file=paste0(path.lib, .Platform$file.sep, "Axiom_GW_Hu_SNP.na34.annot.csv"),
#' output.name="ax",
#' output.type="penncnv"
#' )
#' }
Ax2Mad <- function( calls.file, summary.file, annotation.file, output.name,
    output.type = "mad", metrix.file = NA, metrix.tm = 97.0,
    verbose = FALSE ) {
    
    # SETUP
    OptionsAffy.Init("Ax2Mad", verbose)
    OptionsAffy.Cores(1)
    # /SETUP

    # SYSTEM
    CheckSystem(bits=FALSE)
    # /SYSTEM


    output.type = value( output.type, c( "mad", "penncnv" ) );
    
    
    # Checking input arguments
    if( metrix.tm > 1 & metrix.tm < 0 ) {
        screen_ax_mad( TRUE, "Argument 'metrix.tm' must be between 0 and 1." )
        stop( "HALT." );
    }
    if( is.na( output.type ) ) {
        screen_ax_mad( TRUE, "Argument 'output.type' must be 'mad' or 'penncnv'." )
        stop( "HALT." );
    }
    if( is.na( output.name ) | output.name == '' ) {
        screen_gw_mad( TRUE, "Argument 'output.name' must be filled." )
        stop( "HALT." );   
    }
    # /Checking input arguments
    
    # Loading the calls file
    ScreenAffy( "Loading calls-file file." );
    calls <- read.table( calls.file, comment.char="#", header=TRUE );

    # Coercing the call file to matrix
    ScreenAffy( "Coercing calls-file to matrix." );
    rownames( calls ) <- calls$probeset_id;
    calls <- as( calls[ , -1 ], "matrix" );
    calls <- calls[ order( rownames( calls ) ), ];

    # Loading the summary file with the normalized intensities
    ScreenAffy( "Loading intensities-file." );
    intensities <- read.table( summary.file, comment.char="#", header=TRUE );

    # Coercing intensities to matrix
    ScreenAffy( "Coercing intensities-file to matrix." );
    rownames( intensities ) <- intensities$probeset_id;
    intensities <- as( intensities[ ,-1 ], "matrix" );

    # Loading the metrix file
    if( !is.na( metrix.file ) ) {
        ScreenAffy( "Loading metrix-file file." );
        metrix <- read.table( metrix.file, comment.char="#", header=TRUE );

        # Selecting only the SNPs with a higher Call Rate than 97%
        # [ This came from the Affymetrix Suplement-Best practices supplement to 
        # Axiom Genotyping Solution Data Analysis User Guide, Pg 22,23 ]
        metrix <- metrix[ metrix[ "CR" ] > metrix.tm, ];
        calls <- calls[ rownames( calls ) %in% metrix$probeset_id, ]
    }
    else {
        ScreenAffy( "No metrix-file file given." );
    }

    
    ScreenAffy( "Loading annotation file." );
    ann <- read.csv ( annotation.file, comment.char="#", header=TRUE );

    ScreenAffy( "Filtering annotations." );
    rownames( ann ) <- ann$Probe.Set.ID;
    ann <- ann[ , c( "Probe.Set.ID", "Chromosome", "Physical.Position" ) ];
    ann <- ann[ order( rownames( ann ) ), ];
    ann <- ann[ rownames( ann ) %in% rownames( calls ), ];


    ScreenAffy( "Filtering intensities." );
    names <- gsub( "\\-B", "", gsub( "\\-A", "", as.character( rownames( intensities ) ) ) );
    intensities <- intensities[ names %in% rownames( calls ), ];


    ScreenAffy( "Getting alleles intensities." );
    alleleA.i <- intensities[ seq( 1, nrow( intensities ), 2 ), ]
    alleleB.i <- intensities[ seq( 2, nrow( intensities ), 2 ), ]
    alleleA.i <- alleleA.i[ order( as.character( rownames( alleleA.i ) ) ), ]
    alleleB.i <- alleleB.i[ order( as.character( rownames( alleleB.i ) ) ), ]


    
    ScreenAffy( "Calculating the genotype centers." );
    centers <- sapply( 1:length( colnames( alleleA.i ) ), CalcCenters, alleleA=alleleA.i, alleleB=alleleB.i, geno=calls );
    colnames( centers ) <- colnames( calls )
    rownames( centers ) <- c( "AA", "AB", "BB" );


    ScreenAffy( "Calculating LRR." );
    # R <- as( alleleA.s + alleleB.s, "matrix" );
    R <- as( alleleA.i + alleleB.i, "matrix" );
    medians <- rowMedians( R );
    LRR <- log2( sweep( R, 1, medians, FUN="/" ) )/4;
    rownames( LRR ) <- rownames( calls );
    # colnames( LRR ) <- colnames( summ );
    colnames( LRR ) <- colnames( intensities );


    ScreenAffy( "Calculating thita." );
    theta <- atan2( alleleB.i,alleleA.i ) * 2 / pi;
    rownames( theta ) <- rownames( calls );

    ScreenAffy( "Calculating BAF." );
    BAF <- sapply( 1:length( colnames( theta ) ), CalcBafCor, theta=theta, centers=centers, geno=calls );
    colnames( BAF ) <- colnames( theta );
    rownames( BAF ) <- rownames( theta );
    

    if( output.type == "mad" ) {
        unlink( output.name, recursive=TRUE, force=TRUE );
        dir.create( output.name, showWarnings=FALSE );
        
        for( ind.i in colnames( intensities ) ) {
            file.name <- file.path( output.name, gsub( ".CEL", "", ind.i ), fsep=.Platform$file.sep );
            #g <- calls[,ind.i]
            #g[ g==-1 ] <- NA;
            #g[ g==0 ] <- "AA";
            #g[ g==1 ] <- "AB";
            #g[ g==2 ] <- "BB";
            total <- cbind( ann, LRR[,ind.i], BAF[,ind.i], updateCalls( calls[ ,ind.i ] ) );
            #total <- cbind( ann, BAF[,ind.i], LRR[,ind.i], g );
            colnames( total ) <- c( "Name", "Chr", "Position", "Log.R.Ratio", "B.Allele.Freq", "GType" );
            
            ScreenAffy( "Writting mad for '", ind.i, "." );
            write.table( total, file=file.name, quote=FALSE, row.names=FALSE, sep="\t", na="NA" );
        }
    }
    else {
        ScreenAffy( "Flatting individuals to PennCNV file." );
        file.name = paste0( output.name, ".penncnv" );
        
        ind.i <- colnames( intensities )[ 1 ];
        
        result <- cbind( ann, BAF[,ind.i], LRR[,ind.i], updateCalls( calls[ ,ind.i ] ) );
        ind.i <- gsub( "CEL", "", ind.i );
        names <- c( "Name", "Chr", "Position", paste0( ind.i, "B Allele Freq" ),
            paste0( ind.i, "Log R Ratio" ), paste0( ind.i, "GType" ) );
        
        for( ind.i in colnames( intensities )[ -1 ] ) {
            result <- cbind( result, BAF[,ind.i], LRR[,ind.i], updateCalls( calls[ ,ind.i ] ) );
            ind.i <- gsub( "CEL", "", ind.i );
            names <- c( names, paste0( ind.i, "B Allele Freq" ),
                paste0( ind.i, "Log R Ratio" ), paste0( ind.i, "GType" ) );
        }
        colnames( result ) <- names;
        
        ScreenAffy( "Writting '", file.name, "'." );
        write.table( result, file=file.name, quote=FALSE, row.names=FALSE, sep="\t", na="NA" );
    }

}

updateCalls <- function( calls ) {
    calls[ calls==-1 ] <- NA;
    calls[ calls==0 ] <- "AA";
    calls[ calls==1 ] <- "AB";
    calls[ calls==2 ] <- "BB";

    return( calls );
}


CalcCenters <- function( ii, alleleA, alleleB, geno ) {
    alleleA <- alleleA[ , ii ];
    alleleB <- alleleB[ , ii ];
    geno <- geno[ , ii ];

    AmAA <- median( alleleA[ geno==0 ] );
    AmAB <- median( alleleA[ geno==1 ] );
    AmBB <- median( alleleA[ geno==2 ] );

    BmAA <- median( alleleB[ geno==0 ] );
    BmAB <- median( alleleB[ geno==1 ] );
    BmBB <- median( alleleB[ geno==2 ] );

    c(  atan2( BmAA, AmAA ) * 2 /pi,
        atan2( BmAB, AmAB ) * 2 /pi,
        atan2( BmBB, AmBB ) * 2 /pi
        );
}



CalcBafCor <- function( ii, theta, centers, geno ) {
    theta <- theta[ , ii ];
    centers <- centers[ , ii ];
    geno <- geno[ , ii ];

    lessAA <- theta < centers[ "AA" ] & geno != 1;
    lessAB <- theta < centers[ "AB" ] & geno != 1;
    lessBB <- theta < centers[ "BB" ] & geno != 1;
    isAB <- geno == 1;

    down  <- lessAA 
    up    <- !lessBB
    mdown <- !lessAA & lessAB
    mup   <- !lessAB & lessBB

    M <- rep( NA, length( theta ) );
    M[ mdown ] <- 0.5 * ( theta[ mdown ] - centers[ "AA" ] ) / ( centers[ "AB" ] - centers[ "AA" ] )
    M[ mup   ] <- 0.5 * ( theta[ mup   ] - centers[ "AB" ] ) / ( centers[ "BB" ] - centers[ "AB" ] ) + 0.5
    M[ down  ] <- 0
    M[ up    ] <- 1
    M[ isAB  ] = theta[ isAB ]; 

    M[ M < 0 ] <- 0
    M[ M > 1 ] <- 1
    M;
}