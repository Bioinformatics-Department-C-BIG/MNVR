\name{CSQfilter}
\alias{CSQfilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CSQfilter(OUTPUTcsq = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OUTPUTcsq}{
%%     ~~Describe \code{OUTPUTcsq} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (OUTPUTcsq = NULL)
{
    bcf <- read.table(OUTPUTcsq)
    colnames(bcf)[-ncol(bcf)] <- c("CHR", "POS", "ID", "REF",
        "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    bcf <- add_column(bcf, transcript_num = str_count(bcf$INFO,
        pattern = ",") + 1, .after = "INFO")
    distance_threshold <- 2
    sn <- str_match(OUTPUTcsq, ".+/(.+)\\..+vcf.gz$")[, 2]
    indm <- which(diff(as.numeric(as.factor(bcf[, 1]))) == 0 &
        diff(as.numeric(bcf$POS)) <= distance_threshold)
    ind <- c(indm, (indm + 1))
    ind <- unique(ind[order(ind)])
    bcf_mnvs <- bcf[ind, ]
    bcf_mnvs_PASS <- bcf_mnvs[bcf_mnvs$FILTER == "PASS", ]
    bcf_mnvs_PASS$INFO <- gsub("BCSQ=", "", bcf_mnvs_PASS$INFO)
    bcf_mnvs_PASS <- bcf_mnvs_PASS[grepl(".>.\\+.", bcf_mnvs_PASS$INFO),
        ]
    bcf_mnvs_tidy <- data.frame()
    for (i in 1:nrow(bcf_mnvs_PASS)) {
        temp <- bcf_mnvs_PASS[i, ]
        if (temp$transcript_num > 1) {
            temp_e <- temp
            temp_e[c(1:temp$transcript_num), ] <- temp[1, ]
            temp_e$INFO <- as.vector(str_split(temp$INFO, ",",
                simplify = T))
            temp <- temp_e
        }
        csq <- data.frame(str_split(temp$INFO, "\\|", simplify = T))
        names(csq) <- c("consequence", "gene", "transcript",
            "biotype", "strand", "aa_change", "dna_change")
        temp <- cbind.data.frame(temp, csq)
        bcf_mnvs_tidy <- rbind.data.frame(bcf_mnvs_tidy, temp)
    }
    bcf_mnvs_tidy <- bcf_mnvs_tidy[grepl("missense$", bcf_mnvs_tidy$consequence) |
        grepl("^synonymous$", bcf_mnvs_tidy$consequence) | grepl("^stop_gained$",
        bcf_mnvs_tidy$consequence), ]
    bcf_mnvs_tidy$consequence[bcf_mnvs_tidy$consequence == "*missense"] <- "missense"
    bcf_mnvs_tidy <- bcf_mnvs_tidy[which(bcf_mnvs_tidy$transcript \%in\%
        hgnc$ENselect), ]
    bcf_mnvs_tidy <- na.omit(bcf_mnvs_tidy)
    write.table(bcf_mnvs_tidy, file = paste0(outpath, sn, "_MNV_results.tsv"),
        quote = F, col.names = T, row.names = F)
    sumtab <- as.data.frame.matrix(table(bcf_mnvs_tidy$CHR, bcf_mnvs_tidy$consequence))
    sumtab <- rbind.data.frame(sumtab, (colSums(sumtab)))
    rownames(sumtab)[nrow(sumtab)] <- "Total"
    write.table(sumtab, file = paste0(outpath, sn, "_MNV_Summary.tsv"),
        quote = F, col.names = T, row.names = T)
    ccod <- list(missense = "orange", synonymous = "royalblue",
        stop_gained = "orangered")
    pdf(file = paste0(outpath, sn, "_KP_plot.pdf"), width = 12,
        height = 8)
    plotKaryotype() \%>\% karyoploteR::kpPoints(chr = paste0("chr",
        bcf_mnvs_tidy$CHR), x = bcf_mnvs_tidy$POS, y = 0.2, col = as.vector(unlist(ccod[bcf_mnvs_tidy$consequence]))) \%>\%
        karyoploteR::kpPoints(chr = paste0("chr", bcf_mnvs_tidy$CHR),
            x = bcf_mnvs_tidy$POS, y = 0.2, col = as.vector(unlist(ccod[bcf_mnvs_tidy$consequence])))
    dev.off()
    bcf_mnvs_tidy$CHR <- ordered(bcf_mnvs_tidy$CHR, levels = c(1:22,
        "X", "Y", "MT", paste0("chr", c(1:22, "X", "Y", "M"))))
    p1 <- ggplot(bcf_mnvs_tidy, aes(x = consequence, fill = consequence)) +
        geom_bar(show.legend = F) + theme_bw()
    p2 <- ggplot(bcf_mnvs_tidy, aes(x = CHR, fill = consequence)) +
        geom_bar() + theme_bw() + facet_grid(. ~ consequence)
    pdf(file = paste0(outpath, sn, "_Consequences_plot.pdf"),
        width = 12, height = 8)
    gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
