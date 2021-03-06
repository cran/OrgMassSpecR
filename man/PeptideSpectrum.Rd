\name{PeptideSpectrum}
\alias{PeptideSpectrum}
\title{Plot Annotated Peptide Fragmentation Mass Spectrum}
\description{Identify the b- and y-ions or c- and z-ions in a peptide fragmentation mass spectrum given the sequence.  Generate a plot and table identifying the fragment ions.
}  
\usage{
PeptideSpectrum(expt, theory, t = 0.4, b = 5, label = "", xlim = c(100, 1500), 
                supress = FALSE)
}
\arguments{
    \item{expt}{data frame containing the experimental MS/MS peak list with the m/z values in the first column and corresponding intensities in the second}
    \item{theory}{data frame of theoretical fragment ions generated by \code{\link{FragmentPeptide}}}
    \item{t}{numeric value specifying the m/z tolerance for matching theoretical fragment ions to experimental peaks}
    \item{b}{numeric value specifying the baseline threshold for peak identification.  Expressed as a percent of the maximum intensity.}
    \item{label}{character string to label the spectrum.}
    \item{xlim}{numeric vector of length 2, defining the beginning and ending values of the x-axis.}
    \item{supress}{logical specifying whether or not to supress the peak annotations.}
}
\details{
A graphics window of width = 10 and height = 5.5 seems to work well.  A crude attempt is made to prevent overlapping peak annotations.  If this is unsuccessful, use \code{supress = TRUE} and then add the annotations using the \code{\link{text}} function.  If two or more peak identifications are made within the set tolerance, all identifications will be printed, overlapping, on the plot.
}
\value{
Generates an annotated spectrum and a data frame showing the matches between the experimental peaks and the theoretical fragment ions.  The column names refer to the experimental \emph{m/z} value (\code{expt_mz}), the experimental peak intensity (\code{expt_int}) and the error between the experimental and theoretical \emph{m/z} values (\code{error}).  The other column names are the same as for \code{\link{FragmentPeptide}}.
}
\author{Nathan G. Dodder}
\seealso{\code{\link{FragmentPeptide}}}
\examples{
## plot spectrum of peptide NIDALSGMEGR
t <- FragmentPeptide("NIDALSGMEGR")   # generate theoretical fragment ions
PeptideSpectrum(example.spectrum.peptide, t, label = "CID", xlim = c(100, 1200))
mtext("Peptide fragmentation mass spectrum", line = 1)
}  
