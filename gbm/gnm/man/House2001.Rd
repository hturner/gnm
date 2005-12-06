\name{House2001}
\alias{House2001}
\docType{data}
\title{ Data on twenty roll calls in the US House of Representatives, 2001 }
\description{
  The voting record of every representative in the 2001 House, on 20
  roll calls selected by \emph{Americans for Democratic Action}.  Each row is
  the record of one representative; the first column records the
  representative's registered party allegiance.
}
\usage{data(House2001)}
\format{
  A data frame with 439 observations on the following 21 variables.
  \describe{
    \item{\code{party}}{a factor with levels \code{D} \code{I} \code{N} \code{R}}
    \item{\code{HR333.BankruptcyOverhaul.Yes}}{a numeric vector}
    \item{\code{SJRes6.ErgonomicsRuleDisapproval.No}}{a numeric vector}
    \item{\code{HR3.IncomeTaxReduction.No}}{a numeric vector}
    \item{\code{HR6.MarriageTaxReduction.Yes}}{a numeric vector}
    \item{\code{HR8.EstateTaxRelief.Yes}}{a numeric vector}
    \item{\code{HR503.FetalProtection.No}}{a numeric vector}
    \item{\code{HR1.SchoolVouchers.No}}{a numeric vector}
    \item{\code{HR1836.TaxCutReconciliationBill.No}}{a numeric vector}
    \item{\code{HR2356.CampaignFinanceReform.No}}{a numeric vector}
    \item{\code{HJRes36.FlagDesecration.No}}{a numeric vector}
    \item{\code{HR7.FaithBasedInitiative.Yes}}{a numeric vector}
    \item{\code{HJRes50.ChinaNormalizedTradeRelations.Yes}}{a numeric vector}
    \item{\code{HR4.ANWRDrillingBan.Yes}}{a numeric vector}
    \item{\code{HR2563.PatientsRightsHMOLiability.No}}{a numeric vector}
    \item{\code{HR2563.PatientsBillOfRights.No}}{a numeric vector}
    \item{\code{HR2944.DomesticPartnerBenefits.No}}{a numeric vector}
    \item{\code{HR2586.USMilitaryPersonnelOverseasAbortions.Yes}}{a numeric vector}
    \item{\code{HR2975.AntiTerrorismAuthority.No}}{a numeric vector}
    \item{\code{HR3090.EconomicStimulus.No}}{a numeric vector}
    \item{\code{HR3000.TradePromotionAuthorityFastTrack.No}}{a numeric vector}
  }
}
\details{
  Coding of the votes is as described in ADA (2002).
}
\source{
  Originally printed in ADA (2002).  Kindly supplied in electronic
  format by Jan deLeeuw, who used the data to illustrate methods
  developed in deLeeuw (2006).
}
\references{
Americans for Democratic Action, ADA (2002).  2001 voting record:
Shattered promise of liberal progress.  \emph{ADA Today} \bold{57}(1),
1--17.

deLeeuw, J (2006).  Principal component analysis of binary data by
iterated singular value decomposition.  \emph{Computational Statistics
  and Data Analysis} \bold{50}, 21--39.
}
\examples{
data(House2001)
summary(House2001)
}
\keyword{datasets}
