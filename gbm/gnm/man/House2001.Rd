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
## Put the votes in a matrix, and discard members with too many NAs etc:
House2001m <- as.matrix(House2001[-1])
informative <- apply(House2001m, 1, function(row){
    valid <- !is.na(row)
    validSum <- if (any(valid)) sum(row[valid]) else 0
    nValid <- sum(valid)
    uninformative <- (validSum == nValid) || (validSum == 0) || (nValid < 10)
    !uninformative})
House2001m <- House2001m[informative, ]
## Make a vector of colours, blue for Republican and red for Democrat:
parties <- House2001$party[informative]
partyColors <- rep("black", length(parties))
partyColors <- ifelse(parties == "D", "red", partyColors)
partyColors <- ifelse(parties == "R", "blue", partyColors)
## Expand the data for statistical modelling:
House2001v <- as.vector(House2001m)
House2001f <- data.frame(member = rownames(House2001m),
                         party = parties,
                         rollCall = factor(rep((1:20),
                                           rep(nrow(House2001m), 20))),
                         vote = House2001v)
## Now fit an "empty" model, in which all members vote identically:
baseModel <- glm(vote ~ -1 + rollCall, family = binomial, data = House2001f)
## From this, get starting values for a one-dimensional multiplicative term:
Start <- residSVD(baseModel, rollCall, member)
House2001model1 <- gnm(0.5 + 0.94*(vote - 0.5) ~ Mult(rollCall - 1, member - 1),
              eliminate = ~ rollCall,
              family = binomial, data = House2001f,
              na.action = na.exclude, trace = TRUE,
              start = -Start)
## Deviance is 2234.847, df = 5574
##
## Plot the members' positions as estimated in the model:
##
plot(coef(House2001model1)[21:321], col = partyColors,
     xlab = "Alphabetical index (Abercrombie 1 to Young 301)",
     ylab = "Member's relative position, one-dimensional model")
## Can do the same thing with two dimensions, but gnm takes around 50
## slow iterations to converge (there are more than 600 parameters): 
\dontrun{
Start <- residSVD(baseModel, rollCall, member)
House2001model2 <- gnm(0.5 + 0.94*(vote - 0.5) ~
                          Mult(rollCall - 1, member - 1, multiplicity = 2),
              eliminate = ~ rollCall,
              family = binomial, data = House2001f,
              na.action = na.exclude, trace = TRUE,
              start = Start2)
## Deviance is 1545.166, df = 5257
##
## The second dimension is mainly due to rollCall 12, which does not
## correlate well with the rest.  (Look at the coefficients of House2001model1)
}
}
\keyword{datasets}
