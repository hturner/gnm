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
\dontrun{
## This example takes some while to run 
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
##
## Now fit the logistic model with one multiplicative term.
## For the response variable, instead of vote=0,1 we use 0.03 and 0.97,
## corresponding approximately to a bias-reducing adjustment of p/(2n),
## where p is the number of parameters and n the number of observations.
##
voteAdj <- 0.5 + 0.94*(House2001f$vote - 0.5)
House2001model1 <- gnm(voteAdj ~ Mult(rollCall - 1, member - 1),
              eliminate = rollCall,
              family = binomial, data = House2001f,
              na.action = na.exclude, trace = TRUE, tolerance = 1e-03,
              start = -Start)
## Deviance is 2234.847, df = 5574
##
## Plot the members' positions as estimated in the model:
##
memberParameters <- pickCoef(House2001model1, "member")
plot(coef(House2001model1)[memberParameters], col = partyColors,
     xlab = "Alphabetical index (Abercrombie 1 to Young 301)",
     ylab = "Member's relative position, one-dimensional model")
## Can do the same thing with two dimensions, but gnm takes around 40
## slow iterations to converge (there are more than 600 parameters): 
Start2 <- residSVD(baseModel, rollCall, member, d = 2)
House2001model2 <- gnm(
              voteAdj ~ Mult(rollCall - 1, member - 1, multiplicity = 2),
              eliminate = rollCall,
              family = binomial, data = House2001f,
              na.action = na.exclude, trace = TRUE, tolerance = 1e-03,
              start = Start2, lsMethod = "qr")
## Deviance is 1545.166, df = 5257
##
memberParameters1 <- pickCoef(House2001model, "Mult1.*member")
memberParameters2 <- pickCoef(House2001model, "Mult2.*member")
plot(coef(House2001model2)[memberParameters1],
     coef(House2001model2)[memberParameters2],
     col = partyColors,
     xlab = "Dimension 1",
     ylab = "Dimension 2",
     main = "House2001 data: Member positions, 2-dimensional model")
##
## The second dimension is mainly due to rollCall 12, which does not
## correlate well with the rest -- look at the coefficients of
## House2001model1, or at
cormat <- cor(na.omit(House2001m))
}
}
\keyword{datasets}
