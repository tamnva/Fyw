#' Find the youngwater fraction (Fyw) and other associtate parameters
#'
#' @param AP A numeric vectors (or a scalar), amplitude of sine wave of
#' @details Fitting non linear sine function of the following form:
#' @return A list of object
#' @seealso
#' findAlphaBeta
#' @examples
#' Fyw(AP = runif(3), phiP = runif(3), AS = runif(4), phiS = runif(4))
#' @export

Fyw <- function(AP = NULL, phiP = NULL, AS = NULL, phiS = NULL){

  # Check null
  if (is.null(AP)) stop("AP cannot be null")
  if (is.null(phiP)) stop("phiP cannot be null")
  if (is.null(AS)) stop("AS cannot be null")
  if (is.null(phiS)) stop("phiS cannot be null")

  # Check missing
  if(all(is.na(AP))) stop("AP cannot contain missing values")
  if(all(is.na(phiP))) stop("phiP cannot contain missing values")
  if(all(is.na(AS))) stop("AS cannot contain missing values")
  if(all(is.na(phiS))) stop("phiS cannot contain missing values")

  # Check length
  if(length(AP)/length(phiP) != 1) stop("AP and phiP must have the same length")
  if(length(AS)/length(phiS) != 1) stop("AS and phiS must have the same length")

  # Combine different AP, phiP, AS, and phiS to get Fyw
  output <- data.frame(AP = NA, phiP = NA,
                       AS = NA, phiS = NA,
                       Fyw = NA, alpha = NA,
                       beta = NA, meanTT_year = NA,
                       P = NA, S = NA)

  for (i in 1:length(AS)){
    for (j in 1:length(AP)){
      # Check if phiS > phiP and AS < AP
      if((phiS[i] > phiP[j]) & (AS[i] < AP[j])){
        # Now find Fyw and other parameters
        Fyw = AS[i]/AP[j]
        temp <- findAlphaBeta(phiS = phiS[i], phiP = phiP[j], Fyw = Fyw)
        output <- rbind(output, c(round(AP[j], 3),round(phiP[j], 3),
                                  round(AS[i], 3),round(phiS[i], 3),
                                  round(Fyw, 3), round(temp$alpha, 3),
                                  round(temp$beta, 3), round(temp$meanTT, 3),
                                  P = j, S = i))
      } else {
        print("Warning: The following input combination was not used for Fyw estimation")
        print(paste("phiS = ", round(phiS[i],3), ", phiP = ", round(phiP[i],3),
                    ", AS = ", round(AS[i],3), ", AP = ", round(AP[i],3)))
      }
    }
  }

  # Print output message
  print("Successful done")

  # Return output
  return(output[-c(1),])
}
