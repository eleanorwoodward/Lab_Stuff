## Takes PoN filtered data, and returns a binary yes or no based on cutoff
run.pon <- function(maf, max){
  theDecider <- function(val, max){
    if (is.nan(val)){
      FALSE
    }else if(val > max){
      TRUE
    }else{
      FALSE
    }
  }
  
  germline <- sapply(maf$pon_loglike, theDecider, max)
  maf$pon_germline <- germline
  maf
}


