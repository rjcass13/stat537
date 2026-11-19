orpol <- function(v,d){
	ll <- v - mean(v)
	x <- cbind(1,ll)
	for (i in 2:d){
		y <- v^i
		ll <- y - x%*%solve(t(x)%*%x)%*%t(x)%*%y
		x <- cbind(x,ll)
		}
		for (i in 1:(d+1)){
                        x[,i] <- x[,i]/sqrt(sum(x[,i]^2))
			}
		return(x)
		}
