
solve_branch <- function(branch, coalescent_times, start, end) {
	# If there's at most one lineage in the branch during the time period, nothing needs to be done.
	if (length(branch$inhabitants) > 1) {
    # If there's at least one coalescent event during the time frame, we record it and call the function recursively.
	  minimum_time <- rexp(1, rate = choose(length(branch$inhabitants), 2)*branch$lambda)
	  if (minimum_time < end - start) {
	  	pair <- sample(branch$inhabitants, 2)
	  	new <- paste(pair[1], pair[2], sep = "")
	  	branch$inhabitants <- setdiff(c(branch$inhabitants, new), pair)
	  	coalescent_times[[new]] <- minimum_time + start
	  	inner <- solve_branch(branch, coalescent_times, minimum_time + start, end)
	  	branch <- inner$branch
	  	coalescent_times <- inner$coalescent_times
	  }
	}
	return(list(branch = branch, coalescent_times = coalescent_times))
}

solve_tree <- function(branches) {
	coalescent_times <- list()
	# Compile a vector of pivotal moments (individual sampling times, population divergence and admixture times).
	pivotal <- numeric(0)
	for (branch in branches) {
    for (sample in branch$samples) {pivotal <- c(pivotal, sample)}
		for (switch in branch$switches) {pivotal <- c(pivotal, switch[1])}
	}
	pivotal <- unique(pivotal)
	pivotal <- pivotal[order(pivotal)]
	# Add infinity as the last pivotal moment.
	pivotal[length(pivotal) + 1] <- Inf
	# For all pivotal moments except the infinity:
	for (t in seq(1, length(pivotal) - 1)) {
		now <- pivotal[t]
	  soon <- pivotal[t + 1]
	  # First perform the special actions of the pivotal moment (samplings, switches).
	  for (branch in names(branches)) {
      for (sample in names(branches[[branch]]$samples)) {
      	if (branches[[branch]]$samples[[sample]] == now) {
      		branches[[branch]]$inhabitants <- c(branches[[branch]]$inhabitants, sample)
      		coalescent_times[[sample]] <- now
      	}
      }
	  	for (switch in names(branches[[branch]]$switches)) {
	  		if (branches[[branch]]$switches[[switch]][1] == now) {
	  		  # Moving inhabitants from a branch to another, possibly.
	  			for (inhabitant in branches[[branch]]$inhabitants) {
	  			  if (runif(1) < branches[[branch]]$switches[[switch]][2]) {
	  			  	branches[[branch]]$inhabitants <- setdiff(branches[[branch]]$inhabitants, inhabitant)
	  			  	branches[[switch]]$inhabitants <- c(branches[[switch]]$inhabitants, inhabitant)
	  			  }
	  			}
	  		}
	  	}
	  }
	  # Then solve each branch until the next pivotal moment.
	  for (branch in names(branches)) {
	  	solved <- solve_branch(branches[[branch]], coalescent_times, now, soon)
	  	branches[[branch]] <- solved$branch
	  	coalescent_times <- solved$coalescent_times
	  }
	}
	return(coalescent_times)
}

subsetstring <- function(x, y) {
	return(all(unlist(strsplit(x, "")) %in% unlist(strsplit(y, ""))))
}

pattern <- function(tree) {
  # Adding one point mutation somewhere in the gene tree and reporting the ABBABABA-pattern.
	# Ignoring the fact that some trees have higher chance of point mutation than others.
	# Ignoring recurrent and convergent mutations as well.
	# Pease and Hahn cite Durand (2011) saying it's no big deal.
	L <- rep(0, length(names(tree)))
	for (n in seq(1, length(names(tree)))) {
	  node <- names(tree)[n]
		candidates <- names(tree)[nchar(names(tree)) > nchar(node)]
		if (length(candidates > 0)) {
      candidates <- candidates[sapply(candidates, subsetstring, x = node)]
			if (length(candidates) > 0) {
			  parent <- candidates[which.min(nchar(candidates))]
				L[n] <- tree[[parent]] - tree[[node]]
			}
	  }
	}
	random <- runif(1)*sum(L)
	lsum <- 0
	for (l in seq(1, length(L))) {
		lsum <- lsum + L[l]
		if (lsum > random) {break()}
	} # names(tree)[l] is the first node affected by the mutation.
	leaves <- max(nchar(names(tree)))/2
	pattern <- rep("A", leaves)
	for (j in seq(1, leaves)) {
		if (grepl(paste(j), names(tree)[l])) {pattern[j] <- "B"}
	}
	return(paste(pattern, collapse = ""))
}

n <- function(X, x) {
	if (x %in% names(X)) {return(X[x])} else {return(0)}
}

DPA <- function(X, strict = FALSE) {
	L <- n(X, "BABA")
	R <- n(X, "ABBA")
	if (strict == FALSE) {
		L <- L + n(X, "ABAB")
		R <- R + n(X, "BAAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- as.numeric(1 - pchisq(abs((L - R)**2/(L + R)), 1))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > 0.01) {
			sign <- "0"
		} else {
			if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}
		}
	} else {sign <- "?"}
  return(list(value, pvalue, sign))
}

DFO <- function(X, strict = FALSE) {
	L <- n(X, "BABAA") + n(X, "BBBAA") + n(X, "ABABA") + n(X, "AAABA")
	R <- n(X, "BAABA") + n(X, "BBABA") + n(X, "ABBAA") + n(X, "AABAA")
	if (strict == FALSE) {
		L <- L + n(X, "ABABB") + n(X, "AAABB") + n(X, "BABAB") + n(X, "BBBAB")
		R <- R + n(X, "ABBAB") + n(X, "AABAB") + n(X, "BAABB") + n(X, "BBABB")
	}
	value <- as.numeric((L - R)/(L + R))
	print(1/sqrt(L + R))
	pvalue <- as.numeric(1 - pchisq(abs((L - R)**2/(L + R)), 1))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > 0.01) {
			sign <- "0"
		} else {
			if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}
		}
	} else {sign <- "?"}
  return(list(value, pvalue, sign))
}

DIL <- function(X, strict = FALSE) {
	L <- n(X, "ABBAA") + n(X, "BBBAA") + n(X, "BAABA") + n(X, "AAABA")
	R <- n(X, "ABABA") + n(X, "BBABA") + n(X, "BABAA") + n(X, "AABAA")
	if (strict == FALSE) {
		L <- L + n(X, "BAABB") + n(X, "AAABB") + n(X, "ABBAB") + n(X, "BBBAB")
		R <- R + n(X, "BABAB") + n(X, "AABAB") + n(X, "ABABB") + n(X, "BBABB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- as.numeric(1 - pchisq(abs((L - R)**2/(L + R)), 1))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > 0.01) {
			sign <- "0"
		} else {
			if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}
		}
	} else {sign <- "?"}
  return(list(value, pvalue, sign))
}

DFI <- function(X, strict = FALSE) {
	L <- n(X, "BABAA") + n(X, "BABBA") + n(X, "ABABA") + n(X, "ABAAA")
	R <- n(X, "ABBAA") + n(X, "ABBBA") + n(X, "BAABA") + n(X, "BAAAA")
	if (strict == FALSE) {
		L <- L + n(X, "ABABB") + n(X, "ABAAB") + n(X, "BABAB") + n(X, "BABBB")
		R <- R + n(X, "BAABB") + n(X, "BAAAB") + n(X, "ABBAB") + n(X, "ABBBB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- as.numeric(1 - pchisq(abs((L - R)**2/(L + R)), 1))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > 0.01) {
			sign <- "0"
		} else {
			if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}
		}
	} else {sign <- "?"}
  return(list(value, pvalue, sign))
}

DOL <- function(X, strict = FALSE) {
	L <- n(X, "BAABA") + n(X, "BABBA") + n(X, "ABBAA") + n(X, "ABAAA")
	R <- n(X, "ABABA") + n(X, "ABBBA") + n(X, "BABAA") + n(X, "BAAAA")
	if (strict == FALSE) {
		L <- L + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BAABB") + n(X, "BABBB")
		R <- R + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABABB") + n(X, "ABBBB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- as.numeric(1 - pchisq(abs((L - R)**2/(L + R)), 1))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > 0.01) {
			sign <- "0"
		} else {
			if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}
		}
	} else {sign <- "?"}
  return(list(value, pvalue, sign))
}

DELTA_S1 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BABAA")
	R <- n(X, "BAABA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABABB")
		R <- R + n(X, "ABBAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S2 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "ABBAA")
	R <- n(X, "ABABA")
	if (ancestral == FALSE) {
		L <- L + n(X, "BAABB")
		R <- R + n(X, "BABAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S3 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BABAA")
	R <- n(X, "ABBAA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABABB")
		R <- R + n(X, "BAABB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S4 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BAABA")
	R <- n(X, "ABABA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABBAB")
		R <- R + n(X, "BABAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S5 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BABAA")
	R <- n(X, "ABABA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABABB")
 		R <- R + n(X, "BABAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S6 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BAABA")
	R <- n(X, "ABBAA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABBAB")
		R <- R + n(X, "BAABB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S7 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "ABBBA")
	R <- n(X, "BABBA")
	if (ancestral == FALSE) {
		L <- L + n(X, "BAAAB")
		R <- R + n(X, "ABAAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S8 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BBABA")
	R <- n(X, "BBBAA")
	if (ancestral == FALSE) {
		L <- L + n(X, "AABAB")
    R <- R + n(X, "AAABB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S9 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BAAAA")
	R <- n(X, "ABAAA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABBBB")
    R <- R + n(X, "BABBB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_S10 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "AABAA")
	R <- n(X, "AAABA")
	if (ancestral == FALSE) {
		L <- L + n(X, "BBABB")
    R <- R + n(X, "BBBAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_A1 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BABAA")
	R <- n(X, "ABBAA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABABB")
    R <- R + n(X, "BAABB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_A2 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BAABA")
	R <- n(X, "ABABA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABBAB")
    R <- R + n(X, "BABAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_A3 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "ABBBA")
	R <- n(X, "BABBA")
	if (ancestral == FALSE) {
		L <- L + n(X, "BAAAB")
    R <- R + n(X, "ABAAB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_A4 <- function(X, alpha = 0.01, ancestral = FALSE) {
	L <- n(X, "BAAAA")
	R <- n(X, "ABAAA")
	if (ancestral == FALSE) {
		L <- L + n(X, "ABBBB")
    R <- R + n(X, "BABBB")
	}
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q1 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB")
	R <- n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q2 <- function(X, alpha = 0.01) {
	L <- n(X, "ABABA") + n(X, "BABAB")
	R <- n(X, "ABAAB") + n(X, "BABBA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q3 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB")
	R <- n(X, "ABABA") + n(X, "BABAB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q4 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAB") + n(X, "ABBBA")
	R <- n(X, "ABAAB") + n(X, "BABBA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q5 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB")
	R <- n(X, "BABBA") + n(X, "ABAAB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q6 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAB") + n(X, "ABBBA")
	R <- n(X, "ABABA") + n(X, "BABAB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q7 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB")
	R <- n(X, "ABBAA") + n(X, "BAABB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q8 <- function(X, alpha = 0.01) {
	L <- n(X, "AABBA") + n(X, "BBAAB")
	R <- n(X, "AABAB") + n(X, "BBABA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q9 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DELTA_Q10 <- function(X, alpha = 0.01) {
	L <- n(X, "AAABA") + n(X, "BBBAB")
	R <- n(X, "AAAAB") + n(X, "BBBBA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DS18 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "AAABB") + n(X, "BBBAA")
	R <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "AABAB") + n(X, "BBABA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DS28 <- function(X, alpha = 0.01) {
	L <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "AAABB") + n(X, "BBBAA")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "AABAB") + n(X, "BBABA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DS37 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "ABAAB") + n(X, "BABBA")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DS47 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BABBA")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DS79 <- function(X, alpha = 0.01) {
	L <- n(X, "BAAAB") + n(X, "ABBBA") + n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABAAB") + n(X, "BABBA") + n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DS80 <- function(X, alpha = 0.01) {
	L <- n(X, "AABAB") + n(X, "BBABA") + n(X, "AABAA") + n(X, "BBABB")
	R <- n(X, "AAABB") + n(X, "BBBAA") + n(X, "AAABA") + n(X, "BBBAB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DS1280 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "ABBAA") + n(X, "BAABB") + n(X, "AAABB") + n(X, "BBBAA") + n(X, "AABAA") + n(X, "BBABB")
	R <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "AABAB") + n(X, "BBABA") + n(X, "AAABA") + n(X, "BBBAB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DS3479 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BABBA") + n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABBBA") + n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DA12 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "ABABA") + n(X, "BABAB")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "BAABA") + n(X, "ABBAB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DA13 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "ABAAB") + n(X, "BABBA")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DA23 <- function(X, alpha = 0.01) {
	L <- n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BABBA")
	R <- n(X, "ABABA") + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABBBA")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}

DA1234 <- function(X, alpha = 0.01) {
	L <- n(X, "BABAA") + n(X, "ABABB") + n(X, "BAABA") + n(X, "ABBAB") + n(X, "ABAAB") + n(X, "BABBA") + n(X, "BAAAA") + n(X, "ABBBB")
	R <- n(X, "ABBAA") + n(X, "BAABB") + n(X, "ABABA") + n(X, "BABAB") + n(X, "BAAAB") + n(X, "ABBBA") + n(X, "ABAAA") + n(X, "BABBB")
	value <- as.numeric((L - R)/(L + R))
	pvalue <- 2*as.numeric(1 - pnorm(abs((L - R)/sqrt(L + R))))
	if (is.numeric(pvalue) == TRUE && is.nan(pvalue) == FALSE) {
		if (pvalue > alpha) {
			sign <- "0"
		} else {
	    if (value > 0) {sign <- "+"}
	    if (value < 0) {sign <- "-"}		  	
		}
	} else {sign <- "?"}
	return(list(value, pvalue, sign))
}
