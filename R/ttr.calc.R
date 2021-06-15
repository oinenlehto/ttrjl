#' Converts INR Values to TTR
#'
#' This function converts INR (International Normalized Ratio) values to TTR (Time in Therapeutic Range) for the desired time and interval. The TTR values are calculated via the Rosendaal method (1). Any time periods exceeding the specified TTR interval from the previous measured INR value return '999'. If the patient is not included in the INR data a value '888' will return.
#' @usage ttr.calc(day, ttr_interval, inr_min, inr_max,
#' patient_id, data, inr_variable, time_variable, id_variable,
#' print_number_of_inr_values = FALSE, print_matrix1 = FALSE,
#' print_matrix2 = FALSE)
#' @param day The time of the desired TTR interval. TTR will be calculated for the period before the specified day.
#' @param ttr_interval The lenght of the desired TTR interval.
#' @param inr_min The lower limit of the INR target.
#' @param inr_max The upper limit of the INR target.
#' @param patient_id The individual patient number matching with the id_variable.
#' @param data The data frame including the INR values, patient numbers and time variable.
#' @param inr_variable The variable name in the data frame including the INR values.
#' @param time_variable The variable name including the time of INR measurements.
#' @param id_variable The variable name including the individual patient numbers.
#' @param print_number_of_inr_values If TRUE, will print a matrix including the TTR value and the number of INR measurements during the examined interval. The default is FALSE.
#' @param print_matrix1 If TRUE, will print a matrix including all the INR values for the specified patients in addition to time periods, days in therapeutic range for each time interval, etc. The default is FALSE.
#' @param print_matrix2 If TRUE, will print the aforementioned matrix for the desired time interval in which the first and last time in therapeutic range values are recalculated for the TTR calculation. The default is FALSE.
#' @import dplyr
#' @references Rosendaal F, Cannegieter S, van der Meer F, et al. A method to determine the optimal intensity of oral anticoagulant therapy. Thromb Haemost 1993; 69(3): 236–239.
#' @examples
#' ##Not run:
#' patient_number <- ('BIO1005')
#' time <- c(3, 12, 70)
#' INR <- c(2.2, 2.6, 3.9)
#' INR_data <- data.frame(patient_number, time, INR)
#' ttr.calc(day = 65, ttr_interval = 60, inr_min = 2, inr_max = 3,
#' patient_id = 'BIO1005', data = INR_data, inr_variable = 'INR',
#' time_variable = 'time', id_variable = 'patient_number',
#' print_number_of_inr_values = TRUE)
#' ## End(**Not run**)
#' @export
ttr.calc = function (day, ttr_interval, inr_min, inr_max, patient_id, data, inr_variable, time_variable, id_variable, print_number_of_inr_values = FALSE, print_matrix1 = FALSE, print_matrix2 = FALSE) { 
  
# Checking the missing arguments
  if(missing(day) | missing(ttr_interval) | missing(inr_min) | missing(inr_max) | missing(patient_id) | missing(data) | missing(inr_variable) | missing(time_variable) | missing(id_variable)) {
    return(print('argument missing'))
  } else {
  
# Creating the matrix 'S' for TTR calculations
S <- subset(data, data[, id_variable] == patient_id)
if((!nrow(S) == 0)) {S <- arrange(S, S[, time_variable])
S <- S[c(inr_variable, time_variable)]

# Creating variable 'rank' containing running number for the INR samples
S$rank <- 1:nrow(S)

# Creating variable 'D' containing the date difference between INR samples
if(!(nrow(S) == 1)) { x <- c(2:nrow(S))
for(i in x) {
S[i, 'D'] <- S[i, time_variable] - S[i-1, time_variable]
}
} else {S[1, 'D'] <- NA}

# Creating variable 'A' containing the higher INR value
if(!(nrow(S) == 1)) { for(i in x) {
  if(S[i, inr_variable] >= S[i-1, inr_variable]) {
    S[i, 'A'] <- S[i, inr_variable]
  } else {S[i, 'A'] <- S[i-1, inr_variable]}
}
}else {S[1, 'A'] <- NA}

# Creating variable 'B' containing lower INR value
  if(!(nrow(S) == 1)) { for(i in x) {
  if(S[i, inr_variable] <= S[i-1, inr_variable]) {
    S[i, 'B'] <- S[i, inr_variable]
  } else {S[i, 'B'] <- S[i-1, inr_variable]}
  }
  } else {S[1, 'B'] <- NA}

# Creating variable 'change_per_day' containing the INR change per day during the interval
S$change_per_day <- (S$A - S$B)/S$D

# Creating variable 'days_out_of_range_on_low_end' containing the number of days below the INR target during the interval
S$days_out_of_range_on_low_end <- (inr_min-S$B)/S$change_per_day
if(!(nrow(S) == 1)) { for(i in x) {
  if(S[i, 'change_per_day'] == 0 & S[i, inr_variable] < inr_min) {
    S[i, 'days_out_of_range_on_low_end'] <- S[i, 'D']
  }
  if(S[i, 'change_per_day'] == 0 & S[i, inr_variable] >= inr_min) {
    S[i, 'days_out_of_range_on_low_end'] <- 0
  }}
} else {S[1, 'days_out_of_range_on_low_end'] <- NA}

# Creating variable 'days_out_of_range_on_high_end' containing the number of days above the INR target during the interval
S$days_out_of_range_on_high_end <- (S$A-inr_max)/S$change_per_day
if(!(nrow(S) == 1)) { for(i in x) {
  if(S[i, 'change_per_day'] == 0 & S[i, inr_variable] > inr_max) {
    S[i, 'days_out_of_range_on_high_end'] <- S[i, 'D']
  }
  if(S[i, 'change_per_day'] == 0 & S[i, inr_variable] <= inr_max) {
      S[i, 'days_out_of_range_on_high_end'] <- 0
  }}
} else {S[1, 'days_out_of_range_on_high_end'] <- NA}

# Creating variable 'days_in_therapeutic_range' containing the number of days in INR target during the interval
S$days_in_therapeutic_range[S$rank == '1'] <- 0
if(nrow(S) == 1) {
    S[1, 'days_in_therapeutic_range'] <- 0
  } else {
for(i in x) {
  if(S[i, 'days_out_of_range_on_low_end'] < 0 & S[i, 'days_out_of_range_on_high_end'] < 0) {
    S[i, 'days_in_therapeutic_range'] <- S[i, 'D']
  } else if(S[i, 'days_out_of_range_on_low_end'] < 0) {
    S[i, 'days_in_therapeutic_range'] <- S[i, 'D'] - S[i, 'days_out_of_range_on_high_end']
  } else if(S[i, 'days_out_of_range_on_high_end'] < 0) {
    S[i, 'days_in_therapeutic_range'] <- S[i, 'D'] - S[i, 'days_out_of_range_on_low_end']
  } else {
    S[i, 'days_in_therapeutic_range'] <- S[i, 'D'] - S[i, 'days_out_of_range_on_low_end'] - S[i, 'days_out_of_range_on_high_end']
  }
  if(S[i, 'days_in_therapeutic_range'] < 0) {
    S[i, 'days_in_therapeutic_range'] <- 0
  }
  if((S[i, time_variable] - S[i-1, time_variable]) > ttr_interval) {
    S[i, 'days_in_therapeutic_range'] <- NA
  }
}}

# Creating the matrix Sd1 containing the INR values during the e.g. 60 d interval before the time 'd' and Sd2 containing also values before and after the interval
Sd1 <- subset(S, S[, time_variable] <= day & S[, time_variable] >= day - ttr_interval)
Sd2 <- subset(S, S$rank <= Sd1[nrow(Sd1), 'rank'] + 1 & S$rank >= Sd1[1, 'rank'] - 1)

# Calculating the 60 d TTR for a specified day
if(day >= S[nrow(S), time_variable] + ttr_interval) {
  TTR <- 999
} else if(day < S[1, time_variable] + ttr_interval) {
  TTR <- 999
} else if(any(is.na(Sd1[, 'days_in_therapeutic_range']))) {
  TTR <- 999
} else if(nrow(Sd1) == 0) {
  TTR <- 999
  # If day - ttr_interval == Sd1[1, time_variable] & day == Sd1[nrow(Sd1), time_variable]
} else if(Sd1[nrow(Sd1), time_variable] == day & Sd1[1, time_variable] == day - ttr_interval) {
  Sd1[1, 'days_in_therapeutic_range'] <- 0
  TTR <- sum(Sd1[, 'days_in_therapeutic_range']) / ttr_interval
  # If day - ttr_interval == Sd1[1, time_variable] & day != Sd1[nrow(Sd1), time_variable]
} else if(Sd1[1, time_variable] == day - ttr_interval & Sd1[nrow(Sd1), time_variable] != day) {
  if(Sd1[nrow(Sd1), 'rank'] == Sd2[nrow(Sd2), 'rank']) {
    if(Sd2[nrow(Sd2), inr_variable] >= inr_min & Sd2[nrow(Sd2), inr_variable] <= inr_max) {
      Sd2[nrow(Sd2) + 1, 'days_in_therapeutic_range'] <- day - Sd2[nrow(Sd2), time_variable]
    } else {Sd2[nrow(Sd2) + 1, 'days_in_therapeutic_range'] <- 0}
  } else if(is.na(Sd2[nrow(Sd2), 'days_in_therapeutic_range'])) {
    if(Sd2[nrow(Sd2) - 1, inr_variable] >= inr_min & Sd2[nrow(Sd2) - 1, inr_variable] <= inr_max) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- day - Sd2[nrow(Sd2) - 1, time_variable]
    } else {Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- 0}
  } else if(Sd2[nrow(Sd2), 'days_in_therapeutic_range'] == Sd2[nrow(Sd2), 'D']) {
    Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- Sd2[nrow(Sd2), 'D'] - (Sd2[nrow(Sd2), time_variable] - day)
  } else if(Sd2[nrow(Sd2) - 1, inr_variable] >= inr_min & Sd2[nrow(Sd2) - 1, inr_variable] <= inr_max) {
    if(Sd2[nrow(Sd2), 'days_in_therapeutic_range'] > day - Sd2[nrow(Sd2) - 1, time_variable]) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- day -  Sd2[nrow(Sd2) - 1, time_variable]
    }
  } else if((Sd2[nrow(Sd2) - 1, inr_variable] < inr_min | Sd2[nrow(Sd2) - 1, inr_variable] > inr_max) & Sd2[nrow(Sd2), inr_variable] >= inr_min & Sd2[nrow(Sd2), inr_variable] <= inr_max) {
    if(Sd2[nrow(Sd2), 'days_in_therapeutic_range'] < Sd2[nrow(Sd2), time_variable] - day) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- 0
    } else {Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- Sd2[nrow(Sd2), 'days_in_therapeutic_range'] - Sd2[nrow(Sd2), time_variable] + day}
  } else if(Sd2[nrow(Sd2) - 1, inr_variable] > inr_max & Sd2[nrow(Sd2), inr_variable] < inr_min) {
    if(Sd2[nrow(Sd2), 'days_out_of_range_on_high_end'] >= day - Sd2[nrow(Sd2) - 1, time_variable]) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- 0
    } else if(Sd2[nrow(Sd2), time_variable] - Sd2[nrow(Sd2), 'days_out_of_range_on_low_end'] > day) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- day - Sd2[nrow(Sd2), 'days_out_of_range_on_high_end'] - Sd2[nrow(Sd2) - 1, time_variable]
    }
  } else if(Sd2[nrow(Sd2) - 1, inr_variable] < inr_min & Sd2[nrow(Sd2), inr_variable] > inr_max) {
    if(Sd2[nrow(Sd2), 'days_out_of_range_on_low_end'] >= day - Sd2[nrow(Sd2) - 1, time_variable]) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- 0
    } else if(Sd2[nrow(Sd2), time_variable] - Sd2[nrow(Sd2), 'days_out_of_range_on_high_end'] > day) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- day - Sd2[nrow(Sd2), 'days_out_of_range_on_low_end'] - Sd2[nrow(Sd2) - 1, time_variable]
    }}
  if(Sd1[1, 'rank'] == Sd2[1, 'rank']) {
    Sd2[1, 'days_in_therapeutic_range'] <- 0
    TTR <- sum(Sd2[1:nrow(Sd2), 'days_in_therapeutic_range']) / ttr_interval
  } else {
    Sd2[2, 'days_in_therapeutic_range'] <- 0
    TTR <- sum(Sd2[2:nrow(Sd2), 'days_in_therapeutic_range']) / ttr_interval}
  # If day - ttr_interval != Sd1[1, time_variable] & day == Sd1[nrow(Sd1), time_variable]
} else if(Sd1[1, time_variable] != day - ttr_interval & Sd1[nrow(Sd1), time_variable] == day) {
  if(Sd1[1, 'rank'] == Sd2[1, 'rank']) {
      TTR <- 999
  } else if(Sd2[2, 'days_in_therapeutic_range'] == Sd2[2, 'D']) {
    Sd2[2, 'days_in_therapeutic_range'] <- Sd2[2, time_variable] - day + ttr_interval
  } else if(Sd2[2, inr_variable] >= inr_min & Sd2[2, inr_variable] <= inr_max) {
    if(Sd2[2, 'days_in_therapeutic_range'] >= Sd2[2, time_variable] - day + ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- Sd2[2, time_variable] - day + ttr_interval
    }
  } else if((Sd2[2, inr_variable] < inr_min | Sd2[2, inr_variable] > inr_max) & Sd2[1, inr_variable] >= inr_min & Sd2[1, inr_variable] <= inr_max) {
    if(Sd2[2, 'days_in_therapeutic_range'] < day - ttr_interval - Sd2[1, time_variable]) {
      Sd2[2, 'days_in_therapeutic_range'] <- 0
    } else {Sd2[2, 'days_in_therapeutic_range'] <- Sd2[1, time_variable] + Sd2[2, 'days_in_therapeutic_range'] - day + ttr_interval}
  } else if(Sd2[1, inr_variable] > inr_max & Sd2[2, inr_variable] < inr_min) {
    if(Sd2[2, 'days_out_of_range_on_low_end'] >= Sd2[2, time_variable] - day + ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- 0
    } else if(Sd2[2, time_variable] - Sd2[2, 'days_in_therapeutic_range'] - Sd2[2, 'days_out_of_range_on_low_end'] < day - ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- Sd2[2, time_variable] - Sd2[2, 'days_out_of_range_on_low_end'] - day + ttr_interval
    }
  } else if(Sd2[1, inr_variable] < inr_min & Sd2[2, inr_variable] > inr_max) {
    if(Sd2[2, 'days_out_of_range_on_high_end'] >= Sd2[2, time_variable] - day + ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- 0
    } else if(Sd2[2, time_variable] - Sd2[2, 'days_in_therapeutic_range'] - Sd2[2, 'days_out_of_range_on_high_end'] < day - ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- Sd2[2, time_variable] - Sd2[2, 'days_out_of_range_on_high_end'] - day + ttr_interval
    }}
  if(Sd1[1, 'rank'] != Sd2[1, 'rank'] & Sd1[nrow(Sd1), 'rank'] == Sd2[nrow(Sd2), 'rank']) {
    TTR <- sum(Sd2[2:nrow(Sd2), 'days_in_therapeutic_range']) / ttr_interval
  } else if(Sd1[1, 'rank'] != Sd2[1, 'rank'] & Sd1[nrow(Sd1), 'rank'] != Sd2[nrow(Sd2), 'rank']){
    TTR <- sum(Sd2[2:(nrow(Sd2) - 1), 'days_in_therapeutic_range']) / ttr_interval
  }
  # If day - ttr_interval != Sd1[1, time_variable] & day != Sd1[nrow(Sd1), time_variable]
} else if(Sd1[1, time_variable] != day - ttr_interval & Sd1[nrow(Sd1), time_variable] != day) {
  if(Sd1[1, 'rank'] == Sd2[1, 'rank']) {
    TTR <- 999
  } else if(Sd2[2, 'days_in_therapeutic_range'] == Sd2[2, 'D']) {
    Sd2[2, 'days_in_therapeutic_range'] <- Sd2[2, time_variable] - day + ttr_interval
  } else if(Sd2[2, inr_variable] >= inr_min & Sd2[2, inr_variable] <= inr_max) {
    if(Sd2[2, 'days_in_therapeutic_range'] >= Sd2[2, time_variable] - day + ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- Sd2[2, time_variable] - day + ttr_interval
    }
  } else if((Sd2[2, inr_variable] < inr_min | Sd2[2, inr_variable] > inr_max) & Sd2[1, inr_variable] >= inr_min & Sd2[1, inr_variable] <= inr_max) {
    if(Sd2[2, 'days_in_therapeutic_range'] < day - ttr_interval - Sd2[1, time_variable]) {
      Sd2[2, 'days_in_therapeutic_range'] <- 0
    } else {Sd2[2, 'days_in_therapeutic_range'] <- Sd2[1, time_variable] + Sd2[2, 'days_in_therapeutic_range'] - day + ttr_interval}
  } else if(Sd2[1, inr_variable] > inr_max & Sd2[2, inr_variable] < inr_min) {
    if(Sd2[2, 'days_out_of_range_on_low_end'] >= Sd2[2, time_variable] - day + ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- 0
    } else if(Sd2[2, time_variable] - Sd2[2, 'days_in_therapeutic_range'] - Sd2[2, 'days_out_of_range_on_low_end'] < day - ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- Sd2[2, time_variable] - Sd2[2, 'days_out_of_range_on_low_end'] - day + ttr_interval
    }
  } else if(Sd2[1, inr_variable] < inr_min & Sd2[2, inr_variable] > inr_max) {
    if(Sd2[2, 'days_out_of_range_on_high_end'] >= Sd2[2, time_variable] - day + ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- 0
    } else if(Sd2[2, time_variable] - Sd2[2, 'days_in_therapeutic_range'] - Sd2[2, 'days_out_of_range_on_high_end'] < day - ttr_interval) {
      Sd2[2, 'days_in_therapeutic_range'] <- Sd2[2, time_variable] - Sd2[2, 'days_out_of_range_on_high_end'] - day + ttr_interval
    }}
  if(Sd1[nrow(Sd1), 'rank'] == Sd2[nrow(Sd2), 'rank']) {
    if(Sd2[nrow(Sd2), inr_variable] >= inr_min & Sd2[nrow(Sd2), inr_variable] <= inr_max) {
      Sd2[nrow(Sd2) + 1, 'days_in_therapeutic_range'] <- day - Sd2[nrow(Sd2), time_variable]
    } else {Sd2[nrow(Sd2) + 1, 'days_in_therapeutic_range'] <- 0}
  } else if(is.na(Sd2[nrow(Sd2), 'days_in_therapeutic_range'])) {
    if(Sd2[nrow(Sd2) - 1, inr_variable] >= inr_min & Sd2[nrow(Sd2) - 1, inr_variable] <= inr_max) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- day - Sd2[nrow(Sd2) - 1, time_variable]
    } else {Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- 0}
  } else if(Sd2[nrow(Sd2), 'days_in_therapeutic_range'] == Sd2[nrow(Sd2), 'D']) {
    Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- Sd2[nrow(Sd2), 'D'] - (Sd2[nrow(Sd2), time_variable] - day)
  } else if(Sd2[nrow(Sd2) - 1, inr_variable] >= inr_min & Sd2[nrow(Sd2) - 1, inr_variable] <= inr_max) {
    if(Sd2[nrow(Sd2), 'days_in_therapeutic_range'] > day - Sd2[nrow(Sd2) - 1, time_variable]) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- day -  Sd2[nrow(Sd2) - 1, time_variable]
    }
  } else if((Sd2[nrow(Sd2) - 1, inr_variable] < inr_min | Sd2[nrow(Sd2) - 1, inr_variable] > inr_max) & Sd2[nrow(Sd2), inr_variable] >= inr_min & Sd2[nrow(Sd2), inr_variable] <= inr_max) {
    if(Sd2[nrow(Sd2), 'days_in_therapeutic_range'] < Sd2[nrow(Sd2), time_variable] - day) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- 0
    } else {Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- Sd2[nrow(Sd2), 'days_in_therapeutic_range'] - Sd2[nrow(Sd2), time_variable] + day}
  } else if(Sd2[nrow(Sd2) - 1, inr_variable] > inr_max & Sd2[nrow(Sd2), inr_variable] < inr_min) {
    if(Sd2[nrow(Sd2), 'days_out_of_range_on_high_end'] >= day - Sd2[nrow(Sd2) - 1, time_variable]) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- 0
    } else if(Sd2[nrow(Sd2), time_variable] - Sd2[nrow(Sd2), 'days_out_of_range_on_low_end'] > day) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- day - Sd2[nrow(Sd2), 'days_out_of_range_on_high_end'] - Sd2[nrow(Sd2) - 1, time_variable]
    }
  } else if(Sd2[nrow(Sd2) - 1, inr_variable] < inr_min & Sd2[nrow(Sd2), inr_variable] > inr_max) {
    if(Sd2[nrow(Sd2), 'days_out_of_range_on_low_end'] >= day - Sd2[nrow(Sd2) - 1, time_variable]) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- 0
    } else if(Sd2[nrow(Sd2), time_variable] - Sd2[nrow(Sd2), 'days_out_of_range_on_high_end'] > day) {
      Sd2[nrow(Sd2), 'days_in_therapeutic_range'] <- day - Sd2[nrow(Sd2), 'days_out_of_range_on_low_end'] - Sd2[nrow(Sd2) - 1, time_variable]
    }}
  if(Sd1[1, 'rank'] != Sd2[1, 'rank']) {
    TTR <- sum(Sd2[2:nrow(Sd2), 'days_in_therapeutic_range']) / ttr_interval
  }}
if(print_matrix1 == TRUE) {
  return(S)
}
if(print_matrix2 == TRUE) {
  return(Sd2)
}
if(print_number_of_inr_values == TRUE) {
number_of_INR_values <- nrow(Sd1)
return(data.frame(TTR, number_of_INR_values))
}
if(print_number_of_inr_values == FALSE) {
  return(TTR)
}} else {
    TTR <- 888
    if(print_matrix1 == TRUE | print_matrix2 == TRUE) {
      return(S)
    }
    if(print_number_of_inr_values == TRUE) {
      number_of_INR_values <- 888
      return(data.frame(TTR, number_of_INR_values))
    }
    if(print_number_of_inr_values == FALSE){
      return(TTR)
}}}}
