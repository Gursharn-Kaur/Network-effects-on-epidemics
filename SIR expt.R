library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(gridExtra)
library(tergm)
library(deSolve)
library(network)
library(EpiModel)
library(patchwork)
set.seed(123456789)




# 1. Function to create Erdos-Renyi Network
er_graph_unif <- function(p, nodes, unif_matrix){
  adj_matrix <- matrix(0, nrow= nodes, ncol=nodes)
  adj_matrix[which(unif_matrix <= p)] <- 1
  er_graph <- graph_from_adjacency_matrix(adj_matrix, mode="undirected")
  return(er_graph)
}



# 2. Function for SIR model
SIR <- function(states, er_network, infection_rate, recovery_day){
  
  #Vector for storing number of newly infected nodes time to time
  new_inf_vector <- c() 
  
  #Starting day of infection for initially infected node is 0 and for all other nodes is taken to be -1 (can be taken anything else <0)
  start_day <- rep(-1,length(states))
  start_day[which(states == "I")] <- 0
  
  t=1
  
  while (sum(states == "I") != 0){
    
    # infection process
    new_infected <- c() #Vector to store newly infected nodes at time point t
    infected_nodes <- which(states == "I") #Infected nodes before starting infection at time point t
    for (i in infected_nodes){
      neighbors_i <- neighbors(er_network, i) #Neighbor nodes of infected node i
      s_neighbors <- neighbors_i[states[neighbors_i] == "S"] #Susceptible Neighbors of infected node i
      
      for (s_node in s_neighbors){
        if (runif(1) < infection_rate){
          states[s_node] <- "I"  #Susceptible neighbor gets infected using Bernoulli(infection_rate)
          new_infected <- c(new_infected, s_node)
        }
      }
    }
    start_day[new_infected] <- t
    new_inf_vector <- c(new_inf_vector, length(new_infected))
    
    # recovery process
    if (t>=recovery_day){
      recovered <- which(start_day == t-recovery_day) #Stores the nodes who gets recovered at time t
      states[recovered] <- "R"
    }
    
    t <- t+1
  }
  
  total_inf <- sum(states == "R") #Total number of people got infected in the epidemic
  max_inf <- max(new_inf_vector) #Maximum number of people got infected on a single day
  peak_day <- min(which(new_inf_vector == max_inf)) #The day when new infection peaked
  end_day <- t #Duration of the epidemic
  
  #Returning Proportions instead of numbers
  return(c(total_inf/length(states), max_inf/length(states), peak_day, end_day))
}



# 3. Function for each simulation
Evaluation <- function(nodes, probs, beta, recovery_day, initial_infected){
  
  #Uniform numbers generation, which is same for 10 SIRs in a simulation
  unif_matrix <- matrix(0, nodes, nodes)
  unif_matrix[lower.tri(unif_matrix, diag = FALSE)] <- runif((nodes*(nodes-1))/2)
  unif_matrix <- unif_matrix +t(unif_matrix)
  unif_matrix[which(unif_matrix==0)] <- 1
  
  all_sir <- list()
  
  for (l in 1:10){
    
    #initial states of nodes, same for all simulation
    states <- rep("S", nodes)
    states[initial_infected] <- "I"
    
    vec_measures <- matrix(nrow=0, ncol=4)
    
    for (b in beta){
      for (p in probs){
        
        #Erdos Renyi graph from the same uniform matrix 
        er_graph <- er_graph_unif(p, nodes, unif_matrix) 
        
        #The four summary measures for each pair (beta, probs)
        values <- SIR(states, er_graph, b, recovery_day)
        
        #A matrix for summary measures for (beta, probs), starting with fixed beta and varying probs
        vec_measures <- rbind(vec_measures, values)
      }
    }
    
    #The list of matrix of summary measures for each of the 10 SIRs in a single simulation
    all_sir[[l]] <- vec_measures
  }
  
  #Returning the list of summary measures and the unique uniform matrix used in a single simulation
  return(list(list = all_sir, matrix = unif_matrix))
}




# 4. Function for generating summary measures
Simulated_df <- function(sim,nodes, probs, beta, recovery_day, initial_infected){
  
  result <- lapply(1:sim, function(x) Evaluation(nodes, probs, beta, recovery_day, initial_infected))
  
  # Extract all the lists of summary measures
  mat_list <- lapply(result, function(x) x$list)
  
  mat_array <- array(unlist(mat_list), dim = c(length(beta)*length(probs), 4, sim*10))
  
  mean_matrix <- apply(mat_array, c(1, 2), mean)
  sd_matrix   <- apply(mat_array, c(1, 2), sd)
  sd_matrix <- sd_matrix*1.96/sqrt(sim)
  median_matrix <- apply(mat_array, c(1, 2), median, na.rm = TRUE) #median of the summary measures for each pair (beta, probs)
  lower_matrix <- apply(mat_array, c(1, 2), function(x) quantile(x, 0.025, na.rm = TRUE)) # 2.5% quantile
  upper_matrix <- apply(mat_array, c(1, 2), function(x) quantile(x, 0.975, na.rm = TRUE)) # 97.5% quantile
  
  num_r0 <- length(beta)
  split_indices <- split(1:(length(probs) * num_r0), rep(1:num_r0, each=length(probs)))
  mean_list <- lapply(split_indices, function(idx) mean_matrix[idx, ])
  sd_list <- lapply(split_indices, function(idx) sd_matrix[idx, ])
  median_list <- lapply(split_indices, function(idx) median_matrix[idx, ])
  lower_list <- lapply(split_indices, function(idx) lower_matrix[idx, ])
  upper_list <- lapply(split_indices, function(idx) upper_matrix[idx, ])
  
  
  #dataframe for total proportion of infected with quantile
  df_total_inf1 <- data.frame()
  for (i in 1:num_r0){
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Median = median_list[[i]][,1],
      Upper = upper_list[[i]][,1],
      Lower = lower_list[[i]][,1]
    )
    df_total_inf1 <- rbind(df_total_inf1, df_temp)
  }
  
  #dataframe for maximum proportion of infected on a single day with quantile
  df_max_inf1 <- data.frame()
  for (i in 1:num_r0){
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Median = median_list[[i]][,2],
      Upper = upper_list[[i]][,2],
      Lower = lower_list[[i]][,2]
    )
    df_max_inf1 <- rbind(df_max_inf1, df_temp)
  }
  
  #dataframe for the day when number of new infection peaked with quantile
  df_peak_day1 <- data.frame()
  for (i in 1:num_r0){
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Median = median_list[[i]][,3],
      Upper = upper_list[[i]][,3],
      Lower = lower_list[[i]][,3]
    )
    df_peak_day1 <- rbind(df_peak_day1, df_temp)
  }
  
  #dataframe for the duration of epidemic with quantile
  df_end_day1 <- data.frame()
  for (i in 1:num_r0){
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Median = median_list[[i]][,4],
      Upper = upper_list[[i]][,4],
      Lower = lower_list[[i]][,4]
    )
    df_end_day1 <- rbind(df_end_day1, df_temp)
  }
  
  
  #Total infected
  df_total_inf <- data.frame()
  for (i in 1:num_r0){
    upper <- mean_list[[i]][,1] + sd_list[[i]][,1]
    upper[which(upper>1)] <- 1
    lower <- mean_list[[i]][,1] - sd_list[[i]][,1]
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Mean = mean_list[[i]][,1],
      Upper = upper,
      Lower = lower
    )
    df_total_inf <- rbind(df_total_inf, df_temp)
  }
  
  #Max infected
  df_max_inf <- data.frame()
  for (i in 1:num_r0){
    upper <- mean_list[[i]][,2] + sd_list[[i]][,2]
    lower <- mean_list[[i]][,2] - sd_list[[i]][,2]
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Mean = mean_list[[i]][,2],
      Upper = upper,
      Lower = lower
    )
    df_max_inf <- rbind(df_max_inf, df_temp)
  }
  
  #Peak day
  df_peak_day <- data.frame()
  for (i in 1:num_r0){
    upper <- mean_list[[i]][,3] + sd_list[[i]][,3]
    lower <- mean_list[[i]][,3] - sd_list[[i]][,3]
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Mean = mean_list[[i]][,3],
      Upper = upper,
      Lower = lower
    )
    df_peak_day <- rbind(df_peak_day, df_temp)
  }
  
  #End day
  df_end_day <- data.frame()
  for (i in 1:num_r0){
    upper <- mean_list[[i]][,4] + sd_list[[i]][,4]
    lower <- mean_list[[i]][,4] - sd_list[[i]][,4]
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Mean = mean_list[[i]][,4],
      Upper = upper,
      Lower = lower
    )
    df_end_day <- rbind(df_end_day, df_temp)
  }
  
  
  
  # Extract the uniform matrices
  matrices <- lapply(result, function(x) x$matrix)
  M <- list(matrices)
  
  #The list of dataframes for summary measures
  D <- list(
    total_inf_q = df_total_inf1,
    max_inf_q = df_max_inf1,
    peak_q = df_peak_day1,
    end_q = df_end_day1,
    total_inf_sd = df_total_inf,
    max_inf_sd = df_max_inf,
    peak_sd = df_peak_day,
    end_sd = df_end_day,
    list = matrices
  )
  
  return(D)
}


# Main code :
nodes <- 1000
probs <- c(0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.008, 0.01, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
           0.12, 0.14, 0.16, 0.18,0.25, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 0.97)

initial_infected <- sample(1:nodes, 1)

sim <- 15 #No of simulations

R0 <- c(1.2,1.4,1.5)
recovery_day <- 4
beta <- as.vector(R0/recovery_day) # Infection rate

sim_df1 <- Simulated_df(sim, nodes, probs, beta, recovery_day, initial_infected)


#Plot for Total infected
ggplot(sim_df1$total_inf_q, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Total Proportion of Infected Curve with CI (using quantiles) for rec_day = 4",
       x = "Connection Probability",
       y = "Total Proportion of Infected")


pset <- c(0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.008)
df1 <- sim_df1$total_inf_q
df1 <- df1[df1$Probs %in% pset,]

#Zoomed Plot for Total infected
ggplot(df1, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Total Proportion (Zoomed) of Infected Curve with CI (using quantiles) for rec_day = 4",
       x = "Connection Probability",
       y = "Total Proportion of Infected")


#Plot for Max infected
ggplot(sim_df1$max_inf_q, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Maximum Proportion of new infection on a day Curve with CI (using quantiles) for rec_day = 4",
       x = "Connection Probability",
       y = "Maximum Proportion of new infection")

df2 <- sim_df1$max_inf_q
df2$deviation <- df2$Upper - df2$Lower

#Plot for 95% quantile deviation of Max infected
ggplot(df2, aes(x = Probs, y = deviation, color = factor(R0), fill = factor(R0))) +
  geom_line(linewidth = 0.8) +  ylim(0,0.4)+
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "95% quantile deviation of Maximum Proportion of new infection on a day Curve for rec_day =4",
       x = "Connection Probability",
       y = "Quantile deviation")

df2a <- df2[df2$R0==1.2,]

#Plot for 95% quantile deviation of Max infected only for R0 = 1.2
ggplot(df2a, aes(x = Probs, y = deviation, color = factor(R0), fill = factor(R0))) +
  geom_line(linewidth = 0.8) +  ylim(0,0.4)+
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "95% quantile deviation of Maximum Proportion of new infection on a day Curve for rec_day=4",
       x = "Connection Probability",
       y = "Quantile deviation")


#Plot for Peak day
ggplot(sim_df1$peak_q, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Day of maximum infection Curve with CI (using quantiles) for rec_day = 4",
       x = "Connection Probability",
       y = "Day of maximum infection")

df3 <- sim_df1$peak_q
df3$deviation <- df3$Upper - df3$Lower

#Plot for 95% quantile deviation of Peak day
ggplot(df3, aes(x = Probs, y = deviation, color = factor(R0), fill = factor(R0))) +
  geom_line(linewidth = 0.8) + 
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "95% quantile deviation of Day of maximum infection Curve for rec_day =4",
       x = "Connection Probability",
       y = "Quantile deviation")

df3a <- df3[df3$R0==1.2,]

#Plot for 95% quantile deviation of Peak day only for R0 = 1.2
ggplot(df3a, aes(x = Probs, y = deviation, color = factor(R0), fill = factor(R0))) +
  geom_line(linewidth = 0.8) + 
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "95% quantile deviation of Day of maximum infection Curve for rec_day =4",
       x = "Connection Probability",
       y = "Quantile deviation")


#Plot for Max infected
ggplot(sim_df1$max_inf_sd, aes(x = Probs, y = Mean, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Maximum Proportion of new infection on a day Curve with CI (sd) for rec_day = 4",
       x = "Connection Probability",
       y = "Maximum Proportion of new infection")


df4 <- sim_df1$max_inf_sd
df4$sd <- (df4$Mean - df4$Lower)*sqrt(sim)/1.96

#Plot for standard deviation of Max infected
ggplot(df4, aes(x = Probs, y = sd, color = factor(R0), fill = factor(R0))) +
  geom_line(linewidth = 0.8) +  ylim(0,0.1)+
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Standard deviation of Maximum Proportion of new infection on a day Curve for rec_day = 4",
       x = "Connection Probability",
       y = "Standard deviation")

df4a <- df4[df4$R0==1.2,]

#Plot for 95% quantile deviation of Peak day only for R0 = 1.2
ggplot(df4a, aes(x = Probs, y = sd, color = factor(R0), fill = factor(R0))) +
  geom_line(linewidth = 0.8) + ylim(0,0.1)+
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Standard deviation of Maximum Proportion of new infection on a day Curve for rec_day = 4",
       x = "Connection Probability",
       y = "Standard deviation")


#Plot for Peak day
ggplot(sim_df1$peak_sd, aes(x = Probs, y = Mean, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Day of maximum infection Curve with CI (sd) for rec_day = 4",
       x = "Connection Probability",
       y = "Day of maximum infection")


df5 <- sim_df1$peak_sd
df5$sd <- (df5$Mean - df5$Lower)*sqrt(sim)/1.96

#Plot for standard deviation of Peak day
ggplot(df5, aes(x = Probs, y = sd, color = factor(R0), fill = factor(R0))) +
  geom_line(linewidth = 0.8) + ylim(0,15)+ 
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Standard deviation of Day of maximum infection Curve for rec_day = 4",
       x = "Connection Probability",
       y = "Standard deviation")

df5a <- df5[df5$R0==1.2,]

#Plot for 95% quantile deviation of Peak day only for R0 = 1.2
ggplot(df5a, aes(x = Probs, y = sd, color = factor(R0), fill = factor(R0))) +
  geom_line(linewidth = 0.8) + ylim(0,15)+
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Standard deviation of Day of maximum infection Curve for rec_day = 4",
       x = "Connection Probability",
       y = "Standard deviation")


#Plot for End day
ggplot(sim_df1$end_q, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "R0", palette = "Set1") +
  scale_fill_brewer(name = "R0", palette = "Set1") +
  theme_minimal() +
  labs(title = "Total time for epidemic with CI(using quantiles) for rec_day = 4",
       x = "Connection Probability",
       y = "Epidemic duration")



#### Extra code to construct epicurve for specific values of beta, recovery_day and connection probability

# 5. Function for SIR model for epicurve
SIR_epi <- function(states, er_network, infection_rate, recovery_day){
  sus_vector <- c(length(which(states == "S")))
  inf_vector <- c(length(which(states == "I")))
  rec_vector <- c(length(which(states == "R")))
  
  start_day <- rep(-1,length(states))
  start_day[which(states == "I")] <- 0
  t=1
  
  while (sum(states == "I") != 0){
    # infection process
    new_infected <- c()
    infected_nodes <- which(states == "I")
    for (i in infected_nodes){
      neighbors_i <- neighbors(er_network, i)
      s_neighbors <- neighbors_i[states[neighbors_i] == "S"]
      for (s_node in s_neighbors){
        if (runif(1) < infection_rate){
          states[s_node] <- "I"
          new_infected <- c(new_infected, s_node)
        }
      }
    }
    start_day[new_infected] <- t
    
    # recovery process
    if (t>=recovery_day){
      recovered <- which(start_day == t-recovery_day)
      states[recovered] <- "R"
    }
    
    sus_vector <- c(sus_vector, length(which(states == "S")))
    inf_vector <- c(inf_vector, length(which(states == "I")))
    rec_vector <- c(rec_vector, length(which(states == "R")))
    
    t <- t+1
  }
  df <- data.frame(nsus= sus_vector/length(states), ninf = inf_vector/length(states), nrec = rec_vector/length(states))
  return(df)
}

# 6. Function to construct epi curve with confidence interval
epi_curve <- function(sim, nodes, p, b, rec_day, initial_infected, matrices){
  sus <- matrix(0, nrow=1)
  inf <- matrix(0, nrow=1)
  rec <- matrix(0, nrow=1)
  
  for (l in 1:sim){
    # initial states of nodes for a simulation
    states <- rep("S", nodes)
    states[initial_infected] <- "I"
    
    er_graph <- er_graph_unif(p, nodes, matrices[[l]]) # Using the same "sim" number of uniform matrices from the previous part
    values <- SIR_epi(states, er_graph, b, rec_day) # vector of proportion of susceptible, infected and recovered at all time points
    
    sus_values <- as.vector(values$nsus)
    inf_values <- as.vector(values$ninf)
    rec_values <- as.vector(values$nrec)
    
    # Find max length
    len <- max(ncol(sus), length(sus_values))
    
    # Pad both with zeros to the same length
    sus_pad <- matrix(c(sus, rep(0, nrow(sus)*(len - ncol(sus)))), nrow=nrow(sus))
    sus_values_pad <- c(sus_values, rep(0, len - length(sus_values)))
    
    # Add elementwise
    sus <- rbind(sus_pad, sus_values_pad)
    
    # Pad both with zeros to the same length
    inf_pad <- matrix(c(inf, rep(0, nrow(inf)*(len - ncol(inf)))), nrow=nrow(inf))
    inf_values_pad <- c(inf_values, rep(0, len - length(inf_values)))
    
    # Add elementwise
    inf <- rbind(inf_pad, inf_values_pad)
    
    # Pad both with zeros to the same length
    rec_pad <- matrix(c(rec, rep(0, nrow(rec)*(len - ncol(rec)))), nrow=nrow(rec))
    rec_values_pad <- c(rec_values, rep(0, len - length(rec_values)))
    
    # Add elementwise
    rec <- rbind(rec_pad, rec_values_pad)
  }
  
  #Removing the first rows as they have all the entries zero
  inf <- inf[-1,]
  sus <- sus[-1,]
  rec <- rec[-1,]
  
  sus_median = apply(sus, 2, median, na.rm = TRUE)
  sus_lower <- apply(sus, 2, quantile, probs = 0.025, na.rm = TRUE)
  sus_upper <- apply(sus, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  inf_median = apply(inf,2, median, na.rm = TRUE)
  inf_lower <- apply(inf, 2, quantile, probs = 0.025, na.rm = TRUE)
  inf_upper <- apply(inf, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  rec_median = apply(rec,2, median, na.rm = TRUE)
  rec_lower <- apply(rec, 2, quantile, probs = 0.025, na.rm = TRUE)
  rec_upper <- apply(rec, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  df <- data.frame(
    Time = c(c(0:(ncol(inf)-1)),c(0:(ncol(inf)-1)), c(0:(ncol(inf)-1))),
    Variable = c(rep("Susceptible", ncol(sus)), rep("Infected", ncol(inf)), rep("Recovered", ncol(rec))),
    Median = c(sus_median, inf_median, rec_median),
    Upper = c( sus_upper, inf_upper, rec_upper), 
    Lower = c(sus_lower, inf_lower, rec_lower)
  )
  
  return(df)
}

Umatrices <- sim_df1$list

R0_epi <- 1.2
rec_day_epi <- 4 
beta_epi <- R0_epi / rec_day_epi 


epi_df1 <- epi_curve(sim, nodes, 0.25, beta_epi, rec_day_epi, initial_infected, Umatrices)

ggplot(epi_df1, aes(x = Time, y = Median, color = Variable, fill = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "Epi Curve with median values and 95% C.I.(quantile) for p= 0.25, R0 = 1.2, rec_day=4",
       x = "Time",
       y = "Proportions",
       color = "Variable",
       fill = "Variable"
  ) +
  theme(legend.position = "right")

epi_df2 <- epi_curve(sim, nodes, 0.09, beta_epi, rec_day_epi, initial_infected, Umatrices)

ggplot(epi_df2, aes(x = Time, y = Median, color = Variable, fill = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "Epi Curve with median values and 95% C.I.(quantile) for p= 0.09, R0 = 1.2, rec_day=4",
       x = "Time",
       y = "Proportions",
       color = "Variable",
       fill = "Variable"
  ) +
  theme(legend.position = "right")

epi_df3 <- epi_curve(sim, nodes, 0.05, beta_epi, rec_day_epi, initial_infected, Umatrices)

ggplot(epi_df3, aes(x = Time, y = Median, color = Variable, fill = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "Epi Curve with median values and 95% C.I.(quantile) for p= 0.05, R0 = 1.2, rec_day=4",
       x = "Time",
       y = "Proportions",
       color = "Variable",
       fill = "Variable"
  ) +
  theme(legend.position = "right")


epi_df4 <- epi_curve(sim, nodes, 0.65, beta_epi, rec_day_epi, initial_infected, Umatrices)

ggplot(epi_df4, aes(x = Time, y = Median, color = Variable, fill = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "Epi Curve with median values and 95% C.I.(quantile) for p= 0.65, R0 = 1.2, rec_day=4",
       x = "Time",
       y = "Proportions",
       color = "Variable",
       fill = "Variable"
  ) +
  theme(legend.position = "right")


epi_df5 <- epi_curve(sim, nodes, 0.85, beta_epi, rec_day_epi, initial_infected, Umatrices)

ggplot(epi_df5, aes(x = Time, y = Median, color = Variable, fill = Variable)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "Epi Curve with median values and 95% C.I.(quantile) for p= 0.85, R0 = 1.2, rec_day=4",
       x = "Time",
       y = "Proportions",
       color = "Variable",
       fill = "Variable"
  ) +
  theme(legend.position = "right")

#### Extra code for Network visualization

clustering <- function(er_graph){
  #Global Clustering Coefficient
  global <- transitivity(er_graph, type="global")
  
  
  #Local Clustering Coefficient
  local <- transitivity(er_graph, type="local", isolates="zero")
  avg_cc <- mean(local, na.rm=TRUE)
  
  
  #Average Path Length
  avg_path <- mean_distance(er_graph, directed = FALSE)
  
  
  return(c(global, avg_cc, avg_path))
}

centrality <- function(er_graph){
  #Degree Centrality
  deg_cent <- degree(er_graph,mode= "all", normalized = TRUE)
  
  
  #Betweenness Centrality
  bet_cent <- betweenness(er_graph, directed = FALSE, normalized = TRUE)
  
  
  #Closeness Centrality
  clo_cent <- closeness(er_graph, normalized = TRUE)
  
  
  #Eigenvector Centrality
  eig_cent <- eigen_centrality(er_graph)$vector
  
  Centrality_df <- data.frame(Node = V(er_graph), Degree = round(deg_cent,3), 
                              Betweenness = round(bet_cent,3), Closeness = round(clo_cent,3),
                              Eigenvector = round(eig_cent,3))
  return(Centrality_df)
}

# initial states of nodes
states <- rep("S", nodes)
states[initial_infected] <- "I"
# Map states to colors
state_colors <- ifelse(states == "I", "tomato", "green")

l=4 # Matrix is taken from 10th simulation, l can vary from 1 to sim, here sim=15


g1 <- er_graph_unif(0.25, nodes, Umatrices[[l]])
val_1 <- clustering(g1)
centrality1 <- centrality(g1)
c1 <- colMeans(centrality1[,-1])

g2 <- er_graph_unif(0.09, nodes, Umatrices[[l]])
val_2 <- clustering(g2)
centrality2 <- centrality(g2)
c2 <- colMeans(centrality2[,-1])

g3 <- er_graph_unif(0.05, nodes, Umatrices[[l]])
val_3 <- clustering(g3)
centrality3 <- centrality(g3)
c3 <- colMeans(centrality3[,-1])

g4 <- er_graph_unif(0.65, nodes, Umatrices[[l]])
val_4 <- clustering(g4)
centrality4 <- centrality(g4)
c4 <- colMeans(centrality4[,-1])

g5 <- er_graph_unif(0.85, nodes, Umatrices[[l]])
val_5 <- clustering(g5)
centrality5 <- centrality(g5)
c5 <- colMeans(centrality5[,-1])




A <- matrix(c(c(0.05,val_3,c3), c(0.09, val_2,c2),
              c(0.25, val_1, c1), c(0.65, val_4, c4), c(0.85, val_5, c5)),nrow=5, byrow=TRUE)
dd1 <- as.data.frame(A)
colnames(dd1) <- c("Probs", "Global Clustering Coeff", "Avg Local Clustering Coeff", "Avg Path Length",
                   "Avg Degree", "Betweeenness", "Avg Closeness", "Avg Eigen Vector Centrality")

dd1 <- round(dd1, 3)
