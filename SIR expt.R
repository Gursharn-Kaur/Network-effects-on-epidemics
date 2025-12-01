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
  
  median_matrix <- apply(mat_array, c(1, 2), median, na.rm = TRUE) #median of the summary measures for each pair (beta, probs)
  lower_matrix <- apply(mat_array, c(1, 2), function(x) quantile(x, 0.025, na.rm = TRUE)) # 2.5% quantile
  upper_matrix <- apply(mat_array, c(1, 2), function(x) quantile(x, 0.975, na.rm = TRUE)) # 97.5% quantile
  
  num_r0 <- length(beta)
  split_indices <- split(1:(length(probs) * num_r0), rep(1:num_r0, each=length(probs)))
  median_list <- lapply(split_indices, function(idx) median_matrix[idx, ])
  lower_list <- lapply(split_indices, function(idx) lower_matrix[idx, ])
  upper_list <- lapply(split_indices, function(idx) upper_matrix[idx, ])
  
  
  #dataframe for total proportion of infected
  df_total_inf <- data.frame()
  for (i in 1:num_r0){
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Median = median_list[[i]][,1],
      Upper = upper_list[[i]][,1],
      Lower = lower_list[[i]][,1]
    )
    df_total_inf <- rbind(df_total_inf, df_temp)
  }
  
  #dataframe for maximum proportion of infected on a single day
  df_max_inf <- data.frame()
  for (i in 1:num_r0){
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Median = median_list[[i]][,2],
      Upper = upper_list[[i]][,2],
      Lower = lower_list[[i]][,2]
    )
    df_max_inf <- rbind(df_max_inf, df_temp)
  }
  
  #dataframe for the day when number of new infection peaked
  df_peak_day <- data.frame()
  for (i in 1:num_r0){
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Median = median_list[[i]][,3],
      Upper = upper_list[[i]][,3],
      Lower = lower_list[[i]][,3]
    )
    df_peak_day <- rbind(df_peak_day, df_temp)
  }
  
  #dataframe for the duration of epidemic
  df_end_day <- data.frame()
  for (i in 1:num_r0){
    df_temp <- data.frame(
      Probs = probs,
      R0 = beta[i]*recovery_day,
      Median = median_list[[i]][,4],
      Upper = upper_list[[i]][,4],
      Lower = lower_list[[i]][,4]
    )
    df_end_day <- rbind(df_end_day, df_temp)
  }
  
  # Extract the uniform matrices
  matrices <- lapply(result, function(x) x$matrix)
  M <- list(matrices)
  
  #The list of dataframes for summary meaures
  D <- list(
    total_inf = df_total_inf,
    max_inf = df_max_inf,
    peak = df_peak_day,
    end = df_end_day,
    list = matrices
  )
  
  return(D)
}


# Main code :
nodes <- 20
probs <- c(0.001, 0.005, 0.008, 0.01, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
           0.12, 0.14, 0.16, 0.18,0.25, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 0.97)

initial_infected <- sample(1:nodes, 1)

sim <- 15 #No of simulations

R0 <- c(1.2,1.4,1.5)
recovery_day <- 4
beta <- as.vector(R0/recovery_day) # Infection rate

sim_df1 <- Simulated_df(sim, nodes, probs, beta, recovery_day, initial_infected)



#Plot for Total infected
ggplot(sim_df1$total_inf, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "Beta", palette = "Set1") +
  scale_fill_brewer(name = "Beta", palette = "Set1") +
  theme_minimal() +
  labs(title = "Total Proportion of Infected Curve with Confidence Interval",
       x = "Connection Probability",
       y = "Total Proportion of Infected")


#Plot for Max infected
ggplot(sim_df1$max_inf, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "Beta", palette = "Set1") +
  scale_fill_brewer(name = "Beta", palette = "Set1") +
  theme_minimal() +
  labs(title = "Maximum Proportion of new infection on a day Curve with Confidence Interval",
       x = "Connection Probability",
       y = "Maximum Proportion of new infection")


#Plot for Peak day
ggplot(sim_df1$peak, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "Beta", palette = "Set1") +
  scale_fill_brewer(name = "Beta", palette = "Set1") +
  theme_minimal() +
  labs(title = "Day of maximum infection Curve with Confidence Interval",
       x = "Connection Probability",
       y = "Day of maximum infection")


#Plot for End day
ggplot(sim_df1$end, aes(x = Probs, y = Median, color = factor(R0), fill = factor(R0))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_brewer(name = "Beta", palette = "Set1") +
  scale_fill_brewer(name = "Beta", palette = "Set1") +
  theme_minimal() +
  labs(title = "Total time for epidemic with Confidence Interval",
       x = "Connection Probability",
       y = "Epidemic duration")


# Threshold Value of probability exceeding which $99%$ of the individuals got infected
maxp1 <- c()
for(i in 1:length(beta)){
  d1 <- sim_df1$total_inf[((i-1)*length(probs)+1):(i*length(probs)),]
  maxp1 <- c(maxp1, d1$Probs[min(which(d1$Median>= 0.99))])
}
threshold_p_99_infection1 <- data.frame(Beta = beta, p = maxp1)
print(threshold_p_99_infection1)


# Threshold Values of probability at which day of maximum infection peaked
maxp1 <- c()
for(i in 1:length(beta)){
  d1 <- sim_df1$peak[((i-1)*length(probs)+1):(i*length(probs)),]
  maxp1 <- c(maxp1, d1$Probs[which(d1$Median== max(d1$Median))])
}

threshold_p_day_max_infection1 <- data.frame(Beta = beta, p = maxp1)
print(threshold_p_day_max_infection1)


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

matrices <- sim_df1$list

R0_epi <- 1.2
rec_day_epi <- 4 
beta_epi <- R0_epi / rec_day_epi 
prob_epi <- 0.25

epi_df <- epi_curve(sim, nodes, prob_epi, beta_epi, rec_day_epi, initial_infected, matrices)

ggplot(epi_df, aes(x = Time, y = Median, color = Variable, fill = Variable)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  labs(x = "Time",
       y = "Proporions",
       color = "Variable",
       fill = "Variable"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")


#### Extra code for Network visualization

# initial states of nodes
states <- rep("S", nodes)
states[initial_infected] <- "I"
# Map states to colors
state_colors <- ifelse(states == "I", "tomato", "green")

l=10 # Matrix is taken from 10th simulation, l can vary from 1 to sim, here sim=15
p=0.7
er_graph <- er_graph_unif(p, nodes, matrices[[l]])
plot(er_graph, vertex.color=state_colors, vertex.size=5, vertex.label.color="black",
     main="Network Diagram with Nodes Colored by State")
