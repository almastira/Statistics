# Installazione dei pacchetti necessari
install.packages("BiocManager")
BiocManager::install("graph")
BiocManager::install("RBGL")
install.packages("pcalg")
install.packages("bnlearn")
BiocManager::install("Rgraphviz")

# Caricamento dei pacchetti
library(pcalg)
library(bnlearn)
library(graph)
library(Rgraphviz)
library(ggplot2)

# Simulazione di un dataset osservazionale
set.seed(123)
n <- 100  # Numero di osservazioni
p <- 10   # Numero di covariate (geni)

# Simulazione delle covariate (espressione genica)
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("Gene", 1:p)

# Simulazione della risposta (produzione di riboflavina)
beta <- runif(p)
Y <- X %*% beta + rnorm(n)

# Creazione del dataframe
data <- data.frame(Y, X)

# Apprendimento della struttura del grafo
dag <- hc(data)  # Hill-Climbing

# Conversione in CPDAG
cpdag <- cpdag(dag)

# Conversione dell'oggetto CPDAG in un oggetto 'graphNEL'
cpdag_graphNEL <- as.graphNEL(cpdag)

# Funzione per il calcolo degli effetti causali usando 'pcalg'
calc_effects <- function(cpdag, data) {
  p <- ncol(data) - 1
  effects_list <- list()
  max_length <- 0
  for (i in 1:p) {
    effect <- ida(i, 1, cov(data), cpdag)
    effects_list[[i]] <- effect
    if (length(effect) > max_length) {
      max_length <- length(effect)
    }
  }
  # Creazione di una matrice per memorizzare gli effetti
  effects <- matrix(NA, nrow = p, ncol = max_length)
  for (i in 1:p) {
    effects[i, 1:length(effects_list[[i]])] <- effects_list[[i]]
  }
  return(effects)
}

effects <- calc_effects(cpdag_graphNEL, data)

# Creazione della matrice Θ
Theta <- effects

# Nomi delle righe e delle colonne
rownames(Theta) <- colnames(data)[-1]
colnames(Theta) <- paste0("Effect", 1:ncol(Theta))

# Funzione per calcolare il multiset degli effetti causali
calc_multiset <- function(Theta) {
  multiset <- list()
  for (i in 1:nrow(Theta)) {
    multiset[[i]] <- sort(unique(Theta[i, ], na.rm = TRUE))
  }
  return(multiset)
}

multiset_Theta <- calc_multiset(Theta)

# Calcolo del valore assoluto minimo per determinare l'importanza causale
importance <- apply(Theta, 1, function(x) min(abs(x), na.rm = TRUE))

# Ordinamento delle variabili per importanza
important_genes <- order(importance, decreasing = TRUE)

# Visualizzazione dei risultati
important_genes
# Plot del modello appreso
plot(cpdag_graphNEL, main = "CPDAG - Struttura del Grafo Appreso", sub = "Modello Hill-Climbing")



# Creazione di un dataframe per il multiset degli effetti causali
multiset_df <- data.frame(
  Gene = rep(rownames(Theta), each = ncol(Theta)),
  Effect = as.vector(Theta)
)

# Rimozione dei valori NA
multiset_df <- multiset_df[!is.na(multiset_df$Effect), ]

# Creazione del boxplot
ggplot(multiset_df, aes(x = Gene, y = Effect)) +
  geom_boxplot() +
  labs(title = "Distribuzione degli Effetti Causali per Gene",
       x = "Gene", y = "Effetto Causale") +
  theme_minimal()




##secondo punto
# Matrice di adiacenza
suffStat <- list(C = cor(data), n = n)

# Applicazione dell'algoritmo PC
pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = 0.05, labels = colnames(data))

# Conversione in CPDAG
cpdag <- as(pc.fit@graph, "graphNEL")
# Visualizzazione del CPDAG
plot(cpdag)


####CALCOLO DEGLI INTERVENTI
# Funzione per calcolare la distribuzione markoviana
markovian_dist <- function(data, dag) {
  p <- ncol(data) - 1
  dist_list <- list()
  for (i in 1:(p+1)) {
    pa_i <- which(as(dag, "matrix")[, i] != 0)
    if (length(pa_i) == 0) {
      dist_list[[i]] <- function(x) {
        return(dnorm(x, mean = mean(data[, i]), sd = sd(data[, i])))
      }
    } else {
      dist_list[[i]] <- function(x, pa) {
        return(dnorm(x, mean = lm(data[, i] ~ data[, pa_i])$fitted.values, sd = sd(residuals(lm(data[, i] ~ data[, pa_i])))))
      }
    }
  }
  return(dist_list)
}

# Calcolo della distribuzione markoviana
markov_dist <- markovian_dist(data, cpdag)

# Funzione per calcolare la distribuzione di Y dopo un intervento
calc_intervention_Y_dist <- function(dag, data, intervention_var, intervention_val) {
  p <- ncol(data) - 1
  pa_i <- which(as(dag, "matrix")[, intervention_var] != 0)
  
  if (length(pa_i) == 0) {
    return(mean(data$Y))
  } else {
    f_pa_i <- function(pa_i) {
      prod(sapply(pa_i, function(pa) dnorm(pa, mean = mean(data[, pa]), sd = sd(data[, pa]))))
    }
    E_Y_given_Xi_pa <- lm(data$Y ~ data[, intervention_var] + data[, pa_i])$fitted.values
    return(integrate(function(pa) {
      E_Y_given_Xi_pa * f_pa_i(pa)
    }, lower = -Inf, upper = Inf)$value)
  }
}

# Esempio di utilizzo
intervention_var <- 2
intervention_val <- 1
Y_dist <- calc_intervention_Y_dist(cpdag, data, intervention_var, intervention_val)

# Funzione per calcolare la media della distribuzione dopo un intervento
calc_intervention_mean <- function(dag, data, intervention_var, intervention_val) {
  pa_i <- which(as(dag, "matrix")[, intervention_var] != 0)
  
  if (length(pa_i) == 0) {
    mean_Y <- mean(data$Y)
  } else {
    f_pa_i <- function(pa_i) {
      prod(sapply(pa_i, function(pa) dnorm(pa, mean = mean(data[, pa]), sd = sd(data[, pa]))))
    }
    E_Y_given_Xi_pa <- lm(data$Y ~ data[, intervention_var] + data[, pa_i])$fitted.values
    mean_Y <- integrate(function(pa) {
      E_Y_given_Xi_pa * f_pa_i(pa)
    }, lower = -Inf, upper = Inf)$value
  }
  return(mean_Y)
}

# Esempio di utilizzo
intervention_var <- 2
intervention_val <- 1
mean_Y <- calc_intervention_mean(cpdag, data, intervention_var, intervention_val)


###CALCOLO DEGLI INTERVENTI VS ASSOCIAZIONE
# Installazione e caricamento del pacchetto 'lm' per la regressione lineare

library(MASS)

# Esempio 2.1
set.seed(123)
n <- 100
epsilon1 <- rnorm(n, mean = 0, sd = sqrt(0.36))
epsilon2 <- rnorm(n, mean = 0, sd = sqrt(1))
epsilon3 <- rnorm(n, mean = 0, sd = sqrt(0.36))
epsilon <- rnorm(n, mean = 0, sd = sqrt(1))

X2 <- epsilon2
X1 <- 0.8 * X2 + epsilon1
X3 <- 0.8 * X2 + epsilon3
Y <- -X1 + 2 * X2 - X3 + epsilon

data1 <- data.frame(X1, X2, X3, Y)

# Esempio 2.2
Y2 <- X1 + X3 + epsilon
data2 <- data.frame(X1, X2, X3, Y = Y2)

# Esempio 2.1: Regressione lineare multipla
fit1 <- lm(Y ~ X1 + X2 + X3, data = data1)
summary(fit1)

# Esempio 2.2: Regressione lineare multipla
fit2 <- lm(Y ~ X1 + X2 + X3, data = data2)
summary(fit2)

#ESMPIO 2.1
# Effetti causali diretti per Esempio 2.1
theta1 <- -1  # Effetto di X1 su Y
theta2 <- 0.4  # Effetto di X2 su Y (quando X2 è considerato separatamente)
theta3 <- -1  # Effetto di X3 su Y

# Effetti causali diretti per Esempio 2.2
theta1_2 <- 1  # Effetto di X1 su Y
theta2_2 <- 0  # Effetto di X2 su Y (quando X2 è considerato separatamente)
theta3_2 <- 1  # Effetto di X3 su Y

# Coefficienti di regressione per Esempio 2.1
coef_regression_1 <- summary(fit1)$coefficients[, "Estimate"]

# Effetti causali calcolati per Esempio 2.1
effects_causal_1 <- c(theta1, theta2, theta3)

# Visualizzazione del confronto
print("Esempio 2.1 - Confronto tra coefficienti di regressione e effetti causali:")
print(data.frame(Variable = c("X1", "X2", "X3"), Regression = coef_regression_1[-1], Causal = effects_causal_1))

# Coefficienti di regressione per Esempio 2.2
coef_regression_2 <- summary(fit2)$coefficients[, "Estimate"]

# Effetti causali calcolati per Esempio 2.2
effects_causal_2 <- c(theta1_2, theta2_2, theta3_2)



library(ggplot2)

# Creazione dei dati per il grafico a barre
coefficients <- data.frame(
  Variable = rep(c("X1", "X2", "X3"), 2),
  Value = c(coef(fit1)[-1], coef(fit2)[-1]),
  Type = rep(c("Regression", "Causal"), each = 3)
)

effects_causal <- data.frame(
  Variable = rep(c("X1", "X2", "X3"), 2),
  Value = c(theta1, theta2, theta3, theta1_2, theta2_2, theta3_2),
  Type = rep(c("Regression", "Causal"), each = 3)
)

# Unione dei dati
data_combined <- rbind(coefficients, effects_causal)

# Grafico a barre
ggplot(data_combined, aes(x = Variable, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Type) +
  labs(title = "Confronto tra Coefficienti di Regressione e Effetti Causali",
       y = "Value", x = "Variable") +
  theme_minimal()



# Visualizzazione del confronto
print("Esempio 2.2 - Confronto tra coefficienti di regressione e effetti causali:")
print(data.frame(Variable = c("X1", "X2", "X3"), Regression = coef_regression_2[-1], Causal = effects_causal_2))

# Esempio 2.1: Definizione del DAG
dag1 <- empty.graph(c("X1", "X2", "X3", "Y"))
modelstring(dag1) <- "[X2][X1|X2][X3|X2][Y|X1:X2:X3]"
graphviz.plot(dag1, main = "Esempio 2.1")

# Esempio 2.2: Definizione del DAG
dag2 <- empty.graph(c("X1", "X2", "X3", "Y"))
modelstring(dag2) <- "[X2][X1][X3][Y|X1:X3]"
graphviz.plot(dag2, main = "Esempio 2.2")

# Esempio 2.1: Diagramma di dispersione
pairs(data1, main = "Esempio 2.1: Relazioni tra le variabili")

# Esempio 2.2: Diagramma di dispersione
pairs(data2, main = "Esempio 2.2: Relazioni tra le variabili")


#PUNTO 3
# Definizione della struttura del DAG noto (popolazione)
adj_matrix <- matrix(c(0, 1, 1, 0, 0,
                       0, 0, 1, 1, 0,
                       0, 0, 0, 1, 1,
                       0, 0, 0, 0, 1,
                       0, 0, 0, 0, 0), 
                     byrow = TRUE, nrow = 5)

# Funzione per generare i dati
generate_data <- function(n, adj_matrix) {
  p <- ncol(adj_matrix)
  X <- matrix(0, n, p)
  for (i in 1:p) {
    parents <- which(adj_matrix[, i] == 1)
    if (length(parents) == 0) {
      X[, i] <- rnorm(n)
    } else if (length(parents) == 1) {
      X[, i] <- X[, parents] + rnorm(n)
    } else {
      X[, i] <- rowSums(X[, parents]) + rnorm(n)
    }
  }
  return(X)
}

# Generazione dei dati con il DAG noto
n <- 100  # numero di campioni
data_pop <- generate_data(n, adj_matrix)
colnames(data_pop) <- paste0("X", 1:5)
data_pop <- data.frame(data_pop)
print(head(data_pop))

# Creazione del grafo noto
nodes <- colnames(data_pop)
edges <- list(
  X1 = c("X2", "X3"),
  X2 = c("X3", "X4"),
  X3 = c("X4", "X5"),
  X4 = c("X5"),
  X5 = c()
)

# Conversione del grafo in formato graphNEL
dag_pop_graphNEL <- new("graphNEL", nodes = nodes, edgeL = edges, edgemode = "directed")

# Conversione del DAG in CPDAG
cpdag_pop <- dag2cpdag(dag_pop_graphNEL)
plot(cpdag_pop)


# Funzione per calcolare gli effetti causali nella popolazione
calc_population_effects <- function(dag, data) {
  p <- ncol(data)
  effects_list <- list()
  max_length <- 0
  for (i in 1:(p-1)) {  # Modificato per considerare tutte le variabili tranne l'ultima (la risposta)
    effect <- ida(i, p, cov(data), dag)
    effects_list[[i]] <- effect
    if (length(effect) > max_length) {
      max_length <- length(effect)
    }
  }
  # Creazione di una matrice per memorizzare gli effetti
  effects <- matrix(NA, nrow = p-1, ncol = max_length)  # Modificato per nrow = p-1
  for (i in 1:(p-1)) {
    effects[i, 1:length(effects_list[[i]])] <- effects_list[[i]]
  }
  return(effects)
}

# Calcolo degli effetti nella popolazione
population_effects <- calc_population_effects(cpdag_pop, data_pop)

# Creazione della matrice Θ per la popolazione
Theta_population <- population_effects

# Nomi delle righe e delle colonne
rownames(Theta_population) <- colnames(data_pop)[-ncol(data_pop)]  # Modificato per escludere l'ultima colonna
colnames(Theta_population) <- paste0("Effect", 1:ncol(Theta_population))

# Visualizzazione della matrice Theta della popolazione
print(Theta_population)


#3.2 algoritmi
#ottieni i fratelli
subsets <- function(set) {
  n <- length(set)
  all_subsets <- list()
  
  for (i in 0:(2^n - 1)) {
    subset <- c()
    for (j in 1:n) {
      if (bitwAnd(i, 2^(j-1)) > 0) {
        subset <- c(subset, set[j])
      }
    }
    all_subsets <- c(all_subsets, list(subset))
  }
  
  return(all_subsets)
}

# Funzione per calcolare gli effetti causali
computeCausalEffect <- function(dag, X, Y, S = NULL) {
  # Questa funzione calcola l'effetto causale
  # Implementazione di esempio (sostituire con la logica effettiva)
  return(runif(1))  # numero casuale per ora
}

# Funzione per ottenere i DAG estendibili
extendableDAGs <- function(cpdag, S) {
  # Questa funzione restituisce i DAG estendibili
  # Implementazione di esempio (necessita di logica effettiva)
  return(list(cpdag))  # restituisce il CPDAG originale per ora
}


isLocallyValid <- function(cpdag, S, node) {
  # Questa funzione verifica se il DAG è localmente valido
  # Per semplicità, supponiamo che sia sempre valido (necessita di logica effettiva)
  return(TRUE)
}
dagEquivalenceClass <- function(cpdag) {
  # Questa funzione restituisce tutti i DAG nella classe di equivalenza
  # Implementazione di esempio (necessita di logica effettiva)
  return(list(cpdag))  # restituisce il CPDAG originale per ora
}


# Funzione per ottenere i fratelli in un CPDAG
sib <- function(cpdag, node) {
  adj <- as(cpdag, "matrix")
  siblings <- which((adj[node, ] == 1 & adj[, node] == 1))
  return(names(siblings))
}



algorithm_1 <- function(cpdag, X, Y) {
  dag_list <- dagEquivalenceClass(cpdag)
  m <- length(dag_list)
  p <- length(X)
  theta <- matrix(0, p, m)
  
  for (j in 1:m) {
    for (i in 1:p) {
      theta[i, j] <- computeCausalEffect(dag_list[[j]], X[i], Y)
    }
  }
  return(theta)
}

# Esempio di utilizzo
theta_matrix <- algorithm_1(cpdag_pop, colnames(data_pop)[-ncol(data_pop)], "X5")
print(theta_matrix)

# Implementazione alternativa dell'Algoritmo 2
algorithm_2 <- function(cpdag, X, Y) {
  p <- length(X)
  theta_list <- vector("list", p)
  names(theta_list) <- X
  
  for (i in 1:p) {
    node <- X[i]
    theta_list[[node]] <- c()
    sibi <- sib(cpdag, node)
    subsets_list <- subsets(sibi)
    
    for (S in subsets_list) {
      extendable <- extendableDAGs(cpdag, S)
      mS <- length(extendable)
      if (mS > 0) {
        effect <- computeCausalEffect(cpdag, node, Y, S)
        theta_list[[node]] <- c(theta_list[[node]], rep(effect, mS))
      }
    }
  }
  return(theta_list)
}

# Esempio di utilizzo
theta_multisets <- algorithm_2(cpdag_pop, colnames(data_pop)[-ncol(data_pop)], "X5")
print(theta_multisets)

# Implementazione dell'Algoritmo 3
algorithm_3 <- function(cpdag, X, Y) {
  p <- length(X)
  theta_local_list <- vector("list", p)
  names(theta_local_list) <- X
  
  for (i in 1:p) {
    node <- X[i]
    theta_local_list[[node]] <- c()
    sibi <- sib(cpdag, node)
    subsets_list <- subsets(sibi)
    
    for (S in subsets_list) {
      if (isLocallyValid(cpdag, S, node)) {
        effect <- computeCausalEffect(cpdag, node, Y, S)
        theta_local_list[[node]] <- c(theta_local_list[[node]], effect)
      }
    }
  }
  return(theta_local_list)
}

# Esempio di utilizzo
theta_local_multisets <- algorithm_3(cpdag_pop, colnames(data_pop)[-ncol(data_pop)], "X5")
print(theta_local_multisets)


####PUNTO 4
# Funzione per determinare tutti i DAG nella classe di equivalenza di G
enumerate_dags <- function(cpdag) {
  return(pdag2allDags(cpdag))
}
# Algoritmo di Base
algorithm_base <- function(cpdag, data) {
  all_dags <- enumerate_dags(cpdag)
  p <- ncol(data)
  effects <- list()
  
  for (j in 1:length(all_dags)) {
    dag <- all_dags[[j]]
    for (i in 1:(p-1)) {
      theta_ij <- ida(i, p, cov(data), dag)
      effects[[paste0("G", j, "_X", i)]] <- theta_ij
    }
  }
  return(effects)
}
# Funzione per generare i sottoinsiemi di un insieme
subsets <- function(s) {
  l <- length(s)
  return(lapply(0:(2^l - 1), function(x) s[which(as.logical(intToBits(x) %/% 2^(0:(l-1)) %% 2))]))
}

# Variazione dell'Algoritmo 1
algorithm_variation <- function(cpdag, data) {
  p <- ncol(data)
  effects <- vector("list", p-1)
  
  for (i in 1:(p-1)) {
    effects[[i]] <- c()
    sib_i <- which(cpdag[i,] != 0)
    
    for (S in subsets(sib_i)) {
      G_S <- cpdag
      G_S[i, S] <- 0  # Rimuove gli archi da X_i a S
      G_S <- as(G_S, "graphNEL")
      
      if (isValidGraph(G_S)) {
        m_S <- length(pdag2allDags(G_S))
        effects[[i]] <- c(effects[[i]], rep(ida(i, p, cov(data), G_S), m_S))
      }
    }
  }
  return(effects)
}
# Funzione di utilità per verificare che non ci siano nuove strutture-v
no_new_v_structures <- function(G_S, i) {
  # Implementa la logica per verificare che non ci siano nuove strutture-v con colliders X_i
  return(TRUE)  # Placeholder
}

# Algoritmo Locale
algorithm_local <- function(cpdag, data) {
  p <- ncol(data)
  effects <- vector("list", p-1)
  
  for (i in 1:(p-1)) {
    effects[[i]] <- c()
    sib_i <- which(cpdag[i,] != 0)
    
    for (S in subsets(sib_i)) {
      G_S <- cpdag
      G_S[i, S] <- 0  # Rimuove gli archi da X_i a S
      G_S <- as(G_S, "graphNEL")
      
      if (isValidGraph(G_S)) {
        if (no_new_v_structures(G_S, i)) {
          effects[[i]] <- c(effects[[i]], ida(i, p, cov(data), G_S))
        }
      }
    }
  }
  return(effects)
}

# Funzione di confronto
compare_effects <- function(effects_sample, expected_effects) {
  comparison <- data.frame(
    Variable = names(effects_sample),
    Calculated = sapply(effects_sample, function(x) min(abs(unlist(x)), na.rm = TRUE)),
    Expected = expected_effects
  )
  return(comparison)
}
# Funzione per rilevare e risolvere strutture-v in conflitto
resolve_conflict_v_structures <- function(cpdag) {
  # Implementazione della risoluzione dei conflitti di strutture-v
  # Placeholder: seguiamo semplicemente l'ordine delle variabili
  return(cpdag)
}
# Funzione per controllare se un CPDAG è estendibile
is_extendible <- function(cpdag) {
  # Implementa la logica per verificare se il CPDAG è estendibile
  return(TRUE)  # Placeholder
}

# Funzione per rimuovere le strutture-v in conflitto
remove_conflict_v_structures <- function(cpdag) {
  # Implementazione della rimozione delle strutture-v in conflitto
  return(cpdag)
}

# Funzione per ottenere un CPDAG estendibile
get_extendible_cpdag <- function(cpdag) {
  cpdag <- resolve_conflict_v_structures(cpdag)
  if (is_extendible(cpdag)) {
    return(cpdag)
  } else {
    cpdag <- remove_conflict_v_structures(cpdag)
    if (is_extendible(cpdag)) {
      return(cpdag)
    } else {
      # Se non riusciamo a risolvere, scegliamo casualmente un DAG nel CPDAG stimato
      all_dags <- enumerate_dags(cpdag)
      chosen_dag <- sample(all_dags, 1)
      cpdag <- dag2cpdag(chosen_dag)
      return(cpdag)
    }
  }
}


#_________________________________________________________________________________
#SIMULAZIONI E ANALISI DI DATI 
#Installa e carica i pacchetti necessari
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RBGL")
BiocManager::install("Rgraphviz")
BiocManager::install("graph")
install.packages("pcalg")
install.packages("huge")
install.packages("ggplot2")
install.packages("microbenchmark")
install.packages("hdi")  

# Carica le librerie necessarie
library(graph)
library(pcalg)
library(Rgraphviz)
library(huge)
library(ggplot2)
library(microbenchmark)
library(hdi)
#______________________________
#FUNZIONI

#Funzione per generare dati DAG casuali
generate_DAG_data <- function(p, en, n, nreps) {
  results <- list()
  
  for (k in 1:nreps) {
    # Genera DAG casuale
    g <- randomDAG(p, prob = en / p, lB = 1, uB = 2)
    
    # Simulare i dati dal DAG
    X <- rmvDAG(n, g, errDist = "normal")
    
    # Scegli casualmente la variabile di risposta e la covariata di interesse
    set.seed(k)
    response_var <- sample(1:p, 1)
    covariate_var <- sample(setdiff(1:p, response_var), 1)
    
    results[[k]] <- list(g = g, X = X, response_var = response_var, covariate_var = covariate_var)
  }
  
  return(results)
}

# Funzione per analizzare i dati generati con l'algoritmo 1 (PC)
analyze_data_algorithm1 <- function(data, alpha = 0.01) {
  results <- list()
  
  for (k in 1:length(data)) {
    g <- data[[k]]$g
    X <- data[[k]]$X
    response_var <- data[[k]]$response_var
    covariate_var <- data[[k]]$covariate_var
    
    # Utilizzare l'algoritmo PC per stimare CPDAG
    suffStat <- list(C = cor(X), n = nrow(X))
    pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(X), alpha = alpha)
    
    # Analizzare il CPDAG stimato
    results[[k]] <- list(pc.fit = pc.fit, response_var = response_var, covariate_var = covariate_var)
  }
  
  return(results)
}

# Funzione per analizzare i dati generati con l'algoritmo 3 (versione PC modificata)
analyze_data_algorithm3 <- function(data, alpha = 0.01) {
  results <- list()
  
  for (k in 1:length(data)) {
    g <- data[[k]]$g
    X <- data[[k]]$X
    response_var <- data[[k]]$response_var
    covariate_var <- data[[k]]$covariate_var
    
    # Utilizzare l'algoritmo del PC locale per stimare CPDAG
    suffStat <- list(C = cor(X), n = nrow(X))
    pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(X), alpha = alpha)
    
    # Analizzare il CPDAG stimato
    results[[k]] <- list(pc.fit = pc.fit, response_var = response_var, covariate_var = covariate_var)
  }
  
  return(results)
}

# Funzione per tracciare un DAG stimato
plot_DAG <- function(dag, title) {
  graph <- as(dag@graph, "graphNEL")
  plot(graph, main = title)
}

# Funzione per generare dati DAG casuali con struttura a blocchi
generate_DAG_data_block <- function(p, en, n, nreps, blocks) {
  results <- list()
  
  for (k in 1:nreps) {
    # Genera DAG casuale con struttura a blocchi
    g <- randomDAG(p, prob = en / (p / blocks), lB = 1, uB = 2)
    X <- rmvDAG(n, g, errDist = "normal")
    
    # Assicurati che la covariata e la risposta siano nello stesso blocco
    block_size <- p / blocks
    block <- sample(1:blocks, 1)
    response_var <- sample(((block - 1) * block_size + 1):(block * block_size), 1)
    covariate_var <- sample(setdiff(((block - 1) * block_size + 1):(block * block_size), response_var), 1)
    
    results[[k]] <- list(g = g, X = X, response_var = response_var, covariate_var = covariate_var)
  }
  
  return(results)
}


#_____________________________________________
#SETTINGS 1 CON ALG 1 E ALG3

# Parametri di configurazione per Setting 1
p1 <- 9  # p + 1 = 10 vértices
en1 <- 4
n1a <- 20
n1b <- 2000
nreps1 <- 1  # Usamos un solo ejemplo para graficar

# Genera dati per Setting 1a e Setting 1b
set.seed(123)
setting1a_data <- generate_DAG_data(p1, en1, n1a, nreps1)
setting1b_data <- generate_DAG_data(p1, en1, n1b, nreps1)

# Analizzare i dati generati con l'algoritmo 1 (PC)
results1a_alg1 <- analyze_data_algorithm1(setting1a_data)
results1b_alg1 <- analyze_data_algorithm1(setting1b_data)

# Analizzare i dati generati con l'algoritmo 3 (versione PC modificata)
results1a_alg3 <- analyze_data_algorithm3(setting1a_data)
results1b_alg3 <- analyze_data_algorithm3(setting1b_data)

# Tracciare i risultati per l'impostazione 1a con l'algoritmo 3 (versione PC modificata)
plot_DAG(results1a_alg3[[1]]$pc.fit, "Setting 1a con Algoritmo 3: p=9, n=20")

# Tracciare i risultati per l'impostazione 1b con l'algoritmo 3 (versione PC modificata)
plot_DAG(results1b_alg3[[1]]$pc.fit, "Setting 1b con Algoritmo 3: p=9, n=2000")

# Traccia i risultati per l'impostazione 1a con l'algoritmo 1
plot_DAG(results1a_alg1[[1]]$pc.fit, "Setting 1a con Algoritmo 1: p=9, n=20")

# Traccia i risultati per l'impostazione 1b con l'algoritmo 1
plot_DAG(results1b_alg1[[1]]$pc.fit, "Setting 1b con Algoritmo 1: p=9, n=2000")

#SETTING 2 CON ALG 3
# Parametri di configurazione per Setting 2
p2 <- 999  # p + 1 = 1000 vertici
en2 <- 4
n2 <- 100
nreps2 <- 1  # Usiamo un singolo esempio per rappresentare graficamente

# Genera dati per Setting 2
set.seed(123)
setting2_data <- generate_DAG_data_block(p2, en2, n2, nreps2, blocks = 100)

# Analizzare i dati generati con l'algoritmo 3
results2_alg3 <- analyze_data_algorithm3(setting2_data)

#___________FUNZIONE MSE_________

# Funzione per calcolare l'errore quadratico medio (MSE)
calculate_mse <- function(true_dag, estimated_dag) {
  true_adj <- as(true_dag, "matrix")
  estimated_adj <- wgtMatrix(estimated_dag@graph, transpose = FALSE)
  mse <- mean((true_adj - estimated_adj)^2)
  return(mse)
}

#____________________FUNZIONE MAE___________
# Funzione per calcolare l'errore medio assoluto (MAE)
calculate_mae <- function(true_dag, estimated_dag) {
  true_adj <- as(true_dag, "matrix")
  estimated_adj <- wgtMatrix(estimated_dag@graph, transpose = FALSE)
  mae <- mean(abs(true_adj - estimated_adj))
  return(mae)
}

#__________________________________________________________________
#PARAMETRI
# Parametri di configurazione per Setting 1
p1 <- 9  # p + 1 = 10 vertici
en1 <- 4
n1a <- 20
n1b <- 2000
nreps1 <- 10  # Per ottenere un migliore confronto statistico

# Genera dati per Setting 1a e Setting 1b
set.seed(123)
setting1a_data <- generate_DAG_data(p1, en1, n1a, nreps1)
setting1b_data <- generate_DAG_data(p1, en1, n1b, nreps1)

# Parametri di configurazione per Setting 2
p2 <- 999  # p + 1 = 1000 vertici
en2 <- 4
n2 <- 100
nreps2 <- 10  # Per ottenere un migliore confronto statistico

# Genera dati per Setting 2
set.seed(123)
setting2_data <- generate_DAG_data_block(p2, en2, n2, nreps2, blocks = 100)

# Inizializza gli elenchi per memorizzare MSE e MAE
mse1a_alg1 <- numeric(nreps1)
mse1a_alg3 <- numeric(nreps1)
mse1b_alg1 <- numeric(nreps1)
mse1b_alg3 <- numeric(nreps1)
mae1a_alg1 <- numeric(nreps1)
mae1a_alg3 <- numeric(nreps1)
mae1b_alg1 <- numeric(nreps1)
mae1b_alg3 <- numeric(nreps1)

mse2_alg3 <- numeric(nreps2)
mae2_alg3 <- numeric(nreps2)


# Calcolare MSE e MAE per Setting 1a e 1b utilizzando l'algoritmo 1 e 3
for (i in 1:nreps1) {
  results1a_alg1 <- analyze_data_algorithm1(list(setting1a_data[[i]]))
  results1a_alg3 <- analyze_data_algorithm3(list(setting1a_data[[i]]))
  mse1a_alg1[i] <- calculate_mse(setting1a_data[[i]]$g, results1a_alg1[[1]]$pc.fit)
  mse1a_alg3[i] <- calculate_mse(setting1a_data[[i]]$g, results1a_alg3[[1]]$pc.fit)
  mae1a_alg1[i] <- calculate_mae(setting1a_data[[i]]$g, results1a_alg3[[1]]$pc.fit)
  mae1a_alg3[i] <- calculate_mae(setting1a_data[[i]]$g, results1a_alg3[[1]]$pc.fit)
  
  results1b_alg1 <- analyze_data_algorithm1(list(setting1b_data[[i]]))
  results1b_alg3 <- analyze_data_algorithm3(list(setting1b_data[[i]]))
  mse1b_alg1[i] <- calculate_mse(setting1b_data[[i]]$g, results1b_alg1[[1]]$pc.fit)
  mse1b_alg3[i] <- calculate_mse(setting1b_data[[i]]$g, results1b_alg3[[1]]$pc.fit)
  mae1b_alg1[i] <- calculate_mae(setting1a_data[[i]]$g, results1a_alg3[[1]]$pc.fit)
  mae1b_alg3[i] <- calculate_mae(setting1a_data[[i]]$g, results1a_alg3[[1]]$pc.fit)
}

# Mostra risultati MSE E MAE
cat("Setting 1a (p=9, n=20):\n")
cat("MSE Algoritmo 1:", mse1a_alg1, "\n\n")
cat("MSE Algoritmo 3:", mse1a_alg3, "\n\n")
cat("MAE Algoritmo 1:", mae1a_alg1, "\n\n")
cat("MAE Algoritmo 3:", mae1a_alg3, "\n\n")

cat("Setting 1b (p=9, n=2000):\n")
cat("MSE Algoritmo 1:", mse1b_alg1, "\n\n")
cat("MSE Algoritmo 3:", mse1b_alg3, "\n\n")
cat("MAE Algoritmo 1:", mae1b_alg1, "\n\n")
cat("MAE Algoritmo 3:", mae1b_alg3, "\n\n")


# Calcola MSE e MAE per l'impostazione 2 utilizzando l'algoritmo 3
for (i in 1:nreps2) {
  results2_alg3 <- analyze_data_algorithm3(list(setting2_data[[i]]))
  mse2_alg3[i] <- calculate_mse(setting2_data[[i]]$g, results2_alg3[[1]]$pc.fit)
  mae2_alg3[i] <- calculate_mae(setting1a_data[[i]]$g, results1a_alg3[[1]]$pc.fit)
}

# Mostra risultati MSE e MAE
cat("Setting 1a (p=999, n=100):\n")
cat("MSE Algoritmo 3:", mse2_alg1, "\n\n")
cat("MAE Algoritmo 3:", mae2_alg3, "\n\n")


#________________MSE E MAE SETTING 2______nrep =1
# Calcolare l'MSE per Setting 2
mse2_alg3 <- calculate_mse(setting2_data[[1]]$g, results2_alg3[[1]]$pc.fit)

# Calcolare l'MAE per Setting 2
mae2_alg3 <- calculate_mae(setting2_data[[1]]$g, results2_alg3[[1]]$pc.fit)

# Mostra i risultati MSE
cat("Setting 2 (p=999, n=100):\n")
cat("MSE Algoritmo 3:", mse2_alg3, "\n")
cat("MAE Algoritmo 3:", mae2_alg3, "\n")

#CONFRONTO RISULTATI MSE E MAE
# Mostra i risultati
results <- data.frame(
  Setting = c("Setting 1a", "Setting 1b" , "Setting 2"),
  Algorithm = c("Algorithm 1", "Algorithm 3", "Algorithm 3"), 
  MSE = c(mse1a_alg1, mse1b_alg3, mse2_alg3),
  MAE = c(mae1a_alg1, mae1b_alg3, mae2_alg3)
)

#TAVOLA MSE E MAE
print(results)

# boxplot MSE
ggplot(results, aes(x = interaction(Setting, Algorithm), y = MSE, fill = Algorithm)) +
  geom_boxplot() +
  labs(title = "Comparación del Error Cuadrático Medio (MSE) entre Settings y Algoritmos",
       x = "Combinación de Setting y Algoritmo",
       y = "Error Cuadrático Medio (MSE)") +
  theme_minimal()

#  boxplot MAE
ggplot(results, aes(x = interaction(Setting, Algorithm), y = MAE, fill = Algorithm)) +
  geom_boxplot() +
  labs(title = "Comparación del Error Absoluto Medio (MAE) entre Settings y Algoritmos",
       x = "Combinación de Setting y Algoritmo",
       y = "Error Absoluto Medio (MAE)") +
  theme_minimal()
#____________________________________________________________________________________________
#MULTISET
# Funzione per generare un DAG casuale con pesi dei bordi
generate_fixed_DAG <- function(p, en) {
  set.seed(123)
  g <- randomDAG(p, prob = en / p, lB = 1, uB = 2)
  return(g)
}

# Funzione per generare più set di dati da un DAG
generate_multiset_data <- function(g, n, nreps) {
  results <- list()
  
  for (k in 1:nreps) {
    # Simulare i dati dal DAG
    X <- rmvDAG(n, g, errDist = "normal")
    
    # Scegli casualmente la variabile di risposta e la covariata di interesse
    set.seed(k)
    response_var <- sample(1:ncol(X), 1)
    covariate_var <- sample(setdiff(1:ncol(X), response_var), 1)
    
    results[[k]] <- list(X = X, response_var = response_var, covariate_var = covariate_var)
  }
  
  return(results)
}


# Funzione per analizzare i dati generati con l'algoritmo 1 (PC)
analyze_data_algorithm1 <- function(data, alpha = 0.01) {
  results <- list()
  
  for (k in 1:length(data)) {
    X <- data[[k]]$X
    response_var <- data[[k]]$response_var
    covariate_var <- data[[k]]$covariate_var
    
    # Utilizzare l'algoritmo PC per stimare CPDAG
    suffStat <- list(C = cor(X), n = nrow(X))
    pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(X), alpha = alpha)
    
    results[[k]] <- list(pc.fit = pc.fit, response_var = response_var, covariate_var = covariate_var)
  }
  
  return(results)
}

# Funzione per analizzare i dati generati con l'algoritmo 3 (versione PC modificata)
analyze_data_algorithm3 <- function(data, alpha = 0.01) {
  results <- list()
  
  for (k in 1:length(data)) {
    X <- data[[k]]$X
    response_var <- data[[k]]$response_var
    covariate_var <- data[[k]]$covariate_var
    
    # Utilizzare l'algoritmo del PC locale per stimare CPDAG
    suffStat <- list(C = cor(X), n = nrow(X))
    pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(X), alpha = alpha)
    
    results[[k]] <- list(pc.fit = pc.fit, response_var = response_var, covariate_var = covariate_var)
  }
  
  return(results)
}

# Genera un DAG fisso con p = 7 e en = 3
p_fixed <- 7
en_fixed <- 3
g_fixed <- generate_fixed_DAG(p_fixed, en_fixed)

# Genera 50 set di dati di dimensione 1000 dal DAG fisso
n_samples <- 1000
nreps <- 50
data_multiset <- generate_multiset_data(g_fixed, n_samples, nreps)

# Inizializza gli elenchi per archiviare gli MSE
mse_alg3 <- numeric(nreps)

# Calcola gli MSE utilizzando l'algoritmo 1 e 3
for (i in 1:nreps) {
  results_alg3 <- analyze_data_algorithm3(list(data_multiset[[i]]))
  mse_alg3[i] <- calculate_mse(g_fixed, results_alg3[[1]]$pc.fit)
}

# Inizializza gli elenchi per archiviare i MAE
mae_alg3 <- numeric(nreps)

# Calcola MAE utilizzando l'algoritmo 1 e 3
for (i in 1:nreps) {
  results_alg3 <- analyze_data_algorithm3(list(data_multiset[[i]]))
  mae_alg3[i] <- calculate_mae(g_fixed, results_alg3[[1]]$pc.fit)
}

# Creare un frame di dati per i risultati del grafico della densità
data_density <- data.frame(
  MSE = c(mse_alg3),
  MAE = c(mae_alg3),
  Algorithm = factor(rep(c("Algoritmo 3"), each = nreps))
)

#VALORI
mse_alg3

mae_alg3

# Creare il grafico della densità MSE
ggplot(data_density, aes(x = MSE, fill = Algorithm)) +
  geom_density(alpha = 0.5) +
  labs(title = "Gráfico de Densidad de MSE para Algoritmos 3",
       x = "MSE",
       y = "Densidad") +
  theme_minimal()

# Creare il grafico della densità MAE
ggplot(data_density, aes(x = MAE, fill = Algorithm)) +
  geom_density(alpha = 0.5) +
  labs(title = "Gráfico de Densidad de MSE para Algoritmos 3",
       x = "MAE",
       y = "Densidad") +
  theme_minimal()

#___________trovare gli orari________
# Funzione per generare un DAG fisso casuale
generate_fixed_DAG <- function(p, en) {
  set.seed(123)
  g <- randomDAG(p, prob = en / p, lB = 1, uB = 2)
  return(g)
}

# Funzione per applicare l'algoritmo 1 (semplice PC)
apply_algorithm1 <- function(X, alpha) {
  suffStat <- list(C = cor(X), n = nrow(X))
  pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(X), alpha = alpha)
  return(pc.fit)
}

# Funzione per applicare l'algoritmo 3 (PC locale)
apply_algorithm3 <- function(X, alpha) {
  suffStat <- list(C = cor(X), n = nrow(X))
  pc.fit <- pc(suffStat, indepTest = gaussCItest, p = ncol(X), alpha = alpha, u2pd = "relaxed")
  return(pc.fit)
}

# Parametri per la generazione di DAG
p <- 7
en <- 3
alpha <- 0.01  # Valore alfa utilizzato

# Genera un DAG fisso con pesi dei bordi
dag <- generate_fixed_DAG(p, en)

# Genera 50 set di dati di dimensione 1000
set.seed(123)
data_list <- replicate(50, rmvDAG(1000, dag, errDist = "normal"), simplify = FALSE)

# Applica gli algoritmi a ciascun set di dati e misura il tempo di esecuzione
results <- data.frame()
for (p in c(4, 9, 14, 29, 49, 99)) {
  dag_p <- generate_fixed_DAG(p, en)
  data_list_p <- replicate(50, rmvDAG(1000, dag_p, errDist = "normal"), simplify = FALSE)
  
  times <- microbenchmark(
    Algorithm1 = lapply(data_list_p, function(data) apply_algorithm1(data, alpha)),
    Algorithm3 = lapply(data_list_p, function(data) apply_algorithm3(data, alpha)),
    times = 10
  )
  
  times_df <- as.data.frame(times)
  times_df$p <- p
  results <- rbind(results, times_df)
}

# Costruisci il grafico della densità di runtime
ggplot(results, aes(x = time / 1e9, fill = expr)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ p, scales = "free") +
  labs(title = "Densidad de tiempos de ejecución para diferentes valores de p",
       x = "Tiempo de ejecución (segundos)",
       y = "Densidad") +
  theme_minimal()

# Mostra la tabella dei tempi di esecuzione
library(dplyr)
execution_times <- results %>%
  group_by(p, expr) %>%
  summarise(mean_time = mean(time / 1e9))

print(execution_times)

