# load the required R scripts
source("ipss_m_evo_calculator.R")
source("utils.R")

# load the data
load("example_data.RData")

# data preprocessing
colnames(example_data) <- tolower(colnames(example_data))
colnames(example_data)[1] <- "ipss_m"
input_data <- build_inputs(example_data)

# compute IPSS-M-Evo
ipss_m_evo <- calculate_ipss_m_evo(input_data)
print(ipss_m_evo)
