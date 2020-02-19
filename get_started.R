source("read_QX.R")

# read files from folder
d.data <- read.QX(directory = "example/triplex_1")

# plot a single reaction
plot(sample(d.data$Ch1[,8]))
source("cloudy_modified.R")
cloudy(d.data$Ch1[,8])

test=na.omit(d.data$Ch1[,1])
density(test)
plot(density(test))
# Example of batch analysis
apply(d.data$Ch1, 2, cloudy)

