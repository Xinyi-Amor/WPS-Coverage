# Read the data from the file
data <- read.delim("output_wps.txt", header = FALSE)

# Select columns 2 to 7
columns <- data[, 2:7]

# Subtract the minimum value of each column from all values in that column
normalized_columns <- apply(columns, 2, function(col) col - min(col))

# Combine the normalized columns with other columns, if needed
normalized_data <- cbind(data[, 1], normalized_columns)

# Write the normalized data to a new text file
write.table(normalized_data, file = "processed_wps.txt", sep = "\t", quote = FALSE, row.names = FALSE)
