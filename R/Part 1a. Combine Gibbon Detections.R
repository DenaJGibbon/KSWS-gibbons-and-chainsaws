# List all files in the target directory where high-confidence detections are stored
Files <- list.files('/Volumes/DJC Files/KSWS Gibbon Occupancy/Deop01/Detections/CrestedGibbons/above_0.96/Positive/')

# Extract the 7th element (usually the plot info) from each filename, using underscores as delimiter
Plot <- str_split_fixed(Files, pattern = '_', n = 8)[,7]

# Extract the 9th element (the date in YYYYMMDD format) from each filename
Date <- str_split_fixed(Files, pattern = '_', n = 10)[,9]

# Convert the extracted date strings to Date objects
Date <- as.Date(Date, format = "%Y%m%d")

# Extract the 10th element from each filename and grab the first two characters as the hour (HH)
Time <- substr(str_split_fixed(Files, pattern = '_', n = 10)[,10], 1, 2)

# Extract the second part of the Plot ID if it includes a dash (e.g., "KSWS-T19" → "T19")
Plot <- str_split_fixed(Plot, pattern = '-', n = 2)[,2]

# Combine extracted fields into a single data frame
CombinedDF <- cbind.data.frame(Plot, Date, Time)

# Add a constant column to label all detections with the species
CombinedDF$Common.Name <- 'CrestedGibbons'

# Redundant, but ensures Date is in Date format
CombinedDF$Date <- as.Date(CombinedDF$Date, format = "%Y-%m-%d")
CombinedDF$Date <- as.Date(as.character(CombinedDF$Date), format = "%Y-%m-%d")

# Optionally write to CSV if needed
# write.csv(CombinedDF, 'data/gibbonverifieddetections0.96.csv', row.names = FALSE)
