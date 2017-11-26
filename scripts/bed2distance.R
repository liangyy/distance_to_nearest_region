library(optparse)

option_list <- list(
    make_option(c("-i", "--intersection"), type="character", default=NULL,
                help="Prefiltered intersection file",
                metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output file name",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

library(dplyr)
inter <- tryCatch({
  read.table(opt$intersection, header = F, sep = '\t')
}, error <- function(err){
  gz <- gzfile(opt$output, "w")
  write.table(data.frame(), gz, col.names = T, row.names = F, quote = F, sep = '\t')
  close(gz)
  quit()
})
inter <- read.table(opt$intersection, header = F, sep = '\t')
inter$junction.pos <- (inter$V5 + inter$V6 - 1) / 2
inter$distance <- abs(inter$V2 - inter$junction.pos)
inter$variant.id <- paste(inter$V1, inter$V2, inter$V3, sep = '\t')
result <- inter %>%
  group_by(variant.id) %>%
  summarize(nearest.distance = min(distance), valid = sum(junction.pos) > 0)

result[!result$valid, 'nearest.distance'] <- NA
gz <- gzfile(opt$output, "w")
write.table(result, gz, col.names = T, row.names = F, quote = F, sep = '\t')
close(gz)
