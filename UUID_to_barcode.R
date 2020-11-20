library(TCGAutils)

df <- read.table('UUID_TNBC.txt', sep = '\t', header = T)
id <- df$GDC.Case.UUID %>% as.data.frame()
names(id) <- 'UUID'

barcode <- UUIDtoBarcode(id$UUID, from_type = "case_id")

barcode <- merge(df, barcode, by = 1)

write.csv(barcode,'TNBC_barcode_uuid.csv')
