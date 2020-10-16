GWAS_cor <- function(Assoc_list, table_name = "_", write_table = TRUE){
  #Takes a vector of association tests (given as the "short_test_title" associated with each test) to read in and a text to name the outputted correlation matrix
  #write_table is a logical, if TRUE, writes the output to a table, if false, returns the table as an object
  Assoc_df = data.frame(matrix(ncol = 0, nrow = 510610))
  for(Assoc_vect in Assoc_list){
    Assoc_df <- cbind(Assoc_df, read.table(paste("~/raw_data_analysis/AssociationTests_ManhattanPlots_qqPlots/",Assoc_vect,".txt", sep = ""))$Score.pval)
  }
  colnames(Assoc_df) <- Assoc_list
  Assoc_cor <<- cor(Assoc_df)
  if(write_table == TRUE){
    write.table(Assoc_cor, paste("~/raw_data_analysis/AssociationTests_ManhattanPlots_qqPlots/",table_name, "_corr.txt", sep = ""), sep = "\t")
  } else{
    return(Assoc_cor)
    }
}
