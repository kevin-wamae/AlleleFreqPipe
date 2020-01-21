#function for merging dataframes
#merge function - usage = new.df <- Reduce(MyMerge, list(df1, df2, df3, df4))
MyMerge <- function(x, y){
  df <- merge(x, y, by= "isolate", all.x= TRUE, all.y= TRUE)
  return(df)
}