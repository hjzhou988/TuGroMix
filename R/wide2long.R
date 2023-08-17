#' Converts tumor growth data from "wide" format to "long" format
#'
#' Most tumor growth data are in "wide" format, in which each row represents a mouse, and each column represents a time point.
#' This format cannot be read by lmmFit and other model fitting functions, and has to be converted to long format.
#' @param tv_data Tumor volume data frame in "wide" format, that is, each row is a mouse and each column is a time point. In addition, need to have a Group column, indicating each mouse's group information.
#' @param id.vars The "ID" columns that describes mouse information, such as Mouse ID, Group information, and everything except Day columns. Can use a vector of strings, or indices of the columns in the data table.
#'
#' @export
#' @return A long data formatted data frame

wide2long = function(tv_data,id.vars){
  data.long = reshape2::melt(tv_data, id.vars = id.vars, variable.name = "Day", value.name = "TV", na.rm = T)
  data.long$Day = base::sub("X","",data.long$Day)
  data.long$Day = base::as.numeric(data.long$Day)
  data.long
}


