mkdir_if_not_exist <- function(fld_name){
  if(!dir.exists(fld_name))
    dir.create(fld_name)
}