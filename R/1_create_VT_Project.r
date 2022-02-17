
# path <- "file:///C:/Users/titi7/Desktop/test_VT_project"
#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
create_VT_Project <- function(path){

  path <- gsub("file:/*", "", path)

  active_VT_project <<- path

  if(!dir.exists(path)){

    dir.create(path)


    dir.create(file.path(path, "1_user_inputs"))
    dir.create(file.path(path, "2_celltheques"))
    dir.create(file.path(path, "2_celltheques", "celltheques"))
    dir.create(file.path(path, "2_celltheques", "celltheque_one_per_profil"))
    dir.create(file.path(path, "3_virtual_tumors"))

    # Config File

    text <- readLines(file.path(  find.package(package = "VirtualTumor"), "config.R"))


fileConn<-file(file.path(path, "1_user_inputs", "1_config.r"))
  writeLines(text, fileConn)
  close(fileConn)

  # verif File

  text <- readLines(file.path(  find.package(package = "VirtualTumor"), "verification.R"))


  fileConn<-file(file.path(path, "1_user_inputs", "2_verif.r"))
  writeLines(text, fileConn)
  close(fileConn)

  # use File

  text <- readLines(file.path(  find.package(package = "VirtualTumor"), "user_file.R"))

  text <- gsub("#path_project", active_VT_project, text)
  text <- gsub("#path_file_config", paste0("shell.exec(\"", file.path(path, "1_user_inputs", "1_config.r"), "\")"), text)
  text <- gsub("#path_file_verif", paste0("shell.exec(\"", file.path(path, "1_user_inputs", "2_verif.r"), "\")"), text)

  fileConn<-file(file.path(path,"user_file.R"))
  writeLines(text, fileConn)
  close(fileConn)



  shell.exec(file.path(active_VT_project, "1_user_inputs", "1_config.R"))
  shell.exec(file.path(active_VT_project, "1_user_inputs", "2_verif.R"))
  shell.exec(file.path(active_VT_project, "user_file.R"))

  }else{

    source(file.path(active_VT_project, "1_user_inputs", "1_config.R"))

  }




}


#' VT_source
#'
#' @return
#' @export
#'
#' @examples
VT_source <- function(){

  source(file.path(active_VT_project, "1_user_inputs", "1_config.r"))

}

# create_VT_Project(path)

# file.create("b")
