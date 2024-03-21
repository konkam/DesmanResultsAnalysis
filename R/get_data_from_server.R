#' Downloads data from server
#' @description download data from the server via scp, by default the file 
#' @param destdirorfile character string, path to a destination directory or file
#' @param dirorfilepathonwillow   character string, path to a  directory or file on willow 
#' @examples 
#' get_data_from_server(server="front.migale.inrae.fr",
#'  dirorfilepathonserver="/work_projet/ala/metachick-fugace/Analyses_Ouléye/DESMAN_2023/DESMAN_iterAddPost/DESMAN_erm_F__3_M17808/freq2Desman.txt")
get_data_from_server<-function(
    dirorfilepathonserver,
    destdirorfile=tempfile(fileext=paste0(".",tools::file_ext(dirorfilepathonserver))),  
    server=c("front.migale.inrae.fr"),
    cred=if(is.element("montruc",installed.packages()[,1])){montruc::inraetruc()
    }else{
      list(username=readline(prompt="Enter user name: "),password=readline(prompt="Password: "))}){
  if(file.exists(dirorfilepathonserver)){destdirorfile<-dirorfilepathonserver}else{
  system(paste0(if(!is.null(cred$password)){if(cred$password!=""){paste0("sshpass -p ",cred$password)}}else{}," scp -r ",cred$username,"@",server,":",dirorfilepathonserver," ",destdirorfile|>gsub(pattern = " ",replacement = "\\ ",fixed = T)|>gsub(pattern = "\\\\ ",replacement = "\\ ",fixed = T)))}
  destdirorfile}

#'@examples
#'"timestampedexports/Cube_20220629.csv"|>get_inst_files()
get_inst_files<-function(path,package="DesmanResultsAnalysis"){
  suppressWarnings(c(file.path("inst",path),                           
                     system.file(path,package=package))|>
                     Filter(f=file.exists)|>
                     (`[`)(1))}

load_data_file<-function(.data,package="DesmanResultsAnalysis"){
  if(file.exists(file.path("data",paste0(.data,".rda")))){
    file.path("data",paste0(.data,".rda"))|>load(envir = parent.frame())
  }else{
    data(list=.data,package=package,envir = parent.frame())
  }}

get_load_data_file<-function(.data,package="DesmanResultsAnalysis"){
  if(exists(file.path("data",paste0(.data,".rda")))){
    file.path("data",paste0(.data,".rda"))|>load()|>get()
  }else{
    data(list=.data,package=package)|>get()
  }}

#'@examples
#'"timestampedexports"|>get_inst_files()
get_inst_dirs<-function(path,packagefolder="inst",package="DesmanResultsAnalysis"){
  suppressWarnings(c(file.path("inst",path),                           
                     system.file(path,package=package))|>
                     Filter(f=dir.exists)|>
                     (`[`)(1))}

