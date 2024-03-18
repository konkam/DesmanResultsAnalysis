#'@examples
#'"/work_projet/ala/metachick-fugace/Analyses_Ouléye/DESMAN_2023/DESMAN_iterAddPost/DESMAN_erm_F__3_M17808/freq2Desman.txt"|>
#'get_data_from_server()|>
#'read_desman_input_files()->X
read_desman_input_files<-
  function(desman_input_file){
    read.csv(desman_input_file)|>tidyr::pivot_longer(cols=!c("Position","X"),
                           names_to=c("s","a"),
                           names_pattern = "([[:alnum:]]+_+[[:alnum:]]+)\\.([[:alnum:]]+)")|>
      plyr::dlply(~X,
                  function(d){
                    d|>dplyr::rename(v="Position")|>
                    plyr::daply(.drop_i=FALSE,
                                .drop_o=TRUE,
                                .variables=~v+s+a,
                                .fun=function(dd){dd$value})})}
