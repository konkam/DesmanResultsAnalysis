#'@examples
#'desman_input_file<-"/work_projet/ala/metachick-fugace/Analyses_Ouléye/DESMAN_2023/DESMAN_iterAddPost/DESMAN_erm_F__3_M17808/freq2Desman.txt"|>get_data_from_server()
read_desman_input_files<-
  function(desman_input_file){
    read.csv(desman_input_file)|>tidyr::pivot_longer(cols=!c("Position","X"),
                           names_to=c("sample","nucleotide"),
                           names_pattern = "([[:alnum:]]+_+[[:alnum:]]+)\\.([[:alnum:]]+)")|>
      plyr::dlply(~X,
                  function(d){
                    plyr::daply(d,
                                .drop_i=FALSE,
                                .drop_o=FALSE,
                    .variables=~sample+Position+nucleotide,
                    .fun=function(dd){dd$value})})}
