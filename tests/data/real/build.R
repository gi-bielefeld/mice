source("./soft_processing.R")
k = 31

for (file_list in list.files("/home/luca/_mice/data/Test/list", pattern="*.txt", full.names=TRUE)) {
    name = noext(basename(file_list))
    fastas = readLines(file_list)
    output_folder = "/home/luca/_mice/mice/tests/data/real/gfa"
    run_bifrost_build(k, fastas, file.path(output_folder, name))
    gfa_file = file.path(output_folder, paste0(name,".gfa"))
    run_print_paths(k, gfa_file, fastas, append=T)
}
