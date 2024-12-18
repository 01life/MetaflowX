//Parse parameter configuration into a channel.
def params2Channel(parameter){
    if(parameter){
        fparam = file(parameter)
        if(fparam.exists()){
            ch_parameter = Channel.value(fparam)
        }else{
            exit 1, "Unable to access '"+ fparam + "': No such file or directory"
        }
    }else{
        ch_parameter = Channel.empty()
    }
    return ch_parameter
}

//Check necessary parameters for module execution.
def checkEssentialParams(list){
    def flag = true
    for (param in list) {
        if (!param) {
            flag = false
            break;
        }
    }
    return flag
}

// Minitools check .fa files.
def checkFaFiles(String folder) {  
    def faFiles = file("${folder}/*.fa").size()  
    if (faFiles == 0) {  
        throw new Exception("Error: The specified bin_genome_folder path (${folder}) does not contain any .fa files.")  
    } else {  
        println "FA files found in ${folder}."  
    }  
}  