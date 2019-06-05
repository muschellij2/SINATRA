#### Aligning and Scaling the meshes ####
library(Morpho)
library(rgl)
library(Rvcg)
dir='Data/doug_new_teeth_raw/'
out_path = 'Data/doug_new_teeth_scaled/'

clean_files=function(input_dir,output_dir){
    files=list.files(path = dir,full.names = TRUE)
    filenames=list.files(path=dir,full.names = FALSE)
    num_files=length(files)
    for (i in 1:num_files){
        print(files[i])
        file=vcgImport(files[i])
        temp_file_name = filenames[i]
        new_file_name = tolower(gsub(pattern = "\\.|\\?|\\!|\\-", replacement = '_',x = temp_file_name))
        new_file_name = tolower(gsub(pattern = "_off", replacement = '.off',x = new_file_name))
        filename=paste(out_path,new_file_name,sep='')
        area=vcgArea(file)
        file2=scalemesh(file,1/sqrt(area),center='mean')
        centroid <- colMeans(vert2points(file2))
        file3=translate3d(file2,-centroid[1],-centroid[2],-centroid[3])
        print(filename)
        vcgOffWrite(file3,filename = filename)
    }
}
#clean_files(dir,out_path)



