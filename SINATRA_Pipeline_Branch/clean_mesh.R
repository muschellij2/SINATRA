#### Aligning and Scaling the meshes ####
library(Morpho)
library(rgl)
library(Rvcg)
dir='Data/all_files_tingran_clean/'
out_path = 'Data/all_files_tingran_scaled/'

clean_files=function(input_dir,output_dir){
    files=list.files(path = dir,full.names = TRUE)
    filenames=list.files(path=dir,full.names = FALSE)
    num_files=length(files)
    for (i in 1:num_files){
        print(files[i])
        file=vcgImport(files[i])
        filename=paste(out_path,filenames[i],sep='')
        area=vcgArea(file)
        file2=scalemesh(file,1/sqrt(area),center='mean')
        centroid <- colMeans(vert2points(file2))
        file3=translate3d(file2,-centroid[1],-centroid[2],-centroid[3])
        print(filename)
        vcgOffWrite(file3,filename = filename)
    }
}
#clean_files(dir)



