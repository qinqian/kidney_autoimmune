library(VoltRon)

# import Xenium 
Xen_R1 <- importXenium("../data/kidney/20240803__182820__BWH_20240803_skin_Shruti_kidney/output-XETG00150__0018462__BS22_12012A1__20240803__183643/", sample_name = "XeniumR1")


# import H&E image and build a VoltRon object
Xen_R1_image <- importImageData("BS22_12012A1.TIF",
                                sample_name = "XeniumR1image", 
                                channel_names = "H&E")

test<-vrImages(Xen_R1_image, as.raster = TRUE)
dev.off()


image.list <- list(Xen_R1, Xen_R1_image)

## # This depend on a shiny browser
## options(shiny.launch.browser = F)
## xen_reg <- registerSpatialData(object_list = image.list)



