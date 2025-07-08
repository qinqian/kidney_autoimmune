import spatialdata
import sopa
import spatialdata
import sopa
import spatialdata
import sopa
import sys

arg = sys.argv

zarr = arg[1]
label = arg[2]

sdata = spatialdata.read_zarr(zarr)
###sopa.segmentation.tissue(sdata)

## run cellpose
### sopa.make_image_patches(sdata, patch_width=1500, patch_overlap=50)
### print(sopa.utils.get_channel_names(sdata))
### sopa.segmentation.cellpose(sdata, channels=["DAPI"], diameter=50)

#sopa.settings.parallelization_backend = "dask"
#sopa.settings.dask_client_kwargs["n_workers"] = 2

# run baysor
sopa.make_transcript_patches(sdata, patch_width=1000, prior_shapes_key="cellpose_boundaries")
sopa.segmentation.baysor(sdata, min_area=10)
sopa.aggregate(sdata)

sopa.io.explorer.write(f"{label}.explorer", sdata)
sopa.io.write_report(f"{label}_report.html", sdata)
sdata.write(f"{label}")
