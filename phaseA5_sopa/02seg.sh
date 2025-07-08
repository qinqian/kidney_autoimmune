
zarr=../data/kidney/20241025__200743__BWH_20241025_SHRUTI_RACHEL/output-XETG00392__0045655__BS21-N65682A2__20241025__201009.zarr
sopa patchify image ${zarr} --patch-width-pixel 1200 --patch-overlap-pixel 50
sopa segmentation cellpose ${zarr} \
    --channels DAPI \
    --diameter 35 \
    --min-area 2000


#below not work yet
#sopa patchify transcripts ${zarr} --patch-width-microns 500 --prior-shapes-key cellpose_boundaries
###export SOPA_PARALLELIZATION_BACKEND=dask && sopa segmentation baysor ${zarr} --config '"config.toml"'
###SOPA_PARALLELIZATION_BACKEND=dask sopa segmentation baysor ${zarr} --config '"config.toml"'
#SOPA_PARALLELIZATION_BACKEND=dask sopa segmentation baysor ${zarr} --config '"config.toml"'

python run_sopa.py
