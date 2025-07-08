import scanpy as sc
import os
import anndata as ad
import spatialdata
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', help='zarr data path', nargs="+", action='store', type=str)
    parser.add_argument('-p', help='output prefix')

    args = parser.parse_args()
    print(args)

    adatas = []
    samples = []
    for a in args.d:
        samples.append(os.path.basename(a))
        sdata = sc.read(a)
        adatas.append(sdata)

    print(adatas)
    sdata = ad.concat(adatas,  merge="same", join='inner', label="sample", keys=samples)
    sdata.write(f"{args.p}_orig_seg_merged.h5ad")


if __name__ == '__main__':
    main()

