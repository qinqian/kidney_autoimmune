import spatialdata
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', help='zarr data path')
    parser.add_argument('-p', help='output prefix')

    args = parser.parse_args()

    sdata = spatialdata.read_zarr(args.d)
    print(sdata.tables)

    sdata.tables['table'].write(f"{args.p}_orig_seg.h5ad")


if __name__ == '__main__':
    main()

