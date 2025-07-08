import tangram as tg
import scanpy as sc

ad_sp = sc.read_h5ad("norm/RA_orig_seg_lognorm.h5ad")
ad_sc = sc.read_h5ad("reference/kpmp_sc_orig_seg_lognorm.h5ad")
tg.pp_adatas(ad_sc, ad_sp, genes=None)
ad_map = tg.map_cells_to_space(
               ad_sc,
               ad_sp,
               mode='clusters',
               cluster_label='celltype')

print(ad_map)
ad_map.write("ad_map.h5ad")
