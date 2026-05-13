# %%
import cns
import cns.data_utils as cdu

# %%
# Load CNS data a display first 5 rows
samples_df, raw_df = cdu.main_load("raw", "TRACERx")
cns.fig_heatmap(cns.cns_head(raw_df, 5), max_cn=6)

# %%
# Add missing segments, display first 5 rows
imp_df = cns.main_impute(raw_df, print_info=True)
cns.fig_heatmap(cns.cns_head(imp_df, 5), max_cn=6)

# %%
# Create 3 mb segments, convert to a 3D feature array
segs = cns.main_segment(split_size=3_000_000)
seg_df = cns.main_aggregate(imp_df, segs, print_info=True)
features, rows, columns = cns.bins_to_features(seg_df)
print("Samples: {0}, Alleles: {1}, Bins: {2}.".format(*features.shape))

# %%
# Group segments by cancer type, sum the CNs and create mean linear profile
type_groups = {c: cns.select_cns_by_type(seg_df, samples_df, c, "type") for c in ["LUAD", "LUSC"]}
groups_df = cns.stack_groups([cns.group_samples(v, group_name=k) for k, v in type_groups.items()])
cns.fig_lines(cns.add_total_cn(groups_df), cn_columns="total_cn")

