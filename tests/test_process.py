import unittest
import numpy as np
import pandas as pd
import io

from cns.process import *
from cns.utils import hg19, hg38, segments_to_cns_df, tuples_to_segments
from cns.pipelines import main_segment

class TestSegments(unittest.TestCase):    
    def setUp(self):
        self.assembly = type('Assembly', (object,), {
            'aut_names': ['chr1', 'chr2', 'chr3'],
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300,
            'sex_names': ['chrX', 'chrY'],
            'chr_names': ['chr1', 'chr2', 'chr3'],
            'chr_x': 'chrX',
            'chr_y': 'chrY'
        })
        self.samples_df = pd.DataFrame({
            'sex': ['xx', 'xy']
        }, index=['s1', 's2'])
    
    
    def test_do_segments_overlap(self):
        segs_a = {1: [(0, 5), (4, 8)], 2: [(10, 15)]}
        segs_b = {1: [(5, 10)], 2: [(14, 20)], 3: [(25, 30)]}
        gaps_hg19_segs = tuples_to_segments(hg19.gaps)
        gaps_hg38_segs = tuples_to_segments(hg38.gaps)
        self.assertTrue(do_segments_overlap(segs_a))
        self.assertFalse(do_segments_overlap(segs_b))
        self.assertFalse(do_segments_overlap(gaps_hg19_segs))
        self.assertFalse(do_segments_overlap(gaps_hg38_segs))

    def test_merge_segments(self):
        segs = {1: [(5, 10), (10, 15)], 2: [(20, 25), (25, 30)], 3: [(35, 40)]}
        exp = {1: [(5, 15)], 2: [(20, 30)], 3: [(35, 40)]}
        self.assertEqual(merge_segments(segs), exp)

    def test_segment_union(self):
        segs_a = {1: [(0, 10)], 2: [(10, 15)]}
        segs_b = {1: [(5, 10)], 2: [(15, 20)], 3: [(25, 30)]}
        exp = {1: [(0, 10)], 2: [(10, 20)], 3: [(25, 30)]}
        self.assertEqual(segment_union(segs_a, segs_b), exp)

    def test_find_overlaps(self):
        segs = {1: [(1, 3), (7, 9), (3, 4), (8, 10)], 2: [(10, 15), (12, 20)], 3: [(20, 25), (22, 30)]}
        exp = {1: [(8, 9)], 2: [(12, 15)], 3: [(22, 25)]}
        self.assertEqual(find_overlaps(segs), exp)

    def test_segment_difference(self):
        segs_a = {1: [(0, 10)], 2: [(15, 25)], 3: [(20, 30)]}
        segs_b = {1: [(3, 5), (7, 8)], 2: [(20, 23)], 3: [(22, 25), (26, 29)]}
        exp = {
            1: [(0, 3), (5, 7), (8, 10)],
            2: [(15, 20), (23, 25)],
            3: [(20, 22), (25, 26), (29, 30)],
        }
        self.assertEqual(segment_difference(segs_a, segs_b), exp)

        segs_a = {1: [(0, 10)]}
        segs_b = {1: [(9, 10)]}
        res = segment_difference(segs_a, segs_b)
        self.assertEqual(res, {1: [(0, 9)]})

    def test_filter_min_size(self):
        segs = {1: [(0, 10)], 2: [(15, 20)], 3: [(20, 30)]}
        min_size = 6
        exp = {1: [(0, 10)], 2: [], 3: [(20, 30)]}
        self.assertEqual(filter_cons_size(segs, min_size), exp)

    def test_split_segment(self):
        actual_output = split_segment(1, 11, None, 2, "scale")
        expected_output = [(1, 3), (3, 5), (5, 7), (7, 9), (9, 11)]
        self.assertEqual(actual_output, expected_output)

        expected_output = [(1, 5), (5, 8), (8, 11)]
        actual_output = split_segment(1, 11, None, 3, "pad")
        self.assertEqual(actual_output, expected_output)

    def test_regions_select(self):        
        res = make_segments("whole")
        self.assertEqual(len(res), 24)
        self.assertEqual(len(res["chr1"]), 1)
        self.assertEqual(res["chr1"][0][0], 0)
        self.assertEqual(res["chr1"][0][2], "chr1")

        res = make_segments("arms")
        self.assertEqual(len(res), 24)
        self.assertEqual(len(res["chr1"]), 2)
        self.assertEqual(res["chr1"][0][2], "chr1p")

        res = make_segments("bands")
        self.assertEqual(len(res), 24)
        self.assertEqual(res["chr1"][0][2], "p36.33")
        self.assertEqual(res["chr10"][-1][1], hg19.chr_lens["chr10"])

        select = make_segments("whole")
        remove = make_segments("gaps")
        self.assertGreater(len(remove), 0)
        segs = main_segment(select, remove, filter_size=0)
        self.assertGreater(len(segs["chr1"]), 0)
        self.assertEqual(remove["chr1"][0][1], segs["chr1"][0][0])  # check if the first segment is a gap

    def test_arms_gaps(self):
        select = make_segments("whole")
        remove = make_segments("gaps")
        segs = main_segment(select, remove, filter_size=1000000)
        segs_df = segments_to_cns_df(segs)
        self.assertEqual(segs_df.query("chrom == 'chr1'").shape[0], 2)
        self.assertEqual(segs_df.query("chrom == 'chr13'").shape[0], 1)

    def test_cent_regions(self):
        regions = make_segments("centromeres")
        self.assertEqual(len(regions), 24)
        self.assertEqual(regions["chr1"][0][0], 121500000)

    def test_get_genome_segments(self):
        select = {1: [(0, 10), (20, 30)], 2: [(0, 5)]}
        remove = {1: [(5, 15)]}

        filter_size = 1
        expected_result = {1: [(15, 20), (20, 30)]}
        result = main_segment(select, remove, filter_size=filter_size)

        filter_size = 6
        expected_result = {1: [(20, 30)], 2: []}
        result = main_segment(select, remove, filter_size=filter_size)

        self.assertEqual(result, expected_result)

    def test_gene_segs(self):
        dummy_file = """gene	chrom	start	end
SKI	chr1	2160134	2241558
TNFRSF14	chr1	2487078	2496821
BIRC6	chr2	32582096	32843966
STRN	chr2	37070783	37193615
EML4	chr2	42396490	42559688
SRGAP3	chr3	9022275	9404737
BCORL1	chrX	129115083	129192058"""
        gene_df = pd.read_csv(io.StringIO(dummy_file), sep="\t")
        gene_segs = {}
        for i, row in gene_df.iterrows():
            if row["chrom"] not in gene_segs:
                gene_segs[row["chrom"]] = []
            gene_segs[row["chrom"]].append((row["start"], row["end"], row["gene"]))
        other_segs = {"chr1": [(0, 100000), (2000000, 2200000)], "chrY": [(0, 5)]}
        res = merge_segments(gene_segs)
        self.assertEqual(gene_segs, res)
        res = segment_union(gene_segs, other_segs)
        self.assertEqual(gene_segs["chrX"], res["chrX"])
        self.assertEqual(other_segs["chrY"], res["chrY"])
        res = segment_difference(gene_segs, other_segs)
        self.assertEqual(res["chr1"][0][0], 2200000) # cut by the other segment
        res = split_segments(gene_segs, 100000)
        self.assertGreater(len(res["chr2"]), len(gene_segs["chr2"]))

    def test_align_segs_to_assembly(self):
        segs = {
            "chr1": [(0, 50, "gene1"), (50, 100, "gene2")],
            "chr2": [(50, 150, "gene3")],
            "chr3": []
        }
        aligned = align_segs_to_assembly(segs, assembly=self.assembly)
        self.assertEqual(len(aligned["chr1"]), 2)
        self.assertEqual(len(aligned["chr2"]), 3)
        self.assertEqual(len(aligned["chr3"]), 1)
        self.assertEqual(aligned["chr1"][0][0], 0)
        self.assertEqual(aligned["chr1"][-1][1], 100)
        self.assertEqual(aligned["chr2"][0][0], 0)
        self.assertEqual(aligned["chr2"][-1][1], 200)
        self.assertEqual(aligned["chr3"][0][0], 0)
        self.assertEqual(aligned["chr3"][-1][1], 300)
        self.assertEqual(len(segs["chr2"]), 1)

# TODO: Add sex chromosome checks
class TestImputation(unittest.TestCase):
    def setUp(self):
        self.cns_df = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's2', 's2', 's2'],
            'chrom': ['chr1', 'chr2', 'chr2', 'chr2', 'chr2', 'chr2', 'chr3'],
            'start': [0, 0, 50, 125, 150, 175, 0],
            'end': [100, 150, 100, 150, 175, 200, 350],
            'major_cn': [1, 2, 3, np.nan, 1, 1, 2],
            'minor_cn': [1, 2, 1, 0, 0, 0, 1]
        }) 
        self.assembly = type('Assembly', (object,), {
            'aut_names': ['chr1', 'chr2', 'chr3'],
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300,
            'sex_names': ['chrX', 'chrY'],
            'chr_names': ['chr1', 'chr2', 'chr3'],
            'chr_x': 'chrX',
            'chr_y': 'chrY'
        })
        self.samples_df = pd.DataFrame({
            'sex': ['xx', 'xy']
        }, index=['s1', 's2'])

    def test_add_tails(self):
        result = add_tails(self.cns_df, self.assembly, print_info=True)
        self.assertEqual(result.shape[0], 9)
        self.assertEqual(result.at[3, "start"], 0)
        self.assertEqual(result.at[3, "end"], 50)

    def test_fill_gaps(self):
        result = fill_gaps(self.cns_df, print_info=True)
        self.assertEqual(result.shape[0], 8)
        self.assertEqual(result.at[3, "start"], 100)
        self.assertEqual(result.at[3, "end"], 125)

    def test_add_missing(self):
        result = add_missing(self.cns_df, self.samples_df, self.assembly, print_info=False)
        self.assertEqual(result.shape[0], 9)
        self.assertEqual(result.at[3, "start"], 0)
        self.assertEqual(result.at[3, "end"], 100)

    def test_merge_neighbours(self):
        result = merge_cns_df(self.cns_df, print_info=False)
        self.assertEqual(result.shape[0], 6)        
        self.assertEqual(result.at[4, "start"], 150)
        self.assertEqual(result.at[4, "end"], 200)

    def test_fill_nans_with_zeros(self):
        result = fill_nans_with_zeros(self.cns_df, print_info=False)
        self.assertEqual(result.major_cn.isnull().sum(), 0)
        self.assertEqual(result.minor_cn.isnull().sum(), 0)

    def test_impute_extend(self):
        result = add_tails(self.cns_df, self.assembly)
        result = fill_gaps(result, print_info=False)    
        result = add_missing(result, self.samples_df, self.assembly, print_info=False)
        result = cns_infer(result, self.samples_df, print_info=False)
        result = merge_cns_df(result, print_info=False)
        result = fill_nans_with_zeros(result, print_info=False)  
        self.assertEqual(result.at[4, "end"], 112)
        self.assertEqual(result.at[5, "end"], 125)
        self.assertEqual(result.at[6, "end"], 200)
        self.assertEqual(result.major_cn.isnull().sum(), 0)
        self.assertEqual(result.query("sample_id == 's1'")["chrom"].unique().shape[0], 3)
        self.assertEqual(result.query("sample_id == 's2'")["chrom"].unique().shape[0], 3)

    def test_impute_diploid(self):
        result = add_tails(self.cns_df, self.assembly)
        result = fill_gaps(result, print_info=False)    
        result = add_missing(result, self.samples_df, self.assembly, print_info=False)
        result = cns_infer(result, self.samples_df, method='diploid', print_info=False)
        result = merge_cns_df(result, print_info=False)
        result = fill_nans_with_zeros(result, print_info=False)  
        self.assertEqual(result.major_cn.isnull().sum(), 0)
        self.assertEqual(result.shape[0], 10)
        self.assertEqual(result.at[3, "minor_cn"], 1)

    def test_infer_segs(self):
        cns_df = pd.DataFrame({
            'sample_id': ['s1', 's1', 's1', 's1'],
            'chrom': ['chr1', 'chr1', 'chr1', 'chr1'],
            'start': [0, 100, 200, 300],
            'end': [100, 150, 300, 400],
            'major_cn': [1, np.nan, 3, 2],
            'minor_cn': [1, 2, 1, 0]
        }) 
        res_df = cns_infer(cns_df, self.samples_df)
        self.assertEqual(res_df.at[1, "major_cn"], 1.0) # should ignore the following 3


class TestBreakpoints(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_arm_breaks(self):
        result = make_breaks("arms")
        self.assertEqual(list(result.keys())[0], 'chr1')
        self.assertEqual(list(result.values())[0], [0, 125000000, 249250621])
        
    def test_cytoband_breaks(self):
        result = make_breaks("cytobands")
        self.assertEqual(list(result.keys())[0], 'chr1')
        sum_of_breaks = sum([len(breaks) - 1 for breaks in result.values()])
        sum_of_bands = len(hg19.cytobands)
        self.assertEqual(sum_of_breaks, sum_of_bands)

    def test_bin_breaks(self):
        result = make_breaks(1000000, "scale")
        self.assertEqual(list(result.keys())[0], 'chr1')
        self.assertEqual(list(result.values())[0][0], 0)
        self.assertEqual(list(result.values())[0][-1], 249250621)

    def test_dist_breaks_scaled(self):
        act = split_into_bins(10, 1, "scale")
        exp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 2, "scale")
        exp = [0, 2, 4, 6, 8, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 2.5, "scale")
        exp = [0, 3, 5, 8, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 3, "scale")
        exp = [0, 3, 7, 10]
        self.assertEqual(act, exp)        
        act = split_into_bins(13, 10, "scale")
        exp = [0, 13]
        self.assertEqual(act, exp)
        act = split_into_bins(17, 10, "scale")
        exp = [0, 9, 17]
        self.assertEqual(act, exp)
        act = split_into_bins(33, 10, "scale")
        exp = [0, 11, 22, 33]
        self.assertEqual(act, exp)
        act = split_into_bins(37, 10, "scale")
        exp = [0, 9, 19, 28, 37]
        self.assertEqual(act, exp)
        
    def test_dist_breaks_padded(self):
        act = split_into_bins(10, 1, "pad")
        exp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 2, "pad")
        exp = [0, 2, 4, 6, 8, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 2.5, "pad")
        exp = [0, 3, 5, 8, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 3, "pad")
        exp = [0, 4, 7, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(13, 10, "pad")
        exp = [0, 13]
        self.assertEqual(act, exp)
        act = split_into_bins(17, 10, "pad")
        exp = [0, 9, 17]
        self.assertEqual(act, exp)
        act = split_into_bins(33, 10, "pad")
        exp = [0, 12, 22, 33]
        self.assertEqual(act, exp)
        act = split_into_bins(37, 10, "pad")
        exp = [0, 9, 19, 29, 37]
        self.assertEqual(act, exp)

    def test_dist_breaks_after(self):
        act = split_into_bins(10, 1, "after")
        exp = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 2, "after")
        exp = [0, 2, 4, 6, 8, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 2.5, "after")
        exp = [0, 3, 5, 8, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(10, 3, "after")
        exp = [0, 3, 6, 10]
        self.assertEqual(act, exp)
        act = split_into_bins(13, 10, "after")
        exp = [0, 13]
        self.assertEqual(act, exp)
        act = split_into_bins(17, 10, "after")
        exp = [0, 10, 17]
        self.assertEqual(act, exp)
        act = split_into_bins(33, 10, "after")
        exp = [0, 10, 20, 33]
        self.assertEqual(act, exp)
        act = split_into_bins(37, 10, "after")
        exp = [0, 10, 20, 30, 37]
        self.assertEqual(act, exp)

    def test_diffs(self):
        seg_breaks = make_breaks(10_000_000)
        for chrom, chrom_breaks in seg_breaks.items():
            diffs = np.diff(np.diff(chrom_breaks))
            self.assertTrue(np.abs(np.sum(diffs)) <= 1)
            self.assertTrue(np.max(np.abs(diffs) <= 1))


class TestAggregation(unittest.TestCase):
    def setUp(self):
        self.cns = pd.DataFrame({
            'sample_id': ['s1', 's1', 's2', 's2', 's2', 's3', 's4', 's4', 's4', 's4', 's4', 's4'],
            'chrom': ['chr1', 'chr1', 'chr2', 'chr2', 'chrY', 'chr3', 'chr1', 'chr1', 'chr1', 'chr2', 'chr2', 'chr2'],
            'start': [0, 50, 100, 200, 300, 400, 0, 50, 99, 50, 100, 120],
            'end': [50, 100, 150, 300, 400, 500, 50, 99, 100, 100, 120, 130],
            'major_cn': [1, 2, 1, 3, 4, 5, 2, 1, 0, 2, 1, 1],
            'minor_cn': [0, 2, np.nan, 0, 4, 3, 1, 0, 0, 1, 0, 1],
        })       
        self.samples = pd.DataFrame({
            'sample_id': ['s1', 's2', 's3', 's4'],
            'sex': ['xx', 'xy', 'xx', 'xy']
        }).set_index('sample_id')
        self.assembly = type('Assembly', (object,), {
            'aut_names': ['chr1', 'chr2', 'chr3'],
            'chr_lens':{'chr1': 100, 'chr2': 200, 'chr3': 300, 'chrX': 100, 'chrY': 100},
            'cum_starts': {'chr1': 0, 'chr2': 100, 'chr3': 300, 'chrX': 600, 'chrY': 700},
            'aut_len': 300,
            'sex_names': ['chrX', 'chrY']
        })

    def test_agg_by_breaks(self):
        segments = { 'chr1': [(0, 100, "chr1_0")], 'chr2': [(100, 200, "chr2_0")] }
        breaks = {'chr1': [0, 100], 'chr2': [100, 200]}
        seg_bin = aggregate_by_segments(self.cns, segments, print_info=False)
        break_bin = aggregate_by_breaks(self.cns, breaks, print_info=False)
        self.assertEqual(seg_bin.shape[0], 4)
        pd.testing.assert_frame_equal(seg_bin, break_bin)
    
    def test_agg_by_segments(self):
        segments = {'chr1': [(0, 100, 0)], 'chr2': [(100, 200, 1)]}
        res = aggregate_by_segments(self.cns, segments)
        self.assertEqual(res.shape[0], 4)
        self.assertEqual(res.at[0, "start"], 0)
        self.assertEqual(res.at[0, "end"], 100)
        self.assertEqual(res.at[0, "major_cn"], 1.5)
        self.assertEqual(res.at[0, "minor_cn"], 1.0)
        self.assertEqual(res.at[1, "start"], 100)
    
    def test_agg_none(self):        
        segments = {'chr1': [(0, 100, 0)], 'chr2': [(100, 200, 1)]}
        res = aggregate_by_segments(self.cns, segments, how="none")
        self.assertEqual(res.shape[0], 8)
        for i in range(res.shape[0]):
            if res.at[i, "chrom"] == "chr1":
                self.assertTrue(res.at[i, "start"] >= 0)
                self.assertTrue(res.at[i, "end"] <= 100)
            elif res.at[i, "chrom"] == "chr2":
                self.assertTrue(res.at[i, "start"] >= 100)
                self.assertTrue(res.at[i, "end"] <= 200)

    def test_agg_round(self):
        segments = {'chr1': [(0, 100, 0)], 'chr2': [(100, 200, 1)]}
        mean_res = aggregate_by_segments(self.cns, segments, how="mean", print_info=False)
        round_res = aggregate_by_segments(self.cns, segments, how="round", print_info=False)
        self.assertEqual(round_res.shape, mean_res.shape)
        # rounding must match the rounded length-weighted mean for every value
        for col in ["major_cn", "minor_cn"]:
            for i in range(round_res.shape[0]):
                mean_val = mean_res.at[i, col]
                round_val = round_res.at[i, col]
                if np.isnan(mean_val):
                    self.assertTrue(np.isnan(round_val))
                else:
                    self.assertEqual(round_val, np.round(mean_val))
        # spot-check a known value: s1 chr1 major_cn weighted mean 1.5 -> 2.0
        self.assertEqual(round_res.at[0, "major_cn"], 2.0)


class TestClustering(unittest.TestCase):
    def test_cluster_breaks(self):
        breaks = {'chr1': [0, 50, 149, 200, 299, 300], 'chr2': [100, 200,300] }
        segs = breaks_to_segments(breaks)
        dist = 100
        res = cluster_segments(segs, dist, True)
        self.assertEqual(len(res), 2)
        self.assertEqual(res['chr1'][0][0], 0)
        self.assertEqual(res['chr1'][2][0], 250)
        self.assertEqual(res['chr2'][0][0], 100)

    def test_cluster_segs(self):
        breakpoints = {"clust": [0, 5, 7, 10, 14, 20]}
        segments = breaks_to_segments(breakpoints)
        clust1 = cluster_segments(segments, 5, True)["clust"]
        # Check bounds preservation
        self.assertEqual(clust1[0][0], 0)
        self.assertEqual(clust1[-1][1], 20)
        clust2 = cluster_segments(segments, 5, False)["clust"]
        self.assertEqual(len(clust2), 3)
        self.assertEqual(clust2[0], (2, 8, 'clust_0'))
        self.assertEqual(clust2[2], (14, 20, 'clust_2'))

    def test_pad_segments(self):
        from types import SimpleNamespace
        from cns.process.segments import pad_segments
        # Mock assembly with chromosome lengths
        assembly = SimpleNamespace(chr_lens={'chr1': 1000, 'chr2': 2000})
        # Input segments: no names
        segs = {
            'chr1': [(100, 200), (300, 400)],
            'chr2': [(0, 100), (1900, 2000)]
        }
        pad_size = 50
        padded = pad_segments(segs, pad_size, assembly=assembly)
        self.assertEqual(padded['chr1'][0], (50, 250))
        self.assertEqual(padded['chr1'][1], (250, 450))
        self.assertEqual(padded['chr2'][0], (0, 150))
        self.assertEqual(padded['chr2'][1], (1850, 2000))

        merged = merge_segments(padded)
        self.assertEqual(merged['chr1'], [(50, 450)])

        # With names
        assembly2 = SimpleNamespace(chr_lens={'chr1': 500})
        segs2 = {'chr1': [(10, 20, 'segA'), (480, 500, 'segB')]
        }
        pad_size2 = 15
        padded2 = pad_segments(segs2, pad_size2, assembly=assembly2)
        self.assertEqual(padded2['chr1'][0], (0, 35, 'segA'))
        self.assertEqual(padded2['chr1'][1], (465, 500, 'segB'))
        # Boundary case
        assembly3 = SimpleNamespace(chr_lens={'chrX': 100})
        segs3 = {'chrX': [(0, 10), (90, 100)]}
        pad_size3 = 20
        padded3 = pad_segments(segs3, pad_size3, assembly=assembly3)
        self.assertEqual(padded3['chrX'][0], (0, 30))
        self.assertEqual(padded3['chrX'][1], (70, 100))


if __name__ == "__main__":
    unittest.main()