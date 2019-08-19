import unittest
import numpy as np
import deeptools.heatmapper as heatmapper
import computeMatrixStranded
import deeptools.computeMatrix
from deeptools import computeMatrix




class TestMethods(unittest.TestCase):

    def test_computeMatrixStranded_singleGroup(self):
        return(1)
        args_stranded = ["reference-point", "-Rp", "test_groupP_plus.bed", "-Rm", "test_groupP_minus.bed", "-Sp", "bw/norm_xPap_Mex67AA_ip_0_1_plus_KevinRoyAfiltered_BGsub.bw", "-Sm", "bw/norm_xPap_Mex67AA_ip_0_1_minus_KevinRoyAfiltered_BGsub.bw", "-b", "100", "-a", "100", \
               "--outFileName", "tmp/stranded_test.mat.gz", "--quiet"]

        args = computeMatrixStranded.process_args(args_stranded)
        hm = computeMatrixStranded.compute_matrix(args)
        self.assertEqual(hm.matrix.group_labels, hm.parameters['group_labels'])
        self.assertEqual(hm.matrix.group_boundaries, hm.parameters['group_boundaries'])
        self.assertEqual(hm.matrix.group_labels, ['genes'])
        self.assertEqual(hm.matrix.group_boundaries, [0, 3957])
        self.assertEqual(20, hm.matrix.matrix.shape[1])
        self.assertEqual(3957, hm.matrix.matrix.shape[0])


    def test_computeMatrixStranded_multiGroup(self):
        return(1)
        args_stranded = ["reference-point", "-Rp", "test_groupP_plus.bed", "test_groupnonP_plus.bed", "-Rm", "test_groupP_minus.bed", "test_groupnonP_minus.bed", "-Sp", "bw/norm_xPap_Mex67AA_ip_0_1_plus_KevinRoyAfiltered_BGsub.bw", "-Sm", "bw/norm_xPap_Mex67AA_ip_0_1_minus_KevinRoyAfiltered_BGsub.bw", "-b", "100", "-a", "100", \
               "--outFileName", "tmp/stranded_test.mat.gz", "--quiet"]

        args = computeMatrixStranded.process_args(args_stranded)
        hm = computeMatrixStranded.compute_matrix(args)
        self.assertEqual(hm.matrix.group_labels, hm.parameters['group_labels'])
        self.assertEqual(hm.matrix.group_boundaries, hm.parameters['group_boundaries'])
        self.assertEqual(hm.matrix.group_labels, ['test_groupP', 'test_groupnonP'])
        self.assertEqual(hm.matrix.group_boundaries, [0, 3957, 6685])
        self.assertEqual(20, hm.matrix.matrix.shape[1])
        self.assertEqual(6685, hm.matrix.matrix.shape[0])

    def test_computeMatrixStranded2(self):
        args_stranded = ["reference-point", "-Rp", "bed/PROMPTs_TSStoMaxIn5kb_hg38_withend_plus.bed", "-Rm", "bed/PROMPTs_TSStoMaxIn5kb_hg38_withend_minus.bed", "-Sp", "bw/siGFP_noPAP_in_plus.bw", "-Sm", "bw/siGFP_noPAP_in_minus.bw", "-b", "100", "-a", "100", \
               "--outFileName", "tmp/stranded_test.mat.gz", "--quiet"]

        args = computeMatrixStranded.process_args(args_stranded)
        hm = computeMatrixStranded.compute_matrix(args)
        self.assertEqual(hm.matrix.group_labels, hm.parameters['group_labels'])
        self.assertEqual(hm.matrix.group_boundaries, hm.parameters['group_boundaries'])
        self.assertEqual(hm.matrix.group_labels, ['genes'])
        self.assertEqual(hm.matrix.group_boundaries, [0, 5247])
        self.assertEqual(hm.matrix.sample_boundaries[-1], hm.matrix.matrix.shape[1])
        self.assertEqual(hm.matrix.group_boundaries[-1], hm.matrix.matrix.shape[0])

if __name__ == '__main__':
    unittest.main()