import unittest
import numpy as np
import deeptools.heatmapper as heatmapper
import computeMatrixStranded
import deeptools.computeMatrix
from deeptools import computeMatrix



class TestMethods(unittest.TestCase):

    def test_computeMatrixStranded_singleGroup(self):
        args_stranded = ["reference-point", "-Rp", "test_groupP_plus.bed", "-Rm", "test_groupP_minus.bed", "-Sp", "bw/norm_xPap_Mex67AA_ip_0_1_plus_KevinRoyAfiltered_BGsub.bw", "-Sm", "bw/norm_xPap_Mex67AA_ip_0_1_minus_KevinRoyAfiltered_BGsub.bw", "-b", "100", "-a", "100", \
               "--outFileName", "tmp/stranded_test.mat.gz"]

        computeMatrixStranded.run(args_stranded)

        hm_stranded = heatmapper.heatmapper()
        hm_stranded.read_matrix_file('tmp/stranded_test.mat.gz')
        self.assertEqual(20, hm_stranded.matrix.matrix.shape[1])
        self.assertEqual(3957, hm_stranded.matrix.matrix.shape[0])


    def test_computeMatrixStranded_multiGroup(self):
        args_stranded = ["reference-point", "-Rp", "test_groupP_plus.bed", "test_groupnonP_plus.bed", "-Rm", "test_groupP_minus.bed", "test_groupnonP_minus.bed", "-Sp", "bw/norm_xPap_Mex67AA_ip_0_1_plus_KevinRoyAfiltered_BGsub.bw", "-Sm", "bw/norm_xPap_Mex67AA_ip_0_1_minus_KevinRoyAfiltered_BGsub.bw", "-b", "100", "-a", "100", \
               "--outFileName", "tmp/stranded_test.mat.gz"]

        computeMatrixStranded.run(args_stranded)

        hm_stranded = heatmapper.heatmapper()
        hm_stranded.read_matrix_file('tmp/stranded_test.mat.gz')
        self.assertEqual(20, hm_stranded.matrix.matrix.shape[1])
        self.assertEqual(6685, hm_stranded.matrix.matrix.shape[0])


if __name__ == '__main__':
    unittest.main()