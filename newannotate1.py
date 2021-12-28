import os
import copy
import logging 
import itertools
from collections import namedtuple
from collections import OrderedDict
from nest.parsers.vcfReader import Reader 
from nest.parsers.vcfwriter import Writer

class Merge:

    def __init__(self, tmp_dir, vcf_dict, ref_path):
         ## Edit: (10/10/19): Added fasta path as requirement for initialization of merge
         ## Reason: Reader was modified to require fasta path, see vcfReader for more info
         self.out_path = tmp_dir
         self.vcf_dict = vcf_dict
         self.ref_path = ref_path

    def merge(self, vcflist):
        print([1])