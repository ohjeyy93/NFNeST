import subprocess
from os import listdir
import pandas as pd
import csv

subprocess.run(["nextflow", "run", "trim1testVCF10-11.nf", "-resume"])
subprocess.run(["python", "errortrack1.py"])