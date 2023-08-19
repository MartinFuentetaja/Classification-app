import os
import subprocess
import sys

import Packages.argparse_decision as argparse

args = argparse.parser_argument()

print(args)

column_name    = ",".join(args.column_names)
dirname        = os.path.dirname(args.config_path)

MachineLearningPath = os.path.join(dirname, "Scripts/Machine_Learning/Machine_Learning.py")
DataSetCreationPath = os.path.join(dirname, "Scripts/Data_Set_Creation_3/Data_Set_Creation_3.py")

if os.path.exists(os.path.join(dirname, f"RESULTS/{args.gen_name}/{args.gen_name}_encoded.xlsx")):
   subprocess.run(["python", MachineLearningPath, "-g", args.gen_name, "-l", args.length, "-o", args.outfile, "-n_splits", args.number_splits, "-n_splits_grid", args.number_splits_grid, "-cpu", args.cpu_number, "-column_names", column_name, "-config_path", args.config_path])
else:
   subprocess.run(["python", DataSetCreationPath, "-g", args.gen_name, "-l", args.length, "-o", args.outfile, "-s", args.status, "-w", args.web_decision, "-cpu", args.cpu_number, "-config_path", args.config_path])
   subprocess.run(["python", MachineLearningPath, "-g", args.gen_name, "-l", args.length, "-o", args.outfile, "-n_splits", args.number_splits, "-n_splits_grid", args.number_splits_grid, "-column_names", column_name, "-cpu", args.cpu_number, "-config_path", args.config_path])
      
