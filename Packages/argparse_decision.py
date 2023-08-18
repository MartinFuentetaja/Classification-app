#####################################################################################################################################################################################################################
"""
Lo primero que vamos a hacer es importar los módulos necesarios:
"""

import argparse
import os
import json

#####################################################################################################################################################################################################################

def parser_argument():
    parser = argparse.ArgumentParser(description='Creation of database')
    database_group = parser.add_argument_group(title = 'Database group', description = 'Parser arguments for database creation')
    machinelearning_group = parser.add_argument_group(title = 'Machine Learning group', description = 'Parser arguments for Machine Learning analysis')    

    #Parser arguments for DataBase creation
    input_group = database_group.add_argument_group(title = 'Input group', description = 'There are two ways of introducing the input info: files or gen_name and sequence length')
    input_group.add_argument('-i', '--infile', dest='filename', nargs=8, help = 'Input files: <AlphaFold_file> <ClinVar_file> <FASTA_file> <Lovd3_file> <GenomAD_V2> <GenomAD_V2_control> <GenomAD_V2_non_neuro> <GenomAD_V3_non_V2>. The name of the <FASTA_file> must be strcutured as: <gen_name.fasta> ')
    input_group.add_argument('-g', '--genname', dest='gen_name', help="Name of the gene for constructing the database: <gen_name>")
    input_group.add_argument('-l', '--length', dest='length', help="Length of the amino acid sequence: <length>.")
    output_group = database_group.add_argument_group('Output group')
    output_group.add_argument('-o', '--outfile', dest='outfile', required = True, help='Output file: <outfile>')
    database_group.add_argument('-s', '--status', dest='status', required = False, default = 'Regular', help = "Status for secondary structure analysis. Regular or Extend. Regular will contain three different types of structures while Extend eight. (default = Regular) <Regular>")
    database_group.add_argument('-w', '--web', dest = 'web_decision', required = False, default = "Internet/DataBase", help = 'Status for file downloading: -w <Internet> or -w <Internet/DataBase>. Default <Internet/DataBase>.')
    database_group.add_argument('-cpu', dest = 'cpu_number', required = True, help = 'Number of CPUs for hhblits. Recommendation: 1-4')
    
    #Parser arguments for Machine Learning analysis
    machinelearning_group.add_argument('-n_splits', dest = 'number_splits', required = False, default = "5", help = "Number of splits for cross validation. <default = 5>")
    machinelearning_group.add_argument('-n_splits_grid', dest = 'number_splits_grid', required = False, default = "5", help = "Number of splits for grid search. <default = 5>")
    machinelearning_group.add_argument('-column_names', dest = 'column_names', required = False, default = 'My_Label,residue_conserv,d_size', help = 'Column names of the database for the Machine Learning. Default: <My_Label> <residue_conserv> <d_size>', type=lambda s: [str(item) for item in s.split(',')])  
    
    #Parser argument for init config.txt
    parser.add_argument('-config_path', required = True, help = "Path for config.txt")
    # Validar los argumentos
    args = parser.parse_args()
    
    if not os.path.isfile(args.config_path):
       parser.error("The path for config.txt is wrong")    


    if args.filename:
       if args.gen_name or args.length:
          parser.error('-i and -g/-l are mutually exclusive.')
       elif args.gen_name and args.length:
          if not args.gen_name or not args.length:
             parser.error('-g and -l must be provided together.')
       else:
          parser.error('Either -i or both -g and -l must be provided.')
	
    if args.filename:
       for filename in args.filename:
          if not os.path.isfile(filename):
             parser.error('Input file %s does not exist.' % filename)
    
    if args.length and not args.length.isdigit():
       parser.error('-l length is not a digit.')

    if args.status != "Regular" and args.status != "Extend":
       parser.error('-s is badly indicated. <Regular> or <Extend>')

    if not args.cpu_number.isdigit():
       parser.error('-cpu is not a number.')

    if "My_Label" not in args.column_names:
       parser.error('You must include "My_Label" column in -column_names')

    if not args.number_splits.isdigit():
       parser.error('-n_splits is not a number!')
    if not args.number_splits_grid.isdigit():
       parser.error('-n_splits_grid is not a number!')

    return args



