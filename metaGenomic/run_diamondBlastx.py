__author__ = 'Roshan Padmanabhan'
__version__ = '0.1'

import os
import sys
from glob import glob
import shlex
import subprocess
import argparse


def diamond_blastx_cdl(fas, fas_aa, db):
    cmd = "diamond blastx --quiet -d " + db + " -q " + fas + " -a " +  fas_aa + " -p 8 "
    c =  shlex.split(cmd)
    return c

def diamond_view_cdl(fas_aa):
    fas_aax = fas_aa + '.daa'
    res = fas_aa + '.m8'
    cmd = "diamond view --quiet -o " + res + " -a " +  fas_aax 
    c =  shlex.split(cmd)
    return c

def run_diamond(args):
    self.args = args
    try :
        s = subprocess.call(args, shell=False)
    except :
        raise
    return s


if __name__ == '__main__':

    ######################          Command Line Argument Parsing       ##########################

    des="""
    \n
    This script runs diamond blsatx on server.\n
    
    Requirements: diamond
    thread is set to 8.
    Database : hard-coded into the script. modify. make it using 'diomand makedb' commad.\n
    Input : Path to the directory in which all the files are located and each fasta file is has _chunk_ in it \n
    DISCLAIMER : This script was made to run the shotgun metagenomic data set against a kegg database.
        Each big fasta file is split into chunks and each chunk was blastx against kegg protein database.\n
    """
    parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', help='path to the fasta reads',nargs='+', action='store',dest='file_path',required=True)
    args = parser.parse_args()
    
    ######################          DataBase Location       ##########################

    # location of diamond Database file without extension
    db_loc='/work/rpadmana/kegg/keggBact'

    ######################          The Variables           ##########################
    
    dir_in = args.file_path
    dir_path = os.path.realpath( dir_in )
    dir_bn = os.path.basename( os.path.realpath(dir_in )) + ".res"
    out_file = os.path.join( dir_path , dir_bn )

    ######################          Start Blastx            ##########################

    if os.path.isdir(dir_path) is False:
        print ("\nError : Dir don't exists\n")
        sys.exit()
    else :
        files = sorted(glob( os.path.join( dir_path, "*_chunk_*")))

    # this loop has to be modified depening on the delimiters
    # its lame ...
    counter = 0
    for i in files:
        ix = i
        i = os.path.basename(i)
        i_bn = i.split(".")[0]
        i_bn2 = i_bn + '_' + counter
        counter += 1
        i_aa = os.path.join( dir_path , i_bn2 )
        print "processing ", i_aa
        run_diamond(diamond_blastx_cdl( ix, i_aa, db_loc ))
        run_diamond(diamond_view_cdl( i_aa ))

    ######################          Results and Clean Up       ##########################

    # Concatenate all the .m8 files into a single result file and delete the rest
    m8s = sorted( glob(os.path.join( dir_path, "*.m8")) )
    for f in m8s:
         os.system("cat "+f+" >> " +  out_file)
    [ os.remove(i) for i in m8s ]
    print "Done...Check the result file : " , out_file
