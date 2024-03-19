#!/data3/Group7/wangjiaxuan/biosoft/miniconda3/bin/python
import os
import sys
import subprocess
import argparse


def print_usage():
    print("Usage:")
    print("[-i]: The input tsv without header including two columns which mean fq dictroys path (column 1) and sample name (column 2)")
    print("[-r]: The refer database, select human or mouse, the default is human")
    print("[--snRNA]: the snRNA mean the sequence use single nuclear ,not cell suppsion.")
    print("[-h]: The help document")

    
def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', '--input', help='The input tsv without header including two columns which mean fq dictroys path (column 1) and sample name (column 2)')
    parser.add_argument('-r', '--refer', default='human', help='The refer database, select human or mouse, the default is human')
    parser.add_argument('--snRNA', action='store_true', help='the snRNA mean the sequence use single nuclear, not cell suppsion.')
    parser.add_argument('-h', '--help', action='store_true', help='The help document')

    args = parser.parse_args()

    if args.help:
        print_usage()
        exit()

    input_file = args.input
    if os.path.exists(input_file):
        input_file = os.path.abspath(input_file)
    else:
        sys.exit("the -i parmeter is file not exist")

    refer = args.refer

    if args.snRNA:
        sn = "true"
    else:
        sn = "false"

    # 检查参考文件的路径--------------  
    if refer == "human":
        refer = "/data3/Group7/wangjiaxuan/refer/cellranger/refdata-gex-GRCh38-2020-A"
    elif refer == "mouse":
        refer = "/data3/Group7/wangjiaxuan/refer/cellranger/refdata-gex-mm10-2020-A"
    elif os.path.exists(refer):
        refer = os.path.abspath(refer)
    else:
        sys.exit("The Genome Refer path"+ refer + "is wrong")

    wkdir = "./cellranger_count_out"
    fq_file = input_file

    if not os.path.exists(wkdir):
        os.mkdir(wkdir)

    os.chdir(wkdir)

    with open(fq_file, 'r') as f:
        for line in f:
            col1, col2 = line.strip().split()

            if not os.path.exists(col2):
                os.mkdir(col2)
            os.chdir(col2)

            result = subprocess.run(['which', 'qsub'], capture_output=True, text=True)
            exit_code = result.returncode
            if exit_code == 0:
                with open(f'qsub_{col2}_cellranger_count.sh', 'w') as script_file:
                    script_file.write(f"#!/bin/bash\n")                    
                    script_file.write(f"/data3/Group7/wangjiaxuan/biosoft/cellranger-7.0.0/cellranger count --id={col2}_result --fastqs={col1} --transcriptome={refer} --include-introns {sn}\n")
                
                subprocess.run(['/data3/Group7/wangjiaxuan/biosoft/miniconda3/bin/python', '/data3/Group7/wangjiaxuan/script/qsub.py','-s','1','-g','30g','-c','8','-l','4','-r' ,f'qsub_{col2}_cellranger_count.sh'])
            
            result = subprocess.run(['which', 'sbatch'], capture_output=True, text=True)
            exit_code = result.returncode
            if exit_code == 0:
                with open(f'sbacth_{col2}_cellranger_count.sh', 'w') as script_file:
                    script_file.write(f"#!/bin/bash\n")
                    script_file.write(f"#SBATCH -J {col2}_cellranger\n")
                    script_file.write(f"#SBATCH -o {col2}_cellranger.out\n")
                    script_file.write(f"#SBATCH -e {col2}_cellranger.err\n")
                    script_file.write(f"#SBATCH -N 1\n")
                    script_file.write(f"#SBATCH --ntasks-per-node=1\n")
                    script_file.write(f"#SBATCH --cpus-per-task=16\n")
                    script_file.write(f"#SBATCH --mem=32G\n")
                    script_file.write(f"/data3/Group7/wangjiaxuan/biosoft/cellranger-7.0.0/cellranger count --id={col2}_result --fastqs={col1} --transcriptome={refer} --include-introns {sn}\n")

                subprocess.run(['sbatch', f'sbacth_{col2}_cellranger_count.sh'])

            os.chdir("../")
if __name__ == "__main__":
    main()
