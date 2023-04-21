#!/data/wangjiaxuan/biosoft/miniconda3/bin/python
import argparse
import os
import subprocess

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
    parser.add_argument('--snRNA', action='store_true',help=' the snRNA mean the sequence use single nuclear ,not cell suppsion.')
    parser.add_argument('-h', '--help', action='store_true', help='The help document')

    args = parser.parse_args()

    if args.help:
        print_usage()
        exit()

    input_file = args.input
    refer = args.refer
    snRNA = args.snRNA

    print(input_file)

    if refer == "human":
        refer = "/data3/wangjiaxuan/refer/cellranger/refdata-gex-GRCh38-2020-A"
    elif refer == "mouse":
        refer = "/data3/wangjiaxuan/refer/cellranger/refdata-gex-mm10-2020-A"

    wkdir = "./cellranger_count_out"
    fq_file = input_file

    if not os.path.exists(wkdir):
        os.mkdir(wkdir)

    os.chdir(wkdir)

    with open("../" + fq_file, 'r') as f:
        for line in f:
            col1, col2 = line.strip().split()

            if not os.path.exists(col2):
                os.mkdir(col2)
            os.chdir(col2)
            
            if  not snRNA:
                with open(f'qsub_{col2}_cellraner_count.sh', 'w') as script_file:
                    script_file.write(f"/data/wangjiaxuan/biosoft/cellranger-7.0.0/cellranger count --id={col2}_result --fastqs={col1} --sample={col2} --transcriptome={refer}")
            else:
                with open(f'qsub_{col2}_cellraner_count.sh', 'w') as script_file:
                    script_file.write(f"/data/wangjiaxuan/biosoft/cellranger-7.0.0/cellranger count --id={col2}_result --fastqs={col1} --sa'mple={col2} --transcriptome={refer}")
            subprocess.call(['/data/wangjiaxuan/biosoft/miniconda3/envs/meta/bin/python', '/data/wangjiaxuan/script/qsub.py',"-s","2", "-g", "1g" ,"-c", "2" ,"-l", "1" ,"-r", f'qsub_{col2}_cellraner_count.sh'])

if __name__ == "__main__":
    main()