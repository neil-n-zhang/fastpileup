import sys,os,subprocess
import argparse
import numpy as np
import multiprocessing as mp
from count_utils import base_counter, pileup2count

def worker(region,out,bam_files,ref,parameters):
    cmd = 'samtools mpileup {p} -f {ref} -r {reg} -o {o} {b}'.format(p=parameters,ref=ref, reg=region, o=out, b=bam_files)
    #print (cmd)
    subprocess.check_output(cmd, shell=True)
    pileup2count(out, out.replace("mpileup","txt"))
    return

def pileup_region(args):
    chr,pos_start,pos_end=args.reg.replace("-",":").split(":")
    pos_start=int(pos_start)
    pos_end=int(pos_end)

    boundary=np.linspace(pos_start,pos_end,args.threads+1,dtype=int)

    blocks=[]
    for i in range(len(boundary)-1):
        if i<(len(boundary)-2):
            blocks.append((boundary[i],boundary[i+1]-1))
        else:
            blocks.append((boundary[i],pos_end))

    pool = mp.Pool(args.threads)    
    outfile_pileups = []    

    for start,end in blocks:        
        print (start, end)
        region = '{c}:{s}-{e}'.format(c=chr,s=start,e=end)
        out = '{sa}_{c}_{s}_{e}.mpileup'.format(sa=args.sample,c=chr,s=start,e=end)
        f = pool.apply_async(worker, [region,out,args.bam,args.ref,args.parameters])
        outfile_pileups.append(out)
    outfiles_txt=[outfile.replace("mpileup","txt") for outfile in outfile_pileups]
    pool.close()
    pool.join()

    output_file = open(args.sample+"_count.txt", "w")
    output_file.write("chr\tpos\tA\tC\tG\tT\tN\ta\tc\tg\tt\tn\tIN\tDEL\tindel\n")
    output_file.close()
    cmd = 'cat {i} >> {o}'.format(i=' '.join(outfiles_txt), o=args.sample+"_count.txt")
    subprocess.check_output(cmd, shell=True)
    for f in outfile_pileups:
        os.remove(f)
    for f in outfiles_txt:
        os.remove(f)



def worker_bed(bed_file,out,bam_files,ref,parameters):
    cmd = 'samtools mpileup {p} -f {ref} -l {bed} -o {o} {b}'.format(p=parameters, ref=ref, bed=bed_file, o=out, b=bam_files)
    #print (cmd)
    subprocess.check_output(cmd, shell=True)
    pileup2count(out, out.replace("mpileup","txt"))
    return


def split_file(samplename, input_file, num_subfiles):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    total_lines = len(lines)
    lines_per_file = total_lines // num_subfiles
    remainder = total_lines % num_subfiles

    start = 0
    for i in range(num_subfiles):
        lines_to_write = lines_per_file + 1 if i < remainder else lines_per_file
        end = start + lines_to_write
        subfile_lines = lines[start:end]
        output_file = f'{samplename}_{i + 1}.bed'

        with open(output_file, 'w') as f_out:
            f_out.writelines(subfile_lines)

        start = end

def pileup_bed(args):
    
    sub_bed_files=[ args.sample+"_"+str(i+1)+".bed" for i in range(args.threads)]
 
    pool = mp.Pool(args.threads)
    outfile_pileups = []    
    
    split_file(args.sample, args.bed, args.threads)

    for sub_bed in sub_bed_files:        
        print (sub_bed)
        out = '{sa}_{s}_out.mpileup'.format(sa=args.sample,s=sub_bed.replace(".bed",""))
        f = pool.apply_async(worker_bed, [sub_bed,out,args.bam,args.ref,args.parameters])
        outfile_pileups.append(out)
    outfiles_txt=[outfile.replace("mpileup","txt") for outfile in outfile_pileups]
    pool.close()
    pool.join()

    output_file = open(args.sample+"_count.txt", "w")
    output_file.write("chr\tpos\tA\tC\tG\tT\tN\ta\tc\tg\tt\tn\tIN\tDEL\tindel\n")
    output_file.close()
    cmd = 'cat {i} >> {o}'.format(i=' '.join(outfiles_txt), o=args.sample+"_count.txt")
    subprocess.check_output(cmd, shell=True)

    for f in outfile_pileups:
        os.remove(f)
    for f in outfiles_txt:
        os.remove(f)
    for f in sub_bed_files:
        os.remove(f)

if __name__ == '__main__':
    
    # create the top-level parser
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(required=True)

    # create the parser for the "foo" command
    parser_reg = subparsers.add_parser('region')
    parser_reg.add_argument('-sample')
    parser_reg.add_argument('-bam')
    parser_reg.add_argument('-reg')
    parser_reg.add_argument('-threads', type=int)
    parser_reg.add_argument('-parameters')
    parser_reg.add_argument('-ref')
    parser_reg.set_defaults(func=pileup_region)

    parser_bed = subparsers.add_parser('bed')
    parser_bed.add_argument('-sample')
    parser_bed.add_argument('-bam')
    parser_bed.add_argument('-bed')
    parser_bed.add_argument('-threads', type=int)
    parser_bed.add_argument('-parameters')
    parser_bed.add_argument('-ref')
    parser_bed.set_defaults(func=pileup_bed)

    args=parser.parse_args()
    args.func(args)

    # python3 fastpileup.py region -sample test -bam test1.bam -reg 1:11931-12090 -threads 4 -parameters "-a -B -q0 -Q10" -ref /hs37d5.fa
    # python3 fastpileup.py bed -sample test -bam test1.bam -bed test.bed -threads 4 -parameters "-a -B -q0 -Q10" -ref /hs37d5.fa