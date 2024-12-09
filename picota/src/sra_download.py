
from tqdm import tqdm
import urllib.request
import os
import shutil
import subprocess

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def run_sra_down(sra_acc, out_path, out_dir, path_of_fastq_dump, keep_sra_file, the_force):

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.exists(f"{out_dir}/{sra_acc}"):
        os.mkdir(f"{out_dir}/{sra_acc}")
    else:
        if the_force == False:
            
            print('\n\033[0;33mWARNING:\033[0m Folder Exist!, Do you want to download again? Delete the folder, or run with force option: --force\n')
            return


    the_url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{sra_acc}/{sra_acc}"
    args2 = f"{path_of_fastq_dump} --split-3 {out_path}/*.sra --outdir {out_path}"
    print(out_path)

    try:
        download_url(the_url, out_path + "/" + sra_acc + ".sra")
    except Exception as e:
        print('url download error: ', e)
        return False

    print('Start for Fastq Dump')
    my_process = subprocess.run(args2, shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("Error in Fastq Dump")
    print('End of Fastq Dump')
    if keep_sra_file == False:
        if os.path.exists(f"{out_path}/{sra_acc}.sra"):
            os.remove(f"{out_path}/{sra_acc}.sra")
            print('SRA Files deleted.')
    print('fastq files were created in:', out_path)
    return True


def parse_sra_acc_file(file_path):
    the_return_list = []
    the_xeno_f = open(file_path,'r')
    the_lines = the_xeno_f.readlines()
    for line in the_lines:
        if line.replace('\n','') != '':
            the_return_list.append(line.replace('\n',''))
    the_xeno_f.close()
    return the_return_list


def sra_download_main(sra_acc_file, out_dir, sra_folder, path_of_fastq_dump, keep_sra_file, the_force):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    the_sra_list = parse_sra_acc_file(sra_acc_file)
    for index, sra_acc in enumerate(the_sra_list, start=0):
        print('Downloading:', sra_acc, ' : ', index ,'/', len(the_sra_list))
        out_dir_sra = os.path.join(out_dir, sra_folder)
        if not os.path.exists(out_dir_sra):
            os.mkdir(out_dir_sra)
        out_path = os.path.join(out_dir, sra_folder, sra_acc)
        run_sra_down(sra_acc, out_path, out_dir_sra, path_of_fastq_dump, keep_sra_file, the_force)

    return out_dir_sra

