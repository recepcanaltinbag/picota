
from tqdm import tqdm
import urllib.request
import os
import shutil
import subprocess
import logging


logger: logging.Logger = None


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def run_sra_down(sra_acc, out_path, out_dir, path_of_fastq_dump, keep_sra_file, the_force, logger_name):
    global logger
    logger = logging.getLogger(logger_name)

    def _fastq_paths():
        return [os.path.join(out_path, f"{sra_acc}_{i}.fastq") for i in (1, 2)]

    def _valid_fastqs():
        return [f for f in _fastq_paths() if os.path.exists(f) and os.path.getsize(f) > 0]

    def _clear_partial():
        for f in _fastq_paths():
            try:
                if os.path.exists(f):
                    os.remove(f)
            except Exception:
                pass
        temp_dir = os.path.join(out_path, f"{sra_acc}.sra")
        if os.path.isdir(temp_dir):
            try:
                shutil.rmtree(temp_dir)
            except Exception:
                pass

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.exists(os.path.join(out_dir, sra_acc)):
        os.mkdir(os.path.join(out_dir, sra_acc))

    current_fastqs = _valid_fastqs()
    if len(current_fastqs) == 2 and not the_force:
        logger.info(f"FASTQ already present for {sra_acc}, skipping download")
        return True

    if os.path.exists(os.path.join(out_dir, sra_acc)) and not the_force and len(current_fastqs) < 2:
        logger.warning("Folder Exist but FASTQ incomplete; will re-download with force")
        _clear_partial()
        the_force = True

    the_url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{sra_acc}/{sra_acc}"
    file_sra = os.path.join(out_path, f"{sra_acc}.sra")

    if not os.path.isfile(file_sra) or the_force:
        try:
            download_url(the_url, file_sra)
        except Exception as e:
            logger.error('url download error: ', e)
            return False

    # run parallel-fastq-dump or fallback to fasterq-dump
    def _run_fastq_dump(cmd):
        try:
            p = subprocess.run(cmd, shell=True, text=True, capture_output=True)
            if p.returncode != 0:
                logger.warning(f"fastq-dump command failed (ret {p.returncode}): {p.stderr}")
                return False
            return True
        except Exception as e:
            logger.warning(f"fastq-dump command exception: {e}")
            return False

    logger.info('Start for Fastq Dump')

    # Try 1: parallel-fastq-dump with local .sra file
    pfd_cmd = f"{path_of_fastq_dump} -s {file_sra} -O {out_path} --split-files -t 12"
    success = _run_fastq_dump(pfd_cmd)

    # fasterq-dump writes large temp files; use out_path as temp dir to avoid
    # filling /tmp on the system disk (fasterq-dump -t flag controls temp location)
    _tmp_dir = out_path

    # Try 2: fasterq-dump with local .sra file, temp on output drive
    if not success:
        logger.warning('parallel-fastq-dump failed, trying fasterq-dump with .sra file')
        _clear_partial()
        if shutil.which('fasterq-dump'):
            fq_cmd = f"fasterq-dump --split-files -O {out_path} -t {_tmp_dir} {file_sra}"
            success = _run_fastq_dump(fq_cmd)
        else:
            logger.warning('fasterq-dump not available')

    # Try 3: fasterq-dump directly by accession (no local .sra needed — SRA toolkit streams)
    if not success:
        logger.warning('fasterq-dump with .sra failed, trying direct accession streaming')
        _clear_partial()
        if shutil.which('fasterq-dump'):
            fq_cmd = f"fasterq-dump --split-files -O {out_path} -t {_tmp_dir} {sra_acc}"
            success = _run_fastq_dump(fq_cmd)
        else:
            logger.warning('fasterq-dump not available for direct streaming')

    # Try 4: prefetch + fasterq-dump (re-download via NCBI toolkit, bypasses S3 URL issues)
    if not success:
        logger.warning('Direct streaming failed, trying prefetch + fasterq-dump')
        _clear_partial()
        if shutil.which('prefetch') and shutil.which('fasterq-dump'):
            prefetch_dir = os.path.join(out_path, sra_acc)
            prefetch_cmd = f"prefetch {sra_acc} -O {prefetch_dir}"
            if _run_fastq_dump(prefetch_cmd):
                fq_cmd = f"fasterq-dump --split-files -O {out_path} -t {_tmp_dir} {prefetch_dir}"
                success = _run_fastq_dump(fq_cmd)
        else:
            logger.warning('prefetch or fasterq-dump not available')

    if not success:
        raise Exception("Error in Fastq Dump")

    logger.info('End of Fastq Dump')

    current_fastqs = _valid_fastqs()
    if len(current_fastqs) != 2:
        logger.warning(f"FASTQ download incomplete for {sra_acc} after dump ({len(current_fastqs)}/2)")
        raise Exception("Incomplete FASTQ files")

    if keep_sra_file == False:
        if os.path.exists(file_sra):
            os.remove(file_sra)
            logger.info('SRA Files deleted.')

    logger.info('fastq files were created in: %s', out_path)
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


def sra_download_main(sra_acc_file, out_dir, sra_folder, path_of_fastq_dump, keep_sra_file, the_force, logger_name):
    
    global logger
    logger = logging.getLogger(logger_name)
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    the_sra_list = parse_sra_acc_file(sra_acc_file)
    for index, sra_acc in enumerate(the_sra_list, start=0):
        logger.info('Downloading:', sra_acc, ' : ', index ,'/', len(the_sra_list))
        out_dir_sra = os.path.join(out_dir, sra_folder)
        if not os.path.exists(out_dir_sra):
            os.mkdir(out_dir_sra)
        out_path = os.path.join(out_dir, sra_folder, sra_acc)
        run_sra_down(sra_acc, out_path, out_dir_sra, path_of_fastq_dump, keep_sra_file, the_force, logger_name)

    return out_dir_sra

