import subprocess
import os
import logging
from subprocess import getstatusoutput
import warnings
import shutil


def executable(compiler, command):
    if getstatusoutput(f'{compiler} {command}')[0] == 0:
        return True
    else:
        return False


def find_script_path(script_name, path_search_list):
    searching_paths = path_search_list + os.environ['PATH'].split(os.pathsep)

    for path in searching_paths:
        script_path = os.path.join(path, script_name)
        if os.path.exists(script_path):
            if executable('python', f'{script_path} -h'):
                return script_path

    return None


def run_python_script(script_name, path_search_list, script_params):
    script_path = find_script_path(script_name, path_search_list)
    if script_path:
        script_command = ' '.join(['python', script_path] + script_params)
        try:
            logging.info(f'Running command: {script_command}')
            script_log = subprocess.Popen(script_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                          shell=True)
            stdout, stderr = script_log.communicate()
            logging.info(stdout.decode())
            logging.info(stderr.decode())
        except Exception as e:
            logging.error(f'Error executing command: {script_command}')
            logging.error(f'Exception: {e}')
            raise
    else:
        raise FileNotFoundError(f'Executable {script_name} not found.')


def check_file_exists(file_names, path='.', success='File found!', fail='File not found!'):
    if isinstance(file_names, str):
        file_names = [file_names]
    file_num = len(file_names)
    file_count = 0
    for file_name in file_names:
        full_path = os.path.join(path, file_name)
        if os.path.isfile(full_path):
            file_count = file_count + 1
        else:
            print(fail)
            raise FileNotFoundError(f'File {file_name} not found in path {path}')
    if file_count == file_num:
        print(success)


def check_bed(bed_path):
    try:
        with open(bed_path, 'r') as file:
            first_line = file.readline().strip()

            elements = first_line.split('\t')

            if len(elements) >= 3:
                return True
            else:
                return False
    except FileNotFoundError:
        raise FileNotFoundError(f'{bed_path} not found.')


def diff_file_line_count(file1_path, file2_path):
    try:
        count_file1 = int(subprocess.check_output(['wc', '-l', file1_path]).split()[0])
        count_file2 = int(subprocess.check_output(['wc', '-l', file2_path]).split()[0])

        return int(count_file1) - int(count_file2)
    except subprocess.CalledProcessError:
        print('Error in counting file line number！')
        return None


def delete_files(file_list):
    for file in file_list:
        try:
            os.remove(file)
        except FileNotFoundError:
            warnings.warn(f'File "{file}" does not exist. Skipping deletion.', Warning)


def get_unique_file_extensions(path):
    extensions = set()
    for root, dirs, files in os.walk(path):
        for file in files:
            _, ext = os.path.splitext(file)
            extensions.add(ext)
    return extensions


def check_col_name_in_df(lst, df):
    if set(lst).intersection(df.columns):
        return True
    return False


def make_temp_dir(tmp_dir):
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    else:
        shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir)
