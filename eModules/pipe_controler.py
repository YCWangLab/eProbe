import subprocess
import os
import logging

# check if a command is executable
def executable(compiler, command):
    try:
        subprocess.run(f"{compiler} {command}", shell=True, check=False)
        return True
    except subprocess.CalledProcessError:
        return False


# find path of a script
def find_script_path(path_list, script):
    which_script = None
    for path in path_list:
        if executable("python", os.path.join(path, script)):
            which_script = path
        elif script in os.environ["PATH"]:
            which_script = "PATH"
    return which_script


def diff_file_line_count(file1_path, file2_path):
    try:
        count_file1 = int(subprocess.check_output(["wc", "-l", file1_path]).split()[0])
        count_file2 = int(subprocess.check_output(["wc", "-l", file2_path]).split()[0])

        return int(count_file1) - int(count_file2)
    except subprocess.CalledProcessError:
        print("Error in counting file line number！")
        return None


def delete_files(file_list):
    for file in file_list:
        try:
            os.remove(file)
        except OSError as e:
            print(f"Error deleting temp file '{file}': {e}")


def run_python_script(script_name, path_search_list, script_params):
    which_script = find_script_path(path_search_list, script_name)
    if which_script:
        if which_script != "PATH":
            script_path = os.path.join(which_script, script_name)
        else:
            script_path = script_name

        script_command = ' '.join(['python', script_path] + script_params)
        try:
            script_log = subprocess.Popen(script_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                                shell=True)
            stdout, stderr = script_log.communicate()
            logging.info(f"Running: {script_log}")
            logging.info(stdout.decode())
            logging.info(stderr.decode())
        except Exception as e:
            logging.error(f"Error executing command: {script_log}")
            logging.error(f"Exception: {e}")
    else:
        raise FileNotFoundError(f"{script_name} not found.")
