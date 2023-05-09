# Only works and necessary for linux distributions
import subprocess


data_dir = "../bin/"       # Output directory. If it is not there, it will be created.
url = 'https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar'


if __name__ == '__main__':
    subprocess.run("wget -P {} {}".format(data_dir, url), shell=True)