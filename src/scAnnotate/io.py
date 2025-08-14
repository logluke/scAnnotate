import os
import sys
import urllib.request
import pathlib as path

def download_gff(gff_path: path.Path | str):
    url = f"https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/{gff_path.name}"
    try:
        urllib.request.urlretrieve(url, gff_path)
        print(f"Download complete: {gff_path}")
    except Exception as e:
        print(f"Download failed: {e}")
        if os.path.exists(gff_path):
            os.remove(gff_path)
        sys.exit(1)

