# utils.py
import os

def ensure_directory(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
        print(f"Directory {path} created.")
    else:
        print(f"Directory {path} already exists.")