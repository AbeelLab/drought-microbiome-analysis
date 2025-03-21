import os

def get_qiime_extract_dir(parent_dir):
    # extracted file is stored in a nested directory created by qiime extract
    # (qiime export would result in a cleaner directory structure,
    # but right now it doesn't work due to permission errors on the cluster)
    subdirs = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
    # Get latest extracted directory if there are more
    nested_id_dir = max(subdirs, key=lambda d: os.path.getctime(os.path.join(parent_dir, d)))
    nested_id_dir = os.path.join(parent_dir, nested_id_dir)
    data_dir = os.path.join(nested_id_dir, "data")

    return data_dir
