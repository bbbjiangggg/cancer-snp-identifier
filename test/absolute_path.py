import os

jar_file = '~/Trimmomatic-0.39/trimmomatic-0.39.jar'
jar_path = os.path.expanduser(jar_file)
print(jar_path)

if os.path.exists(jar_path):
    trim_path = jar_path
    replace_in_untrimmed_bash_srr('trim_path', trim_path)
    print(f"{jar_path} is the absolute path.")
else:
    print(f"NOTE: {jar_path} is not an absolute path.")
    trim_path = input('Copy and paste the absolute path to your trimmomatic-0.39.jar file: ')
    replace_in_untrimmed_bash_srr('trim_path', trim_path)


