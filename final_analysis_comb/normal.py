
import requests
import io

# Create a dictionary of Google Drive links
link_list = {
    "1": "https://drive.google.com/file/d/1-09an3LHEJYnf-0s4S70ATis81dZPV1G/view?usp=sharing",
    "2": "https://drive.google.com/file/d/1xP8ELyW5Kr7fay0NulbgiPQOAx5Eq9E3/view?usp=sharing",
    "3": "https://drive.google.com/file/d/1-JNXwA4VvwcSOoulcmRmkLq8qYueltyQ/view?usp=sharing",
    "4": "https://drive.google.com/file/d/1n2T5xzMOqOobm5w97uthlvW5d7e3nV9M/view?usp=sharing",
    "5": "https://drive.google.com/file/d/1s8OPif3dKAesCgoOTIQGdAr9CdcUhe2q/view?usp=sharing",
    "6": "https://drive.google.com/file/d/1IzA-rUe8ELPyFtCePiCQgQeWDnOEPJxN/view?usp=sharing",
    "7": "https://drive.google.com/file/d/16eGfwwe1yjvyY6TO-hWj1Rx7YMvXuSbO/view?usp=sharing",
    "8": "https://drive.google.com/file/d/1uSqwMSvQpyBxPwiHhgUcNVkL7YdkcOHg/view?usp=sharing",
    "9": "https://drive.google.com/file/d/1chh47ez6Vqe1W0bdoKqxPHO8iOgMMiJe/view?usp=sharing",
    "10": "https://drive.google.com/file/d/12NXqBfLFJZwifp7_LLAzZFsfpE60czgN/view?usp=sharing",
    "11": "https://drive.google.com/file/d/1nsicoTeVg3t4AL592QGLa0QCiOmrhZz4/view?usp=sharing",
    "12": "https://drive.google.com/file/d/1hhX-NaqjzkpsA9cvFMj4F-1cEPCl6JmK/view?usp=sharing",
    "13": "https://drive.google.com/file/d/112aOwZiP4kOJhrBsBockjM0q2i1a1lx-/view?usp=sharing",
    "14": "https://",
    "15": "https://drive.google.com/file/d/1-umaSF-rAFqqSwVeQI6djyW7QKxBm5A7/view?usp=sharing",
    "16": "https://drive.google.com/file/d/1cOFSi6ujbJV_I144je6PnBUmiAZL1-J0/view?usp=share_link",
    "17": "https://drive.google.com/file/d/1yhKZrhhw2YlBMySWa4L-DYaJ4UVYH-y8/view?usp=share_link",
    "18": "https://drive.google.com/file/d/1yvj7_SQE93sJ2kwceFr5ikYsFstOZ36I/view?usp=sharing",
    "19": "https://drive.google.com/file/d/1hRKEmwtw7uBkx8n_5wBWvdY2WqQjuOZz/view?usp=sharing",
    "20": "http",
    "21": "http",
    "22": "https://drive.google.com/file/d/1zqOFpGV8HEVBo_uDIqVi0_YVQOCL6LZZ/view?usp=share_link",
    "X": "https",
    "Y": "https",
    # and so on
}

while True:
    normal_file = input('\033[1;45mWould you like to download a normal combined VCF file (yes or no)?\033[0;0;0m ')
    print('\n')
    if normal_file.lower() == 'yes':
        # Get chromosome number from user input
        chromosome = input('\033[1;45mEnter the chromosome number of interest (1-22, X or Y): \033[0;0;0m ')

        # Check if chromosome number is in link list
        if chromosome in link_list:
            # Extract the file ID from the Google Drive link
            file_id = link_list[chromosome].split("/")[-2]

            # Download the VCF file as bytes
            print("Downloading chromosome " + chromosome + " normal combined VCF file...")
            url = f"https://drive.google.com/uc?id={file_id}"
            response = requests.get(url)
            vcf_bytes = io.BytesIO(response.content).read()

            # Write the VCF bytes to a file
            filename = f"ch{chromosome}_norm_comb.vcf"
            with open(filename, "wb") as f:
                f.write(vcf_bytes)
            print(f"{filename} downloaded successfully!")
            break
        else:
            print("Chromosome number not found in link list.")
    elif normal_file.lower() == 'no':
        break
    else:
        print("Invalid input. Please enter 'yes' or 'no'.")

