import requests

# define the URL of The Human Protein Atlas API
hpa_api_url = 'https://www.proteinatlas.org'

# ask the user for the gene name to check
gene_name = input("Please enter the name of the gene to check: ")

# build the URL for The Human Protein Atlas API request to check if the gene is cancer-related
hpa_api_request_url = f"{hpa_api_url}/api/search?search={gene_name}&format=json"

# make The Human Protein Atlas API request
response = requests.get(hpa_api_request_url, headers={'Content-Type': 'application/json'})

# check if the API request was successful
if response.ok:
    # extract the response data as a JSON object
    data = response.json()

    # check if the gene is cancer-related
    if 'rna_tissue_consensus' in data and data['rna_tissue_consensus']['blood'] == 'not detected':
        # if it is, print a message indicating that the gene is cancer-related
        print(f"{gene_name} is a cancer-related gene.")
    else:
        # if it is not, print a message indicating that the gene is not cancer-related
        print(f"{gene_name} is not a cancer-related gene.")
else:
    # if the API request failed, print an error message with the HTTP status code
    print(f"Error: API request failed with status code {response.status_code}")

