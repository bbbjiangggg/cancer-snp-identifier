import re
from colorama import Fore, Style, init

# Initialize colorama
init(autoreset=True)

# Function to get a number from the user, supporting scientific notation
def get_numeric_input(prompt):
    input_str = input(Fore.CYAN + prompt + Fore.RESET + ' (e.g., 2000000, 2e7, or 1.63x10^6 cells/ml): ')
    
    # Attempt to convert the input to a float, handling scientific notation
    try:
        if 'x10^' in input_str:
            # Split the input into base and exponent
            base, exponent = input_str.split('x10^')
            return float(base) * 10 ** float(exponent)
        else:
            return float(input_str)
    except ValueError:
        print(Fore.RED + "Invalid input format. Please enter the number (e.g., 2000000, 2e7, or 1.63x10^6)." + Fore.RESET)
        return get_numeric_input(prompt)

# Function to convert milliliters to microliters
def ml_to_ul(volume_ml):
    return volume_ml * 1000

# Brief instructions for the user
print(Fore.GREEN + "Please enter the cell count and concentration as numbers (e.g., 2000000, 2e7, or 1.63x10^6 cells/ml)." + Fore.RESET)

# Get the live cell count from the user
live_count = get_numeric_input(Fore.YELLOW + 'Enter the live count' + Fore.RESET)

# Get the desired cell concentration from the user
desired_conc = get_numeric_input(Fore.YELLOW + 'Enter the desired cell concentration' + Fore.RESET)

while True:
    # Get the desired volume from the user and remove the 'ml' suffix
    volume_input = input(Fore.CYAN + 'Enter your desired volume (e.g., 0.5): ' + Fore.RESET)
    try:
        volume_ml = float(volume_input)
        break  # Exit the loop if the input is valid
    except ValueError:
        print(Fore.RED + "Invalid input format. Please enter the volume as a number (e.g., 0.5)." + Fore.RESET)

# Calculate the amount of live cells needed in milliliters and microliters
live_cells_needed_ml = desired_conc * volume_ml / live_count
live_cells_needed_ul = ml_to_ul(live_cells_needed_ml)

# Calculate the volume of PBS or complete medium needed in milliliters and microliters
medium_vol_ml = volume_ml - live_cells_needed_ml
medium_vol_ul = ml_to_ul(medium_vol_ml)

# Display the final mixing instructions
print(Fore.GREEN + f'For a volume of {volume_ml} ml, mix {live_cells_needed_ml:.3f} ml ({live_cells_needed_ul:.2f} µl) of your cells with')
print(f'{medium_vol_ml:.3f} ml ({medium_vol_ul:.2f} µl) of PBS for injection or of complete medium for well-plates.' + Fore.RESET)
print('\nMake sure you pipette up and down to mix the cells!')

