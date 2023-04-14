# Get the live cell count
count1 = float(input('Enter the live count without the base or exponent: '))
power1 = int(input('Enter the exponent: '))
live_count = count1 * 10 ** power1

# Display the live cell count in scientific notation
print(f'Your live cell count is: {live_count:.2e}\n')

# Get the desired cell concentration
count2 = float(input('Enter the desired cell concentration without the base or exponent: '))
power2 = int(input('Enter the exponent: '))
desired_conc = count2 * 10 ** power2

# Display the desired cell concentration in scientific notation
print(f'Your desired cell concentration is: {desired_conc:.2e}\n')

# Get the desired volume and calculate the amount of live cells needed
volume = float(input('Enter your desired volume in ml: '))
live_cells_needed = desired_conc * volume / live_count

# Convert the live cells needed to microliters
live_cells_needed_ul = live_cells_needed * 1000

# Calculate the volume of PBS or complete medium needed
medium_vol = volume - live_cells_needed

# Convert the medium volume to microliters
medium_vol_ul = medium_vol * 1000

# Display the final mixing instructions
print(f'For a volume of {volume} ml, mix {live_cells_needed:.3f} ml ({live_cells_needed_ul:.2f} µl) of your cells with')
print(f'{medium_vol:.3f} ml ({medium_vol_ul:.2f} µl) of PBS for injection or of complete medium for well-plates.')
print('\nMake sure you pipette up and down to mix the cells!')
