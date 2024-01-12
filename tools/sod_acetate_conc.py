#!/usr/bin/env python3

def main():
    # Ask the user for the total volume of the plasmid in microliters
    plasmid_volume_ul = float(input("Enter the total volume of your plasmid solution in microliters: "))

    # Define the initial and final concentrations
    initial_concentration = 3  # 3M sodium acetate
    final_concentration = 0.3  # Desired concentration

    # Calculate the final volume (V2) in mL, as the input is in microliters
    final_volume_ml = plasmid_volume_ul / 1000

    # Using the dilution formula C1V1 = C2V2, calculate V1 (the volume of 3M sodium acetate to be added)
    # V1 = (C2 * V2) / C1
    volume_to_add_ml = (final_concentration * final_volume_ml) / initial_concentration

    # Convert the volume to add from mL to microliters for consistency
    volume_to_add_ul = volume_to_add_ml * 1000

    # Output the result
    print(f"Add {volume_to_add_ul:.2f} microliters of 3M Sodium Acetate.")

if __name__ == "__main__":
    main()

