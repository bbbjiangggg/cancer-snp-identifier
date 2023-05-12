import numpy as np

# Accept user input for protein standard absorbance values
protein_std_absorbance = []
for i in [2, 1.0, 0.5, 0.25, 0.125]:
    absorbance = float(input(f"Enter protein absorbance value for {i} mg/ml standard: "))
    protein_std_absorbance.append(absorbance)

# Calculate slope and y-intercept
linear_fit = np.polyfit([2, 1.0, 0.5, 0.25, 0.125], protein_std_absorbance, 1)
slope = linear_fit[0]
y_intercept = linear_fit[1]

# Accept user input for unknown protein absorbance value
unknown_absorbance = float(input("Enter the absorbance value of the unknown protein: "))

# Calculate unknown protein concentration
unknown_concentration = (unknown_absorbance - y_intercept) / slope

print(f"The concentration of the unknown protein is {unknown_concentration} mg/ml.")

