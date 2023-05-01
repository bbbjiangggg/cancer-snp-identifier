def calculate_line(x_vals, y_vals):
    """
    Given x and y coordinate values, calculate the slope and y-intercept of the line that
    passes through these points.
    """
    n = len(x_vals)
    assert n == len(y_vals) == 5, "Expected 5 x and y coordinate values"

    # Calculate the slope
    delta_y = y_vals[-1] - y_vals[0]
    delta_x = x_vals[-1] - x_vals[0]
    slope = delta_y / delta_x

    # Calculate the y-intercept
    y_intercept = y_vals[0] - slope * x_vals[0]

    return slope, y_intercept


def calculate_x(slope, y_intercept, y_vals):
    """
    Given the slope, y-intercept, and a list of y values, calculate the corresponding
    x value for each y.
    """
    x_vals = []
    for y in y_vals:
        x = (y - y_intercept) / slope
        x_vals.append(x)
    return x_vals


# Ask the user for 5 x and y values
x_vals = []
y_vals = []
for i in range(5):
    x = float(input("Enter x{}: ".format(i+1)))
    y = float(input("Enter y{}: ".format(i+1)))
    x_vals.append(x)
    y_vals.append(y)

# Calculate the slope and y-intercept
slope, y_intercept = calculate_line(x_vals, y_vals)
print("Slope: {}".format(slope))
print("Y-Intercept: {}".format(y_intercept))

# Ask the user for 4 y values and calculate the corresponding x values
y_vals = []
for i in range(4):
    y = float(input("Enter y{}: ".format(i+1)))
    y_vals.append(y)

x_vals = calculate_x(slope, y_intercept, y_vals)
for i, x in enumerate(x_vals):
    print("X{}: {}".format(i+1, x))

