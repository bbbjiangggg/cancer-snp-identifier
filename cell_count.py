#!/usr/bin/env python3


plate = input('\033[1;45m Enter the plate you wish to use (96-wells, 24-wells, 6-wells): \033[0;0;0m')

#number of wells
wells = input('\033[1;45m Enter the number of wells from the ' + plate + ' you wish to use: \033[0;0;0m ')
print('\n')

inwellvol = input('\033[1;45m Enter the intended volume (ml) for each well of the ' + plate + ' (number only): \033[0;0;0m ')
if inwellvol == int:
    inwellvol = float(inwellvol)
else:
    print('\033[1;31m Please enter a number only \033[0;0;0m') 
print('\n')

#well wellvol calculation

wellvol = float(inwellvol) * float(wells) + 1

print('\033[1;45m You will now enter the total number of live cell count in the sample. \033[0;0;0m ')
print('\033[1;45m You will enter the number using the following format: \033[0;0;0m 7.67e6  ')
print('\033[1;45m This means 7.67 million cells \033[0;0;0m ')
input('\033[1;45m Press enter to continue \033[0;0;0m ')
print('\n')

#number of live cells in the sample
cells = float(input('\033[1;45m Enter the number of live cell count: \033[0;0;0m '))
print('\n')
density = float(input('\033[1;45m Enter the seeding density you wish to use in your experiment (e.g. 1.0e5): \033[0;0;0m '))
print('\n')


final = density * wellvol/cells
final2 = round(final, 2)
micro = final2 * 1000
new = wellvol - final2

if final > 1:
    print('\033[1;45m ERROR:\033[0;0;0m You do not have enough cells to seed the wells you have chosen. Please choose a lower seeding density or a lower number of wells.')
    exit()
else:
    print('\033[1;45m The number of CLC you will need are: \033[0;0;0m ' + str(final2) + ' ml')

print('\n')

print('\033[1;45m Therefore, mix: \033[0;0;0m ')
print(str(micro) + ' \u03BCl of the live cell count \033[0;0;0m ')
print('with ' + str(new) + ' ml of new medium ')
print('to get a seeding density of ' + str(density) + ' cells/ml')


