#!/usr/bin/env python3

import os
import sys
import pip 

#print("\033[1;45m This text is cool  \033[0;0;0m")

#Check Python version
if sys.version_info.major == 3:
    print('\033[1;45m Python3 is installed. \033[0;0;0m')
else:
    print('\033[0;101m You need to install a current version of Python3 \033[0;0;0m')

#Check Pip version
pip_versn = pip.__version__
print('\033[1;45m Your Pip version is ' + pip_versn + '\033[0;0;0m')
print('\033[1;45m Note, if you do not have Pip installed, install it using the command: sudo apt install pip \033[0;0;0m')




