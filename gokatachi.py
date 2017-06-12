#!/usr/bin/env python3

import os
import time
lists = os.listdir()

for item in lists:
    if item.find('mm')==0:
        name = item 
        print(name)


os.system('katachi.py ihf_result_'+name+' 0 calcall 10')
os.system('fkatachi.py ihf_result_'+name+' 0 10')


