#!/bin/bash
#find E/ -name E_\*.txt ! -name E_0.txt -print;
find V/ -name V_\*.txt ! -name V_0.txt -exec rm {} \;
find E/ -name E_\*.txt ! -name E_0.txt -exec rm {} \;
find E/ -name E_\*.bin ! -name E_0.bin -exec rm {} \;
find G/ -name G_\*.txt -exec rm {} \;
find C/ -name C_\*.txt -exec rm {} \;
find N/ -name N_\*.txt -exec rm {} \;
find P/ -name P_\*.txt -exec rm {} \;
find T/ -name T_\*.txt -exec rm {} \;
find tga/ -name image_\*.tga -exec rm {} \;
find jpg/ -name image_\*.jpg -exec rm {} \;
