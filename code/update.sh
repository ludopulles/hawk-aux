#!/bin/bash

sage --preparse _attack_Zn.sage
mv _attack_Zn.sage.py _attack_Zn.py

sage --preparse hawk.sage
mv hawk.sage.py hawk.py
