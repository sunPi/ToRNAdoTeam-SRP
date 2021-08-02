#!/bin/bash

# This is a comment!
echo Opening a tunnel with LAMP through spectre2...	# This is a comment, too!
#ssh -X jr429@spectre2.le.ac.uk
#ssh -X jr429@spectre2.le.ac.uk 'ssh -X jr429@login406.lamp.le.ac.uk '
ssh -L 65001:login406.lamp.le.ac.uk:22 saj39@uol-jh.le.ac.uk
